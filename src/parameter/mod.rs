//! Structures and traits that can be used to build model parameters for equations of state.

use indexmap::IndexSet;
use ndarray::Array2;
use serde::de::DeserializeOwned;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::BufReader;
use std::path::Path;
use thiserror::Error;

mod chemical_record;
mod identifier;
mod model_record;
mod segment;

pub use chemical_record::ChemicalRecord;
pub use identifier::{Identifier, IdentifierOption};
pub use model_record::{BinaryRecord, FromSegments, FromSegmentsBinary, PureRecord};
pub use segment::SegmentRecord;

/// Constructor methods for parameters.
///
/// By implementing `Parameter` for a type, you define how parameters
/// of an equation of state can be constructed from a sequence of
/// single substance records and possibly binary interaction parameters.
pub trait Parameter
where
    Self: Sized,
{
    type Pure: Clone + DeserializeOwned + Default;
    type IdealGas: Clone + DeserializeOwned + Default;
    type Binary: Clone + DeserializeOwned + Default;

    /// Creates parameters from records for pure substances and possibly binary parameters.
    fn from_records(
        pure_records: Vec<PureRecord<Self::Pure, Self::IdealGas>>,
        binary_records: Array2<Self::Binary>,
    ) -> Self;

    /// Creates parameters for a pure component from a pure record.
    fn new_pure(pure_record: PureRecord<Self::Pure, Self::IdealGas>) -> Self {
        let binary_record = Array2::from_elem([1, 1], Self::Binary::default());
        Self::from_records(vec![pure_record], binary_record)
    }

    /// Creates parameters for a binary system from pure records and an optional
    /// binary interaction parameter.
    fn new_binary(
        pure_records: Vec<PureRecord<Self::Pure, Self::IdealGas>>,
        binary_record: Option<Self::Binary>,
    ) -> Self {
        let binary_record = Array2::from_shape_fn([2, 2], |(i, j)| {
            if i == j {
                Self::Binary::default()
            } else {
                binary_record.clone().unwrap_or_default()
            }
        });
        Self::from_records(pure_records, binary_record)
    }

    /// Return the original pure and binary records that were used to construct the parameters.
    #[allow(clippy::type_complexity)]
    fn records(
        &self,
    ) -> (
        &[PureRecord<Self::Pure, Self::IdealGas>],
        &Array2<Self::Binary>,
    );

    /// Creates parameters from substance information stored in json files.
    fn from_json<P>(
        substances: Vec<&str>,
        file_pure: P,
        file_binary: Option<P>,
        search_option: IdentifierOption,
    ) -> Result<Self, ParameterError>
    where
        P: AsRef<Path>,
    {
        Self::from_multiple_json(&[(substances, file_pure)], file_binary, search_option)
    }

    /// Creates parameters from substance information stored in multiple json files.
    fn from_multiple_json<P>(
        input: &[(Vec<&str>, P)],
        file_binary: Option<P>,
        search_option: IdentifierOption,
    ) -> Result<Self, ParameterError>
    where
        P: AsRef<Path>,
    {
        let mut queried: IndexSet<String> = IndexSet::new();
        let mut record_map: HashMap<String, PureRecord<Self::Pure, Self::IdealGas>> =
            HashMap::new();

        for (substances, file) in input {
            substances.iter().try_for_each(|identifier| {
                match queried.insert(identifier.to_string()) {
                    true => Ok(()),
                    false => Err(ParameterError::IncompatibleParameters(format!(
                        "tried to add substance '{}' to system but it is already present.",
                        identifier
                    ))),
                }
            })?;
            let f = File::open(file)?;
            let reader = BufReader::new(f);

            let pure_records: Vec<PureRecord<Self::Pure, Self::IdealGas>> =
                serde_json::from_reader(reader)?;

            pure_records
                .into_iter()
                .filter_map(|record| {
                    record
                        .identifier
                        .as_string(search_option)
                        .map(|i| (i, record))
                })
                .for_each(|(i, r)| {
                    let _ = record_map.insert(i, r);
                });
        }

        // Compare queried components and available components
        let available: IndexSet<String> = record_map
            .keys()
            .map(|identifier| identifier.to_string())
            .collect();
        if !queried.is_subset(&available) {
            let missing: Vec<String> = queried.difference(&available).cloned().collect();
            let msg = format!("{:?}", missing);
            return Err(ParameterError::ComponentsNotFound(msg));
        };
        let p: Vec<_> = queried
            .iter()
            .filter_map(|identifier| record_map.remove(&identifier.clone()))
            .collect();

        // Read binary records from file if provided
        let binary_map = if let Some(path) = file_binary {
            let file = File::open(path)?;
            let reader = BufReader::new(file);
            let binary_records: Vec<BinaryRecord<Identifier, Self::Binary>> =
                serde_json::from_reader(reader)?;
            binary_records
                .into_iter()
                .filter_map(|br| {
                    let id1 = br.id1.as_string(search_option);
                    let id2 = br.id2.as_string(search_option);
                    id1.and_then(|id1| id2.map(|id2| ((id1, id2), br.model_record)))
                })
                .collect()
        } else {
            HashMap::with_capacity(0)
        };

        let n = p.len();
        let br = Array2::from_shape_fn([n, n], |(i, j)| {
            let id1 = p[i].identifier.as_string(search_option).unwrap();
            let id2 = p[j].identifier.as_string(search_option).unwrap();
            binary_map
                .get(&(id1.clone(), id2.clone()))
                .or_else(|| binary_map.get(&(id2, id1)))
                .cloned()
                .unwrap_or_default()
        });

        Ok(Self::from_records(p, br))
    }

    /// Creates parameters from the molecular structure and segment information.
    ///
    /// The [FromSegments] trait needs to be implemented for both the model record
    /// and the ideal gas record.
    fn from_segments(
        chemical_records: Vec<ChemicalRecord>,
        segment_records: Vec<SegmentRecord<Self::Pure, Self::IdealGas>>,
        binary_segment_records: Option<Vec<BinaryRecord<String, Self::Binary>>>,
    ) -> Result<Self, ParameterError>
    where
        Self::Pure: FromSegments,
        Self::IdealGas: FromSegments,
        Self::Binary: FromSegmentsBinary + Default,
    {
        // update the pure records with model and ideal gas records
        // calculated from the gc method
        let pure_records = chemical_records
            .iter()
            .map(|cr| {
                cr.segment_count(&segment_records)
                    .map(|segments| PureRecord::from_segments(cr.identifier().clone(), segments))
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Map: (id1, id2) -> model_record
        // empty, if no binary segment records are provided
        let binary_map: HashMap<_, _> = binary_segment_records
            .into_iter()
            .flat_map(|seg| seg.into_iter())
            .map(|br| ((br.id1, br.id2), br.model_record))
            .collect();

        // For every component:  map: id -> count
        let segment_id_counts: Vec<_> = chemical_records
            .iter()
            .map(|cr| cr.segment_id_count())
            .collect();

        // full matrix of binary records from the gc method.
        // If a specific segment-segment interaction is not in the binary map,
        // the default value is used.
        let n = pure_records.len();
        let binary_records = Array2::from_shape_fn([n, n], |(i, j)| {
            let mut vec = Vec::new();
            for (id1, &n1) in segment_id_counts[i].iter() {
                for (id2, &n2) in segment_id_counts[j].iter() {
                    let binary = binary_map
                        .get(&(id1.clone(), id2.clone()))
                        .or_else(|| binary_map.get(&(id2.clone(), id1.clone())))
                        .cloned()
                        .unwrap_or_default();
                    vec.push((binary, n1, n2));
                }
            }
            Self::Binary::from_segments_binary(&vec)
        });

        Ok(Self::from_records(pure_records, binary_records))
    }

    /// Creates parameters from segment information stored in json files.
    ///
    /// The [FromSegments] trait needs to be implemented for both the model record
    /// and the ideal gas record.
    fn from_json_segments<P>(
        substances: &[&str],
        file_pure: P,
        file_segments: P,
        file_binary: Option<P>,
        search_option: IdentifierOption,
    ) -> Result<Self, ParameterError>
    where
        P: AsRef<Path>,
        Self::Pure: FromSegments,
        Self::IdealGas: FromSegments,
        Self::Binary: FromSegmentsBinary,
    {
        let queried: IndexSet<String> = substances
            .iter()
            .map(|identifier| identifier.to_string())
            .collect();

        let file = File::open(file_pure)?;
        let reader = BufReader::new(file);
        let chemical_records: Vec<ChemicalRecord> = serde_json::from_reader(reader)?;
        let mut record_map: HashMap<_, _> = chemical_records
            .into_iter()
            .filter_map(|record| {
                record
                    .identifier()
                    .as_string(search_option)
                    .map(|i| (i, record))
            })
            .collect();

        // Compare queried components and available components
        let available: IndexSet<String> = record_map
            .keys()
            .map(|identifier| identifier.to_string())
            .collect();
        if !queried.is_subset(&available) {
            let missing: Vec<String> = queried.difference(&available).cloned().collect();
            let msg = format!("{:?}", missing);
            return Err(ParameterError::ComponentsNotFound(msg));
        };

        // collect all pure records that were queried
        let chemical_records: Vec<_> = queried
            .iter()
            .filter_map(|identifier| record_map.remove(&identifier.clone()))
            .collect();

        // Read segment records
        let segment_records: Vec<SegmentRecord<Self::Pure, Self::IdealGas>> =
            SegmentRecord::from_json(file_segments)?;

        // Read binary records
        let binary_records = file_binary
            .map(|file_binary| {
                let reader = BufReader::new(File::open(file_binary)?);
                let binary_records: Result<
                    Vec<BinaryRecord<String, Self::Binary>>,
                    ParameterError,
                > = Ok(serde_json::from_reader(reader)?);
                binary_records
            })
            .transpose()?;

        Self::from_segments(chemical_records, segment_records, binary_records)
    }

    fn subset(&self, component_list: &[usize]) -> Self {
        let (pure_records, binary_records) = self.records();
        let pure_records = component_list
            .iter()
            .map(|&i| pure_records[i].clone())
            .collect();
        let n = component_list.len();
        let binary_records = Array2::from_shape_fn([n, n], |(i, j)| {
            binary_records[(component_list[i], component_list[j])].clone()
        });

        Self::from_records(pure_records, binary_records)
    }
}

/// Error type for incomplete parameter information and IO problems.
#[derive(Error, Debug)]
pub enum ParameterError {
    #[error(transparent)]
    FileIO(#[from] io::Error),
    #[error(transparent)]
    Serde(#[from] serde_json::Error),
    #[error("The following component(s) were not found: {0}")]
    ComponentsNotFound(String),
    #[error("The identifier '{0}' is not known. ['cas', 'name', 'iupacname', 'smiles', inchi', 'formula']")]
    IdentifierNotFound(String),
    #[error("Information missing.")]
    InsufficientInformation,
    #[error("Incompatible parameters: {0}")]
    IncompatibleParameters(String),
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::joback::JobackRecord;
    use serde::{Deserialize, Serialize};

    #[derive(Debug, Clone, Serialize, Deserialize, Default)]
    struct MyPureModel {
        a: f64,
    }

    struct MyParameter {
        pure_records: Vec<PureRecord<MyPureModel, JobackRecord>>,
        binary_records: Array2<f64>,
    }

    impl Parameter for MyParameter {
        type Pure = MyPureModel;
        type IdealGas = JobackRecord;
        type Binary = f64;
        fn from_records(
            pure_records: Vec<PureRecord<MyPureModel, JobackRecord>>,
            binary_records: Array2<f64>,
        ) -> Self {
            Self {
                pure_records,
                binary_records,
            }
        }

        fn records(&self) -> (&[PureRecord<MyPureModel, JobackRecord>], &Array2<f64>) {
            (&self.pure_records, &self.binary_records)
        }
    }

    #[test]
    fn from_records() {
        let r = r#"
        [
            {
                "identifier": {
                    "cas": "123-4-5"
                },
                "molarweight": 16.0426,
                "model_record": {
                    "a": 0.1
                }
            }
        ]
        "#;
        let records: Vec<PureRecord<MyPureModel, JobackRecord>> =
            serde_json::from_str(r).expect("Unable to parse json.");
        let p = MyParameter::from_records(records, Array2::zeros((1, 1)));
        assert_eq!(p.pure_records[0].identifier.cas, "123-4-5")
    }
}
