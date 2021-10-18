//! Structures and traits that can be used to build model parameters for equations of state.

use indexmap::IndexSet;
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
pub use model_record::{BinaryRecord, FromSegments, GroupContributionRecord, NoRecord, PureRecord};
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
    type Binary: DeserializeOwned + Default;

    /// Creates parameters from records for pure substances and possibly binary parameters.
    fn from_records(
        pure_records: Vec<PureRecord<Self::Pure, Self::IdealGas>>,
        binary_records: Option<Vec<BinaryRecord<Identifier, Self::Binary>>>,
    ) -> Result<Self, ParameterError>;

    /// Creates parameters from substance information stored in json files.
    fn from_json<P>(
        substances: &[&str],
        file_pure: P,
        file_binary: Option<P>,
        search_option: IdentifierOption,
    ) -> Result<Self, ParameterError>
    where
        P: AsRef<Path>,
    {
        let queried: IndexSet<String> = substances
            .iter()
            .map(|identifier| identifier.to_string())
            .collect();
        let file = File::open(file_pure)?;
        let reader = BufReader::new(file);

        let pure_records: Vec<PureRecord<Self::Pure, Self::IdealGas>> =
            serde_json::from_reader(reader)?;
        let mut record_map: HashMap<_, _> = pure_records
            .into_iter()
            .filter_map(|record| {
                record
                    .identifier
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
        let p = queried
            .iter()
            .filter_map(|identifier| record_map.remove(&identifier.clone()))
            .collect();

        // Todo: maybe change into buffer
        let bp = if let Some(path) = file_binary {
            let file = File::open(path)?;
            let reader = BufReader::new(file);
            let binary_records: Vec<BinaryRecord<Identifier, Self::Binary>> =
                serde_json::from_reader(reader)?;
            let brs: Vec<BinaryRecord<Identifier, Self::Binary>> = binary_records
                .into_iter()
                .filter(|record| {
                    let id1 = record.id1.as_string(search_option);
                    let id2 = record.id2.as_string(search_option);
                    if let (Some(i1), Some(i2)) = (id1, id2) {
                        let mut s = IndexSet::with_capacity(2);
                        s.insert(i1);
                        s.insert(i2);
                        s.is_subset(&queried)
                    } else {
                        false
                    }
                })
                .collect();
            Some(brs)
        } else {
            None
        };
        Self::from_records(p, bp)
    }

    /// Creates parameters from the molecular structure and segment information.
    ///
    /// The [FromSegments] trait needs to be implemented for both the model record
    /// and the ideal gas record.
    fn from_segments(
        pure_records: Vec<PureRecord<Self::Pure, Self::IdealGas>>,
        segment_records: Vec<SegmentRecord<Self::Pure, Self::IdealGas>>,
        binary_segment_records: Option<Vec<BinaryRecord<String, Self::Binary>>>,
    ) -> Result<Self, ParameterError>
    where
        Self::Pure: FromSegments<Binary = Self::Binary>,
        Self::IdealGas: FromSegments,
    {
        let mut records: Vec<_> = Vec::with_capacity(pure_records.len());
        for pr in pure_records {
            let chemical_record = pr
                .chemical_record
                .clone()
                .ok_or(ParameterError::InsufficientInformation)?;
            let molarweight = chemical_record.molarweight_from_segments(&segment_records)?;
            let segment_count = chemical_record.segment_count(&segment_records)?;

            let model_segments: Vec<_> = segment_count
                .iter()
                .map(|(s, &n)| (s.model_record.clone(), n))
                .collect();
            let model_record =
                Self::Pure::from_segments(&model_segments, binary_segment_records.as_deref())?;

            let ideal_gas_segments: Option<Vec<_>> = segment_count
                .iter()
                .map(|(s, &n)| s.ideal_gas_record.clone().map(|ig| (ig, n)))
                .collect();
            let ideal_gas_record = ideal_gas_segments
                .as_ref()
                .map(|s| Self::IdealGas::from_segments(s, None))
                .transpose()?;

            records.push(PureRecord::new(
                pr.identifier.clone(),
                molarweight,
                Some(chemical_record),
                Some(model_record),
                ideal_gas_record,
            ))
        }
        Self::from_records(records, None)
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
        Self::Pure: FromSegments<Binary = Self::Binary>,
        Self::IdealGas: FromSegments,
    {
        let queried: IndexSet<String> = substances
            .iter()
            .map(|identifier| identifier.to_string())
            .collect();

        let file = File::open(file_pure)?;
        let reader = BufReader::new(file);
        let pure_records: Vec<PureRecord<Self::Pure, Self::IdealGas>> =
            serde_json::from_reader(reader)?;
        let mut record_map: HashMap<_, _> = pure_records
            .into_iter()
            .filter_map(|record| {
                record
                    .identifier
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
        let pure_records: Vec<_> = queried
            .iter()
            .filter_map(|identifier| record_map.remove(&identifier.clone()))
            .collect();

        // Read segment records
        let segment_records: Vec<SegmentRecord<Self::Pure, Self::IdealGas>> =
            serde_json::from_reader(BufReader::new(File::open(file_segments)?))?;

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

        Self::from_segments(pure_records, segment_records, binary_records)
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
    #[error("Building model parameter from homo segments failed: {0}")]
    HomoGc(String),
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

    #[derive(Debug, Clone, Serialize, Deserialize, Default)]
    struct MyBinaryModel {
        a: f64,
    }

    struct MyParameter {
        pure: Vec<PureRecord<MyPureModel, JobackRecord>>,
    }

    impl Parameter for MyParameter {
        type Pure = MyPureModel;
        type IdealGas = JobackRecord;
        type Binary = MyBinaryModel;
        fn from_records(
            pure_records: Vec<PureRecord<MyPureModel, JobackRecord>>,
            _: Option<Vec<BinaryRecord<Identifier, MyBinaryModel>>>,
        ) -> Result<Self, ParameterError>
where {
            Ok(Self {
                pure: pure_records.to_vec(),
            })
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
        let p = MyParameter::from_records(records, None).unwrap();
        assert_eq!(p.pure[0].identifier.cas, "123-4-5")
    }
}
