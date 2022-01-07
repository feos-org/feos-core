use super::identifier::Identifier;
use super::segment::SegmentRecord;
use super::ParameterError;
use either::Either;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

// Auxiliary structure used to deserialize chemical records without bond information.
#[derive(Deserialize)]
#[serde(untagged)]
enum ChemicalRecordJSON {
    List {
        identifier: Identifier,
        segments: Vec<String>,
        bonds: Option<Vec<[usize; 2]>>,
    },
    Count {
        identifier: Identifier,
        segments: HashMap<String, f64>,
        bonds: Option<HashMap<[String; 2], f64>>,
    },
}

/// Chemical information of a substance.
#[derive(Deserialize, Serialize, Debug, Clone)]
#[serde(from = "ChemicalRecordJSON")]
#[serde(untagged)]
pub enum ChemicalRecord {
    List {
        identifier: Identifier,
        segments: Vec<String>,
        bonds: Vec<[usize; 2]>,
    },
    Count {
        identifier: Identifier,
        segments: HashMap<String, f64>,
        bonds: HashMap<[String; 2], f64>,
    },
}

impl From<ChemicalRecordJSON> for ChemicalRecord {
    fn from(record: ChemicalRecordJSON) -> Self {
        match record {
            ChemicalRecordJSON::List {
                identifier,
                segments,
                bonds,
            } => Self::new(identifier, segments, bonds),
            ChemicalRecordJSON::Count {
                identifier,
                segments,
                bonds,
            } => Self::new_count(identifier, segments, bonds),
        }
    }
}

impl ChemicalRecord {
    /// Create a new `ChemicalRecord`.
    ///
    /// If no bonds are given, the molecule is assumed to be linear.
    pub fn new(
        identifier: Identifier,
        segments: Vec<String>,
        bonds: Option<Vec<[usize; 2]>>,
    ) -> ChemicalRecord {
        let bonds = bonds.unwrap_or_else(|| {
            (0..segments.len() - 1)
                .zip(1..segments.len())
                .map(|x| [x.0, x.1])
                .collect()
        });
        Self::List {
            identifier,
            segments,
            bonds,
        }
    }

    /// Create a new `ChemicalRecord` from a segment count.
    pub fn new_count(
        identifier: Identifier,
        segments: HashMap<String, f64>,
        bonds: Option<HashMap<[String; 2], f64>>,
    ) -> ChemicalRecord {
        let bonds = bonds.unwrap_or_default();
        Self::Count {
            identifier,
            segments,
            bonds,
        }
    }

    pub fn segments(&self) -> Either<&Vec<String>, &HashMap<String, f64>> {
        match self {
            Self::List {
                identifier: _,
                segments,
                bonds: _,
            } => Either::Left(segments),
            Self::Count {
                identifier: _,
                segments,
                bonds: _,
            } => Either::Right(segments),
        }
    }

    pub fn identifier(&self) -> &Identifier {
        match self {
            Self::List {
                identifier,
                segments: _,
                bonds: _,
            } => identifier,
            Self::Count {
                identifier,
                segments: _,
                bonds: _,
            } => identifier,
        }
    }

    /// Count the number of occurences of each individual segment in the
    /// chemical record.
    ///
    /// The map contains the segment record as key and the count as (float) value.
    pub fn segment_count<M: Clone, I: Clone>(
        &self,
        segment_records: &[SegmentRecord<M, I>],
    ) -> Result<HashMap<SegmentRecord<M, I>, f64>, ParameterError> {
        let segment_map = self.segment_map(segment_records)?;
        Ok(match self.segments() {
            Either::Left(segments) => {
                let mut counts = HashMap::with_capacity(segments.len());
                for si in segments {
                    let key = segment_map.get(si).unwrap().clone();
                    let entry = counts.entry(key).or_insert(0.0);
                    *entry += 1.0;
                }
                counts
            }
            Either::Right(segments) => segments
                .iter()
                .map(|(id, &c)| (segment_map[id].clone(), c))
                .collect(),
        })
    }

    /// Count the number of occurences of each individual segment identifier in the
    /// chemical record.
    ///
    /// The map contains the segment identifier as key and the count as (float) value.
    pub fn segment_id_count(&self) -> Cow<HashMap<String, f64>> {
        match self.segments() {
            Either::Left(segments) => {
                let mut counts = HashMap::with_capacity(segments.len());
                for si in segments {
                    let entry = counts.entry(si.clone()).or_insert(0.0);
                    *entry += 1.0;
                }
                Cow::Owned(counts)
            }
            Either::Right(segments) => Cow::Borrowed(segments),
        }
    }

    /// Count the number of occurences of each individual segment identifier and each
    /// pair of identifiers for every bond in the chemical record.
    #[allow(clippy::type_complexity)]
    pub fn segment_and_bond_count(
        &self,
    ) -> (Cow<HashMap<String, f64>>, Cow<HashMap<[String; 2], f64>>) {
        match self {
            Self::List {
                identifier: _,
                segments,
                bonds,
            } => {
                let mut segment_counts = HashMap::with_capacity(segments.len());
                for si in segments {
                    let entry = segment_counts.entry(si.clone()).or_insert(0.0);
                    *entry += 1.0;
                }

                let mut bond_counts = HashMap::new();
                for b in bonds {
                    let s1 = segments[b[0]].clone();
                    let s2 = segments[b[1]].clone();
                    let indices = if s1 > s2 { [s2, s1] } else { [s1, s2] };
                    let entry = bond_counts.entry(indices).or_insert(0.0);
                    *entry += 1.0;
                }
                (Cow::Owned(segment_counts), Cow::Owned(bond_counts))
            }
            Self::Count {
                identifier: _,
                segments,
                bonds,
            } => (Cow::Borrowed(segments), Cow::Borrowed(bonds)),
        }
    }

    /// Return the full segment and bond information for the molecule
    /// if possible.
    #[allow(clippy::type_complexity)]
    pub fn segment_and_bond_list(
        &self,
    ) -> Result<(&Vec<String>, &Vec<[usize; 2]>), ParameterError> {
        match self {
            Self::List {
                identifier: _,
                segments,
                bonds,
            } => Ok((segments, bonds)),
            Self::Count {
                identifier: _,
                segments: _,
                bonds: _,
            } => Err(ParameterError::IncompatibleParameters(
                "No detailed structural information available.".into(),
            )),
        }
    }

    /// Build a HashMap from SegmentRecords for the segments.
    ///
    /// The map contains the segment identifier (String) as key
    /// and the SegmentRecord as value.
    pub fn segment_map<M: Clone, I: Clone>(
        &self,
        segment_records: &[SegmentRecord<M, I>],
    ) -> Result<HashMap<String, SegmentRecord<M, I>>, ParameterError> {
        let queried: HashSet<_> = match self.segments() {
            Either::Left(segments) => segments.iter().cloned().collect(),
            Either::Right(segments) => segments.keys().cloned().collect(),
        };
        let mut segments: HashMap<String, SegmentRecord<M, I>> = segment_records
            .iter()
            .map(|r| (r.identifier.clone(), r.clone()))
            .collect();
        let available = segments.keys().cloned().collect();
        if !queried.is_subset(&available) {
            let missing: Vec<String> = queried.difference(&available).cloned().collect();
            let msg = format!("{:?}", missing);
            return Err(ParameterError::ComponentsNotFound(msg));
        };
        Ok(queried
            .iter()
            .map(|s| segments.remove_entry(s).unwrap())
            .collect())
    }
}

impl std::fmt::Display for ChemicalRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::List {
                identifier,
                segments,
                bonds,
            } => {
                write!(f, "ChemicalRecord(")?;
                write!(f, "\n\tidentifier={},", identifier)?;
                write!(f, "\n\tsegments={:?},", segments)?;
                write!(f, "\n\tbonds={:?}\n)", bonds)
            }
            Self::Count {
                identifier,
                segments,
                bonds,
            } => {
                write!(f, "ChemicalRecord(")?;
                write!(f, "\n\tidentifier={},", identifier)?;
                write!(f, "\n\tsegments={:?},", segments)?;
                write!(f, "\n\tbonds={:?}\n)", bonds)
            }
        }
    }
}
