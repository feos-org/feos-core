use super::segment::SegmentRecord;
use super::ParameterError;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

// Auxiliary structure used to deserialize chemical records without bond information.
#[derive(Deserialize)]
struct ChemicalRecordJSON {
    pub segments: Vec<String>,
    pub bonds: Option<Vec<[usize; 2]>>,
}

/// Chemical information of a substance.
#[derive(Deserialize, Serialize, Debug, Clone)]
#[serde(from = "ChemicalRecordJSON")]
pub struct ChemicalRecord {
    pub segments: Vec<String>,
    pub bonds: Vec<[usize; 2]>,
}

impl From<ChemicalRecordJSON> for ChemicalRecord {
    fn from(record: ChemicalRecordJSON) -> Self {
        Self::new(record.segments, record.bonds)
    }
}

impl ChemicalRecord {
    /// Create a new `ChemicalRecord`.
    ///
    /// If no bonds are given, the molecule is assumed to be linear.
    pub fn new(segments: Vec<String>, bonds: Option<Vec<[usize; 2]>>) -> ChemicalRecord {
        let bonds = bonds.unwrap_or_else(|| {
            (0..segments.len() - 1)
                .zip(1..segments.len())
                .map(|x| [x.0, x.1])
                .collect()
        });
        Self { segments, bonds }
    }

    /// Count the number of occurences of each individual segment in the
    /// chemical record.
    ///
    /// The map contains the segment record as key and the count as (float) value.
    pub fn segment_count<M: Clone, I: Clone>(
        &self,
        segment_records: &[SegmentRecord<M, I>],
    ) -> Result<HashMap<SegmentRecord<M, I>, f64>, ParameterError> {
        let mut counts = HashMap::with_capacity(self.segments.len());
        let segments = self.segment_map(segment_records)?;
        for si in &self.segments {
            let key = segments.get(si).unwrap().clone();
            let entry = counts.entry(key).or_insert(0.0);
            *entry += 1.0;
        }
        Ok(counts)
    }

    /// Count the number of occurences of each individual segment identifier in the
    /// chemical record.
    ///
    /// The map contains the segment identifier as key and the count as (float) value.
    pub fn segment_id_count(&self) -> HashMap<String, f64> {
        let mut counts = HashMap::with_capacity(self.segments.len());
        for si in &self.segments {
            let entry = counts.entry(si.clone()).or_insert(0.0);
            *entry += 1.0;
        }
        counts
    }

    /// Build a HashMap from SegmentRecords for the segments.
    ///
    /// The map contains the segment identifier (String) as key
    /// and the SegmentRecord as value.
    pub fn segment_map<M: Clone, I: Clone>(
        &self,
        segment_records: &[SegmentRecord<M, I>],
    ) -> Result<HashMap<String, SegmentRecord<M, I>>, ParameterError> {
        let queried: HashSet<String> = self.segments.iter().cloned().collect();
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
            .intersection(&available)
            .map(|s| segments.remove_entry(s).unwrap())
            .collect())
    }
}

impl std::fmt::Display for ChemicalRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ChemicalRecord(segments={:?}", self.segments)?;
        write!(f, ", bonds={:?})", self.bonds)
    }
}
