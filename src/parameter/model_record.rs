use super::chemical_record::ChemicalRecord;
use super::identifier::Identifier;
use super::ParameterError;
use serde::{Deserialize, Serialize};

/// A collection of parameters of a pure substance.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct PureRecord<M, I> {
    pub identifier: Identifier,
    pub molarweight: f64,
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chemical_record: Option<ChemicalRecord>,
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub model_record: Option<M>,
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ideal_gas_record: Option<I>,
}

impl<M, I> PureRecord<M, I> {
    /// Create a new `PureRecord`.
    pub fn new(
        identifier: Identifier,
        molarweight: f64,
        chemical_record: Option<ChemicalRecord>,
        model_record: Option<M>,
        ideal_gas_record: Option<I>,
    ) -> Self {
        Self {
            identifier,
            molarweight,
            chemical_record,
            model_record,
            ideal_gas_record,
        }
    }
}

impl<M, I> std::fmt::Display for PureRecord<M, I>
where
    M: std::fmt::Display,
    I: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PureRecord(")?;
        write!(f, "\n\tidentifier={},", self.identifier)?;
        write!(f, "\n\tmolarweight={},", self.molarweight)?;
        if let Some(m) = self.chemical_record.as_ref() {
            write!(f, "\n\tchemical_record={},", m)?;
        }
        if let Some(m) = self.model_record.as_ref() {
            write!(f, "\n\tmodel_record={},", m)?;
        }
        if let Some(i) = self.ideal_gas_record.as_ref() {
            write!(f, "\n\tideal_gas_record={},", i)?;
        }
        write!(f, "\n)")
    }
}

/// Empty record type to be used in models (e.g. group contribution
/// methods) that do not require a model or ideal gas record.
#[derive(Serialize, Deserialize, Clone, Debug, Default)]
pub struct NoRecord;

impl std::fmt::Display for NoRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "")
    }
}

/// Convenience type for models in which molecules are purely describes
/// by their chemical record.
pub type GroupContributionRecord = PureRecord<NoRecord, NoRecord>;

/// Trait for models that implement a homosegmented group contribution
/// method
pub trait FromSegments: Clone {
    /// Type of the binary interaction parameter(s) used in this model.
    type Binary;

    /// Constructs the record from a list of segment records with their
    /// number of occurences and possibly binary interaction parameters.
    fn from_segments(
        segments: &[(Self, f64)],
        binary_records: Option<&[BinaryRecord<String, Self::Binary>]>,
    ) -> Result<Self, ParameterError>;
}

/// A collection of parameters that model interactions between two
/// substances or segments.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct BinaryRecord<I, B> {
    /// Identifier of the first component
    pub id1: I,
    /// Identifier of the second component
    pub id2: I,
    /// Binary interaction parameter(s)
    pub model_record: B,
}

impl<I, B> BinaryRecord<I, B> {
    /// Crates a new `BinaryRecord`.
    pub fn new(id1: I, id2: I, model_record: B) -> Self {
        Self {
            id1,
            id2,
            model_record,
        }
    }
}

impl<I, B> std::fmt::Display for BinaryRecord<I, B>
where
    I: std::fmt::Display,
    B: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BinaryRecord(")?;
        write!(f, "\n\tid1={},", self.id1)?;
        write!(f, "\n\tid2={},", self.id2)?;
        write!(f, "\n\tmodel_record={},", self.model_record)?;
        write!(f, "\n)")
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::joback::JobackRecord;

    #[derive(Serialize, Deserialize, Debug, Default, Clone)]
    struct TestModelRecordSegments {
        a: f64,
    }

    #[test]
    fn deserialize() {
        let r = r#"
        {
            "identifier": {
                "cas": "123-4-5"
            },
            "molarweight": 16.0426,
            "model_record": {
                "a": 0.1
            }
        }
        "#;
        let record: PureRecord<TestModelRecordSegments, JobackRecord> =
            serde_json::from_str(r).expect("Unable to parse json.");
        assert_eq!(record.identifier.cas, "123-4-5")
    }

    #[test]
    fn deserialize_list() {
        let r = r#"
        [
            {
                "identifier": {
                    "cas": "1"
                },
                "molarweight": 1.0,
                "model_record": {
                    "a": 1.0
                }
            },
            {
                "identifier": {
                    "cas": "2"
                },
                "molarweight": 2.0,
                "model_record": {
                    "a": 2.0
                }
            }
        ]"#;
        let records: Vec<PureRecord<TestModelRecordSegments, JobackRecord>> =
            serde_json::from_str(r).expect("Unable to parse json.");
        assert_eq!(records[0].identifier.cas, "1");
        assert_eq!(records[1].identifier.cas, "2")
    }
}
