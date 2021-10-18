use super::ParameterError;
use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::hash::{Hash, Hasher};

/// Possible variants to identify a substance.
#[repr(C)]
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum IdentifierOption {
    Cas,
    Name,
    IupacName,
    Smiles,
    Inchi,
    Formula,
}

impl TryFrom<&str> for IdentifierOption {
    type Error = ParameterError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let lower = s.to_lowercase();
        match lower.as_str() {
            "cas" => Ok(Self::Cas),
            "name" => Ok(Self::Name),
            "iupacname" => Ok(Self::IupacName),
            "smiles" => Ok(Self::Smiles),
            "inchi" => Ok(Self::Inchi),
            "formula" => Ok(Self::Formula),
            _ => Err(ParameterError::IdentifierNotFound(s.to_owned())),
        }
    }
}

/// A collection of identifiers for a chemical structure or substance.
#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct Identifier {
    /// CAS number
    pub cas: String,
    /// Commonly used english name
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    /// IUPAC name
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub iupac_name: Option<String>,
    /// SMILES key
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub smiles: Option<String>,
    /// InchI key
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inchi: Option<String>,
    /// Chemical formula
    #[serde(default)]
    #[serde(skip_serializing_if = "Option::is_none")]
    pub formula: Option<String>,
}

impl Identifier {
    /// Create a new identifier.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use feos_core::parameter::Identifier;
    /// let methanol = Identifier::new(
    ///     "67-56-1",
    ///     Some("methanol"),
    ///     Some("methanol"),
    ///     Some("CO"),
    ///     Some("InChI=1S/CH4O/c1-2/h2H,1H3"),
    ///     Some("CH4O")
    /// );
    pub fn new(
        cas: &str,
        name: Option<&str>,
        iupac_name: Option<&str>,
        smiles: Option<&str>,
        inchi: Option<&str>,
        formula: Option<&str>,
    ) -> Identifier {
        Identifier {
            cas: cas.to_owned(),
            name: name.map(Into::into),
            iupac_name: iupac_name.map(Into::into),
            smiles: smiles.map(Into::into),
            inchi: inchi.map(Into::into),
            formula: formula.map(Into::into),
        }
    }

    pub fn as_string(&self, option: IdentifierOption) -> Option<String> {
        match option {
            IdentifierOption::Cas => Some(self.cas.clone()),
            IdentifierOption::Name => self.name.clone(),
            IdentifierOption::IupacName => self.iupac_name.clone(),
            IdentifierOption::Smiles => self.smiles.clone(),
            IdentifierOption::Inchi => self.inchi.clone(),
            IdentifierOption::Formula => self.formula.clone(),
        }
    }
}

impl std::fmt::Display for Identifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Identifier(cas={}", self.cas)?;
        if let Some(n) = &self.name {
            write!(f, ", name={}", n)?;
        }
        if let Some(n) = &self.iupac_name {
            write!(f, ", iupac_name={}", n)?;
        }
        if let Some(n) = &self.smiles {
            write!(f, ", smiles={}", n)?;
        }
        if let Some(n) = &self.inchi {
            write!(f, ", inchi={}", n)?;
        }
        if let Some(n) = &self.formula {
            write!(f, ", formula={}", n)?;
        }
        write!(f, ")")
    }
}

impl PartialEq for Identifier {
    fn eq(&self, other: &Self) -> bool {
        self.cas == other.cas
    }
}
impl Eq for Identifier {}

impl Hash for Identifier {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cas.hash(state);
    }
}
