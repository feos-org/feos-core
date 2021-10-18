use crate::impl_json_handling;
use crate::parameter::{BinaryRecord, ChemicalRecord, Identifier, NoRecord, ParameterError};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

impl From<ParameterError> for PyErr {
    fn from(e: ParameterError) -> PyErr {
        PyRuntimeError::new_err(e.to_string())
    }
}

/// Create an identifier for a pure substance.
///
/// Parameters
/// ----------
/// cas : str
///     CAS number.
/// name : str, optional
///     name
/// iupac_name : str, optional
///     Iupac name.
/// smiles : str, optional
///     Canonical SMILES
/// inchi : str, optional
///     Inchi number
/// formula : str, optional
///     Molecular formula.
///
/// Returns
/// -------
/// Identifier
#[pyclass(name = "Identifier")]
#[derive(Clone)]
#[pyo3(text_signature = "(cas, name=None, iupac_name=None, smiles=None, inchi=None, formula=None)")]
pub struct PyIdentifier(pub Identifier);

#[pymethods]
impl PyIdentifier {
    #[new]
    fn new(
        cas: &str,
        name: Option<&str>,
        iupac_name: Option<&str>,
        smiles: Option<&str>,
        inchi: Option<&str>,
        formula: Option<&str>,
    ) -> Self {
        Self(Identifier::new(
            cas, name, iupac_name, smiles, inchi, formula,
        ))
    }

    #[getter]
    fn get_cas(&self) -> String {
        self.0.cas.clone()
    }

    #[setter]
    fn set_cas(&mut self, cas: &str) {
        self.0.cas = cas.to_string();
    }

    #[getter]
    fn get_name(&self) -> Option<String> {
        self.0.name.clone()
    }

    #[setter]
    fn set_name(&mut self, name: &str) {
        self.0.name = Some(name.to_string());
    }

    #[getter]
    fn get_iupac_name(&self) -> Option<String> {
        self.0.iupac_name.clone()
    }

    #[setter]
    fn set_iupac_name(&mut self, iupac_name: &str) {
        self.0.iupac_name = Some(iupac_name.to_string());
    }

    #[getter]
    fn get_smiles(&self) -> Option<String> {
        self.0.smiles.clone()
    }

    #[setter]
    fn set_smiles(&mut self, smiles: &str) {
        self.0.smiles = Some(smiles.to_string());
    }

    #[getter]
    fn get_inchi(&self) -> Option<String> {
        self.0.inchi.clone()
    }

    #[setter]
    fn set_inchi(&mut self, inchi: &str) {
        self.0.inchi = Some(inchi.to_string());
    }

    #[getter]
    fn get_formula(&self) -> Option<String> {
        self.0.formula.clone()
    }

    #[setter]
    fn set_formula(&mut self, formula: &str) {
        self.0.formula = Some(formula.to_string());
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyIdentifier {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

impl_json_handling!(PyIdentifier);

/// Create a chemical record for a pure substance.
///
/// Parameters
/// ----------
/// segments : [str]
///     List of segments, that the molecule consists of.
/// bonds : [[int, int]], optional
///     List of bonds with the indices starting at 0 and
///     referring to the list of segments passed as first
///     argument. If no bonds are specified, the molecule
///     is assumed to be linear.
///
/// Returns
/// -------
/// ChemicalRecord
#[pyclass(name = "ChemicalRecord")]
#[derive(Clone)]
#[pyo3(text_signature = "(segments, bonds=None)")]
pub struct PyChemicalRecord(pub ChemicalRecord);

#[pymethods]
impl PyChemicalRecord {
    #[new]
    fn new(segments: Vec<String>, bonds: Option<Vec<[usize; 2]>>) -> Self {
        Self(ChemicalRecord::new(segments, bonds))
    }

    #[getter]
    fn get_segments(&self) -> Vec<String> {
        self.0.segments.clone()
    }

    #[setter]
    fn set_segments(&mut self, segments: Vec<String>) {
        self.0.segments = segments
    }

    #[getter]
    fn get_bonds(&self) -> Vec<[usize; 2]> {
        self.0.bonds.clone()
    }

    #[setter]
    fn set_bonds(&mut self, bonds: Vec<[usize; 2]>) {
        self.0.bonds = bonds
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyChemicalRecord {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

impl_json_handling!(PyChemicalRecord);

/// Create a record for a binary interaction parameter.
///
/// Parameters
/// ----------
/// id1 : Identifier
///     The identifier of the first component.
/// id2 : Identifier
///     The identifier of the second component.
/// model_record : float
///     The binary interaction parameter.
///
/// Returns
/// -------
/// BinaryRecord
#[pyclass(name = "BinaryRecord")]
#[pyo3(text_signature = "(id1, id2, model_record)")]
#[derive(Clone)]
pub struct PyBinaryRecord(pub BinaryRecord<Identifier, f64>);

#[pymethods]
impl PyBinaryRecord {
    #[new]
    fn new(id1: PyIdentifier, id2: PyIdentifier, model_record: f64) -> PyResult<Self> {
        Ok(Self(BinaryRecord::new(id1.0, id2.0, model_record)))
    }

    #[getter]
    fn get_id1(&self) -> PyIdentifier {
        PyIdentifier(self.0.id1.clone())
    }

    #[setter]
    fn set_id1(&mut self, id1: PyIdentifier) {
        self.0.id1 = id1.0;
    }

    #[getter]
    fn get_id2(&self) -> PyIdentifier {
        PyIdentifier(self.0.id2.clone())
    }

    #[setter]
    fn set_id2(&mut self, id2: PyIdentifier) {
        self.0.id2 = id2.0;
    }

    #[getter]
    fn get_model_record(&self) -> f64 {
        self.0.model_record
    }

    #[setter]
    fn set_model_record(&mut self, model_record: f64) {
        self.0.model_record = model_record;
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyBinaryRecord {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

impl_json_handling!(PyBinaryRecord);

/// Create a record for a binary segment interaction parameter.
///
/// Parameters
/// ----------
/// id1 : str
///     The identifier of the first segment.
/// id2 : str
///     The identifier of the second segment.
/// model_record : float
///     The binary segment interaction parameter.
///
/// Returns
/// -------
/// BinarySegmentRecord
#[pyclass(name = "BinarySegmentRecord")]
#[pyo3(text_signature = "(id1, id2, model_record)")]
#[derive(Clone)]
pub struct PyBinarySegmentRecord(pub BinaryRecord<String, f64>);

#[pymethods]
impl PyBinarySegmentRecord {
    #[new]
    fn new(id1: String, id2: String, model_record: f64) -> PyResult<Self> {
        Ok(Self(BinaryRecord::new(id1, id2, model_record)))
    }

    #[getter]
    fn get_id1(&self) -> String {
        self.0.id1.clone()
    }

    #[setter]
    fn set_id1(&mut self, id1: String) {
        self.0.id1 = id1;
    }

    #[getter]
    fn get_id2(&self) -> String {
        self.0.id2.clone()
    }

    #[setter]
    fn set_id2(&mut self, id2: String) {
        self.0.id2 = id2;
    }

    #[getter]
    fn get_model_record(&self) -> f64 {
        self.0.model_record
    }

    #[setter]
    fn set_model_record(&mut self, model_record: f64) {
        self.0.model_record = model_record;
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyBinarySegmentRecord {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

impl_json_handling!(PyBinarySegmentRecord);

#[pyclass(name = "NoRecord")]
#[derive(Clone)]
pub struct PyNoRecord(pub NoRecord);

#[macro_export]
macro_rules! impl_pure_record {
    ($model_record:ident, $py_model_record:ident, $ideal_gas_record:ident, $py_ideal_gas_record:ident) => {
        /// All information required to characterize a pure component.
        ///
        /// Parameters
        /// ----------
        /// identifier : Identifier
        ///     The identifier of the pure component.
        /// molarweight : float
        ///     The molar weight (in g/mol) of the pure component.
        /// model_record : ModelRecord, optional
        ///     The pure component model parameters.
        /// ideal_gas_record: IdealGasRecord, optional
        ///     The pure component parameters for the ideal gas model.
        /// chemical_record: ChemicalRecord, optional
        ///     The chemical record of the pure component.
        ///
        /// Returns
        /// -------
        /// PureRecord
        #[pyclass(name = "PureRecord")]
        #[pyo3(text_signature = "(identifier, molarweight, model_record=None, ideal_gas_record=None, chemical_record=None)")]
        #[derive(Clone)]
        pub struct PyPureRecord(PureRecord<$model_record, $ideal_gas_record>);

        #[pymethods]
        impl PyPureRecord {
            #[new]
            fn new(
                identifier: PyIdentifier,
                molarweight: f64,
                model_record: Option<$py_model_record>,
                ideal_gas_record: Option<$py_ideal_gas_record>,
                chemical_record: Option<PyChemicalRecord>,
            ) -> PyResult<Self> {
                Ok(Self(PureRecord::new(
                    identifier.0,
                    molarweight,
                    chemical_record.map(|cr| cr.0),
                    model_record.map(|mr| mr.0),
                    ideal_gas_record.map(|ig| ig.0),
                )))
            }

            #[getter]
            fn get_identifier(&self) -> PyIdentifier {
                PyIdentifier(self.0.identifier.clone())
            }

            #[setter]
            fn set_identifier(&mut self, identifier: PyIdentifier) {
                self.0.identifier = identifier.0;
            }

            #[getter]
            fn get_molarweight(&self) -> f64 {
                self.0.molarweight
            }

            #[setter]
            fn set_molarweight(&mut self, molarweight: f64) {
                self.0.molarweight = molarweight;
            }

            #[getter]
            fn get_model_record(&self) -> Option<$py_model_record> {
                self.0.model_record.clone().map($py_model_record)
            }

            #[setter]
            fn set_model_record(&mut self, model_record: $py_model_record) {
                self.0.model_record = Some(model_record.0);
            }

            #[getter]
            fn get_ideal_gas_record(&self) -> Option<$py_ideal_gas_record> {
                self.0.ideal_gas_record.clone().map($py_ideal_gas_record)
            }

            #[setter]
            fn set_ideal_gas_record(&mut self, ideal_gas_record: $py_ideal_gas_record) {
                self.0.ideal_gas_record = Some(ideal_gas_record.0);
            }

            #[getter]
            fn get_chemical_record(&self) -> Option<PyChemicalRecord> {
                self.0.chemical_record.clone().map(PyChemicalRecord)
            }

            #[setter]
            fn set_chemical_record(&mut self, chemical_record: PyChemicalRecord) {
                self.0.chemical_record = Some(chemical_record.0);
            }
        }

        #[pyproto]
        impl pyo3::class::basic::PyObjectProtocol for PyPureRecord {
            fn __repr__(&self) -> PyResult<String> {
                Ok(self.0.to_string())
            }
        }

        impl_json_handling!(PyPureRecord);
    };
}

#[macro_export]
macro_rules! impl_segment_record {
    ($model_record:ident, $py_model_record:ident, $ideal_gas_record:ident, $py_ideal_gas_record:ident) => {
        /// All information required to characterize a single segment.
        ///
        /// Parameters
        /// ----------
        /// identifier : str
        ///     The identifier of the segment.
        /// molarweight : float
        ///     The molar weight (in g/mol) of the segment.
        /// model_record : ModelRecord
        ///     The segment model parameters.
        /// ideal_gas_record: IdealGasRecord, optional
        ///     The seg,ment ideal gas parameters.
        ///
        /// Returns
        /// -------
        /// SegmentRecord
        #[pyclass(name = "SegmentRecord")]
        #[pyo3(text_signature = "(identifier, molarweight, model_record, ideal_gas_record=None)")]
        #[derive(Clone)]
        pub struct PySegmentRecord(SegmentRecord<$model_record, $ideal_gas_record>);

        #[pymethods]
        impl PySegmentRecord {
            #[new]
            fn new(
                identifier: String,
                molarweight: f64,
                model_record: $py_model_record,
                ideal_gas_record: Option<$py_ideal_gas_record>,
            ) -> PyResult<Self> {
                Ok(Self(SegmentRecord::new(
                    identifier,
                    molarweight,
                    model_record.0,
                    ideal_gas_record.map(|ig| ig.0),
                )))
            }

            #[getter]
            fn get_identifier(&self) -> String {
                self.0.identifier.clone()
            }

            #[setter]
            fn set_identifier(&mut self, identifier: String) {
                self.0.identifier = identifier;
            }

            #[getter]
            fn get_molarweight(&self) -> f64 {
                self.0.molarweight
            }

            #[setter]
            fn set_molarweight(&mut self, molarweight: f64) {
                self.0.molarweight = molarweight;
            }

            #[getter]
            fn get_model_record(&self) -> $py_model_record {
                $py_model_record(self.0.model_record.clone())
            }

            #[setter]
            fn set_model_record(&mut self, model_record: $py_model_record) {
                self.0.model_record = model_record.0;
            }

            #[getter]
            fn get_ideal_gas_record(&self) -> Option<$py_ideal_gas_record> {
                self.0.ideal_gas_record.clone().map($py_ideal_gas_record)
            }

            #[setter]
            fn set_ideal_gas_record(&mut self, ideal_gas_record: $py_ideal_gas_record) {
                self.0.ideal_gas_record = Some(ideal_gas_record.0);
            }
        }

        #[pyproto]
        impl pyo3::class::basic::PyObjectProtocol for PySegmentRecord {
            fn __repr__(&self) -> PyResult<String> {
                Ok(self.0.to_string())
            }
        }

        impl_json_handling!(PySegmentRecord);
    };
}

#[macro_export]
macro_rules! impl_parameter {
    ($parameter:ty, $py_parameter:ty) => {
        #[pymethods]
        impl $py_parameter {
            /// Creates parameters from records.
            ///
            /// Parameters
            /// ----------
            /// pure_records : [PureRecord]
            ///     A list of pure component parameters.
            /// binary_records : [BinaryRecord], optional
            ///     A list of binary interaction parameters.
            #[staticmethod]
            #[pyo3(text_signature = "(pure_records, binary_records=None)")]
            fn from_records(
                pure_records: Vec<PyPureRecord>,
                binary_records: Option<Vec<PyBinaryRecord>>,
            ) -> Result<Self, ParameterError> {
                Ok(Self(<$parameter>::from_records(
                    pure_records.into_iter().map(|pr| pr.0).collect(),
                    binary_records.map(|br| br.into_iter().map(|br| br.0).collect()),
                )?))
            }

            /// Creates parameters from json files.
            ///
            /// Parameters
            /// ----------
            /// substances : List[str]
            ///     The substances to search.
            /// pure_path : str
            ///     Path to file containing pure substance parameters.
            /// binary_path : str, optional
            ///     Path to file containing binary substance parameters.
            /// search_option : str, optional, defaults to "Name"
            ///     Identifier that is used to search substance.
            ///     One of 'Name', 'Cas', 'Inchi', 'IupacName', 'Formula', 'Smiles'
            #[staticmethod]
            #[pyo3(
                text_signature = "(substances, pure_path, binary_path=None, search_option='Name')"
            )]
            fn from_json(
                substances: Vec<&str>,
                pure_path: String,
                binary_path: Option<String>,
                search_option: Option<&str>,
            ) -> Result<Self, ParameterError> {
                let io = match search_option {
                    Some(o) => IdentifierOption::try_from(o)?,
                    None => IdentifierOption::Name,
                };
                Ok(Self(<$parameter>::from_json(
                    &substances,
                    pure_path,
                    binary_path,
                    io,
                )?))
            }
        }
    };
}

#[macro_export]
macro_rules! impl_parameter_from_segments {
    ($parameter:ty, $py_parameter:ty) => {
        #[pymethods]
        impl $py_parameter {
            /// Creates parameters from segment records.
            ///
            /// Parameters
            /// ----------
            /// pure_records : [PureRecord]
            ///     A list of pure component parameters.
            /// segment_records : [SegmentRecord]
            ///     A list of records containing the parameters of
            ///     all individual segments.
            /// binary_segment_records : [BinarySegmentRecord], optional
            ///     A list of binary segment-segment parameters.
            #[staticmethod]
            #[pyo3(text_signature = "(pure_records, segment_records, binary_segment_records=None)")]
            fn from_segments(
                pure_records: Vec<PyPureRecord>,
                segment_records: Vec<PySegmentRecord>,
                binary_segment_records: Option<Vec<PyBinarySegmentRecord>>,
            ) -> Result<Self, ParameterError> {
                Ok(Self(<$parameter>::from_segments(
                    pure_records.into_iter().map(|pr| pr.0).collect(),
                    segment_records.into_iter().map(|sr| sr.0).collect(),
                    binary_segment_records.map(|r| r.into_iter().map(|r| r.0).collect()),
                )?))
            }

            /// Creates parameters using segments from json file.
            ///
            /// Parameters
            /// ----------
            /// substances : List[str]
            ///     The substances to search.
            /// pure_path : str
            ///     Path to file containing pure substance parameters.
            /// segments_path : str
            ///     Path to file containing segment parameters.
            /// binary_path : str, optional
            ///     Path to file containing binary segment-segment parameters.
            /// search_option : str, optional, defaults to "Name"
            ///     Identifier that is used to search substance.
            ///     One of 'Name', 'Cas', 'Inchi', 'IupacName', 'Formula', 'Smiles'
            #[staticmethod]
            #[pyo3(
                text_signature = "(substances, pure_path, segments_path, binary_path=None, search_option='Name')"
            )]
            fn from_json_segments(
                substances: Vec<&str>,
                pure_path: String,
                segments_path: String,
                binary_path: Option<String>,
                search_option: Option<&str>,
            ) -> Result<Self, ParameterError> {
                let io = match search_option {
                    Some(o) => IdentifierOption::try_from(o)?,
                    None => IdentifierOption::Name,
                };
                Ok(Self(<$parameter>::from_json_segments(
                    &substances,
                    pure_path,
                    segments_path,
                    binary_path,
                    io,
                )?))
            }
        }
    };
}

#[macro_export]
macro_rules! impl_json_handling {
    ($py_parameter:ty) => {
        #[pymethods]
        impl $py_parameter {
            /// Creates record from json string.
            #[staticmethod]
            #[pyo3(text_signature = "(json)")]
            fn from_json_str(json: &str) -> Result<Self, ParameterError> {
                Ok(Self(serde_json::from_str(json)?))
            }

            /// Creates a json string from record.
            fn to_json_str(&self) -> Result<String, ParameterError> {
                Ok(serde_json::to_string(&self.0)?)
            }
        }
    };
}
