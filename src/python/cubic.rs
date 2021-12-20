use crate::cubic::{PengRobinson, PengRobinsonParameters, PengRobinsonRecord};
use crate::joback::JobackRecord;
use crate::parameter::{IdentifierOption, Parameter, ParameterError, PureRecord};
use crate::python::joback::PyJobackRecord;
use crate::python::parameter::{PyBinaryRecord, PyChemicalRecord, PyIdentifier};
use crate::python::{PyContributions, PyVerbosity};
use crate::*;
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use quantity::python::*;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::rc::Rc;

/// A pure substance parameter for the Peng-Robinson equation of state.
#[pyclass(name = "PengRobinsonRecord", unsendable)]
#[derive(Clone)]
pub struct PyPengRobinsonRecord(PengRobinsonRecord);

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyPengRobinsonRecord {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

impl_json_handling!(PyPengRobinsonRecord);

impl_pure_record!(
    PengRobinsonRecord,
    PyPengRobinsonRecord,
    JobackRecord,
    PyJobackRecord
);

/// Create a set of Peng-Robinson parameters from records.
///
/// Parameters
/// ----------
/// pure_records : List[PureRecord]
///     pure substance records.
/// binary_records : List[BinaryRecord], optional
///     binary parameter records
/// substances : List[str], optional
///     The substances to use. Filters substances from `pure_records` according to
///     `search_option`.
///     When not provided, all entries of `pure_records` are used.
/// search_option : {'Name', 'Cas', 'Inchi', 'IupacName', 'Formula', 'Smiles'}, optional, defaults to 'Name'.
///     Identifier that is used to search substance.
///
/// Returns
/// -------
/// PengRobinsonParameters
#[pyclass(name = "PengRobinsonParameters", unsendable)]
#[pyo3(
    text_signature = "(pure_records, binary_records=None, substances=None, search_option='Name')"
)]
#[derive(Clone)]
pub struct PyPengRobinsonParameters(pub Rc<PengRobinsonParameters>);

impl_parameter!(PengRobinsonParameters, PyPengRobinsonParameters);

#[pymethods]
impl PyPengRobinsonParameters {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

/* EQUATION OF STATE */

/// A simple version of the Peng-Robinson equation of state.
///
/// Parameters
/// ----------
/// parameters : PengRobinsonParameters
///     The parameters of the Peng-Robinson equation of state to use.
///
/// Returns
/// -------
/// PengRobinson
#[pyclass(name = "PengRobinson", unsendable)]
#[pyo3(text_signature = "(parameters)")]
#[derive(Clone)]
pub struct PyPengRobinson(pub Rc<PengRobinson>);

#[pymethods]
impl PyPengRobinson {
    #[new]
    fn new(parameters: PyPengRobinsonParameters) -> Self {
        Self(Rc::new(PengRobinson::new(parameters.0.clone())))
    }
}

impl_equation_of_state!(PyPengRobinson);
impl_virial_coefficients!(PyPengRobinson);

impl_state!(PengRobinson, PyPengRobinson);
impl_state_molarweight!(PengRobinson, PyPengRobinson);
impl_vle_state!(PengRobinson, PyPengRobinson);
// impl_estimator!(PengRobinson, PyPengRobinson);

#[pymodule]
pub fn cubic(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyIdentifier>()?;
    m.add_class::<PyVerbosity>()?;
    m.add_class::<PyContributions>()?;
    m.add_class::<PyChemicalRecord>()?;
    m.add_class::<PyJobackRecord>()?;

    m.add_class::<PyPengRobinsonRecord>()?;
    m.add_class::<PyPureRecord>()?;
    m.add_class::<PyBinaryRecord>()?;
    m.add_class::<PyPengRobinsonParameters>()?;

    m.add_class::<PyPengRobinson>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagramPure>()?;
    m.add_class::<PyPhaseDiagramBinary>()?;
    m.add_class::<PyPhaseDiagramHetero>()?;
    m.add_class::<PyPhaseEquilibrium>()?;

    // let utils = PyModule::new(py, "utils")?;
    // utils.add_class::<PyDataSet>()?;
    // utils.add_class::<PyEstimator>()?;
    // m.add_submodule(utils)?;
    Ok(())
}
