use crate::impl_json_handling;
use crate::joback::JobackRecord;
use crate::parameter::ParameterError;
use pyo3::prelude::*;

/// Create a set of Joback ideal gas heat capacity parameters
/// for a segment or a pure component.
///
/// The fourth order coefficient `e` is not present in the
/// orginial publication by Joback and Reid but is required
/// for correlations for some pure components that are modeled
/// using the samen polynomial approach.
///
/// Parameters
/// ----------
/// a : float
///     zeroth order coefficint
/// b : float
///     first order coefficint
/// c : float
///     second order coefficint
/// d : float
///     third order coefficint
/// e : float
///     fourth order coefficint
///
/// Returns
/// -------
/// JobackRecord
#[pyclass(name = "JobackRecord")]
#[derive(Clone)]
#[pyo3(text_signature = "(a, b, c, d, e)")]
pub struct PyJobackRecord(pub JobackRecord);

#[pymethods]
impl PyJobackRecord {
    #[new]
    fn new(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        Self(JobackRecord::new(a, b, c, d, e))
    }
}

#[pyproto]
impl pyo3::class::basic::PyObjectProtocol for PyJobackRecord {
    fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }
}

impl_json_handling!(PyJobackRecord);
