use feos_core::cubic::PengRobinson;
use feos_core::python::cubic::{
    PyBinaryRecord, PyPengRobinsonParameters, PyPengRobinsonRecord, PyPureRecord,
};
use feos_core::*;
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use quantity::python::*;
use quantity::si::*;
use std::collections::HashMap;
use std::rc::Rc;

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
impl_phase_equilibrium!(PengRobinson, PyPengRobinson);

#[pymodule]
pub fn cubic(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPengRobinson>()?;
    m.add_class::<PyPengRobinsonParameters>()?;
    m.add_class::<PyPengRobinsonRecord>()?;
    m.add_class::<PyPureRecord>()?;
    m.add_class::<PyBinaryRecord>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagram>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    Ok(())
}
