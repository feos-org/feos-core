use feos_core::python::statehd::*;
use feos_core::python::user_defined::*;
use feos_core::*;
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use quantity::python::*;
use quantity::si::*;
use std::collections::HashMap;
use std::rc::Rc;

/// Equation of state implemented as python class.
///
/// Parameters
/// ----------
/// obj : python class object
///     The class that implements the equation of state.
///
/// Returns
/// -------
/// UserDefinedEos
///
/// Raises
/// ------
/// RunTimeError
///     If the class does not implement all necessary methods.
#[pyclass(name = "Python", unsendable)]
#[derive(Clone)]
#[pyo3(text_signature = "(obj)")]
pub struct PyUserDefinedEos(Rc<PyEoSObj>);

#[pymethods]
impl PyUserDefinedEos {
    #[new]
    fn new(obj: Py<PyAny>) -> PyResult<Self> {
        Ok(Self(Rc::new(PyEoSObj::new(obj)?)))
    }
}

impl_equation_of_state!(PyUserDefinedEos);
impl_virial_coefficients!(PyUserDefinedEos);
impl_state!(PyEoSObj, PyUserDefinedEos);
impl_state_molarweight!(PyEoSObj, PyUserDefinedEos);
impl_phase_equilibrium!(PyEoSObj, PyUserDefinedEos);

#[pymodule]
pub fn user_defined(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyStateHD>()?;
    m.add_class::<PyStateD>()?;
    m.add_class::<PyStateDDV3>()?;
    m.add_class::<PyStateD3>()?;
    m.add_class::<PyStateF>()?;
    m.add_class::<PyStateHDD>()?;
    m.add_class::<PyStateHDDV2>()?;
    m.add_class::<PyStateHDDV3>()?;
    m.add_class::<PyStateD3D>()?;
    m.add_class::<PyStateD3DV2>()?;
    m.add_class::<PyStateD3DV3>()?;
    m.add_class::<PyUserDefinedEos>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagram>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    Ok(())
}
