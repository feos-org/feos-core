use crate::python::{statehd::*, PyContributions, PyVerbosity};
use crate::*;
use ndarray::prelude::*;
use num_dual::python::{
    PyDual3Dual64, PyDual3_64, PyDual64, PyHyperDual64, PyHyperDualDual64, PyInit_num_dual,
};
use num_dual::{Dual3, Dual3_64, Dual64, HyperDual, HyperDual64};
use numpy::convert::{IntoPyArray, ToPyArray};
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::{PySIArray1, PySIArray2, PySINumber};
use std::collections::HashMap;
use std::fmt;
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
/// RunTimeError : if the class does not implement all necessary methods.
#[pyclass(name = "UserDefinedEos", unsendable)]
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

struct PyHelmholtzEnergy(Py<PyAny>);

pub struct PyEoSObj {
    obj: Py<PyAny>,
    contributions: Vec<Box<dyn HelmholtzEnergy>>,
}

impl PyEoSObj {
    pub fn new(obj: Py<PyAny>) -> PyResult<Self> {
        Python::with_gil(|py| {
            let attr = obj.as_ref(py).hasattr("components")?;
            if !attr {
                panic!("Python Class has to have a method 'components' with signature:\n\tdef signature(self) -> int")
            }
            let attr = obj.as_ref(py).hasattr("subset")?;
            if !attr {
                panic!("Python Class has to have a method 'subset' with signature:\n\tdef subset(self, component_list: List[int]) -> Self")
            }
            let attr = obj.as_ref(py).hasattr("molar_weight")?;
            if !attr {
                panic!("Python Class has to have a method 'molar_weight' with signature:\n\tdef molar_weight(self) -> SIArray1\nwhere the size of the returned array has to be 'components'.")
            }
            let attr = obj.as_ref(py).hasattr("max_density")?;
            if !attr {
                panic!("Python Class has to have a method 'max_density' with signature:\n\tdef max_density(self, moles: numpy.ndarray[float]) -> float\nwhere the size of the input array has to be 'components'.")
            }
            let attr = obj.as_ref(py).hasattr("helmholtz_energy")?;
            if !attr {
                panic!("{}", "Python Class has to have a method 'helmholtz_energy' with signature:\n\tdef helmholtz_energy(self, state: StateHD) -> HD\nwhere 'HD' has to be any of {{float, Dual64, HyperDual64, HyperDualDual64, Dual3Dual64, Dual3_64}}.")
            }
            Ok(Self {
                obj: obj.clone(),
                contributions: vec![Box::new(PyHelmholtzEnergy(obj))],
            })
        })
    }
}

impl MolarWeight<SIUnit> for PyEoSObj {
    fn molar_weight(&self) -> SIArray1 {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let py_result = self.obj.as_ref(py).call_method0("molar_weight").unwrap();
        if py_result.get_type().name().unwrap() != "SIArray1" {
            panic!(
                "Expected an 'SIArray1' for the 'molar_weight' method return type, got {}",
                py_result.get_type().name().unwrap()
            );
        }
        py_result.extract::<PySIArray1>().unwrap().into()
    }
}

impl EquationOfState for PyEoSObj {
    fn components(&self) -> usize {
        Python::with_gil(|py| {
            let py_result = self.obj.as_ref(py).call_method0("components").unwrap();
            if py_result.get_type().name().unwrap() != "int" {
                panic!(
                    "Expected an integer for the components() method signature, got {}",
                    py_result.get_type().name().unwrap()
                );
            }
            py_result.extract().unwrap()
        })
    }

    fn subset(&self, component_list: &[usize]) -> Self {
        Python::with_gil(|py| {
            let py_result = self
                .obj
                .as_ref(py)
                .call_method1("subset", (component_list.to_vec(),))
                .unwrap();
            Self::new(py_result.extract().unwrap()).unwrap()
        })
    }

    fn compute_max_density(&self, moles: &Array1<f64>) -> f64 {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let py_result = self
            .obj
            .as_ref(py)
            .call_method1("max_density", (moles.to_owned().into_pyarray(py),))
            .unwrap();
        // if py_result.get_type().name().unwrap() != "numpy.float64" {
        //     panic!(
        //         "Expected an 'numpy.float64' for the 'compute_max_density' method return type, got {}",
        //         py_result.get_type().name().unwrap()
        //     );
        // }
        py_result.extract().unwrap()
    }

    fn residual(&self) -> &[Box<dyn HelmholtzEnergy>] {
        &self.contributions
    }
}

macro_rules! impl_helmholtz_energy {
    ($pystate:ty, $pyhd:ty, $hd:ty) => {
        impl HelmholtzEnergyDual<$hd> for PyHelmholtzEnergy {
            fn helmholtz_energy(&self, state: &StateHD<$hd>) -> $hd {
                let gil = Python::acquire_gil();
                let py = gil.python();
                let py_result = self
                    .0
                    .as_ref(py)
                    .call_method1("helmholtz_energy", (<$pystate>::from(state.clone()),))
                    .unwrap();
                // if py_result.get_type().name() != stringify!($hd) {
                //     panic!(
                //         "Expected an {} for the 'helmholtz_energy' method signature, got {}",
                //         stringify!($hd),
                //         py_result.get_type().name()
                //     );
                // }
                <$hd>::from(py_result.extract::<$pyhd>().unwrap())
            }

            // fn identifier(&self) -> String {
            //     let gil = Python::acquire_gil();
            //     let py = gil.python();
            //     let py_result = self.obj.as_ref(py).call_method0("identifier").unwrap();
            //     if py_result.get_type().name() != "str" {
            //         panic!(
            //             "Expected an 'str' for the 'identifier' method signature, got {}",
            //             py_result.get_type().name()
            //         );
            //     }
            //     py_result.extract().unwrap()
            // }
        }
    };
}

impl fmt::Display for PyHelmholtzEnergy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Custom")
    }
}

impl_helmholtz_energy!(PyStateD, PyDual64, Dual64);
impl_helmholtz_energy!(PyStateHD, PyHyperDual64, HyperDual64);
impl_helmholtz_energy!(PyStateHDD, PyHyperDualDual64, HyperDual<Dual64, f64>);
impl_helmholtz_energy!(PyStateD3, PyDual3_64, Dual3_64);
impl_helmholtz_energy!(PyStateD3D, PyDual3Dual64, Dual3<Dual64, f64>);
impl_helmholtz_energy!(PyStateF, f64, f64);

impl_equation_of_state!(PyUserDefinedEos);
impl_virial_coefficients!(PyUserDefinedEos);
impl_state!(PyEoSObj, PyUserDefinedEos);
impl_state_molarweight!(PyEoSObj, PyUserDefinedEos);
impl_vle_state!(PyEoSObj, PyUserDefinedEos);

#[pymodule]
pub fn user_defined(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(num_dual))?;

    m.add_class::<PyStateHD>()?;
    m.add_class::<PyStateD>()?;
    m.add_class::<PyStateD3>()?;
    m.add_class::<PyStateF>()?;
    m.add_class::<PyStateHDD>()?;
    m.add_class::<PyStateD3D>()?;

    m.add_class::<PyUserDefinedEos>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    m.add_class::<PyPhaseDiagramPure>()?;
    m.add_class::<PyPhaseDiagramBinary>()?;
    m.add_class::<PyPhaseDiagramHetero>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    Ok(())
}
