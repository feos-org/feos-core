use feos_core::python::feos_core;
use pyo3::prelude::*;

#[pymodule]
pub fn build_wheel(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    feos_core(py, m)
}
