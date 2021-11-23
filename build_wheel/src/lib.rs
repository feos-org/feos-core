use feos_core::python::PyInit_feos_core;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;

#[pymodule]
pub fn build_wheel(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(feos_core))
}
