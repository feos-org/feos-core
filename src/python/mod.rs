use crate::EosError;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::{wrap_pymodule, PyErr};
use quantity::python::__PYO3_PYMODULE_DEF_QUANTITY;

mod cubic;
mod equation_of_state;
pub mod joback;
pub mod parameter;
mod phase_equilibria;
mod state;
mod statehd;
mod user_defined;
mod utils;

pub use cubic::__PYO3_PYMODULE_DEF_CUBIC;
pub use user_defined::__PYO3_PYMODULE_DEF_USER_DEFINED;

impl From<EosError> for PyErr {
    fn from(e: EosError) -> PyErr {
        PyRuntimeError::new_err(e.to_string())
    }
}

#[pymodule]
pub fn feos_core(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(quantity))?;
    m.add_wrapped(wrap_pymodule!(user_defined))?;
    m.add_wrapped(wrap_pymodule!(cubic))?;
    py.run(
        "\
import sys
quantity.SINumber.__module__ = 'feos_core.si'
quantity.SIArray1.__module__ = 'feos_core.si'
quantity.SIArray2.__module__ = 'feos_core.si'
quantity.SIArray3.__module__ = 'feos_core.si'
quantity.SIArray4.__module__ = 'feos_core.si'
sys.modules['feos_core.si'] = quantity
sys.modules['feos_core.user_defined'] = user_defined
sys.modules['feos_core.cubic'] = cubic
    ",
        None,
        Some(m.dict()),
    )?;
    Ok(())
}
