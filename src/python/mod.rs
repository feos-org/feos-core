use crate::{Contributions, EosError, Verbosity};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::{wrap_pymodule, PyErr};
use quantity::python::PyInit_quantity;

mod cubic;
mod equation_of_state;
pub mod joback;
pub mod parameter;
mod phase_equilibria;
mod state;
mod statehd;
mod user_defined;
mod utils;

pub use cubic::PyInit_cubic;
pub use user_defined::PyInit_user_defined;

/// Helmholtz energy contributions to consider
/// when computing a property.
#[pyclass(name = "Contributions")]
#[derive(Copy, Clone)]
pub struct PyContributions(pub Contributions);

#[pymethods]
impl PyContributions {
    /// Only compute ideal gas contribution.
    #[classattr]
    #[allow(non_snake_case)]
    pub fn IdealGas() -> Self {
        Self(Contributions::IdealGas)
    }

    /// Only compute residual contribution with respect
    /// to an ideal gas contribution which is defined at
    /// T, V, {n}.
    ///
    /// See also
    /// --------
    /// ResidualP: to use an ideal gas reference defined at T, p, {n}
    #[classattr]
    #[allow(non_snake_case)]
    pub fn Residual() -> Self {
        Self(Contributions::Residual)
    }

    /// Only compute residual contribution with respect
    /// to an ideal gas contribution which is defined at
    /// T, p, {n}.
    ///
    /// See also
    /// --------
    /// Residual: to use an ideal gas reference defined at T, V, {n}
    #[classattr]
    #[allow(non_snake_case)]
    pub fn ResidualP() -> Self {
        Self(Contributions::ResidualP)
    }

    /// Compute all contributions
    ///
    /// Note
    /// ----
    /// This is the default for most properties.
    #[classattr]
    #[allow(non_snake_case)]
    pub fn Total() -> Self {
        Self(Contributions::Total)
    }
}

/// Verbosity levels for iterative solvers.
#[pyclass(name = "Verbosity")]
#[derive(Copy, Clone)]
pub struct PyVerbosity(pub Verbosity);

#[pymethods]
impl PyVerbosity {
    /// Print a status message at the end of the iteration.
    #[classattr]
    #[allow(non_snake_case)]
    pub fn Result() -> Self {
        Self(Verbosity::Result)
    }

    /// Print a detailed progress of the iteration.
    #[classattr]
    #[allow(non_snake_case)]
    pub fn Iter() -> Self {
        Self(Verbosity::Iter)
    }
}

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
