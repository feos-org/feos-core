use feos_core::python::joback::PyJobackRecord;
use feos_core::python::parameter::*;
use feos_core::{Contributions, Verbosity};
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::__PYO3_PYMODULE_DEF_QUANTITY;

mod cubic;
mod user_defined;
use cubic::__PYO3_PYMODULE_DEF_CUBIC;
use user_defined::__PYO3_PYMODULE_DEF_USER_DEFINED;

#[pymodule]
pub fn feos_core(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyIdentifier>()?;
    m.add_class::<Verbosity>()?;
    m.add_class::<Contributions>()?;
    m.add_class::<PyChemicalRecord>()?;
    m.add_class::<PyJobackRecord>()?;

    m.add_wrapped(wrap_pymodule!(user_defined))?;
    m.add_wrapped(wrap_pymodule!(cubic))?;
    m.add_wrapped(wrap_pymodule!(quantity))?;

    py.run(
        "\
import sys
sys.modules['feos_core.cubic'] = cubic
sys.modules['feos_core.user_defined'] = user_defined
quantity.SINumber.__module__ = 'feos_core.si'
quantity.SIArray1.__module__ = 'feos_core.si'
quantity.SIArray2.__module__ = 'feos_core.si'
quantity.SIArray3.__module__ = 'feos_core.si'
quantity.SIArray4.__module__ = 'feos_core.si'
sys.modules['feos_core.si'] = quantity
    ",
        None,
        Some(m.dict()),
    )?;
    Ok(())
}
