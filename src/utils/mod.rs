//! Collection of utility functionalities.
//!
//! # Utilities for parameter estimation
//! - [`DataSet`]: stores experimental data and defines a cost function that can be used to compare to equation of state calculation.
//! - [`Estimator`]: stores multiple `DataSet`
mod dataset;
mod estimator;
mod loss;
pub use dataset::*;
pub use estimator::{Estimator, FitError};
pub use loss::*;
