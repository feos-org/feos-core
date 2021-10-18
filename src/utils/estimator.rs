//! The [`Estimator`] struct can be used to store multiple [`DataSet`]s for convenient parameter
//! optimization.
use super::dataset::*;
use crate::equation_of_state::EquationOfState;
use crate::EosUnit;
use ndarray::{arr1, concatenate, Array1, ArrayView1, Axis};
use quantity::{QuantityError, QuantityScalar};
use std::fmt;
use std::fmt::Write;
use std::num::ParseFloatError;
use std::rc::Rc;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum FitError {
    #[error("Missing input. Need '{needed}' to evaluate '{to_evaluate}'.")]
    MissingInput { needed: String, to_evaluate: String },
    #[error("Input has not the same amount of data as the target.")]
    IncompatibleInput,
    #[error("Keyword '{0}' unknown. Try: 'liquid density', 'vapor pressure', 'equilibrium liquid density'")]
    KeyError(String),
    #[error(transparent)]
    ShapeError(#[from] ndarray::ShapeError),
    #[error(transparent)]
    ParseError(#[from] ParseFloatError),
    #[error(transparent)]
    QuantityError(#[from] QuantityError),
}

/// A collection of [`DataSet`]s and weights that can be used to
/// evaluate an equation of state versus experimental data.
pub struct Estimator<U: EosUnit, E: EquationOfState> {
    data: Vec<Rc<dyn DataSet<U, E>>>,
    weights: Vec<f64>,
}

impl<U: EosUnit, E: EquationOfState> Estimator<U, E>
where
    QuantityScalar<U>: std::fmt::Display + std::fmt::LowerExp,
{
    /// Create a new `Estimator` given `DataSet`s and weights.
    ///
    /// The weights are normalized and used as multiplicator when the
    /// cost function across all `DataSet`s is evaluated.
    pub fn new(data: Vec<Rc<dyn DataSet<U, E>>>, weights: Vec<f64>) -> Self {
        Self { data, weights }
    }

    /// Add a `DataSet` and its weight.
    pub fn add_data(&mut self, data: &Rc<dyn DataSet<U, E>>, weight: f64) {
        self.data.push(data.clone());
        self.weights.push(weight);
    }

    /// Returns the cost of each `DataSet`.
    ///
    /// Each cost contains the inverse weight.
    pub fn cost(&self, eos: &Rc<E>) -> Result<Array1<f64>, FitError> {
        let predictions: Result<Vec<Array1<f64>>, FitError> = self
            .data
            .iter()
            .enumerate()
            .map(|(i, d)| {
                let w_sum = self.weights.iter().sum::<f64>();
                let w = arr1(&self.weights) / w_sum;
                Ok(d.cost(eos)? * w[i])
            })
            .collect();
        if let Ok(p) = predictions {
            let aview: Vec<ArrayView1<f64>> = p.iter().map(|pi| pi.view()).collect();
            Ok(concatenate(Axis(0), &aview)?)
        } else {
            Err(FitError::IncompatibleInput)
        }
    }

    /// Returns the relative difference for each `DataSet`.
    pub fn relative_difference(&self, eos: &Rc<E>) -> Result<Vec<Array1<f64>>, FitError> {
        self.data
            .iter()
            .map(|d| d.relative_difference(eos))
            .collect()
    }

    /// Returns the mean absolute relative difference for each `DataSet`.
    pub fn mean_absolute_relative_difference(&self, eos: &Rc<E>) -> Result<Array1<f64>, FitError> {
        self.data
            .iter()
            .map(|d| d.mean_absolute_relative_difference(eos))
            .collect()
    }

    /// Returns the stored `DataSet`s.
    pub fn datasets(&self) -> Vec<Rc<dyn DataSet<U, E>>> {
        self.data.to_vec()
    }

    pub fn _repr_markdownn_(&self) -> String {
        let mut f = String::new();
        write!(f, "| target | input | datapoints |\n|:-|:-|:-|").unwrap();
        for d in self.data.iter() {
            write!(
                f,
                "\n|{}|{}|{}|",
                d.target_str(),
                d.input_str().join(", "),
                d.datapoints()
            )
            .unwrap();
        }
        f
    }
}

impl<U: EosUnit, E: EquationOfState> fmt::Display for Estimator<U, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for d in self.data.iter() {
            writeln!(f, "{}", d.to_string())?;
        }
        Ok(())
    }
}
