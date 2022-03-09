use crate::parameter::ParameterError;
use argmin::core::Error as ArgminError;
use num_dual::linalg::LinAlgError;
use quantity::QuantityError;
use thiserror::Error;

/// Error type for improperly defined states and convergence problems.
#[derive(Error, Debug)]
pub enum EosError {
    #[error("`{0}` did not converge within the maximum number of iterations.")]
    NotConverged(String),
    #[error("`{0}` encountered illegal values during the iteration.")]
    IterationFailed(String),
    #[error("Iteration resulted in trivial solution.")]
    TrivialSolution,
    #[error("Equation of state is initialized for {0} components while the input specifies {1} components.")]
    IncompatibleComponents(usize, usize),
    #[error("Invalid state in {0}: {1} = {2}.")]
    InvalidState(String, String, f64),
    #[error("Undetermined state: {0}.")]
    UndeterminedState(String),
    #[error("System is supercritical.")]
    SuperCritical(),
    #[error("No phase split according to stability analysis.")]
    NoPhaseSplit,
    #[error(transparent)]
    QuantityError(#[from] QuantityError),
    #[error(transparent)]
    ParameterError(#[from] ParameterError),
    #[error(transparent)]
    ArgminError(#[from] ArgminError),
    #[error(transparent)]
    LinAlgError(#[from] LinAlgError),
}

/// Convenience type for `Result<T, EosError>`.
pub type EosResult<T> = Result<T, EosError>;
