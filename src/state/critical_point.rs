use super::{Contributions, State, StateHD, TPSpec};
use crate::equation_of_state::EquationOfState;
use crate::errors::{EosError, EosResult};
use crate::phase_equilibria::{SolverOptions, Verbosity};
use crate::EosUnit;
use argmin::prelude::{ArgminOp, Error, Executor};
use argmin::solver::brent::Brent;
use ndarray::{arr1, arr2, Array1, Array2};
use num_dual::linalg::{norm, smallest_ev, LU};
use num_dual::{Dual3, Dual64, DualNum, HyperDual};
use num_traits::{One, Zero};
use quantity::{QuantityArray1, QuantityScalar};
use std::rc::Rc;

const MAX_ITER_CRIT_POINT: usize = 50;
const TOL_CRIT_POINT: f64 = 1e-8;

/// # Critical points
impl<U: EosUnit, E: EquationOfState> State<U, E> {
    /// Calculate the pure component critical point of all components.
    pub fn critical_point_pure(
        eos: &Rc<E>,
        initial_temperature: Option<QuantityScalar<U>>,
        options: SolverOptions,
    ) -> EosResult<Vec<Self>>
    where
        QuantityScalar<U>: std::fmt::Display,
    {
        (0..eos.components())
            .map(|i| {
                Self::critical_point(
                    &Rc::new(eos.subset(&[i])),
                    None,
                    initial_temperature,
                    options,
                )
            })
            .collect()
    }

    /// Calculate the critical point of a binary system for given temperature.
    pub fn critical_point_binary_t(
        eos: &Rc<E>,
        temperature: QuantityScalar<U>,
        options: SolverOptions,
    ) -> EosResult<Self>
    where
        QuantityScalar<U>: std::fmt::Display,
    {
        Self::critical_point_binary(eos, TPSpec::Temperature(temperature), options)
    }

    /// Calculate the critical point of a binary system for given pressure.
    pub fn critical_point_binary_p(
        eos: &Rc<E>,
        pressure: QuantityScalar<U>,
        options: SolverOptions,
    ) -> EosResult<Self>
    where
        QuantityScalar<U>: std::fmt::Display,
    {
        Self::critical_point_binary(eos, TPSpec::Pressure(pressure), options)
    }

    pub(crate) fn critical_point_binary(
        eos: &Rc<E>,
        tp: TPSpec<U>,
        options: SolverOptions,
    ) -> EosResult<Self>
    where
        QuantityScalar<U>: std::fmt::Display,
    {
        let solver = Brent::new(1e-10, 1.0 - 1e-10, options.tol.unwrap_or(TOL_CRIT_POINT));
        let cost = CritOp::new(eos, tp);
        let x = Executor::new(cost, solver, 0.5)
            .max_iters(options.max_iter.unwrap_or(MAX_ITER_CRIT_POINT) as u64)
            .run()?
            .state
            .best_param;
        let moles = arr1(&[x, 1.0 - x]) * U::reference_moles();
        State::critical_point(eos, Some(&moles), None, SolverOptions::default())
    }

    /// Calculate the critical point of a system for given moles.
    pub fn critical_point(
        eos: &Rc<E>,
        moles: Option<&QuantityArray1<U>>,
        initial_temperature: Option<QuantityScalar<U>>,
        options: SolverOptions,
    ) -> EosResult<Self>
    where
        QuantityScalar<U>: std::fmt::Display,
    {
        let moles = eos.validate_moles(moles)?;
        let trial_temperatures = [
            300.0 * U::reference_temperature(),
            700.0 * U::reference_temperature(),
            500.0 * U::reference_temperature(),
        ];
        if let Some(t) = initial_temperature {
            return Self::critical_point_hkm(eos, &moles, t, options);
        }
        for &t in trial_temperatures.iter() {
            let s = Self::critical_point_hkm(eos, &moles, t, options);
            if s.is_ok() {
                return s;
            }
        }
        Err(EosError::NotConverged(String::from("Critical point")))
    }

    fn critical_point_hkm(
        eos: &Rc<E>,
        moles: &QuantityArray1<U>,
        initial_temperature: QuantityScalar<U>,
        options: SolverOptions,
    ) -> EosResult<Self>
    where
        QuantityScalar<U>: std::fmt::Display,
    {
        let (max_iter, tol, verbosity) = options.unwrap_or(MAX_ITER_CRIT_POINT, TOL_CRIT_POINT);

        let mut t = initial_temperature.to_reduced(U::reference_temperature())?;
        let max_density = eos
            .max_density(Some(moles))?
            .to_reduced(U::reference_density())?;
        let mut rho = 0.3 * max_density;
        let n = moles.to_reduced(U::reference_moles())?;

        log_iter!(
            verbosity,
            " iter |    residual    |   temperature   |       density        "
        );
        log_iter!(verbosity, "{:-<64}", "");
        log_iter!(
            verbosity,
            " {:4} |                | {:13.8} | {:12.8}",
            0,
            t * U::reference_temperature(),
            rho * U::reference_density(),
        );

        for i in 1..=max_iter {
            // calculate residuals and derivatives w.r.t. temperature and density
            let res_t =
                critical_point_objective(eos, Dual64::from(t).derive(), Dual64::from(rho), &n)?;
            let res_r =
                critical_point_objective(eos, Dual64::from(t), Dual64::from(rho).derive(), &n)?;
            let res = res_t.map(Dual64::re);

            // calculate Newton step
            let h = arr2(&[
                [res_t[0].eps[0], res_r[0].eps[0]],
                [res_t[1].eps[0], res_r[1].eps[0]],
            ]);
            let mut delta = LU::new(h)?.solve(&res);

            // reduce step if necessary
            if delta[0].abs() > 0.25 * t {
                delta *= 0.25 * t / delta[0].abs()
            }
            if delta[1].abs() > 0.03 * max_density {
                delta *= 0.03 * max_density / delta[1].abs()
            }

            // apply step
            t -= delta[0];
            rho -= delta[1];
            rho = f64::max(rho, 1e-4 * max_density);

            log_iter!(
                verbosity,
                " {:4} | {:14.8e} | {:13.8} | {:12.8}",
                i,
                norm(&res),
                t * U::reference_temperature(),
                rho * U::reference_density(),
            );

            // check convergence
            if norm(&res) < tol {
                log_result!(
                    verbosity,
                    "Critical point calculation converged in {} step(s)\n",
                    i
                );
                return State::new_nvt(
                    eos,
                    t * U::reference_temperature(),
                    moles.sum() / (rho * U::reference_density()),
                    moles,
                );
            }
        }
        Err(EosError::NotConverged(String::from("Critical point")))
    }
}

pub fn critical_point_objective<E: EquationOfState>(
    eos: &Rc<E>,
    temperature: Dual64,
    density: Dual64,
    moles: &Array1<f64>,
) -> EosResult<Array1<Dual64>> {
    // calculate second partial derivatives w.r.t. moles
    let t = HyperDual::from_re(temperature);
    let v = HyperDual::from_re(density.recip() * moles.sum());
    let qij = Array2::from_shape_fn((eos.components(), eos.components()), |(i, j)| {
        let mut m = moles.mapv(HyperDual::from);
        m[i].eps1[0] = Dual64::one();
        m[j].eps2[0] = Dual64::one();
        let state = StateHD::new(t, v, m);
        (eos.evaluate_residual(&state).eps1eps2[(0, 0)]
            + eos.ideal_gas().evaluate(&state).eps1eps2[(0, 0)])
            * (moles[i] * moles[j]).sqrt()
    });

    // calculate smallest eigenvalue and corresponding eigenvector of q
    let (eval, evec) = smallest_ev(qij);

    // evaluate third partial derivative w.r.t. s
    let moles_hd = Array1::from_shape_fn(eos.components(), |i| {
        Dual3::new(
            Dual64::from(moles[i]),
            evec[i] * moles[i].sqrt(),
            Dual64::zero(),
            Dual64::zero(),
        )
    });
    let state_s = StateHD::new(
        Dual3::from_re(temperature),
        Dual3::from_re(density.recip() * moles.sum()),
        moles_hd,
    );
    let res = eos.evaluate_residual(&state_s) + eos.ideal_gas().evaluate(&state_s);
    Ok(arr1(&[eval, res.v3]))
}

struct CritOp<U: EosUnit, E: EquationOfState> {
    eos: Rc<E>,
    tp: TPSpec<U>,
}

impl<U: EosUnit, E: EquationOfState> CritOp<U, E> {
    fn new(eos: &Rc<E>, tp: TPSpec<U>) -> Self {
        Self {
            eos: eos.clone(),
            tp,
        }
    }
}

impl<U: EosUnit, E: EquationOfState> ArgminOp for CritOp<U, E>
where
    QuantityScalar<U>: std::fmt::Display,
{
    type Param = f64;
    type Output = f64;
    type Jacobian = ();
    type Hessian = ();
    type Float = f64;

    fn apply(&self, p: &Self::Param) -> Result<Self::Output, Error> {
        let moles = arr1(&[*p, 1.0 - *p]) * U::reference_moles();
        let state = State::critical_point(&self.eos, Some(&moles), None, SolverOptions::default())?;
        match self.tp {
            TPSpec::Pressure(p) => Ok(
                (state.pressure(Contributions::Total) - p).to_reduced(U::reference_pressure())?
            ),
            TPSpec::Temperature(t) => {
                Ok((state.temperature - t).to_reduced(U::reference_temperature())?)
            }
        }
    }
}
