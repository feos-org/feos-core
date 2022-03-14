use super::{PhaseEquilibrium, SolverOptions};
use crate::equation_of_state::EquationOfState;
use crate::errors::EosResult;
use crate::state::{Contributions, State};
use crate::EosUnit;
use ndarray::prelude::*;
use quantity::{QuantityArray1, QuantityScalar};
use std::fmt;
use std::rc::Rc;

/// Pure component phase diagram.
#[derive(Debug)]
pub struct PhaseDiagramPure<U: EosUnit, E: EquationOfState> {
    pub states: Vec<PhaseEquilibrium<U, E, 2>>,
}

impl<U: EosUnit, E: EquationOfState> PhaseDiagramPure<U, E> {
    pub fn new(
        eos: &Rc<E>,
        min_temperature: QuantityScalar<U>,
        npoints: usize,
        critical_temperature: Option<QuantityScalar<U>>,
        options: SolverOptions,
    ) -> EosResult<Self>
    where
        QuantityScalar<U>: std::fmt::Display + std::fmt::LowerExp,
    {
        let mut states = Vec::with_capacity(npoints);

        let sc = State::critical_point(eos, None, critical_temperature, SolverOptions::default())?;

        let max_temperature = min_temperature
            + (sc.temperature - min_temperature) * ((npoints - 2) as f64 / (npoints - 1) as f64);
        let temperatures = Array::linspace(0.0, 1.0, npoints - 1)
            .map(|&i| min_temperature + (max_temperature - min_temperature) * i);

        let mut vle = None;
        for &ti in temperatures.iter() {
            vle = PhaseEquilibrium::pure_t(eos, ti, vle.as_ref(), options).ok();
            if let Some(vle) = vle.as_ref() {
                states.push(vle.clone());
            }
        }
        states.push(PhaseEquilibrium::from_states(sc.clone(), sc));

        Ok(PhaseDiagramPure { states })
    }

    pub fn temperature(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| self.states[i].vapor().temperature)
    }

    pub fn pressure(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].vapor().pressure(Contributions::Total)
        })
    }

    pub fn density_vapor(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| self.states[i].vapor().density)
    }

    pub fn density_liquid(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| self.states[i].liquid().density)
    }

    pub fn molar_enthalpy_vapor(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].vapor().molar_enthalpy(Contributions::Total)
        })
    }

    pub fn molar_enthalpy_liquid(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].liquid().molar_enthalpy(Contributions::Total)
        })
    }

    pub fn molar_entropy_vapor(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].vapor().molar_entropy(Contributions::Total)
        })
    }

    pub fn molar_entropy_liquid(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].liquid().molar_entropy(Contributions::Total)
        })
    }
}

impl<U, E> fmt::Display for PhaseDiagramPure<U, E>
where
    U: EosUnit,
    QuantityScalar<U>: fmt::Display,
    E: EquationOfState,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let temperature = self.temperature();
        let pressure = self.pressure();
        let density_vapor = self.density_vapor();
        let density_liquid = self.density_liquid();
        let molar_enthalpy_vapor = self.molar_enthalpy_vapor();
        let molar_enthalpy_liquid = self.molar_enthalpy_liquid();
        let molar_entropy_vapor = self.molar_entropy_vapor();
        let molar_entropy_liquid = self.molar_entropy_liquid();

        for i in 0..temperature.len() {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                temperature.get(i),
                pressure.get(i),
                density_vapor.get(i),
                density_liquid.get(i),
                molar_enthalpy_vapor.get(i),
                molar_enthalpy_liquid.get(i),
                molar_entropy_vapor.get(i),
                molar_entropy_liquid.get(i)
            )?;
        }
        Ok(())
    }
}
