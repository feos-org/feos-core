use super::{PhaseEquilibrium, SolverOptions};
use crate::equation_of_state::{EquationOfState, MolarWeight};
use crate::errors::EosResult;
use crate::state::{Contributions, State};
use crate::EosUnit;
use ndarray::prelude::*;
use quantity::{QuantityArray1, QuantityScalar};
use std::rc::Rc;

/// Pure component and binary mixture phase diagrams.
pub struct PhaseDiagram<U, E> {
    pub states: Vec<PhaseEquilibrium<U, E, 2>>,
}

impl<U: Clone, E> Clone for PhaseDiagram<U, E> {
    fn clone(&self) -> Self {
        Self {
            states: self.states.clone(),
        }
    }
}

impl<U: EosUnit, E: EquationOfState> PhaseDiagram<U, E> {
    /// Calculate a phase diagram for a pure component.
    pub fn pure(
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
            vle = PhaseEquilibrium::pure(eos, ti, vle.as_ref(), options).ok();
            if let Some(vle) = vle.as_ref() {
                states.push(vle.clone());
            }
        }
        states.push(PhaseEquilibrium::from_states(sc.clone(), sc));

        Ok(PhaseDiagram { states })
    }

    /// Return the vapor states of the diagram.
    pub fn vapor(&self) -> StateVec<'_, U, E> {
        StateVec {
            states: self.states.iter().map(|s| s.vapor()).collect(),
        }
    }

    /// Return the liquid states of the diagram.
    pub fn liquid(&self) -> StateVec<'_, U, E> {
        StateVec {
            states: self.states.iter().map(|s| s.liquid()).collect(),
        }
    }
}

/// A list of states for a simple access to properties
/// of multiple states.
pub struct StateVec<'a, U, E> {
    pub states: Vec<&'a State<U, E>>,
}

impl<'a, U: EosUnit, E: EquationOfState> StateVec<'a, U, E> {
    pub fn temperature(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| self.states[i].temperature)
    }

    pub fn pressure(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].pressure(Contributions::Total)
        })
    }

    pub fn compressibility(&self) -> Array1<f64> {
        Array1::from_shape_fn(self.states.len(), |i| {
            self.states[i].compressibility(Contributions::Total)
        })
    }

    pub fn density(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| self.states[i].density)
    }

    pub fn molefracs(&self) -> Array2<f64> {
        Array2::from_shape_fn(
            (self.states.len(), self.states[0].eos.components()),
            |(i, j)| self.states[i].molefracs[j],
        )
    }

    pub fn molar_enthalpy(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].molar_enthalpy(Contributions::Total)
        })
    }

    pub fn molar_entropy(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].molar_entropy(Contributions::Total)
        })
    }
}

impl<'a, U: EosUnit, E: EquationOfState + MolarWeight<U>> StateVec<'a, U, E> {
    pub fn mass_density(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| self.states[i].mass_density())
    }

    pub fn massfracs(&self) -> Array2<f64> {
        Array2::from_shape_fn(
            (self.states.len(), self.states[0].eos.components()),
            |(i, j)| self.states[i].massfracs()[j],
        )
    }

    pub fn specific_enthalpy(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].specific_enthalpy(Contributions::Total)
        })
    }

    pub fn specific_entropy(&self) -> QuantityArray1<U> {
        QuantityArray1::from_shape_fn(self.states.len(), |i| {
            self.states[i].specific_entropy(Contributions::Total)
        })
    }
}
