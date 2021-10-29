# Welcome to `feos-core`!

> `FeOs` is in active development and so is this documentation. If you have questions or want to report bugs, please [file an issue](https://github.com/feos-org/feos-core/issues) or [discuss with us](https://github.com/feos-org/feos-core/discussions).

`feos-core` contains the basic data types and interfaces used to implement and work with equations of state.
It is used in all projects of the `FeOs` framework.

On these pages you will find information about both the Rust library, in which all data structures and algorithms are implemented,
and the Python interfaces.
If you are interested in the general concepts of `FeOs`, you are at the right place.
If you are interested in specific implementations of equations of state or the concepts of density functional theory in `FeOs`, see here:

- `feos-dft` (TODO: add link): basic data types and interfaces to implement Helmholtz energy functionals
- `feos-pcsaft` (TODO: add link): implementation of the PCP-SAFT equation of state and Helmholtz energy functionals
- `feos-gcpcsaft` (TODO: add link): implementation of the group contribution method for PCP-SAFT and Helmholtu energy functionals

## Features

### Concepts

- `feos-core` defines the interfaces for equations of state and provides objects that can be used to compute thermodynamic properties and phase equilibria
- interfaces use dimensioned properties (in SI units) which makes working with the code less error-prone
- equations of state can be implemented in Rust (more robust, faster at runtime) or as Python `class` (faster prototyping, slower at runtime, no Rust experience needed)
- Python interfaces are written in Rust. To expose an equation of state implemented in Rust, only very little code is needed.
- implementing an equation of state is done by implementing one or more contributions to the Helmholtz energy
  - derivatives are not needed - we use generalized dual numbers to evaluate the Helmholtz energy and its partial derivatives

### Methods

- thermodynamic properties as (partial) derivatives of the Helmholtz energy are computed using generalized dual numbers
- critical point calculations for pure substances and mixtures
- phase equilibrium calculations for pure substances and mixtures
- utilities to construct phase diagrams for pure substances and binary mixtures
- stability analysis
- dynamic properties via entropy scaling
- utilities for parameter I/O (json format)
- example implementation: Peng-Robinson equation of state

## Getting started: Python

- [Browse the examples](examples/index.rst) or download the Jupyter notebooks provided in the github repository,
- or take a look at [the Python API](api.rst).

## Getting started: Rust

- If you want to learn how the Rust library is structured, [take a look at the Rust guide](devguide/equation_of_state/index.rst)
- or check out the Rust API.

