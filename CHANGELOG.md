# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

## [0.2.0] - 2022-03-09
### Added
- Added conversions between `ParameterError` -> `EosError` and `EosError` -> `FitError` to improve the error messages in some cases. [#40](https://github.com/feos-org/feos-core/pull/40)
- Added new struct `StateVec`, that gives easy access to properties of lists of states, e.g. in phase diagrams. [#48](https://github.com/feos-org/feos-core/pull/48)
- Added `ln_symmetric_activity_coefficient` and `ln_phi_pure` to the list of state properties that can be calculated. [#50](https://github.com/feos-org/feos-core/pull/50)

### Changed
- Removed `State` from `EntropyScaling` trait and adjusted associated methods to use temperature, volume and moles instead of state. [#36](https://github.com/feos-org/feos-core/pull/36)
- Replaced the outer loop iterations for the critical point of binary systems with dedicated algorithms. [#34](https://github.com/feos-org/feos-core/pull/34)
- Renamed `VLEOptions` to `SolverOptions`. [#38](https://github.com/feos-org/feos-core/pull/38)
- Renamed methods of `StateBuilder` and the parameters in the `State` constructor in python to `molar_enthalpy`, `molar_entropy`, and `molar_internal_energy`.  [#35](https://github.com/feos-org/feos-core/pull/35)
- Removed `PyContributions` and `PyVerbosity` in favor of a simpler implementation using `PyO3`'s new `#[pyclass]` for fieldless enums feature. [#41](https://github.com/feos-org/feos-core/pull/41)
- Renamed `Contributions::Residual` to `Contributions::ResidualNvt` and `Contributions::ResidualP` to `Contributions::ResidualNpt`. [#43](https://github.com/feos-org/feos-core/pull/43)
- Renamed macro `impl_vle_state!` to `impl_phase_equilibrium`. [#48](https://github.com/feos-org/feos-core/pull/48)
- Removed `_t` and `_p` functions in favor of simpler interfaces. The kind of specification (temperature or pressure) is determined from the unit of the argument. Properties of the phase diagram are available from the `vapor` and `liquid` getters, that return `StateVec`s. [#48](https://github.com/feos-org/feos-core/pull/48)
  - `PhaseEquilibrium::pure_t`, `PhaseEquilibrium::pure_p` -> `PhaseEquilibrium::Pure`
  - `PhaseEquilibrium::vle_pure_comps_t`, `PhaseEquilibrium::vle_pure_comps_p` -> `PhaseEquilibrium::vle_pure_comps`\
  The `PhaseEquilibria` returned by this function now have the same number of components as the (mixture) eos, that it is called with.
  - `PhaseEquilibrium::bubble_point_tx`, `PhaseEquilibrium::bubble_point_px` -> `PhaseEquilibrium::bubble_point`
  - `PhaseEquilibrium::dew_point_tx`, `PhaseEquilibrium::dew_point_px` -> `PhaseEquilibrium::dew_point`
  - `PhaseEquilibrium::heteroazeotrope_t`, `PhaseEquilibrium::heteroazeotrope_p` -> `PhaseEquilibrium::heteroazeotrope`
  - `State::critical_point_binary_t`, `State::critical_point_binary_p` -> `State::crititcal_point_binary`
- Combined `PhaseDiagramPure` and `PhaseDiagramBinary` into a single struct `PhaseDiagram` and renamed its constructors. [#48](https://github.com/feos-org/feos-core/pull/48)
  - `PhaseDiagramPure::new` -> `PhaseDiagram::pure`
  - `PhaseDiagramBinary::new_txy`, `PhaseDiagramBinary::new_pxy` -> `PhaseDiagram::binary_vle`
  - `PhaseDiagramBinary::new_txy_lle`, `PhaseDiagramBinary::new_pxy_lle` -> `PhaseDiagram::lle`
  - `PhaseDiagramHetero::new_txy`, `PhaseDiagramHetero::new_pxy` -> `PhaseDiagram::binary_vlle`\
  which still returns an instance of `PhaseDiagramHetero`
- Changed the internal implementation of the Peng-Robinson equation of state to use contributions like the more complex equations of state and removed the suggestion to overwrite the `evaluate_residual` function of `EquationOfState`. [#51](https://github.com/feos-org/feos-core/pull/51)

### Packaging
- Updated `pyo3` and `numpy` dependencies to 0.16.
- Updated `num-dual` dependency to 0.5.
- Updated `quantity` dependency to 0.5.

## [0.1.5] - 2022-02-21
### Fixed
- Fixed bug in `predict` of `Estimator`. [#30](https://github.com/feos-org/feos-core/pull/30)

### Added
- Add `pyproject.toml`. [#29](https://github.com/feos-org/feos-core/pull/29)

## [0.1.4] - 2022-02-18
### Fixed
- Fix state constructor for `T`, `p`, `V`, `x_i` specification. [#26](https://github.com/feos-org/feos-core/pull/26)

### Added
- Added method `predict` to `Estimator`. [#27](https://github.com/feos-org/feos-core/pull/27)

### Changed
- Changed method for vapor pressure in `DataSet` to `vapor_pressure` (was `pressure` of VLE liquid phase). [#27](https://github.com/feos-org/feos-core/pull/27)

## [0.1.3] - 2022-01-21
### Added
- Added the following properties to `State`: [#21](https://github.com/feos-org/feos-core/pull/21)
  - `dp_drho` partial derivative of pressure w.r.t. density
  - `d2p_drho2` second partial derivative of pressure w.r.t. density
  - `isothermal_compressibility` the isothermal compressibility
- Read a list of segment records directly from a JSON file. [#22](https://github.com/feos-org/feos-core/pull/22)


## [0.1.2] - 2022-01-10
### Changed
- Changed `ChemicalRecord` to an enum that can hold either the full structural information of a molecule or only segment and bond counts and added an `Identifier`. [#19](https://github.com/feos-org/feos-core/pull/19)
- Removed the `chemical_record` field from `PureRecord` and made `model_record` non-optional. [#19](https://github.com/feos-org/feos-core/pull/19)

## [0.1.1] - 2021-12-22
### Added
- Added `from_multiple_json` function to `Parameter` trait that is able to read parameters from separate JSON files. [#15](https://github.com/feos-org/feos-core/pull/15)

### Packaging
- Updated `pyo3` and `numpy` dependencies to 0.15.
- Updated `quantity` dependency to 0.4.
- Updated `num-dual` dependency to 0.4.
- Removed `ndarray-linalg` and `ndarray-stats` dependencies.
- Removed obsolete features for the selection of the BLAS/LAPACK library.

## [0.1.0] - 2021-12-02
### Added
- Initial release
