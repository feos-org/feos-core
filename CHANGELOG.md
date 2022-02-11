# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
- Fix state constructor for $T$, $p$, $V$, $x_i$ specification. [#26](https://github.com/feos-org/feos-core/pull/26)

## [0.1.3] - 2022-01-21
### Added
- Added the following properties to `State`: [#21](https://github.com/feos-org/feos-core/pull/21)
  - `dp_drho` partial derivative of pressure w.r.t. pressure
  - `d2p_drho2` second partial derivative of pressure w.r.t. pressure
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
