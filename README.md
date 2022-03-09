# FeOs - A Framework for Equations of State

[![crate](https://img.shields.io/crates/v/feos-core.svg)](https://crates.io/crates/feos-core)
[![documentation](https://docs.rs/feos-core/badge.svg)](https://docs.rs/feos-core)
[![minimum rustc 1.51](https://img.shields.io/badge/rustc-1.51+-red.svg)](https://rust-lang.github.io/rfcs/2495-min-rust-version.html)

Core traits and functionalities for the `feos` project.

## Installation

Add this to your `Cargo.toml`

```toml
[dependencies]
feos-core = "0.2"
```

## Test building python wheel

From within a Python virtual environment with `maturin` installed, type

```
maturin build --release --out dist --no-sdist -m build_wheel/Cargo.toml
```
