# Repository Guidelines

## Project Structure & Module Organization
This repository is an R package with a Rust backend. Key paths:
- `R/`: R-facing functions and wrappers.
- `man/`: R documentation (Rd files).
- `src/`: native build outputs and entrypoints (see `src/entrypoint.c`).
- `src/rust/`: Rust crate (`Cargo.toml`, `src/`).
- `ref/`: reference R/C/C++ implementation used for parity checks.
- `tools/`: any auxiliary scripts (if added in future).

## Build, Test, and Development Commands
- `R CMD build .`: build the R package tarball from the repository root.
- `R CMD check .`: run standard R package checks locally.
- `cargo build --manifest-path src/rust/Cargo.toml`: compile the Rust static library.
- `cargo test --manifest-path src/rust/Cargo.toml`: run Rust unit tests (if present).
Use `ref/` for validation against legacy behavior, not for building the new crate.

## Coding Style & Naming Conventions
- R: follow base R style (2-space indent, `snake_case` for functions/objects).
- Rust: standard `rustfmt` style (4-space indent, `snake_case` for functions, `CamelCase` for types).
- Keep Rust FFI symbols stable and descriptive; mirror R names where feasible.

## Testing Guidelines
No dedicated test directory is present yet. If adding tests:
- R: use `tests/testthat/` with `testthat` naming (`test-*.R`).
- Rust: add unit tests in `src/rust/src/` with `#[cfg(test)]`.
- Aim for parity tests comparing `ref/` outputs and the Rust backend.

## Commit & Pull Request Guidelines
This directory is not currently a git checkout, so no historical conventions are visible.
When contributing, use clear, imperative commit summaries (e.g., “Add Rust FFI wrapper”).
Pull requests should include:
- A concise description of behavior changes.
- Any affected R function names and matching Rust modules.
- Test commands run (e.g., `R CMD check .`, `cargo test ...`).

## Reference & Parity Notes
The `ref/` folder is the authoritative baseline for algorithm behavior and R API
shape. Keep a mapping between `ref/R` exports and the Rust FFI entrypoints as you
port functionality to avoid regressions.
