# mcmcCheckConvergence

This repository is a full AI-assisted rewrite of the original R/C/C++ convergence checker [`convenience`](https://github.com/lfabreti/convenience), preserving the R API while moving core logic into Rust.

## Goals

- Match the R-facing function signatures from the original package (e.g., `checkConvergence`, `essTracerC`).
- Remove external R package dependencies by moving core logic into Rust or base R.
- Keep test fixtures and outputs comparable to the original package.
- Boost performance for large MCMC traces by parallelizing key computations.

## Project Structure

Top-level layout with key files:

- `R/` (R layer): public API wrappers and base-R orchestration.
  - `R/convergence.R`: main R entrypoints and S3 print helpers.
- `src/rust/` (Rust layer): core implementations and data assets.
  - `src/rust/src/core.rs`: convergence logic, burn-in, ESS/KS, formatting.
  - `src/rust/src/lib.rs`: extendr exports and tree utilities.
  - `src/rust/data/`: expected-diff threshold tables used by split tests.
- `src/` (C layer): minimal entrypoint and build glue for R.
  - `src/entrypoint.c`: registers Rust symbols for R.
  - `src/Makevars*`: generated build configuration.
- `tests/`: RevBayes fixtures and expected outputs (each `tests/test_*` directory).
  - `tests/test_*/posterior_run_*.log|.trees`: input data.
  - `tests/test_*/convergence_*.txt`: expected outputs for parity.
- `tools/`: helper scripts for build and testing.

## Build Script

Use the helper script for common build flows:

```sh
./build.sh                # build CLI + R tarball
./build.sh --cli-only      # build CLI only
./build.sh --r-only        # build R tarball only
./build.sh --check         # run R CMD check --no-manual after build
```

## Install Notes

- R package: install the tarball produced by `./build.sh` with `R CMD INSTALL` or `install.packages(..., repos = NULL, type = "source")`.
- CLI: `./build.sh --cli-only` produces `src/rust/target/release/convergence_cli`.

## Dependencies

Runtime:

- R (>= 4.2) for the R package interface.

Build-time:

- Rust toolchain: `cargo`, `rustc` (>= 1.65) for compiling the Rust core.
- A C compiler toolchain for building the R shared library.

Rust crate deps:

- `extendr-api` (R <-> Rust bindings).
- `regex` (header filtering).

## Usage

R API (same names as the original package):

```r
checkConvergence(list_files = c(
  "tests/test_1/posterior_run_1.log",
  "tests/test_1/posterior_run_2.log",
  "tests/test_1/posterior_run_1.trees",
  "tests/test_1/posterior_run_2.trees"
), format = "revbayes")
```

CLI (no R runtime needed):

```sh
./src/rust/target/release/convergence_cli --path tests/test_1 --format revbayes
```

## Detailed Docs

- R package guide: `inst/doc/convergence_guide.md`
- Rust CLI guide: `src/rust/CLI.md`

## Tooling Scripts

- `tools/compare_tests.R`: compares current R output vs `tools/convergence_check.R`, and checks Rust CLI parity.
- `tools/bench_convergence.R`: lightweight benchmark runner for `checkConvergence` on a given `tests/` dataset.
- `tools/convergence_check.R`: reference-style R script used for parity comparisons.
- `tools/convergence_cli.R`: R helper for invoking the Rust CLI in scripts.
- `tools/config.R`: generates `src/Makevars` during `R CMD INSTALL`.
- `tools/msrv.R`: records the minimum supported Rust version used by the crate.
- `tools/test-ess_tracer.R`: run ESS and convergence checks on a test directory.
- `tools/test-ess_tracer.sh`: shell wrapper for the Rust CLI test.

## Parity Rules

- `tools/compare_tests.R` runs current `checkConvergence` and the legacy-style `tools/convergence_check.R` for each dataset under `tests/`.
- Reference files live under `tests/test_*/output/` and include `convergence_assessment.txt`, `convergence_burnin.txt`, and `convergence_failedNames.txt`.
- If reference files are missing, `compare_tests.R` infers values from `convergence_check.R` stdout.
- CLI parity uses `convergence_cli --tsv` and compares `converged`, `burnin`, and normalized `message_complete` against the R API output.

## Convergence Test Script

The script runs ESS and a full convergence check using the Rust-backed API and prints the full message to stdout:

```sh
Rscript tools/test-ess_tracer.R tests/test_1
```

Pure Rust CLI equivalent:

```sh
tools/test-ess_tracer.sh tests/test_1
```

It expects `posterior_run_1.log`, `posterior_run_2.log`, `posterior_run_1.trees`, and `posterior_run_2.trees` in the directory you pass.
