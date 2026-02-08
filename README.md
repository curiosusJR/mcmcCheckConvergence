# mcmcCheckConvergence

This repository is a full AI-assisted rewrite of the original R/C++ convergence checker [`convenience`](https://github.com/lfabreti/convenience), preserving the R API while moving core logic into Rust.

## Goals

- Match the R-facing function signatures from the original package (e.g., `checkConvergence`, `essTracerC`).
- Remove external R package dependencies by moving core logic into Rust or base R.
- Keep test fixtures and outputs comparable to the original package.
- Boost performance for large MCMC traces by parallelizing key computations.

## Project Status

- Large-scale smoke testing: run against ~9,000 empirical MCMC trace datasets.
- Current results: it always shows the same value in all tests except only 2 edge cases, both due to differing variable column counts between runs.
- Code review: not yet reviewed by a human.

### Known Differences (CLI vs Rscript)

- Burn-in semantics: `--burnin 0` is a fixed burn-in value. Automatic burn-in only runs when `--burnin` is omitted or set to `--burnin auto` / `--burnin -1`. In the R API, `makeControl(burnin = NULL)` uses automatic burn-in.
- The log-file-only test suite compares CLI output to `tools/convergence_check.R` and normalizes noisy fragments (e.g., `list(...)`, `c(...)`, repeated `=`) in the Rscript message before comparing.
- Merged RevBayes logs must have the same set of variable columns per run. If a column is constant in one run and variable in another, the CLI errors with a message that names the constant column(s).

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
- `tests/`: format fixtures and expected outputs (each `tests/test_*` directory).
  - `tests/test_*/posterior_run_*.log|.trees`: RevBayes input data.
  - `tests/test_format_merge/posterior.{log,trees}`: RevBayes merged traces with `Replicate_ID`.
  - `tests/test_format_mb/*.p|*.t`: MrBayes input data.
  - `tests/test_*/convergence_*.txt`: expected outputs for parity.
- `tools/`: helper scripts for build and testing.

## Build Script

Use the helper script for common build flows:

```sh
./build.sh                # build CLI (no R bindings) + R tarball
./build.sh --cli-only      # build CLI only (without R bindings)
./build.sh --r-only        # build R tarball only
./build.sh --check         # run R CMD check --no-manual after build
```

## Install Notes

- R package: install the tarball produced by `./build.sh` with `R CMD INSTALL` or `install.packages(..., repos = NULL, type = "source")`.
- CLI: `./build.sh --cli-only` produces `src/rust/target/release/convergence_cli` without R bindings.
- CLI-only (no R bindings): `cargo build --manifest-path src/rust/Cargo.toml --bin convergence_cli --release --no-default-features`.
- R-only builds do not produce the CLI binary; they only build the Rust static library for the R package.

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
), format = "revbayes", control = makeControl(threads = 4, fastSplits = TRUE))
```

CLI (no R runtime needed):

```sh
./src/rust/target/release/convergence_cli --path tests/test_1 --format revbayes -j 4 --fast-splits
```

MrBayes example:

```sh
./src/rust/target/release/convergence_cli --path tests/test_format_mb --format mrbayes
```

RevBayes merged-trace example:

```sh
./src/rust/target/release/convergence_cli --path tests/test_format_merge --format revbayes
```

Threading defaults to a fixed `N` captured at compile time (based on available cores during compilation). Use `-j/--threads` or `threads` in `makeControl()` to override. Run loading and burn-in scanning parallelize automatically when logs are quiet.

Rust crate reuse (git dependency):

```toml
[dependencies]
mcmcCheckConvergence = { git = "https://github.com/curiosusjr/mcmcCheckConvergence", package = "mcmcCheckConvergence" }
```

```rust
use mcmcCheckConvergence::core::{check_convergence, Control};
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
- `tools/test-suite.py`: CLI-based checker that compares outputs against `output/` reference files for multiple datasets.

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

## Test Suites

Quick ways to exercise the fixtures under `tests/`:

```sh
R CMD check . --no-manual
Rscript tools/compare_tests.R
Rscript tools/test-ess_tracer.R tests/test_1
tools/test-ess_tracer.sh tests/test_1
python3 tools/test-suite.py tools/test-suites.txt
```

The list file passed to `tools/test-suite.py` should contain one dataset
directory per line (each with an `output/` subdirectory). Lines beginning with
`#` are ignored.

`tools/test-suite.py` prints a `RESULTS` summary (PASSED/FAILED/SKIP) and writes
full failure details to `test-suite.failures.txt`.

Rust-only checks:

```sh
cargo test --manifest-path src/rust/Cargo.toml
./src/rust/target/release/convergence_cli --path tests/test_format_mb --format mrbayes
./src/rust/target/release/convergence_cli --path tests/test_format_merge --format revbayes
```
