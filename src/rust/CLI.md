# mcmcCheckConvergence Rust CLI Guide

This document describes the standalone Rust command line tool
`convergence_cli`. It runs convergence checks without an R runtime.

## Build

From the repository root:

```sh
cargo build --manifest-path src/rust/Cargo.toml --bin convergence_cli
```

The binary is created at:
`src/rust/target/debug/convergence_cli`

Release build:
```sh
cargo build --manifest-path src/rust/Cargo.toml --bin convergence_cli --release
```
Release binary:
`src/rust/target/release/convergence_cli`

## Usage

```sh
convergence_cli --files f1,f2[,f3,f4] [options]
convergence_cli --path /path/to/output [options]
```

One of `--files` or `--path` is required.

### Examples

Explicit file list:

```sh
./src/rust/target/debug/convergence_cli \
  --files tests/test_1/posterior_run_1.log,tests/test_1/posterior_run_2.log,tests/test_1/posterior_run_1.trees,tests/test_1/posterior_run_2.trees
```

Directory scan:

```sh
./src/rust/target/debug/convergence_cli --path tests/test_1
```

Continuous-only (log files only):

```sh
./src/rust/target/debug/convergence_cli --path tests/test_1 --continuous-only
```

## Options

- `--format <revbayes>`: input format (currently only `revbayes`).
- `--burnin <fraction|percent>`: burn-in fraction (0-1) or percent (> 1).
- `--precision <float>`: ESS precision threshold.
- `--tracer <true|false>`: enable ESS tracer output.
- `--namesToExclude <regex>`: regex for column names to ignore.
- `--continuous-only`: ignore tree files and analyze logs only.
- `--json`: emit a single JSON object.
- `--tsv`: emit key/value rows.
- `--message-only`: print only the full convergence message.
- `--quiet`: print only converged and burnin.
- `--help` / `-h`: show usage.

## Input Rules Reminder

When using `--path` with `--format revbayes`, the CLI expects:
- `_run_1` and `_run_2` file stems.
- Two `.log` and two `.trees` files (unless `--continuous-only`).

The CLI sorts files lexicographically before running checks.
When using `--files`, provide a comma-separated list in any order; the CLI
sorts them before processing.

## Output Formats

Default output:
```
converged: true|false
burnin: <float>
message:
<full message>
```

JSON output (one line):
```
{"converged":true,"burnin":0.25,"message":"...","message_complete":"..."}
```

TSV output:
```
converged    true
burnin       0.25
message      ...
message_complete     ...
```

### Output Fields

- `converged`: `true`/`false`.
- `burnin`: numeric burn-in used for the run.
- `message`: brief summary.
- `message_complete`: full message including excluded parameters/splits.

## Exit Codes

- `0`: success.
- `1`: errors such as missing inputs or parse failures.

## Environment Controls

- `MCMC_CONVERGENCE_THREADS=1` forces single-thread execution in the Rust core.
- `RAYON_NUM_THREADS=1` also forces Rayon to one thread if Rayon is enabled.

## Troubleshooting

- "expected _run_1/_run_2": RevBayes inputs need paired run files.
- "no .log files provided": use `--continuous-only` with `.log` inputs.
- Missing output: check file paths and extension case sensitivity.
