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

CLI-only build without R bindings:
```sh
cargo build --manifest-path src/rust/Cargo.toml --bin convergence_cli --release --no-default-features
```


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

- `--format <name>`: input format (`revbayes`, `mrbayes`/`mb`, `beast`, `*beast`,
  `phylobayes`, `pyrate`).
- `--burnin <fraction|percent>`: fixed burn-in fraction (0-1) or percent (> 1).
- `--burnin auto` or `--burnin -1`: enable automatic burn-in estimation.
- `--precision <float>`: ESS precision threshold.
- `--tracer <true|false>`: enable ESS tracer output.
- `--namesToExclude <regex>`: regex for column names to ignore.
- `-j`, `--threads <n>`: number of threads for Rust parallelism. Defaults to a fixed `N` captured at compile time (based on available cores during compilation).
- `--fast-splits`: use a faster bitset-based split backend.
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

When using `--path` with `--format mrbayes`/`mb`, the CLI expects:
- `.p` and `.t` files (e.g., `*.run1.p`, `*.run1.t`, `*.run2.p`, `*.run2.t`).

The CLI sorts files lexicographically before running checks.
When using `--files`, provide a comma-separated list in any order; the CLI
sorts them before processing.

## Test Fixtures

The repository includes format-specific fixtures under `tests/`:
- `tests/test_1` and `tests/test_2`: RevBayes `_run_1/_run_2` examples.
- `tests/test_format_merge`: RevBayes merged traces with `Replicate_ID`.
- `tests/test_format_mb`: MrBayes `.p/.t` examples.

Example runs:

```sh
./src/rust/target/debug/convergence_cli --path tests/test_format_mb --format mrbayes
./src/rust/target/debug/convergence_cli --path tests/test_format_merge --format revbayes
```

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
{"converged":true,"burnin":0.25,"message":"...","message_complete":"...","failed_names":"..."}
```

TSV output:
```
converged    true
burnin       0.25
message      ...
message_complete     ...
failed_names ...
```

### Output Fields

- `converged`: `true`/`false`.
- `burnin`: numeric burn-in used for the run.
- `message`: brief summary.
- `message_complete`: full message including excluded parameters/splits.
- `failed_names`: comma-separated list of failed check labels (may be empty).

## Exit Codes

- `0`: success.
- `1`: errors such as missing inputs or parse failures.

## Troubleshooting

- "expected _run_1/_run_2": RevBayes inputs need paired run files.
- "no .log files provided": use `--continuous-only` with `.log` inputs.
- "Filtered continuous parameter columns differ between runs": for merged logs, a column can be constant in one run and variable in another; the error lists the constant column(s).
- Missing output: check file paths and extension case sensitivity.
