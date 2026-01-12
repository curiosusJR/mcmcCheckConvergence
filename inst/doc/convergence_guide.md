# mcmcCheckConvergence R Package Guide

This guide explains the R-facing API, inputs, outputs, and common workflows for
`mcmcCheckConvergence`. The package wraps a Rust backend but keeps the original
R API names for compatibility.

## Supported Formats

Current format strings accepted by `checkConvergence`:
- `revbayes`
- `mb` / `mrbayes`
- `beast`
- `*beast`
- `phylobayes`
- `pyrate`

Formats differ by expected file extensions and headers; RevBayes is the default.

## Quickstart

```r
library(mcmcCheckConvergence)

output <- checkConvergence(
  list_files = c(
    "tests/test_1/posterior_run_1.log",
    "tests/test_1/posterior_run_2.log",
    "tests/test_1/posterior_run_1.trees",
    "tests/test_1/posterior_run_2.trees"
  ),
  format = "revbayes",
  control = makeControl(burnin = 0.25, precision = 0.01, tracer = TRUE)
)

print(output$message)
output$converged
output$burnin
```

## Inputs

The main entrypoint accepts either a directory path or an explicit file list:
- `path`: a directory containing `.log` and `.trees` files.
- `list_files`: explicit file paths, often two runs of `.log` and `.trees`.

RevBayes format expects paired run files with `_run_1` and `_run_2` stems.
If `path` is used, all files in the directory are scanned and filtered by
extension.

## Core API

### `checkConvergence(path = NULL, list_files = NULL, format = "revbayes", control = makeControl())`

Runs the convergence check and returns a structured result. Most users should
call this directly.

### `makeControl(tracer = NULL, burnin = NULL, precision = NULL, namesToExclude = NULL, emitLogs = NULL)`

Builds a list of control settings:
- `tracer`: `TRUE` to compute ESS tracer output for continuous parameters.
- `burnin`: fraction (0-1) or percent (> 1, treated as percent).
- `precision`: numeric threshold for ESS checks.
- `namesToExclude`: regex of column names to ignore.
- `emitLogs`: `TRUE` to keep log output from the Rust backend.

Default behaviors:
- `burnin = 0.0` triggers automatic burn-in estimation (in 10% steps up to 50%).
- `precision = 0.01` sets the ESS threshold used across tests.
- `namesToExclude` defaults to a RevBayes-friendly regex that drops likelihood,
  posterior, and bookkeeping columns.

## Output Structure

The return value is a list with these fields:
- `converged`: `TRUE`/`FALSE`.
- `burnin`: numeric burn-in value used for the analysis.
- `message`: summary message (class `list.fails`).
- `message_complete`: full multi-line message (class `list.fails`).
- `tree_parameters`: list of tree split stats.
- `continuous_parameters`: list of continuous parameter stats.
- `failed`: optional list of failed checks.
- `failed_names`: optional vector of parameter names that failed.

### `tree_parameters`

Fields inside `tree_parameters`:
- `frequencies`: split frequencies (per run or aggregated depending on input).
- `ess`: ESS table per split.
- `compare_runs`: comparisons between runs (RevBayes usually two runs).
- `freq_per_run`: per-run split frequencies.

### `continuous_parameters`

Fields inside `continuous_parameters`:
- `means`: per-parameter means.
- `ess`: per-parameter ESS.
- `compare_runs`: per-parameter comparisons between runs.

## Helper Output Functions

These functions return data frames suitable for inspection or saving:
- `printTableSplits(output, splits_per_run = FALSE, filename = NULL)`
- `printTableContinuous(output, filename = NULL)`
- `printConvergenceTable(output)`
- `printConvergenceDiag(output)`

## Plotting Helpers

Plot helpers are kept for compatibility with the original package:
`plotKS`, `plotKSPooled`, `plotEssSplits`, `plotEssContinuous`,
`plotDiffSplits`.

## Common Workflows

### Directory-based run

```r
output <- checkConvergence(
  path = "tests/test_1",
  format = "revbayes"
)
```

### Excluding parameters

```r
control <- makeControl(namesToExclude = "alpha|beta|rate")
output <- checkConvergence(list_files = files, control = control)
```

### Saving tables

```r
printTableSplits(output, filename = "splits.csv")
printTableContinuous(output, filename = "continuous.csv")
```

### Inspecting failures

```r
if (!output$converged) {
  output$failed_names
  cat(output$message_complete)
}
```

### Using explicit burn-in

```r
control <- makeControl(burnin = 0.25)
output <- checkConvergence(list_files = files, control = control)
```

## Performance and Determinism

- The Rust core uses all available CPU cores by default.
- Set `MCMC_CONVERGENCE_THREADS=1` to force single-thread execution.
- For fully deterministic profiling, also set `RAYON_NUM_THREADS=1`.

## Troubleshooting

- "Provide path or list_files": pass either `path` or `list_files`.
- Missing run files: RevBayes expects `_run_1` and `_run_2` for both `.log`
  and `.trees`.
- Non-numeric values: input tables must be numeric after header rows; `NA` or
  `NaN` values will stop the run.
- Large datasets: set `emitLogs = FALSE` to reduce output verbosity.
