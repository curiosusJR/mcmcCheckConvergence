#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) > 0) args[1] else "tests/test_1"
iters <- if (length(args) > 1) as.integer(args[2]) else 3L

if (!requireNamespace("mcmcCheckConvergence", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    cmd <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd, value = TRUE)
    script_path <- if (length(file_arg) == 0) normalizePath(getwd()) else normalizePath(sub("^--file=", "", file_arg[1]))
    repo_root <- normalizePath(file.path(dirname(script_path), ".."))
    suppressMessages(devtools::load_all(repo_root))
  } else {
    stop("Package not installed. Run `R CMD INSTALL .` or install devtools.")
  }
} else {
  suppressPackageStartupMessages(library(mcmcCheckConvergence))
}

name_files <- c(
  file.path(data_dir, "posterior_run_1.log"),
  file.path(data_dir, "posterior_run_2.log"),
  file.path(data_dir, "posterior_run_1.trees"),
  file.path(data_dir, "posterior_run_2.trees")
)

if (!all(file.exists(name_files))) {
  stop(paste("Missing required files in", data_dir))
}

control <- makeControl(emitLogs = FALSE)

invisible(checkConvergence(list_files = name_files, format = "revbayes", control = control))

times <- numeric()
for (i in seq_len(iters)) {
  t <- system.time(checkConvergence(list_files = name_files, format = "revbayes", control = control))
  times <- c(times, t["elapsed"])
}

cat("Benchmark runs:", iters, "\n")
cat("Elapsed (s):", paste(sprintf("%.4f", times), collapse = ", "), "\n")
cat("Mean elapsed (s):", sprintf("%.4f", mean(times)), "\n")
