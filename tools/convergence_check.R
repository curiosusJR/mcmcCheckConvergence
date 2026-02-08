#!/usr/local/bin/Rscript --

### Check for convergence ###
suppressPackageStartupMessages(library("convenience"))

args <- commandArgs(trailingOnly = TRUE)
flags <- args[grepl("^--", args)]
data_args <- args[!grepl("^--", args)]
data_dir <- if (length(data_args) > 0) data_args[1] else "output"

name_files <- c(
  file.path(data_dir, "posterior_run_1.log"),
  file.path(data_dir, "posterior_run_2.log")
  # file.path(data_dir, "posterior_run_1.trees"),
  # file.path(data_dir, "posterior_run_2.trees")
)


check_convergence <- checkConvergence(list_files = name_files)

cat(as.character(check_convergence$converged))
cat(as.character(check_convergence$message_complete))
cat(as.character(check_convergence$burnin))
cat(as.character(check_convergence$failed_names))
