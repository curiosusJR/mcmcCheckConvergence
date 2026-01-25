#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }
  normalizePath(sub("^--file=", "", file_arg[1]))
}

normalize_text <- function(x) gsub("\\s+", " ", trimws(x))
normalize_bool <- function(x) toupper(trimws(x))

compare_numeric <- function(a, b, tol = 1e-8) {
  na <- suppressWarnings(as.numeric(a))
  nb <- suppressWarnings(as.numeric(b))
  if (is.na(na) || is.na(nb)) return(FALSE)
  abs(na - nb) <= tol
}

load_current_pkg <- function(repo_root) {
  if (!requireNamespace("mcmcCheckConvergence", quietly = TRUE)) {
    if (requireNamespace("devtools", quietly = TRUE)) {
      suppressMessages(devtools::load_all(repo_root))
    } else {
      stop("Package not installed. Run `R CMD INSTALL .` or install devtools.")
    }
  } else {
    suppressPackageStartupMessages(library(mcmcCheckConvergence))
  }
}

detect_format <- function(output_dir) {
  files <- list.files(output_dir, full.names = FALSE)
  has_p <- any(grepl("\\.p$", files))
  has_t <- any(grepl("\\.t$", files))
  has_log <- any(grepl("\\.log$", files))
  has_trees <- any(grepl("\\.trees$", files))
  has_trace <- any(grepl("\\.trace$", files))
  has_treelist <- any(grepl("\\.treelist$", files))

  if (has_p || has_t) {
    return("mrbayes")
  }
  if (has_trace || has_treelist) {
    return("phylobayes")
  }
  if (has_log || has_trees) {
    return("revbayes")
  }
  NULL
}

read_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  paste(readLines(path, warn = FALSE), collapse = " ")
}

flatten_failed <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.list(x)) return(unlist(x, recursive = TRUE, use.names = FALSE))
  as.character(x)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript tools/test-suite.R <dir> [<dir> ...]\n")
  cat("Each <dir> should contain an output/ subdirectory with trace files.\n")
  quit(status = 1)
}

script_path <- get_script_dir()
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
load_current_pkg(repo_root)

any_fail <- FALSE

for (dir_path in args) {
  dir_path <- normalizePath(dir_path, mustWork = FALSE)
  output_dir <- file.path(dir_path, "output")
  if (!dir.exists(output_dir)) {
    cat("SKIP:", dir_path, "(missing output/)\n")
    next
  }

  format <- detect_format(output_dir)
  if (is.null(format)) {
    cat("SKIP:", output_dir, "(no recognized trace files)\n")
    next
  }

  cat("RUN:", output_dir, "format=", format, "\n", sep = " ")

  control <- makeControl(emitLogs = FALSE)
  res <- tryCatch(
    checkConvergence(path = output_dir, format = format, control = control),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    cat("  ERROR:", res$message, "\n")
    any_fail <- TRUE
    next
  }

  ref_assess <- read_if_exists(file.path(output_dir, "convergence_assessment.txt"))
  ref_burnin <- read_if_exists(file.path(output_dir, "convergence_burnin.txt"))
  ref_detail <- read_if_exists(file.path(output_dir, "convergence_detail.txt"))
  ref_failed <- read_if_exists(file.path(output_dir, "convergence_failedNames.txt"))

  got_assess <- if (isTRUE(res$converged)) "TRUE" else "FALSE"
  got_burnin <- as.character(res$burnin %||% NA_real_)
  got_detail <- as.character(res$message_complete %||% "")
  got_failed <- paste(flatten_failed(res$failed_names), collapse = " ")

  if (!is.null(ref_assess)) {
    ok <- normalize_bool(got_assess) == normalize_bool(ref_assess)
    cat("  convergence_assessment:", if (ok) "OK" else "DIFF", "\n")
    if (!ok) {
      cat("    expected:", normalize_text(ref_assess), "\n")
      cat("    got:     ", normalize_text(got_assess), "\n")
      any_fail <- TRUE
    }
  }

  if (!is.null(ref_burnin)) {
    ok <- compare_numeric(got_burnin, ref_burnin)
    cat("  convergence_burnin:", if (ok) "OK" else "DIFF", "\n")
    if (!ok) {
      cat("    expected:", normalize_text(ref_burnin), "\n")
      cat("    got:     ", normalize_text(got_burnin), "\n")
      any_fail <- TRUE
    }
  }

  if (!is.null(ref_detail)) {
    ok <- normalize_text(got_detail) == normalize_text(ref_detail)
    cat("  convergence_detail:", if (ok) "OK" else "DIFF", "\n")
    if (!ok) {
      any_fail <- TRUE
    }
  }

  if (!is.null(ref_failed)) {
    ok <- normalize_text(got_failed) == normalize_text(ref_failed)
    cat("  convergence_failedNames:", if (ok) "OK" else "DIFF", "\n")
    if (!ok) {
      any_fail <- TRUE
    }
  }
}

if (any_fail) {
  quit(status = 1)
}
