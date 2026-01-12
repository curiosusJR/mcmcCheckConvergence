#!/usr/bin/env Rscript

get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }
  normalizePath(sub("^--file=", "", file_arg[1]))
}

script_path <- get_script_dir()
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
tests_root <- file.path(repo_root, "tests")

if (!dir.exists(tests_root)) {
  stop("Missing tests/ directory.")
}

test_dirs <- list.dirs(tests_root, full.names = TRUE, recursive = FALSE)
test_dirs <- c(tests_root, test_dirs)

required_files <- c(
  "posterior_run_1.log",
  "posterior_run_2.log",
  "posterior_run_1.trees",
  "posterior_run_2.trees"
)

run_cmd <- function(cmd, args, wd) {
  out <- tryCatch(
    system2(cmd, args = args, stdout = TRUE, stderr = TRUE, wait = TRUE),
    error = function(e) paste("ERROR:", e$message)
  )
  paste(out, collapse = "\n")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

find_cli <- function(repo_root) {
  candidates <- c(
    file.path(repo_root, "src", "rust", "target", "release", "convergence_cli"),
    file.path(repo_root, "src", "rust", "target", "debug", "convergence_cli")
  )
  for (path in candidates) {
    if (file.exists(path)) return(path)
  }
  NULL
}

parse_tsv_kv <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  out <- list()
  for (line in lines) {
    if (!nzchar(trimws(line))) next
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    key <- parts[1]
    value <- paste(parts[-1], collapse = "\t")
    out[[key]] <- value
  }
  out
}

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

summary_current <- function(dir_path) {
  name_files <- c(
    file.path(dir_path, "posterior_run_1.log"),
    file.path(dir_path, "posterior_run_2.log"),
    file.path(dir_path, "posterior_run_1.trees"),
    file.path(dir_path, "posterior_run_2.trees")
  )
  control <- makeControl(emitLogs = FALSE)
  res <- checkConvergence(list_files = name_files, format = "revbayes", control = control)

  failed_names <- res$failed_names
  failed_dump <- if (is.null(failed_names)) "" else paste(capture.output(dput(failed_names)), collapse = " ")
  failed_count <- 0
  count_failed <- function(x) {
    if (is.null(x)) return(0)
    if (is.list(x)) return(sum(vapply(x, count_failed, numeric(1))))
    length(x)
  }
  failed_count <- count_failed(failed_names)

  list(
    converged = as.character(res$converged),
    burnin = as.character(res$burnin),
    message_complete = as.character(res$message_complete),
    failed_names = failed_dump,
    failed_count = failed_count
  )
}

summary_reference <- function(output_dir, origin_out) {
  read_one <- function(path) {
    if (!file.exists(path)) return("")
    paste(readLines(path, warn = FALSE), collapse = " ")
  }
  converged <- if (nzchar(output_dir)) {
    read_one(file.path(output_dir, "convergence_assessment.txt"))
  } else {
    ""
  }
  burnin <- if (nzchar(output_dir)) {
    read_one(file.path(output_dir, "convergence_burnin.txt"))
  } else {
    ""
  }
  failed <- if (nzchar(output_dir)) {
    read_one(file.path(output_dir, "convergence_failedNames.txt"))
  } else {
    ""
  }

  inferred <- FALSE
  if (nchar(trimws(converged)) == 0) {
    if (grepl("ACHIEVED CONVERGENCE", origin_out, fixed = TRUE)) {
      converged <- "TRUE"
    } else if (grepl("FAILED CONVERGENCE", origin_out, fixed = TRUE)) {
      converged <- "FALSE"
    }
    inferred <- TRUE
  }
  if (nchar(trimws(burnin)) == 0) {
    burn_match <- regexpr("BURN-IN SET AT[ ]+([0-9.]+)", origin_out)
    if (burn_match[1] != -1) {
      burnin <- sub(".*BURN-IN SET AT[ ]+([0-9.]+).*", "\\1", origin_out)
    }
    inferred <- TRUE
  }

  failed_norm <- gsub("\\s+", " ", trimws(failed))
  failed_count <- if (nchar(failed_norm) == 0) {
    if (inferred) NA_integer_ else 0L
  } else {
    length(strsplit(failed_norm, "\\s+")[[1]])
  }
  list(
    converged = converged,
    burnin = burnin,
    failed_names = failed,
    failed_count = failed_count
  )
}

normalize_text <- function(x) gsub("\\s+", " ", trimws(x))

cat("Running stdout for datasets in tests/ ...\n")

load_current_pkg(repo_root)
cli_path <- find_cli(repo_root)

for (dir_path in test_dirs) {
  name <- basename(dir_path)
  if (name == "") name <- "tests"

  has_all <- all(file.exists(file.path(dir_path, required_files)))
  if (!has_all) {
    cat("Skipping:", dir_path, "(missing required files)\n")
    next
  }

  cat("Running:", dir_path, "\n")
  test_script <- file.path(repo_root, "tools", "test-ess_tracer.R")
  orig_script <- file.path(repo_root, "tools", "convergence_check.R")

  current_out <- run_cmd("Rscript", c(test_script, dir_path), dir_path)
  origin_out <- run_cmd("Rscript", c(orig_script, dir_path), dir_path)

  cat("  CURRENT (tools/test-ess_tracer.R):\n")
  cat(current_out, "\n")
  cat("  REFERENCE (tools/convergence_check.R):\n")
  cat(origin_out, "\n")

  current_summary <- summary_current(dir_path)
  ref_summary <- summary_reference("", origin_out)
  cat("  SUMMARY DIFF:\n")
  cat("    converged: current=", current_summary$converged, " ref=", ref_summary$converged, "\n", sep = "")
  cat("    burnin: current=", current_summary$burnin, " ref=", ref_summary$burnin, "\n", sep = "")
  cat("    failed_count: current=", current_summary$failed_count, " ref=", ref_summary$failed_count, "\n", sep = "")
  failed_cur <- normalize_text(current_summary$failed_names)
  failed_ref <- normalize_text(ref_summary$failed_names)
  if (is.na(ref_summary$failed_count) || nchar(failed_ref) == 0) {
    cat("    failed_names: SKIP\n")
  } else if (failed_cur != failed_ref) {
    cat("    failed_names: DIFF\n")
    cat("      current=", failed_cur, "\n", sep = "")
    cat("      ref=", failed_ref, "\n", sep = "")
  } else {
    cat("    failed_names: OK\n")
  }

  if (!is.null(cli_path)) {
    cli_out <- run_cmd(cli_path, c("--path", dir_path, "--format", "revbayes", "--tsv"), dir_path)
    cli_summary <- parse_tsv_kv(cli_out)
    cat("  CLI PARITY (Rust CLI vs R):\n")
    cli_converged <- cli_summary$converged %||% "MISSING"
    cli_burnin <- cli_summary$burnin %||% "MISSING"
    converged_ok <- normalize_bool(cli_converged) == normalize_bool(current_summary$converged)
    burnin_ok <- compare_numeric(cli_burnin, current_summary$burnin)
    cat(
      "    converged: cli=",
      cli_converged,
      " r=",
      current_summary$converged,
      " ",
      if (converged_ok) "OK" else "DIFF",
      "\n",
      sep = ""
    )
    cat(
      "    burnin: cli=",
      cli_burnin,
      " r=",
      current_summary$burnin,
      " ",
      if (burnin_ok) "OK" else "DIFF",
      "\n",
      sep = ""
    )
    cli_message <- normalize_text(cli_summary$message_complete %||% "")
    r_message <- normalize_text(current_summary$message_complete)
    if (!nzchar(cli_message)) {
      cat("    message_complete: SKIP (missing CLI output)\n")
    } else if (cli_message != r_message) {
      cat("    message_complete: DIFF\n")
    } else {
      cat("    message_complete: OK\n")
    }
  } else {
    cat("  CLI PARITY: SKIP (missing convergence_cli binary)\n")
  }
}
