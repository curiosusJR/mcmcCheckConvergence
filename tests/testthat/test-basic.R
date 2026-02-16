test_that("makeControl returns a list with expected names", {
  skip_if_not(requireNamespace("mcmcCheckConvergence", quietly = TRUE))
  control <- mcmcCheckConvergence::makeControl()
  expect_true(is.list(control))
  expect_equal(names(control), c("tracer", "burnin", "precision", "namesToExclude", "emitLogs", "threads", "fastSplits"))
})

test_that("checkConvergence errors without inputs", {
  skip_if_not(requireNamespace("mcmcCheckConvergence", quietly = TRUE))
  expect_error(mcmcCheckConvergence::checkConvergence(), "Provide path or list_files")
})

test_that("checkConvergence errors on identical replicate chains", {
  skip_if_not(requireNamespace("mcmcCheckConvergence", quietly = TRUE))
  log_path <- tempfile("identical_reps_", fileext = ".log")
  data <- paste(
    "Iteration\tReplicate_ID\tx",
    "1\t0\t0.1",
    "2\t0\t0.2",
    "1\t1\t0.1",
    "2\t1\t0.2",
    sep = "\n"
  )
  writeLines(data, con = log_path)
  on.exit(unlink(log_path), add = TRUE)

  expect_error(
    mcmcCheckConvergence::checkConvergence(
      list_files = log_path,
      format = "revbayes",
      control = mcmcCheckConvergence::makeControl(emitLogs = FALSE)
    ),
    "Detected identical MCMC runs"
  )
})
