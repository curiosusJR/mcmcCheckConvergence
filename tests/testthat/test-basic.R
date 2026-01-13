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
