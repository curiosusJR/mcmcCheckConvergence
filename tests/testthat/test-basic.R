test_that("makeControl returns a list with expected names", {
  control <- makeControl()
  expect_true(is.list(control))
  expect_equal(names(control), c("tracer", "burnin", "precision", "namesToExclude", "emitLogs"))
})

test_that("checkConvergence errors without inputs", {
  expect_error(checkConvergence(), "Provide path or list_files")
})
