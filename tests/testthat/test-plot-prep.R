test_that("plot prep helpers return expected shapes", {
  skip_if_not(requireNamespace("mcmcCheckConvergence", quietly = TRUE))
  data_dir <- testthat::test_path("..", "test_1")
  skip_if_not(file.exists(data_dir))

  list_files <- list.files(data_dir, full.names = TRUE)
  output <- mcmcCheckConvergence::checkConvergence(
    list_files = list_files,
    format = "revbayes",
    control = mcmcCheckConvergence::makeControl(emitLogs = FALSE, fastSplits = TRUE)
  )

  ks_prep <- mcmcCheckConvergence:::plot_ks_data(output, 0.01)
  expect_true(is.list(ks_prep))
  expect_true(length(ks_prep$ks_values) > 0)
  expect_true(length(ks_prep$breaks) > 1)
  expect_length(ks_prep$xlim, 2)
  expect_length(ks_prep$ylim, 2)

  ks_pool <- mcmcCheckConvergence:::plot_ks_pooled_data(output, 0.01)
  expect_true(is.list(ks_pool))
  expect_true(length(ks_pool$ks_values) >= 1)
  expect_true(length(ks_pool$breaks) > 1)
  expect_length(ks_pool$xlim, 2)
  expect_length(ks_pool$ylim, 2)

  ess_split_all <- mcmcCheckConvergence:::plot_ess_splits_data(output, FALSE, 0.01, NULL)
  expect_true(is.list(ess_split_all))
  expect_true(length(ess_split_all$values) > 0)
  expect_true(length(ess_split_all$breaks) > 1)

  ess_split_runs <- mcmcCheckConvergence:::plot_ess_splits_data(output, TRUE, 0.01, NULL)
  expect_true(is.list(ess_split_runs))
  expect_true(length(ess_split_runs$runs) >= 1)

  ess_cont_all <- mcmcCheckConvergence:::plot_ess_continuous_data(output, FALSE, 0.01, NULL)
  expect_true(is.list(ess_cont_all))
  expect_true(length(ess_cont_all$values) > 0)
  expect_true(length(ess_cont_all$breaks) > 1)

  ess_cont_runs <- mcmcCheckConvergence:::plot_ess_continuous_data(output, TRUE, 0.01, NULL)
  expect_true(is.list(ess_cont_runs))
  expect_true(length(ess_cont_runs$runs) >= 1)

  diff_prep <- mcmcCheckConvergence:::plot_diff_splits_data(output, 0.01)
  expect_true(is.list(diff_prep))
  expect_true(length(diff_prep$points) >= 1)
  expect_true(is.matrix(diff_prep$expected))
  expect_equal(nrow(diff_prep$expected), 2)
})
