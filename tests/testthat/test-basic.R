test_that("exported functions exist and are functions", {
  expect_true(is.function(dif_statistic))
  expect_true(is.function(fit_statistic_pcm))
  expect_true(is.function(fit_statistic_rm))
  expect_true(is.function(plot_residual_pca))
  expect_true(is.function(q3_statistic))
  expect_true(is.function(plot_ipf))
  expect_true(is.function(RMUreliability))
})
