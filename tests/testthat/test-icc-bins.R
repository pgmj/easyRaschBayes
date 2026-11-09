# Class intervals for plot_icc(), on the pcmdat2 sum-score distribution
# (non-extreme scores 1-7, maximum score 8).
score <- rep(1:7, times = c(4, 23, 35, 80, 35, 51, 28))

test_that("quantile grouping reports fewer groups when scores tie", {
  b4 <- .icc_class_bins(score, "quantile", 4L)
  expect_equal(b4$n_groups, 3L)
  expect_equal(as.vector(table(b4$interval)), c(142L, 86L, 28L))
  expect_false(b4$fallback)
  expect_match(.icc_interval_caption("quantile", 4L, NULL, b4, 5),
               "3 groups.*4 were requested")
  b5 <- .icc_class_bins(score, "quantile", 5L)
  expect_equal(b5$n_groups, 5L)
  expect_match(.icc_interval_caption("quantile", 5L, NULL, b5, 5),
               "5 groups of about equal size \\(quantile\\)")
})

test_that("width, score and manual grouping", {
  w <- .icc_class_bins(score, "width", 3L)
  expect_equal(w$n_groups, 3L)
  # Edges 0.5, 2.83, 5.17, 7.5: scores 1-2, 3-5, 6-7
  expect_equal(as.vector(table(w$interval)), c(27L, 150L, 79L))
  s <- .icc_class_bins(score, "score", 4L)
  expect_equal(s$n_groups, 7L)
  expect_equal(s$interval, score)
  m <- .icc_class_bins(score, "manual", NA, c(4, 5), max_score = 8)
  expect_equal(m$n_groups, 3L)
  expect_equal(as.vector(table(m$interval)), c(62L, 80L, 114L))
  expect_match(.icc_interval_caption("manual", NA, c(4, 5), m, 5),
               "starting at sum scores 4, 5 \\(3 groups\\)")
})

test_that("empty groups and fall back are reported", {
  m <- .icc_class_bins(c(1, 1, 2, 2), "manual", NA, c(3, 5), max_score = 8)
  expect_equal(m$n_groups, 3L)
  expect_equal(m$n_used, 1L)
  expect_match(.icc_interval_caption("manual", NA, c(3, 5), m, 5),
               "1 of them contain persons")
  f <- .icc_class_bins(rep(3, 10), "quantile", 4L)
  expect_true(f$fallback)
  expect_equal(f$n_groups, 1L)
  expect_match(.icc_interval_caption("quantile", 4L, NULL, f, 5),
               "too few distinct sum scores")
})

test_that(".hpcm_invert_curve finds theta where the curve equals the target", {
  theta <- seq(-4, 4, length.out = 401)
  e <- 1 + 2 * stats::plogis(theta)
  target <- 1 + 2 * stats::plogis(c(-1.5, 0, 0.7))
  expect_equal(.hpcm_invert_curve(target, theta, e), c(-1.5, 0, 0.7),
               tolerance = 1e-3)
  # outside the curve's range: nearest end; NA stays NA
  expect_equal(.hpcm_invert_curve(c(0.5, 3.5, NA), theta, e), c(-4, 4, NA))
})
