# Tests for the plot_targeting() helpers in R/targeting_helpers.R. The
# band geometry tests follow easyRasch2's tests for RMtargeting().

test_that("modal boundaries reduce to the thresholds when they are ordered", {
  tau <- c(-1.5, 0.2, 1.4)
  mb <- .modal_boundaries(tau)
  expect_equal(mb$bounds, tau)
  expect_equal(mb$cats, 0:3)
  expect_length(mb$dropped, 0L)
})

test_that("a disordered threshold pair drops the category between them", {
  mb <- .modal_boundaries(c(-1.2, -1.6, 0.4))
  expect_equal(mb$dropped, 1L)
  expect_equal(mb$cats, c(0L, 2L, 3L))
  expect_equal(mb$bounds[1], mean(c(-1.2, -1.6)))
})

test_that("a run of disordered thresholds drops two categories", {
  mb <- .modal_boundaries(c(-2.2, 0.6, 0.0, -0.6, 1.2, 2.2))
  expect_equal(mb$dropped, c(2L, 3L))
})

test_that("modal boundaries agree with the most likely category on a grid", {
  set.seed(2)
  for (rep in 1:20) {
    tau <- sort(rnorm(4)) + rnorm(4, 0, 0.8)
    mb <- .modal_boundaries(tau)
    theta <- seq(-6, 6, by = 0.001)
    cum <- outer(theta, 0:4) - matrix(c(0, cumsum(tau)), length(theta), 5,
                                      byrow = TRUE)
    modal <- max.col(cum, ties.method = "first") - 1L
    expect_setequal(unique(modal), mb$cats)
    # Changes of the modal category sit at the boundaries
    switch_at <- theta[which(diff(modal) != 0) + 1L]
    expect_length(switch_at, length(mb$bounds))
    # Absolute tolerance of two grid steps
    expect_true(all(abs(switch_at - mb$bounds) < 0.002))
  }
})

test_that("threshold reversals are found per item and NULL when ordered", {
  thr <- data.frame(Item = rep(c("a", "b"), each = 3), k = rep(1:3, 2),
                    Location = c(-1, 0, 1, 0.5, -0.5, 1))
  r <- .threshold_reversals(thr, c("a", "b"))
  expect_equal(nrow(r), 1L)
  expect_equal(r$Item, "b")
  expect_equal(r$gap, 1)
  thr$Location[5] <- 0.7
  expect_null(.threshold_reversals(thr, c("a", "b")))
})

test_that("posterior category structure counts draws", {
  # Item x: 3 of 4 draws reversed (category 1 never modal), item y ordered
  tau <- list(
    x = rbind(c(0, -1), c(0.5, -0.5), c(1, 0), c(-1, 1)),
    y = rbind(c(-1, 1), c(-1, 1), c(-1, 1), c(-1, 1))
  )
  cp <- .category_posterior(tau, prob = 0.5)
  cx <- cp$categories[cp$categories$item == "x", ]
  expect_equal(cx$p_never_modal, c(0, 0.75, 0))
  expect_true(is.na(cx$lower_bound[1]))
  # Category 1 is modal only in draw 4, starting at -1
  expect_equal(cx$lower_bound[2], -1)
  expect_equal(cp$disordered$p_disordered, c(0.75, 0))
  expect_equal(cp$disordered$thresholds, c("1-2", "1-2"))
})

test_that("probabilities are formatted without a leading zero", {
  expect_equal(.fmt_prob(c(0.834, 1, 0.001, 0)), c(".83", "1.00", "<.01", ".00"))
})

test_that("the category panel builds with red marks and probability labels", {
  thr <- data.frame(Item = rep(c("a", "b"), each = 2), k = rep(1:2, 2),
                    Location = c(-1, 1, 0.4, -0.4),
                    CI_low = c(-1.3, 0.7, 0.1, -0.7),
                    CI_high = c(-0.7, 1.3, 0.7, -0.1))
  drop_prob <- data.frame(item = "b", category = 1L, p_never_modal = 0.9)
  rev_prob <- data.frame(item = "b", thresholds = "1-2", p_disordered = 0.9)
  p <- .targeting_category_panel(
    thr, c("b", "a"), c(-4, 4), c("0", "1", "2"),
    c("#DEF5E5", "#3497A9", "#382A54"), TRUE, 1.44, "Note. x", "Location",
    drop_prob, rev_prob
  )
  expect_s3_class(p, "ggplot")
  labels <- unlist(lapply(p$layers, function(l) {
    if (inherits(l$geom, "GeomText") && "label" %in% names(l$data)) l$data$label
  }))
  # The collapse of category 1 comes from the reversal of thresholds 1-2,
  # so both are labelled together on the arrow and the tick has no text
  expect_true("1: 0.80 (P = .90)" %in% labels)
  expect_false("1 (P = .90)" %in% labels)
  expect_false("0.80 (P = .90)" %in% labels)
})
