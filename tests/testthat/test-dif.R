# Tests for the DIF extractors in R/brms_dif.R. The model fits themselves
# are too slow for the test suite, so the extractors get synthetic draws.

# Synthetic threshold draws for a non-uniform DIF model: thresholds
# tau[i, k] for the reference group, the focal group's thresholds shifted
# by -impact (the person distribution has mean 0 in both groups) plus
# dif_raw[i, k] (positive = harder for the focal group).
make_threshold_draws <- function(tau, impact, dif_raw, n_draws = 50) {
  out <- list()
  for (i in seq_len(nrow(tau))) {
    for (k in seq_len(ncol(tau))) {
      out[[sprintf("b_Intercept[i%dg1,%d]", i, k)]] <- rep(tau[i, k], n_draws)
      out[[sprintf("b_Intercept[i%dg2,%d]", i, k)]] <-
        rep(tau[i, k] - impact + dif_raw[i, k], n_draws)
    }
  }
  as.data.frame(out, check.names = FALSE)
}

test_that("threshold DIF separates impact, centres on items and flips sign", {
  tau <- matrix(c(-1, 0.5, -0.5, 1, 0, 1.5), ncol = 2, byrow = TRUE)
  dif_raw <- matrix(0, 3, 2)
  dif_raw[2, ] <- 0.9          # item 2 harder for the focal group
  d <- make_threshold_draws(tau, impact = 0.6, dif_raw = dif_raw)
  res <- .extract_dif_threshold(d, c("a", "b", "c"))
  # Mean item DIF is 0.3, so impact is 0.6 - 0.3 and DIF is relative to it
  expect_equal(unique(res$impact), 0.3)
  est <- colMeans(res$draws_mat)
  expect_equal(unname(est[res$labels$item == "b"]), c(-0.6, -0.6))
  expect_equal(unname(est[res$labels$item != "b"]), rep(0.3, 4))
  # DIF averages zero over items
  expect_equal(mean(tapply(est, res$labels$item, mean)), 0)
  expect_equal(res$labels$threshold, rep(c("1", "2"), 3))
})

test_that("threshold DIF with no DIF returns the impact and zeros", {
  tau <- matrix(c(-1, 0.5, -0.5, 1), ncol = 2, byrow = TRUE)
  d <- make_threshold_draws(tau, impact = 0.8, dif_raw = matrix(0, 2, 2))
  res <- .extract_dif_threshold(d, c("x", "y"))
  expect_equal(unique(res$impact), 0.8)
  expect_true(all(abs(res$draws_mat) < 1e-12))
})

test_that("threshold DIF warns on uneven categories and stops on a missing group", {
  tau <- matrix(c(-1, 0.5, -0.5, 1), ncol = 2, byrow = TRUE)
  d <- make_threshold_draws(tau, impact = 0, dif_raw = matrix(0, 2, 2))
  d[["b_Intercept[i2g2,2]"]] <- NULL
  expect_warning(.extract_dif_threshold(d, c("x", "y")), "same response categories")
  d2 <- make_threshold_draws(tau, impact = 0, dif_raw = matrix(0, 2, 2))
  d2 <- d2[, !grepl("i2g2", names(d2))]
  expect_error(.extract_dif_threshold(d2, c("x", "y")), "no responses in one")
})

test_that("slope DIF uses ranef labels, the fixed group effect and stops on gaps", {
  n_draws <- 20
  re <- array(NA_real_, dim = c(n_draws, 3, 1),
              dimnames = list(NULL, c("Item 1", "Item 2", "Item 3"), "difgroup"))
  re[, , 1] <- rep(c(0.2, -0.5, 0.3), each = n_draws)
  fe <- cbind("Intercept[1]" = 0, "Intercept[2]" = 1, difgroup = rep(0.7, n_draws))
  local_mocked_bindings(
    ranef = function(object, summary = TRUE, ...) list(item = re),
    fixef = function(object, summary = TRUE, ...) fe,
    .package = "brms"
  )
  res <- .extract_dif_slopes(NULL, c("Item 2", "Item 1", "Item 3"), "item",
                             fixed_coef = "difgroup", exact = TRUE)
  expect_equal(colnames(res$draws_mat), c("Item 2", "Item 1", "Item 3"))
  expect_equal(unname(colMeans(res$draws_mat)), c(-0.5, 0.2, 0.3))
  expect_equal(unique(res$impact), 0.7)
  expect_error(.extract_dif_slopes(NULL, c("Item 1", "Item 9"), "item",
                                   fixed_coef = "difgroup", exact = TRUE),
               "not found for item")
})

test_that("dichotomous slope DIF finds the group coefficient by prefix", {
  n_draws <- 10
  re <- array(0, dim = c(n_draws, 2, 2),
              dimnames = list(NULL, c("I1", "I2"), c("Intercept", "genderM")))
  re[, , "genderM"] <- rep(c(0.4, -0.4), each = n_draws)
  fe <- cbind(Intercept = 0, genderM = rep(-0.3, n_draws))
  local_mocked_bindings(
    ranef = function(object, summary = TRUE, ...) list(item = re),
    fixef = function(object, summary = TRUE, ...) fe,
    .package = "brms"
  )
  res <- .extract_dif_slopes(NULL, c("I1", "I2"), "item",
                             fixed_coef = "gender", exact = FALSE)
  expect_equal(unname(colMeans(res$draws_mat)), c(0.4, -0.4))
  expect_equal(unique(res$impact), -0.3)
})
