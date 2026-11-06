# Tests for the shared item threshold extractor (R/item_draws.R), the
# eigenvector sign alignment and the per-pattern WLE.

n <- 4  # draws

test_that("grouped thresholds are matched by position, renamed labels allowed", {
  d <- data.frame(
    "b_Intercept[I1,1]" = rep(-1, n), "b_Intercept[I1,2]" = rep(1, n),
    "b_Intercept[I2,1]" = rep(0.5, n),
    "r_id[1,Intercept]" = 0, check.names = FALSE
  )
  # brms writes "I 1" as "I1" in parameter names
  tau <- .item_tau_draws(d, c("I 1", "I2"), "item", ordinal = TRUE)
  expect_named(tau, c("I 1", "I2"))
  expect_equal(unname(tau[["I 1"]][1, ]), c(-1, 1))
  expect_equal(ncol(tau[["I2"]]), 1)
})

test_that("a positional mismatch of plain labels is an error", {
  d <- data.frame("b_Intercept[B,1]" = rep(0, n), "b_Intercept[A,1]" = rep(1, n),
                  check.names = FALSE)
  expect_error(.item_tau_draws(d, c("A", "B"), "item", ordinal = TRUE),
               "by position")
})

test_that("rating scale model: common thresholds minus item random effects", {
  d <- data.frame(
    "b_Intercept[1]" = rep(-1, n), "b_Intercept[2]" = rep(1, n),
    "r_item[I1,Intercept]" = rep(0.5, n), "r_item[I2,Intercept]" = rep(-0.5, n),
    check.names = FALSE
  )
  tau <- .item_tau_draws(d, c("I1", "I2"), "item", ordinal = TRUE)
  expect_equal(unname(tau$I1[1, ]), c(-1.5, 0.5))
  expect_equal(unname(tau$I2[1, ]), c(-0.5, 1.5))
})

test_that("common thresholds without item effects are an error", {
  d <- data.frame("b_Intercept[1]" = rep(0, n), check.names = FALSE)
  expect_error(.item_tau_draws(d, c("I1", "I2"), "item", ordinal = TRUE),
               "all items are identical")
})

test_that("dichotomous random items: difficulty = -(intercept + r_item)", {
  d <- data.frame(
    b_Intercept = rep(0.2, n),
    "r_item[I1,Intercept]" = rep(0.3, n), "r_item[I2,Intercept]" = rep(-0.4, n),
    "r_id[1,Intercept]" = 0, check.names = FALSE
  )
  tau <- .item_tau_draws(d, c("I1", "I2"), "item", ordinal = FALSE)
  expect_equal(c(tau$I1[1, 1], tau$I2[1, 1]), c(-0.5, 0.2))
})

test_that("dichotomous fixed items, with and without an intercept", {
  d0 <- data.frame(b_itemI1 = rep(0.3, n), b_itemI2 = rep(-0.4, n))
  tau0 <- .item_tau_draws(d0, c("I1", "I2"), "item", ordinal = FALSE)
  expect_equal(c(tau0$I1[1, 1], tau0$I2[1, 1]), c(-0.3, 0.4))
  d1 <- data.frame(b_Intercept = rep(0.3, n), b_itemI2 = rep(-0.7, n))
  tau1 <- .item_tau_draws(d1, c("I1", "I2"), "item", ordinal = FALSE)
  expect_equal(c(tau1$I1[1, 1], tau1$I2[1, 1]), c(-0.3, 0.4))
})

test_that("eigenvector alignment follows strong loadings, not sign counts", {
  set.seed(4)
  k <- 12
  truth <- c(rep(0.577, 3), rep(0, 9))
  V <- t(replicate(300, {
    x <- truth + rnorm(k, 0, 0.18)
    x <- x / sqrt(sum(x^2))
    if (runif(1) < 0.5) -x else x
  }))
  A <- .align_eigenvectors(V)
  agree <- sign((A %*% truth)[, 1])
  expect_true(all(agree == agree[1]))
  # Rows with NA stay untouched
  V[5, ] <- NA
  expect_true(all(is.na(.align_eigenvectors(V)[5, ])))
})

test_that("WLE by answered-item pattern uses only the answered items", {
  diffs <- c(-1, -0.5, 0, 0.5, 1, 1.5)
  table_fun <- function(items) .compute_wle_dichotomous(diffs[items], c(-7, 7))
  answered <- rbind(rep(TRUE, 6), c(rep(TRUE, 3), rep(FALSE, 3)))
  res <- .wle_by_pattern(answered, score = c(0, 0), table_fun = table_fun)
  full <- table_fun(1:6)
  part <- table_fun(1:3)
  expect_equal(res$wle[1], full$wle[1])
  expect_equal(res$wle[2], part$wle[1])
  expect_false(isTRUE(all.equal(res$wle[1], res$wle[2])))
})

test_that("gk_gamma_vec drops pairs with missing values", {
  expect_equal(gk_gamma_vec(c(0, 1, NA, 2, 1), c(1, 3, 2, NA, 2)),
               gk_gamma_vec(c(0, 1, 1), c(1, 3, 2)))
})

test_that(".compute_ess uses the chains", {
  set.seed(1)
  x <- c(rnorm(500), rnorm(500), rnorm(500), rnorm(500, 0.5))
  ch <- rep(1:4, each = 500)
  expect_equal(.compute_ess(x, ch), posterior::ess_bulk(matrix(x, ncol = 4)))
  expect_true(is.na(.compute_ess(c(1, 2, 3))))
})
