# Tests for the total-score conditioning machinery in R/conditional_score.R
# and the summaries used by item_restscore_statistic().

# Category probabilities for a partial credit model. tau is a list of
# threshold vectors, one per item. Returns an n_obs x C matrix in long
# format (persons within items) plus the index vectors.
pcm_long <- function(theta, tau, n_cat = max(lengths(tau)) + 1L) {
  n <- length(theta)
  k <- length(tau)
  pid <- rep(seq_len(n), times = k)
  iid <- rep(seq_len(k), each = n)
  prob <- matrix(0, nrow = n * k, ncol = n_cat)
  for (row in seq_along(pid)) {
    ct <- c(0, cumsum(tau[[iid[row]]]))
    num <- exp(seq_along(ct) * 0 + (seq_along(ct) - 1) * theta[pid[row]] - ct)
    prob[row, seq_along(ct)] <- num / sum(num)
  }
  list(prob = prob, pid = pid, iid = iid)
}

# Brute-force P(X_i = x | total r) by enumerating all response patterns.
brute_cond <- function(tau, i, r) {
  grid <- expand.grid(lapply(tau, function(t) 0:length(t)))
  logw <- rowSums(sapply(seq_along(tau), function(j) {
    c(0, -cumsum(tau[[j]]))[grid[[j]] + 1]
  }))
  sel <- rowSums(grid) == r
  w <- exp(logw[sel])
  vapply(0:length(tau[[i]]),
         function(x) sum(w[grid[sel, i] == x]) / sum(w), numeric(1))
}

tau_poly <- list(c(-1, 0.5), c(0.2), c(-0.5, 0.3, 1.2), c(0.8, -0.2))

test_that("conditional item probabilities match brute-force enumeration", {
  set.seed(10)
  pl <- pcm_long(rnorm(5), tau_poly)
  cw <- .cond_item_weights(pl$prob, pl$pid, pl$iid, n_persons = 5, k = 4)
  expect_equal(max(cw$key), 1L)
  W <- matrix(cw$W[1, , ], nrow = 4)
  esf <- .cond_esf(W)
  M <- sum(lengths(tau_poly))
  probs <- .cond_expected_probs(W, esf, 0:M)
  for (i in seq_along(tau_poly)) {
    for (r in 0:M) {
      expected <- brute_cond(tau_poly, i, r)
      expect_equal(probs[[i]][r + 1, seq_along(expected)], expected,
                   tolerance = 1e-10)
      # Categories the item does not have get probability 0
      expect_true(all(probs[[i]][r + 1, -seq_along(expected)] == 0))
    }
  }
})

test_that("conditional quantities do not depend on the person locations", {
  X <- matrix(c(0, 1, 2, 1,
                1, 0, 3, 2,
                2, 1, 1, 0,
                0, 0, 0, 1,
                2, 1, 3, 2,
                1, 1, 2, 0), ncol = 4, byrow = TRUE)
  run <- function(theta) {
    pl <- pcm_long(theta, tau_poly)
    .cond_restscore_draw(pl$prob, pl$pid, pl$iid, X)$gamma_exp
  }
  expect_equal(run(rep(0, 6)), run(c(-3, 2, 0.5, 4, -1, 1)),
               tolerance = 1e-10)
})

test_that("expected gamma matches a brute-force expected table", {
  X <- matrix(c(0, 1, 2, 1,
                1, 0, 3, 2,
                2, 1, 1, 0,
                0, 0, 0, 1,
                2, 1, 3, 2,
                1, 1, 2, 0,
                1, 0, 1, 1), ncol = 4, byrow = TRUE)
  pl <- pcm_long(seq(-1, 1, length.out = 7), tau_poly)
  got <- .cond_restscore_draw(pl$prob, pl$pid, pl$iid, X)$gamma_exp
  r <- rowSums(X)
  M <- sum(lengths(tau_poly))
  for (i in seq_along(tau_poly)) {
    tab <- matrix(0, nrow = M + 1, ncol = 4)
    for (v in seq_len(nrow(X))) {
      p <- brute_cond(tau_poly, i, r[v])
      for (x in seq_along(p) - 1) {
        if (r[v] - x >= 0) {
          tab[r[v] - x + 1, x + 1] <- tab[r[v] - x + 1, x + 1] + p[x + 1]
        }
      }
    }
    expect_equal(got[i], gk_gamma(tab), tolerance = 1e-10)
  }
})

test_that("expected gamma matches iarm for a dichotomous Rasch model", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  X <- as.matrix(eRm::raschdat1[, 1:10])
  fit <- eRm::RM(X)
  beta <- -fit$betapar
  n <- nrow(X)
  k <- ncol(X)
  pid <- rep(seq_len(n), k)
  iid <- rep(seq_len(k), each = n)
  set.seed(3)
  p1 <- stats::plogis(rnorm(n)[pid] - beta[iid])
  got <- .cond_restscore_draw(cbind(1 - p1, p1), pid, iid, X)$gamma_exp
  ia <- iarm::item_restscore(fit)
  expect_equal(got, as.numeric(ia[, 2]), tolerance = 1e-3)
})

test_that("conditional sampling keeps totals and matches exact probabilities", {
  tau <- list(c(-0.5), c(0.3, 0.6), c(1))
  n <- 6000
  set.seed(4)
  pl <- pcm_long(rnorm(n), tau)
  cw <- .cond_item_weights(pl$prob, pl$pid, pl$iid, n_persons = n, k = 3)
  W1 <- matrix(cw$W[1, , ], nrow = 3)
  esf <- .cond_esf(W1)
  suffix <- array(esf$suffix, dim = c(1, dim(esf$suffix)))
  r <- rep(2L, n)
  x <- .cond_sample(cw$W, suffix, cw$key, r)
  expect_true(all(rowSums(x) == 2L))
  for (i in 1:3) {
    exact <- brute_cond(tau, i, 2)
    observed <- tabulate(x[, i] + 1, nbins = length(exact)) / n
    # Binomial SE at n = 6000 is at most 0.0065
    expect_lt(max(abs(observed - exact)), 0.03)
  }
})

test_that("missing responses stay at zero and keep their own group", {
  tau <- list(c(0), c(0.5), c(-0.5))
  pl <- pcm_long(c(0, 1, -1), tau)
  # Person 3 did not answer item 2
  keep <- !(pl$pid == 3 & pl$iid == 2)
  cw <- .cond_item_weights(pl$prob[keep, ], pl$pid[keep], pl$iid[keep],
                           n_persons = 3, k = 3)
  expect_equal(cw$key[1], cw$key[2])
  expect_false(cw$key[3] == cw$key[1])
  U <- dim(cw$W)[1]
  suffix <- array(0, dim = c(U, 4, 4))
  for (u in seq_len(U)) {
    suffix[u, , ] <- .cond_esf(matrix(cw$W[u, , ], nrow = 3))$suffix
  }
  set.seed(5)
  x <- .cond_sample(cw$W, suffix, cw$key, r = c(2L, 1L, 2L))
  expect_equal(x[3, 2], 0L)
  expect_equal(rowSums(x), c(2, 1, 2))
})

test_that("persons with different item parameters form separate groups", {
  tau_a <- list(c(-1), c(0), c(1))
  tau_b <- list(c(-1), c(0.7), c(1))   # DIF in item 2
  pa <- pcm_long(c(0, 1), tau_a)
  pb <- pcm_long(c(-0.5, 2), tau_b)
  prob <- rbind(pa$prob, pb$prob)
  pid <- c(pa$pid, pb$pid + 2L)
  iid <- c(pa$iid, pb$iid)
  cw <- .cond_item_weights(prob, pid, iid, n_persons = 4, k = 3)
  expect_equal(cw$key, c(1L, 1L, 2L, 2L))
})

test_that(".ppp_mid counts ties as one half", {
  rep_draws <- matrix(c(0.1, 0.2, 0.2, 0.3), ncol = 1)
  expect_equal(.ppp_mid(0.2, rep_draws), 0.25 + 0.5 * 0.5)
  expect_equal(.ppp_mid(NA_real_, rep_draws), NA_real_)
})

test_that(".ppp_maxt is two-sided, bounded and at least the per-item value", {
  set.seed(6)
  rep_draws <- matrix(rnorm(4000 * 5), ncol = 5)
  obs <- c(0, 1, -2.5, 3, 0.5)
  adj <- .ppp_maxt(obs, rep_draws)
  expect_true(all(adj >= 0 & adj <= 1))
  m <- colMeans(rep_draws)
  s <- apply(rep_draws, 2, sd)
  z <- abs(obs - m) / s
  unadj <- vapply(1:5, function(i) {
    mean(abs(rep_draws[, i] - m[i]) / s[i] >= z[i])
  }, numeric(1))
  expect_true(all(adj >= unadj))
  # Same deviation in either direction gives the same value
  expect_equal(.ppp_maxt(c(2, 0), rep_draws[, 1:2]),
               .ppp_maxt(c(-2 + 2 * mean(rep_draws[, 1]), 0),
                         rep_draws[, 1:2]))
  # Items with no replicated variation get NA
  expect_true(is.na(.ppp_maxt(c(0, 1), cbind(rep(0, 10), rnorm(10)))[1]))
})

test_that("the sufficiency check accepts Rasch models and rejects others", {
  ok <- function(...) {
    args <- list(family = "acat", link = "logit", dpar_names = character(0),
                 disc_values = rep(1, 3), sd_names = "sd_id__Intercept",
                 person_name = "id", fun_name = "f")
    args[names(list(...))] <- list(...)
    do.call(.check_sufficient_score_model, args)
  }
  expect_true(ok())
  expect_true(ok(family = "bernoulli",
                 sd_names = c("sd_id__Intercept", "sd_item__Intercept")))
  expect_true(ok(sd_names = c("sd_id__Intercept:groupA",
                              "sd_id__Intercept:groupB")))
  expect_error(ok(family = "cumulative"), "Rasch-family")
  expect_error(ok(link = "probit"), "Rasch-family")
  expect_error(ok(dpar_names = "disc"), "discrimination")
  expect_error(ok(disc_values = c(1, 1.3)), "discrimination")
  expect_error(ok(sd_names = c("sd_id__dimA", "sd_id__dimB")),
               "single person location")
  expect_error(ok(sd_names = "sd_item__Intercept"), "no person random")
})

test_that("item_restscore_post rejects output from the old version", {
  old <- list(
    result = data.frame(item = "a", gamma_obs = 0.5, gamma_rep = 0.4,
                        gamma_diff = 0.1, ppp = 0.9),
    draws = data.frame(item = "a", gamma_rep = 0.4, gamma = 0.5)
  )
  expect_error(item_restscore_post(old), "0.3.0.9000 or later")
})

test_that("item_restscore_post uses the item column name it is given", {
  set.seed(7)
  out <- list(
    result = tibble::tibble(question = c("q1", "q2"),
                            gamma_obs = c(0.5, 0.4), gamma_exp = c(0.45, 0.42),
                            gamma_diff = c(0.05, -0.02), ppp = c(0.8, 0.4),
                            p_adj = c(0.6, 1)),
    draws = tibble::tibble(question = rep(c("q1", "q2"), each = 50),
                           draw = rep(1:50, 2),
                           gamma_rep = rnorm(100, 0.44, 0.03),
                           gamma_exp = 0.44)
  )
  res <- item_restscore_post(out)
  expect_named(res$summary, c("question", "gamma_obs", "gamma_exp",
                              "gamma_diff", "ppp", "p_adj"))
  expect_s3_class(res$plot, "ggplot")
})

# --- Conditional infit ------------------------------------------------------

test_that("conditional infit and outfit match brute-force moments", {
  X <- matrix(c(0, 1, 2, 1,
                1, 0, 3, 2,
                2, 1, 1, 0,
                0, 0, 0, 1,
                2, 1, 3, 2,
                1, 1, 2, 0,
                1, 0, 1, 1), ncol = 4, byrow = TRUE)
  pl <- pcm_long(seq(-1, 1, length.out = 7), tau_poly)
  got <- .cond_infit_draw(pl$prob, pl$pid, pl$iid, X)
  r <- rowSums(X)
  E <- V <- matrix(NA_real_, nrow(X), 4)
  for (v in seq_len(nrow(X))) {
    for (i in 1:4) {
      p <- brute_cond(tau_poly, i, r[v])
      x <- seq_along(p) - 1
      E[v, i] <- sum(x * p)
      V[v, i] <- sum(x^2 * p) - E[v, i]^2
    }
  }
  use <- V > 1e-12
  infit <- colSums(ifelse(use, (X - E)^2, 0)) / colSums(ifelse(use, V, 0))
  outfit <- colSums(ifelse(use, (X - E)^2 / V, 0)) / colSums(use)
  expect_equal(got$infit, infit, tolerance = 1e-10)
  expect_equal(got$outfit, outfit, tolerance = 1e-10)
})

test_that("conditional infit matches iarm::out_infit for RM and PCM", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(8)
  # Dichotomous
  X <- as.matrix(eRm::raschdat1[, 1:10])
  fit <- eRm::RM(X)
  beta <- -fit$betapar
  n <- nrow(X)
  k <- ncol(X)
  pid <- rep(seq_len(n), k)
  iid <- rep(seq_len(k), each = n)
  p1 <- stats::plogis(rnorm(n, sd = 2)[pid] - beta[iid])
  got <- .cond_infit_draw(cbind(1 - p1, p1), pid, iid, X)
  ia <- iarm::out_infit(fit)
  expect_equal(unname(got$infit), unname(ia$Infit), tolerance = 1e-8)
  expect_equal(unname(got$outfit), unname(ia$Outfit), tolerance = 1e-8)

  # Partial credit
  Xp <- as.matrix(eRm::pcmdat2)
  pfit <- eRm::PCM(Xp)
  thr <- eRm::thresholds(pfit)[[3]][[1]][, -1, drop = FALSE]
  tau <- lapply(seq_len(nrow(thr)), function(i) stats::na.omit(thr[i, ]))
  n <- nrow(Xp)
  k <- ncol(Xp)
  pl <- pcm_long(rnorm(n), tau)
  got <- .cond_infit_draw(pl$prob, pl$pid, pl$iid, Xp)
  ia <- iarm::out_infit(pfit)
  expect_equal(unname(got$infit), unname(ia$Infit), tolerance = 1e-8)
  expect_equal(unname(got$outfit), unname(ia$Outfit), tolerance = 1e-8)
})

test_that("conditional infit ignores person locations and extreme scores", {
  X <- matrix(c(0, 1, 2, 1,
                1, 0, 3, 2,
                2, 1, 1, 0,
                0, 0, 0, 1,
                1, 1, 2, 0), ncol = 4, byrow = TRUE)
  run <- function(X, theta) {
    pl <- pcm_long(theta, tau_poly)
    set.seed(9)
    .cond_infit_draw(pl$prob, pl$pid, pl$iid, X)$infit
  }
  base <- run(X, rep(0, 5))
  expect_equal(base, run(X, c(2, -1, 0.3, -2, 1)), tolerance = 1e-10)
  # Adding persons with the minimum and maximum total changes nothing
  X_ext <- rbind(X, c(0, 0, 0, 0), c(2, 1, 3, 2))
  expect_equal(run(X_ext, rep(0, 7)), base, tolerance = 1e-10)
})

test_that("replicated infit keeps missing responses out", {
  X <- matrix(c(0, 1, 2, 1,
                1, 0, 3, 2,
                2, 1, 1, 0,
                0, 0, 0, 1,
                2, 1, 3, 2,
                1, 1, 2, 0), ncol = 4, byrow = TRUE)
  pl <- pcm_long(rep(0, 6), tau_poly)
  drop <- pl$pid == 2 & pl$iid == 3
  X_na <- X
  X_na[2, 3] <- NA
  set.seed(11)
  got <- .cond_infit_draw(pl$prob[!drop, ], pl$pid[!drop], pl$iid[!drop], X_na)
  expect_true(all(is.finite(got$infit)))
  expect_true(all(is.finite(got$infit_rep)))
})

# --- infit_post -------------------------------------------------------------

test_that("infit_post ppp is P(observed > replicated) with mid-p ties", {
  draws <- tibble::tibble(
    item = rep(c("a", "b"), each = 4),
    draw = rep(1:4, 2),
    infit = c(1.2, 1.2, 1.2, 1.2, 0.9, 0.9, 0.9, 0.9),
    infit_rep = c(1.0, 1.1, 1.2, 1.3, 1.0, 1.0, 0.8, 1.1)
  )
  s <- infit_post(draws)$summary
  expect_equal(s$infit_ppp[s$item == "a"], 0.5 + 0.5 * 0.25)
  expect_equal(s$infit_ppp[s$item == "b"], 0.25)
})

test_that("infit_post uses the item column name it is given", {
  set.seed(12)
  draws <- tibble::tibble(
    question = rep(c("q1", "q2"), each = 50),
    draw = rep(1:50, 2),
    infit = rnorm(100, 1, 0.02),
    infit_rep = rnorm(100, 1, 0.05)
  )
  res <- infit_post(dplyr::group_by(draws, question))
  expect_named(res$summary, c("question", "infit_obs", "infit_rep",
                              "infit_ppp"))
  expect_named(res$hdi, c("question", "underfit", "overfit"))
  expect_s3_class(res$plot, "ggplot")
})

test_that("q3_post ppp and diff use the per-draw replicated values", {
  q <- tibble::tibble(
    draw = rep(1:4, 2),
    item_pair = rep(c("x", "y"), each = 4), item_1 = "i", item_2 = "j",
    q3 = c(0.1, 0.1, 0.1, 0.1, -0.1, -0.1, -0.1, -0.1),
    q3_rep = c(0, 0.2, 0.05, 0.15, 0, -0.2, 0.1, -0.15)
  )
  s <- q3_post(q)$summary
  expect_named(s, c("item_pair", "item_1", "item_2", "q3_obs", "q3_rep",
                    "q3_diff", "q3_ppp", "q3_p_adj"))
  expect_equal(s$q3_ppp[s$item_pair == "x"], 0.5)
  expect_equal(s$q3_ppp[s$item_pair == "y"], 0.5)
  expect_equal(s$q3_diff[s$item_pair == "y"], round(-0.1 - mean(c(0, -0.2, 0.1, -0.15)), 3))
})

# --- conditional Q3 ---------------------------------------------------------

test_that("conditional Q3 ignores person locations and extreme scores", {
  X <- matrix(c(0, 1, 2, 1,
                1, 0, 3, 2,
                2, 1, 1, 0,
                0, 0, 0, 1,
                1, 1, 2, 0,
                2, 0, 2, 2), ncol = 4, byrow = TRUE)
  run <- function(X, theta) {
    pl <- pcm_long(theta, tau_poly)
    set.seed(9)
    .cond_q3_draw(pl$prob, pl$pid, pl$iid, X)
  }
  base <- run(X, rep(0, 6))
  expect_length(base$obs, 6)
  expect_equal(base$obs, run(X, c(2, -1, 0.3, -2, 1, 0.5))$obs,
               tolerance = 1e-10)
  X_ext <- rbind(X, c(0, 0, 0, 0), c(2, 1, 3, 2))
  expect_equal(run(X_ext, rep(0, 8))$obs, base$obs, tolerance = 1e-10)
})

test_that("conditional Q3 uses residuals from E(X | total) in pair order", {
  set.seed(3)
  n <- 40
  pl <- pcm_long(rnorm(n), tau_poly)
  X <- matrix(NA_real_, n, 4)
  X[cbind(pl$pid, pl$iid)] <- apply(pl$prob, 1, function(p) {
    sample(seq_along(p) - 1, 1, prob = p)
  })
  got <- .cond_q3_draw(pl$prob, pl$pid, pl$iid, X)
  # brute-force conditional expectations
  r <- rowSums(X)
  E <- sapply(1:4, function(i) vapply(r, function(rr) {
    p <- brute_cond(tau_poly, i, rr)
    sum((seq_along(p) - 1) * p)
  }, numeric(1)))
  M <- sum(lengths(tau_poly))
  keep <- r > 0 & r < M
  R <- (X - E)[keep, ]
  C <- cor(R)
  expect_equal(got$obs, c(C[1, 2], C[1, 3], C[1, 4], C[2, 3], C[2, 4],
                          C[3, 4]), tolerance = 1e-8)
})

test_that("conditional Q3 centres on -1/(k-1) under the model", {
  set.seed(21)
  n <- 3000
  b <- seq(-1.5, 1.5, length.out = 6)
  theta <- rnorm(n)
  pid <- rep(seq_len(n), times = 6)
  iid <- rep(seq_len(6), each = n)
  p1 <- plogis(theta[pid] - b[iid])
  prob <- cbind(1 - p1, p1)
  X <- matrix(rbinom(n * 6, 1, p1), n, 6)
  got <- .cond_q3_draw(prob, pid, iid, X)
  expect_equal(mean(got$obs), -1 / 5, tolerance = 0.01)
  expect_equal(mean(got$rep), -1 / 5, tolerance = 0.01)
})

test_that(".ppp_maxt_upper is one-sided and at most the per-pair value", {
  set.seed(6)
  rep_draws <- matrix(rnorm(4000 * 5), ncol = 5)
  obs <- c(0, 1, -2.5, 3, 0.5)
  adj <- .ppp_maxt_upper(obs, rep_draws)
  expect_true(all(adj >= 0 & adj <= 1))
  m <- colMeans(rep_draws)
  s <- apply(rep_draws, 2, sd)
  unadj <- vapply(1:5, function(i) {
    mean((rep_draws[, i] - m[i]) / s[i] >= (obs[i] - m[i]) / s[i])
  }, numeric(1))
  expect_true(all(adj >= unadj))
  # a large negative deviation is not flagged
  expect_gt(adj[3], 0.99)
  expect_lt(adj[4], 0.05)
})

test_that("q3_post adds q3_p_adj and requires equal draws per pair", {
  set.seed(4)
  q <- tibble::tibble(
    draw = rep(1:200, 3),
    item_pair = rep(c("a : b", "a : c", "b : c"), each = 200),
    item_1 = "a", item_2 = "b",
    q3 = c(rnorm(200, 0.3, 0.01), rnorm(200, 0, 0.01), rnorm(200, 0, 0.01)),
    q3_rep = rnorm(600, 0, 0.05)
  )
  s <- q3_post(q)$summary
  expect_lt(s$q3_p_adj[s$item_pair == "a : b"], 0.05)
  expect_gt(s$q3_p_adj[s$item_pair == "a : c"], 0.05)
  expect_error(q3_post(q[-1, ]), "same number of draws")
})
