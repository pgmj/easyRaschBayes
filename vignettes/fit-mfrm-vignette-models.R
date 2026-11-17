# Simulated data and model fits for the many-facet Rasch model vignette
# (easyRaschBayes 0.4.0.9001, MFRM step 5, see dev/mfrm-design-plan.md).
#
# Writes to vignettes/fits/:
#   mfrm_data.rds          simulated ratings (polytomous and dichotomous)
#   mfrm_truth.rds         generating values, to check each diagnostic
#   fit_mfrm_rsm.rds       rating scale MFRM, no halo term
#   fit_mfrm_rsm_halo.rds  rating scale MFRM with (1 | rater:id)
#   fit_mfrm_rsm_scales.rds  rater-specific rating scales + halo term
#   fit_mfrm_rm.rds        dichotomous MFRM, no halo term
#   fit_mfrm_rm_halo.rds   dichotomous MFRM with (1 | rater:id)
#   prior_contr_sum.rds    prior-only fit, items coded with contr.sum
#   prior_contr_ozs.rds    prior-only fit, items coded with contr_sum_zero
#
# Run on the host from the package root:
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 nohup Rscript vignettes/fit-mfrm-vignette-models.R > vignettes/fits/fit_mfrm.log 2>&1 &
# FIT=data | poly | dich | prior fits only that part (data is always
# simulated first, with a fixed seed, so every part sees the same data).
# Uses 4 cores per fit, one fit at a time.

stopifnot(file.exists("DESCRIPTION"),
          read.dcf("DESCRIPTION")[, "Package"] == "easyRaschBayes")
suppressMessages({
  library(brms)
  devtools::load_all(".", quiet = TRUE)
})

which_fits <- Sys.getenv("FIT", "all")
dir.create("vignettes/fits", showWarnings = FALSE)
CORES <- 4L

# ---------------------------------------------------------------------------
# Simulated rating design
# ---------------------------------------------------------------------------
# 150 persons rated on 4 criteria by 3 of 6 raters each. Rater sets follow
# the cyclic offsets {0, 1, 3} (mod 6), so every pair of raters shares
# persons and the design is connected, but each person has only 12 of the
# 24 rater x criterion cells (sparse). 1,800 ratings.
#
# Built-in departures from the rating scale MFRM, each with a diagnostic:
#   rater severity differs               -> facet_parameters() locations
#   R6 is inconsistent (discrimination 0.4) -> underfit for R6 (infit by
#                                            rater, restscore for R6 cells)
#   R3 uses the middle categories more    -> wider thresholds for R3 in the
#     (central tendency, wider thresholds)   rater-specific scales model
#   halo SD 0.6 for all raters            -> positive same-rater Q3 without
#                                            the halo term, halo SD with it
# Ratings have 4 categories, coded 1-4 for acat. The dichotomous version
# uses the same persons, design, severity, halo and inconsistency.

simulate_mfrm <- function(seed = 2026) {
  set.seed(seed)
  n <- 150L
  items <- paste0("crit", 1:4)
  raters <- paste0("R", 1:6)
  truth <- list(
    theta_mean = 0, theta_sd = 1,
    item_loc = setNames(c(-0.5, -0.1, 0.2, 0.4), items),
    severity = setNames(c(-0.6, -0.2, 0, 0.1, 0.3, 0.4), raters),
    disc     = setNames(c(1, 1, 1, 1, 1, 0.4), raters),
    thresholds = c(-1.5, 0, 1.5),
    thresholds_R3 = c(-2.6, 0, 2.6),
    halo_sd = 0.6,
    dich_intercept = 0.3
  )
  theta <- rnorm(n, truth$theta_mean, truth$theta_sd)
  d <- do.call(rbind, lapply(seq_len(n), function(p) {
    r <- raters[((p - 1L) + c(0L, 1L, 3L)) %% 6L + 1L]
    expand.grid(id = p, item = items, rater = r, stringsAsFactors = FALSE)
  }))
  halo <- rnorm(n * 6L, 0, truth$halo_sd)
  halo_idx <- (d$id - 1L) * 6L + match(d$rater, raters)
  loc <- theta[d$id] + halo[halo_idx] - truth$item_loc[d$item] -
    truth$severity[d$rater]
  a <- truth$disc[d$rater]

  # Polytomous: adjacent-category logits a * (loc - tau_k)
  rating <- integer(nrow(d))
  for (row in seq_len(nrow(d))) {
    tau <- if (d$rater[row] == "R3") truth$thresholds_R3 else truth$thresholds
    lp <- c(0, cumsum(a[row] * (loc[row] - tau)))
    p <- exp(lp - max(lp))
    rating[row] <- sample.int(4L, 1L, prob = p)
  }
  d$rating <- rating
  # Dichotomous: one threshold at -dich_intercept
  d$pass <- rbinom(nrow(d), 1L,
                   plogis(a * (loc + truth$dich_intercept)))
  d$id <- factor(d$id)
  d$item <- factor(d$item, levels = items)
  d$rater <- factor(d$rater, levels = raters)
  truth$theta <- theta
  list(data = d, truth = truth)
}

sim <- simulate_mfrm()
d <- sim$data
saveRDS(d, "vignettes/fits/mfrm_data.rds")
saveRDS(sim$truth, "vignettes/fits/mfrm_truth.rds")
message("Simulated ", nrow(d), " ratings of ", nlevels(d$id), " persons.")

# Four criteria as fixed effects with orthonormal sum-to-zero contrasts. A
# normal(0, s) prior on the coefficients gives each criterion effect prior
# SD s * sqrt(3/4), so s = sqrt(4/3) gives SD 1 per criterion.
contrasts(d$item) <- contr_sum_zero(levels(d$item))
s_item <- sqrt(4 / 3)

pri_poly <- prior(normal(0, 3), class = "Intercept") +
  prior_string(paste0("normal(0, ", round(s_item, 4), ")"), class = "b") +
  prior(exponential(1), class = "sd")
pri_dich <- prior(normal(0, 3), class = "Intercept") +
  prior_string(paste0("normal(0, ", round(s_item, 4), ")"), class = "b") +
  prior(exponential(1), class = "sd")

fit_and_save <- function(name, ...) {
  t0 <- Sys.time()
  fit <- brm(..., data = d, chains = 4, cores = CORES, iter = 2000,
             seed = 2026, refresh = 0)
  saveRDS(fit, file.path("vignettes/fits", paste0(name, ".rds")))
  message(name, " done in ", format(Sys.time() - t0),
          ", max Rhat ", round(max(rhat(fit), na.rm = TRUE), 3),
          ", divergent transitions ",
          sum(subset(nuts_params(fit), Parameter == "divergent__")$Value))
  invisible(fit)
}

# ---------------------------------------------------------------------------
# Polytomous fits
# ---------------------------------------------------------------------------
if (which_fits %in% c("all", "poly")) {
  fit_and_save("fit_mfrm_rsm",
               rating ~ 1 + item + (1 | rater) + (1 | id),
               family = acat(), prior = pri_poly)
  fit_and_save("fit_mfrm_rsm_halo",
               rating ~ 1 + item + (1 | rater) + (1 | id) + (1 | rater:id),
               family = acat(), prior = pri_poly)
  # Rater-specific scales: the rater location is in the thresholds, so no
  # (1 | rater) term. Six threshold sets of 3.
  fit_and_save("fit_mfrm_rsm_scales",
               rating | thres(gr = rater) ~ 1 + item + (1 | id) +
                 (1 | rater:id),
               family = acat(), prior = pri_poly)
}

# ---------------------------------------------------------------------------
# Dichotomous fits
# ---------------------------------------------------------------------------
if (which_fits %in% c("all", "dich")) {
  fit_and_save("fit_mfrm_rm",
               pass ~ 1 + item + (1 | rater) + (1 | id),
               family = bernoulli(), prior = pri_dich)
  fit_and_save("fit_mfrm_rm_halo",
               pass ~ 1 + item + (1 | rater) + (1 | id) + (1 | rater:id),
               family = bernoulli(), prior = pri_dich)
}

# ---------------------------------------------------------------------------
# Prior-only fits: contr.sum against contr_sum_zero for the four criteria
# ---------------------------------------------------------------------------
# Same coefficient prior normal(0, s_item) in both. The implied prior SD of
# each criterion effect is s_item * sqrt(3/4) = 1 for all four with
# contr_sum_zero, but s_item = 1.15 for the first three and
# s_item * sqrt(3) = 2 for the last with contr.sum. (Corrected 2026-11-16,
# the comment said sqrt(3/4) for the first three; checked on the fits.)
if (which_fits %in% c("all", "prior")) {
  d_sum <- d
  contrasts(d_sum$item) <- contr.sum(4)
  for (v in c("sum", "ozs")) {
    dd <- if (v == "sum") d_sum else d
    t0 <- Sys.time()
    fit <- brm(pass ~ 1 + item + (1 | rater) + (1 | id), data = dd,
               family = bernoulli(), prior = pri_dich,
               sample_prior = "only", chains = 4, cores = CORES,
               iter = 2000, seed = 2026, refresh = 0)
    saveRDS(fit, file.path("vignettes/fits",
                           paste0("prior_contr_", v, ".rds")))
    message("prior_contr_", v, " done in ", format(Sys.time() - t0))
  }
}

message("All requested fits done. Knit the vignette with vignettes/precompile.R.")
