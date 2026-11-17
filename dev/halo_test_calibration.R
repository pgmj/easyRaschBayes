# A pooled halo test for many-facet Rasch models: null calibration,
# specificity and power (easyRaschBayes 0.4.0.9001, see "Halo test" in
# dev/mfrm-design-plan.md). DRAFT, NOT RUN.
#
# Background. The MFRM null run (2026-11-15) found that a halo in the data
# but not in the model raises the mean Q3 of same-rater cell pairs clearly
# (-0.07 RM, 0.00 RSM, against -0.11 replicated), but the pair-level
# max-T (q3_p_adj) detects it in only 4-64% of datasets, because the
# signal is spread over 80 pairs. This script tests statistics that pool
# the pairs first. All are computed from the conditional q3_statistic()
# draws of a model WITHOUT a halo term, so no new package code is needed
# yet.
#
# Statistics (one-sided, small p = more same-rater dependence than
# replicated; mid-p, p = P(rep > obs) + 0.5 * P(rep == obs) over draws):
#   mean_same  mean Q3 over all same-rater pairs, observed vs replicated
#              in each draw.
#   contrast   mean same-rater Q3 minus mean Q3 of all other pairs. Since
#              residuals given the total score sum to zero, a halo pushes
#              the other pairs down (seen in the null run), so the contrast
#              should carry about twice the signal.
#   per rater  mean Q3 over each rater's same-rater pairs, with p per rater
#              and a one-sided max-T p_adj over raters (.ppp_maxt_upper(),
#              as q3_p_adj in q3_post()). Identifies which raters have a
#              halo.
#   pair_maxt  the existing q3_p_adj from q3_post(): smallest value over
#              same-rater pairs, for comparison.
#
# Questions:
#   1. Null: are mean_same and contrast uniform, and does per-rater p_adj
#      keep the familywise rate at or below .05?
#   2. Specificity: do rater misfits other than a halo (an inconsistent
#      rater, and in the RSM a rater with central tendency) trigger the
#      halo test?
#   3. Power: halo SD 0.3 and 0.6 for all raters, and SD 0.6 for two raters
#      only. Pooled vs pair-level detection. Per-rater identification.
#   4. With the halo term in the model (halo SD 0.6 in the data), is the
#      test calibrated for the halo left after the term? The conditional
#      engine treats the halo effects as known, so it may flag less often
#      than nominal, as infit did.
#
# Design: 2 models x 6 conditions x 2 sample sizes, REPS replications.
# A condensed version (n = 300, 4 conditions, 1,000 fits) is
# dev/halo_test_quick.R.
#   Same rating design as dev/null_mfrm_calibration.R: 5 items, 8 raters,
#   each person rated by 2 adjacent raters (cyclic), 10 of 40 cells. 80
#   same-rater pairs, 40 same-item pairs, 160 neither pairs.
#   model  rm   dichotomous; rsm  rating scale, acat, 4 categories
#   cond   null          no halo, no rater misfit, model without halo term
#          misfit        no halo. R2 inconsistent (discrimination 0.4),
#                        and in the RSM R5 with central tendency
#                        (thresholds x 1.8). Model without halo term.
#          halo03        halo SD 0.3, all raters, not modelled
#          halo06        halo SD 0.6, all raters, not modelled
#          halo_two      halo SD 0.6 for R3 and R7 only, not modelled
#          halo06_model  halo SD 0.6, all raters, + (1 | rater:id) in model
#   n      300 and 600 persons, theta ~ N(0, 1)
#
# Run from the package root on the host:
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 nohup Rscript dev/halo_test_calibration.R > dev/results/halo_test.log 2>&1 &
# Smoke test (2 reps, n = 150, short chains, 2 workers):
#   SMOKE=1 Rscript dev/halo_test_calibration.R
# Summarise with dev/halo_test_calibration_summary.R.
#
# Rough cost (from the MFRM null run, not measured for this script): 4,800
# fits. Fits take 10-30 s (RM) and 80-185 s (RSM). Only Q3 is computed,
# which took a fifth or less of the 30-124 s of statistics per replication
# without a halo term. halo06_model is slow: every person is its own
# weight group, about 2-5 min of Q3 per replication on top of the fit.
# Expect very roughly 15-20 hours on 10 workers. Extrapolate from the
# smoke test's fit_sec and stat_sec.
#
# Each replication is saved on its own, so an interrupted run resumes.

stopifnot(file.exists("DESCRIPTION"),
          read.dcf("DESCRIPTION")[, "Package"] == "easyRaschBayes")
smoke <- Sys.getenv("SMOKE") == "1"

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(brms)
})

N_WORKERS <- if (smoke) 2L else 10L   # max 10 cores on the M4 Max
ITER      <- if (smoke) 400L else 1500L
WARMUP    <- if (smoke) 200L else 500L
CHAINS    <- 2L
NDRAWS    <- if (smoke) 200L else 1000L   # draws used by q3_statistic()
SIZES     <- if (smoke) c(150L) else c(300L, 600L)

out_dir <- file.path("dev", "results",
                     if (smoke) "halo_test_smoke" else "halo_test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Generating values, fixed across replications ---
items  <- sprintf("I%d", 1:5)
raters <- sprintf("R%d", 1:8)
b_item   <- setNames(seq(-0.8, 0.8, length.out = 5), items)
severity <- setNames(seq(-0.5, 0.5, length.out = 8), raters)
tau_rsm  <- c(-1.2, 0, 1.2)

# Rater-level settings per condition: halo SD, discrimination, threshold
# stretch (central tendency).
rater_setup <- function(cond, model_type) {
  halo_sd <- setNames(rep(0, 8), raters)
  disc    <- setNames(rep(1, 8), raters)
  stretch <- setNames(rep(1, 8), raters)
  if (cond == "misfit") {
    disc["R2"] <- 0.4
    if (model_type == "rsm") stretch["R5"] <- 1.8
  }
  if (cond == "halo03") halo_sd[] <- 0.3
  if (cond %in% c("halo06", "halo06_model")) halo_sd[] <- 0.6
  if (cond == "halo_two") halo_sd[c("R3", "R7")] <- 0.6
  list(halo_sd = halo_sd, disc = disc, stretch = stretch)
}

simulate_data <- function(model_type, cond, n) {
  rs <- rater_setup(cond, model_type)
  theta <- rnorm(n)
  d <- do.call(rbind, lapply(seq_len(n), function(p) {
    r <- raters[((p - 1L) + 0:1) %% 8L + 1L]
    expand.grid(id = p, item = items, rater = r, stringsAsFactors = FALSE)
  }))
  # One halo effect per person x rater, scaled by that rater's halo SD
  h <- matrix(rnorm(n * 8L), n, 8L) %*% diag(rs$halo_sd)
  ri <- match(d$rater, raters)
  loc <- theta[d$id] + h[cbind(d$id, ri)] - b_item[d$item] - severity[ri]
  a <- rs$disc[ri]
  if (model_type == "rsm") {
    tau <- outer(rs$stretch[ri], tau_rsm)          # rows x 3
    lp <- cbind(0, t(apply(a * (loc - tau), 1L, cumsum)))
    p <- exp(lp - apply(lp, 1L, max))
    p <- p / rowSums(p)
    u <- runif(nrow(d))
    d$response <- 1L + rowSums(u > t(apply(p, 1L, cumsum))[, -4L,
                                                           drop = FALSE])
  } else {
    d$response <- rbinom(nrow(d), 1L, plogis(a * loc))
  }
  d$id <- factor(d$id)
  d$item <- factor(d$item, levels = items)
  d$rater <- factor(d$rater, levels = raters)
  d
}

halo_in_model <- function(cond) cond == "halo06_model"
model_key <- function(model_type, cond) {
  paste0(model_type, if (halo_in_model(cond)) "_halo" else "")
}
model_formula <- function(model_type, cond) {
  f <- response ~ 1 + item + (1 | rater) + (1 | id)
  if (halo_in_model(cond)) f <- update(f, . ~ . + (1 | rater:id))
  bf(f)
}
model_family <- function(model_type) {
  if (model_type == "rsm") acat() else bernoulli()
}
model_prior <- function() {
  prior(normal(0, 3), class = "Intercept") +
    prior(normal(0, 3), class = "b") +
    prior(exponential(1), class = "sd")
}

conditions <- expand.grid(
  model_type = c("rm", "rsm"),
  cond = c("null", "misfit", "halo03", "halo06", "halo_two", "halo06_model"),
  n = SIZES, stringsAsFactors = FALSE
)
conditions$reps <- if (smoke) 2L else 200L
conditions$id <- paste0(conditions$model_type, "_", conditions$cond, "_n",
                        conditions$n)

# --- Compile each model once, then reuse with update(recompile = FALSE) ---
# chains = 0 compiles without sampling. rstan then prints "the number of
# chains is less than 1; sampling not done", which is expected, so the
# output is captured and replaced by one line per model.
base_fits <- list()
for (i in seq_len(nrow(conditions))) {
  key <- model_key(conditions$model_type[i], conditions$cond[i])
  if (!is.null(base_fits[[key]])) next
  set.seed(1)
  t0 <- Sys.time()
  invisible(utils::capture.output(suppressMessages(
    base_fits[[key]] <- brm(
      model_formula(conditions$model_type[i], conditions$cond[i]),
      data = simulate_data(conditions$model_type[i], conditions$cond[i], 40),
      family = model_family(conditions$model_type[i]),
      prior = model_prior(), chains = 0, refresh = 0
    )
  )))
  cat(format(Sys.time()), "compiled", key, "in",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}

# --- The halo statistics, from q3_statistic(facet_var = rater) draws ---
# Mid-p one-sided p-value: small when obs exceeds rep.
midp_upper <- function(obs, rep) {
  ok <- !is.na(obs) & !is.na(rep)
  mean(rep[ok] > obs[ok]) + 0.5 * mean(rep[ok] == obs[ok])
}

halo_stats <- function(q3) {
  q3 <- as.data.frame(q3)
  draws <- sort(unique(q3$draw))
  same <- q3$pair_type == "same rater"
  # Both cells of a same-rater pair have the same rater: "item:rater"
  q3$rater_1 <- sub("^.*:", "", as.character(q3$item_1))
  per_draw <- function(sel, col) {
    v <- tapply(q3[[col]][sel], q3$draw[sel], mean, na.rm = TRUE)
    as.numeric(v[as.character(draws)])
  }
  obs_same <- per_draw(same, "q3")
  rep_same <- per_draw(same, "q3_rep")
  obs_other <- per_draw(!same, "q3")
  rep_other <- per_draw(!same, "q3_rep")
  # obs_mean / rep_mean, not obs / rep: the saved tables also carry the
  # replication number as `rep`
  pooled <- data.frame(
    statistic = c("mean_same", "contrast"),
    obs_mean = c(mean(obs_same), mean(obs_same - obs_other)),
    rep_mean = c(mean(rep_same), mean(rep_same - rep_other)),
    p = c(midp_upper(obs_same, rep_same),
          midp_upper(obs_same - obs_other, rep_same - rep_other))
  )
  rr <- sort(unique(q3$rater_1[same]))
  obs_r <- sapply(rr, function(r) per_draw(same & q3$rater_1 == r, "q3"))
  rep_r <- sapply(rr, function(r) per_draw(same & q3$rater_1 == r,
                                           "q3_rep"))
  obs_r <- matrix(obs_r, ncol = length(rr))
  rep_r <- matrix(rep_r, ncol = length(rr))
  per_rater <- data.frame(
    rater = rr,
    obs_mean = colMeans(obs_r, na.rm = TRUE),
    rep_mean = colMeans(rep_r, na.rm = TRUE),
    p = vapply(seq_along(rr), function(j) midp_upper(obs_r[, j], rep_r[, j]),
               numeric(1)),
    p_adj = .ppp_maxt_upper(colMeans(obs_r, na.rm = TRUE), rep_r)
  )
  list(pooled = pooled, per_rater = per_rater)
}

run_one <- function(cond_row, rep) {
  f <- file.path(out_dir, cond_row$id, sprintf("rep_%03d.rds", rep))
  if (file.exists(f)) return(invisible(NULL))
  seed <- 1e5 * match(cond_row$id, conditions$id) + rep
  set.seed(seed)
  dat <- simulate_data(cond_row$model_type, cond_row$cond, cond_row$n)
  t0 <- Sys.time()
  fit <- suppressMessages(update(
    base_fits[[model_key(cond_row$model_type, cond_row$cond)]],
    newdata = dat, recompile = FALSE, chains = CHAINS, cores = 1,
    iter = ITER, warmup = WARMUP, refresh = 0, seed = seed
  ))
  fit_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rhat_max <- max(brms::rhat(fit), na.rm = TRUE)
  n_div <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)

  t1 <- Sys.time()
  info <- data.frame(id = cond_row$id, model_type = cond_row$model_type,
                     cond = cond_row$cond, n = cond_row$n, rep = rep)
  set.seed(seed + 1)
  q3 <- suppressMessages(q3_statistic(fit, facet_var = rater,
                                      ndraws_use = NDRAWS))
  hs <- halo_stats(q3)
  # Pair-level max-T as in the package: smallest q3_p_adj over same-rater
  # pairs and over all pairs
  q3s <- q3_post(q3)$summary
  pair <- data.frame(
    info,
    min_p_adj_same = min(q3s$q3_p_adj[q3s$pair_type == "same rater"],
                         na.rm = TRUE),
    min_p_adj_all = min(q3s$q3_p_adj, na.rm = TRUE)
  )
  stat_sec <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  rs <- rater_setup(cond_row$cond, cond_row$model_type)
  per_rater <- data.frame(info, hs$per_rater,
                          halo_sd_true = unname(rs$halo_sd[hs$per_rater$rater]),
                          disc_true = unname(rs$disc[hs$per_rater$rater]),
                          stretch_true = unname(rs$stretch[hs$per_rater$rater]))
  res <- list(
    info = data.frame(info, rhat_max = rhat_max, n_div = n_div,
                      fit_sec = fit_sec, stat_sec = stat_sec),
    pooled = data.frame(info, hs$pooled),
    per_rater = per_rater,
    pair = pair
  )
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, f)
  invisible(NULL)
}

jobs <- do.call(rbind, lapply(seq_len(nrow(conditions)), function(i) {
  data.frame(ci = i, rep = seq_len(conditions$reps[i]))
}))
# Interleave conditions, so a stopped run has some of each
jobs <- jobs[order(jobs$rep, jobs$ci), ]
cat(format(Sys.time()), "starting", nrow(jobs), "fits on", N_WORKERS,
    "workers\n")
status <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  id <- conditions$id[jobs$ci[j]]
  t0 <- Sys.time()
  out <- tryCatch({
    run_one(conditions[jobs$ci[j], ], jobs$rep[j])
    "ok"
  }, error = function(e) paste("error:", conditionMessage(e)))
  cat(format(Sys.time()), id, "rep", jobs$rep[j], out, "in",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      "min\n")
  out
}, mc.cores = N_WORKERS, mc.preschedule = FALSE)
status <- unlist(status)
cat(format(Sys.time()), "finished.", sum(status == "ok"), "ok,",
    sum(status != "ok"), "errors\n")
if (any(status != "ok")) print(table(status[status != "ok"]))

# --- Combine (dated and undated copies) ---
files <- list.files(out_dir, pattern = "^rep_.*\\.rds$", recursive = TRUE,
                    full.names = TRUE)
reps <- lapply(files, readRDS)
parts <- c("info", "pooled", "per_rater", "pair")
all_res <- lapply(parts, function(el) dplyr::bind_rows(lapply(reps, `[[`, el)))
names(all_res) <- parts
name <- if (smoke) "halo_test_smoke" else "halo_test"
saveRDS(all_res, file.path("dev", "results", paste0(name, ".rds")))
saveRDS(all_res, file.path("dev", "results",
                           paste0(name, Sys.Date(), ".rds")))
cat("saved", length(reps), "replications\n")
