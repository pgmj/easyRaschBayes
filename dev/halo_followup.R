# Follow-up to dev/halo_test_quick.R (2026-11-16): can rater infit and
# rater-specific thresholds tell a halo apart from other rater problems?
# See dev/mfrm-design-plan.md, step 6, quick run results.
#
# Background. In the quick run the pooled same-rater Q3 detected halo well
# (RSM .78-.86) but was not specific. Rater misfit without a halo raised
# it (RM .12, RSM .25). The inconsistent rater R2 was not flagged, but its
# neighbours R1 and R3 (who rate the same persons) were, and the
# central-tendency rater R5 was flagged directly.
#
# Questions:
#   (a) With rater-specific thresholds in the model (thres(gr = rater)),
#       does the R5 flag go away, does the halo still show, and is the
#       pooled test still calibrated under the null?
#   (b) Does conditional rater infit (by = "facet") flag R2 in the same
#       datasets, and does it stay quiet for halo raters? If so, infit plus
#       the per-rater Q3 together separate an inconsistent rater, a
#       central-tendency rater and a halo.
#
# Datasets. The same generator, design and seeds as dev/halo_test_quick.R,
# so replication r of a condition is the same dataset as in the quick run.
# The arm "common" refits the quick run's model with the same seed, so its
# Q3 results should reproduce the quick run (checked in the summary).
#
# Arms and conditions (n = 300, 100 replications each, 700 fits):
#   common  model of the quick run: response ~ 1 + item + (1 | rater) +
#           (1 | id), common thresholds. RM and RSM, conditions misfit and
#           halo_two. Adds rater infit to the quick run's Q3 statistics.
#   thres   RSM with rater-specific thresholds:
#           response | thres(gr = rater) ~ 1 + item + (1 | id).
#           Conditions null, misfit, halo_two. Q3 statistics, rater
#           infit, and each rater's threshold estimates (spread).
#   Conditions as in the quick run: misfit = R2 discrimination 0.4 (and in
#   the RSM R5 thresholds x 1.8), halo_two = halo SD 0.6 for R3 and R7.
#   Each person is rated by two adjacent raters (cyclic), so R1 and R3
#   share persons with R2, R4 and R6 with R5, and R2, R4, R6, R8 with the
#   halo raters.
#   The null of the thres arm has 100 replications (MC SE about .022 at a
#   rate of .05). Raise REPS_THRES_NULL and rerun to add more, since
#   finished replications are skipped.
#
# Statistics per replication:
#   Q3 (seed + 1, as in the quick run): pooled mean_same and contrast,
#     per-rater mean same-rater Q3 with one-sided max-T p_adj, smallest
#     pair-level q3_p_adj.
#   Rater infit (seed + 2): conditional infit_statistic(by = "facet"),
#     per rater mid-p ppp (P(obs > rep), near 1 = underfit) and two-sided
#     max-T p_adj over raters (.ppp_maxt(), as p_adj in
#     item_restscore_statistic()).
#   Thresholds (thres arm): facet_parameters()$thresholds, posterior mean
#     of each rater's centred thresholds.
#
# Run from the package root on the host (10 workers):
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 nohup Rscript dev/halo_followup.R > dev/results/halo_followup.log 2>&1 &
# Smoke test (2 reps, n = 150, short chains, 2 workers):
#   SMOKE=1 Rscript dev/halo_followup.R
# Summarise with
#   Rscript dev/halo_followup_summary.R        (SMOKE=1 for the smoke run)
#
# Rough cost: quick run fits took 11 s (RM) and 80 s (RSM), Q3 5-12 s. In
# the smoke test the thres fits took 0.75 times the common RSM fits, and
# infit 0.75 times Q3. Expect about 1.5 hours on 10 workers. Each
# replication is saved on its own, so an interrupted run resumes.
#
# Smoke test 2026-11-16: 14 fits, 0 errors. The common arm reproduced a
# fresh smoke run of dev/halo_test_quick.R exactly (pooled p and per-rater
# p_adj identical), so the full run should reproduce the quick run.

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
NDRAWS    <- if (smoke) 200L else 1000L
N_PERSONS <- if (smoke) 150L else 300L
REPS      <- if (smoke) 2L else 100L
REPS_THRES_NULL <- if (smoke) 2L else 100L

name <- if (smoke) "halo_followup_smoke" else "halo_followup"
out_dir <- file.path("dev", "results", name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Generator, copied unchanged from dev/halo_test_quick.R ---
items  <- sprintf("I%d", 1:5)
raters <- sprintf("R%d", 1:8)
b_item   <- setNames(seq(-0.8, 0.8, length.out = 5), items)
severity <- setNames(seq(-0.5, 0.5, length.out = 8), raters)
tau_rsm  <- c(-1.2, 0, 1.2)

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
  h <- matrix(rnorm(n * 8L), n, 8L) %*% diag(rs$halo_sd)
  ri <- match(d$rater, raters)
  loc <- theta[d$id] + h[cbind(d$id, ri)] - b_item[d$item] - severity[ri]
  a <- rs$disc[ri]
  if (model_type == "rsm") {
    tau <- outer(rs$stretch[ri], tau_rsm)
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

# Seeds of the quick run: 1e5 * (position of the condition id in its
# conditions table) + rep. Its table was expand.grid(model_type =
# c("rm", "rsm"), cond = c("null", "misfit", "halo03", "halo_two")).
quick_ids <- with(expand.grid(model_type = c("rm", "rsm"),
                              cond = c("null", "misfit", "halo03",
                                       "halo_two"),
                              stringsAsFactors = FALSE),
                  paste0(model_type, "_", cond, "_n", N_PERSONS))
quick_seed <- function(model_type, cond, rep) {
  1e5 * match(paste0(model_type, "_", cond, "_n", N_PERSONS), quick_ids) + rep
}

# --- Models ---
model_formula <- function(arm) {
  if (arm == "thres") {
    bf(response | thres(gr = rater) ~ 1 + item + (1 | id))
  } else {
    bf(response ~ 1 + item + (1 | rater) + (1 | id))
  }
}
model_family <- function(model_type) {
  if (model_type == "rsm") acat() else bernoulli()
}
model_prior <- function() {
  prior(normal(0, 3), class = "Intercept") +
    prior(normal(0, 3), class = "b") +
    prior(exponential(1), class = "sd")
}

conditions <- rbind(
  expand.grid(arm = "common", model_type = c("rm", "rsm"),
              cond = c("misfit", "halo_two"), stringsAsFactors = FALSE),
  data.frame(arm = "thres", model_type = "rsm",
             cond = c("null", "misfit", "halo_two"))
)
conditions$n <- N_PERSONS
conditions$reps <- ifelse(conditions$arm == "thres" &
                            conditions$cond == "null", REPS_THRES_NULL, REPS)
conditions$key <- paste0(conditions$arm, "_", conditions$model_type)
conditions$id <- paste0(conditions$key, "_", conditions$cond, "_n",
                        conditions$n)

# Compile each model once (chains = 0), reuse with update(recompile = FALSE)
base_fits <- list()
for (i in seq_len(nrow(conditions))) {
  key <- conditions$key[i]
  if (!is.null(base_fits[[key]])) next
  set.seed(1)
  t0 <- Sys.time()
  invisible(utils::capture.output(suppressMessages(
    base_fits[[key]] <- brm(
      model_formula(conditions$arm[i]),
      data = simulate_data(conditions$model_type[i], "null", 40),
      family = model_family(conditions$model_type[i]),
      prior = model_prior(), chains = 0, refresh = 0
    )
  )))
  cat(format(Sys.time()), "compiled", key, "in",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}

# --- Q3 halo statistics, copied unchanged from dev/halo_test_quick.R ---
midp_upper <- function(obs, rep) {
  ok <- !is.na(obs) & !is.na(rep)
  mean(rep[ok] > obs[ok]) + 0.5 * mean(rep[ok] == obs[ok])
}

halo_stats <- function(q3) {
  q3 <- as.data.frame(q3)
  draws <- sort(unique(q3$draw))
  same <- q3$pair_type == "same rater"
  q3$rater_1 <- sub("^.*:", "", as.character(q3$item_1))
  per_draw <- function(sel, col) {
    v <- tapply(q3[[col]][sel], q3$draw[sel], mean, na.rm = TRUE)
    as.numeric(v[as.character(draws)])
  }
  obs_same <- per_draw(same, "q3")
  rep_same <- per_draw(same, "q3_rep")
  obs_other <- per_draw(!same, "q3")
  rep_other <- per_draw(!same, "q3_rep")
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

# --- Rater infit: mid-p ppp and two-sided max-T over raters ---
# As ppp_by() and infit_maxt() in dev/null_mfrm_calibration.R.
infit_by_rater <- function(inf) {
  d <- as.data.frame(inf)
  keys <- sort(unique(as.character(d$rater)))
  rep_mat <- sapply(keys, function(k) {
    x <- d[d$rater == k, ]
    x$infit_rep[order(x$draw)]
  })
  obs <- vapply(keys, function(k) mean(d$infit[d$rater == k]), numeric(1))
  data.frame(
    rater = keys,
    infit_obs = obs,
    infit_rep = colMeans(rep_mat),
    infit_ppp = vapply(keys, function(k) {
      x <- d[d$rater == k, ]
      mean(x$infit > x$infit_rep) + 0.5 * mean(x$infit == x$infit_rep)
    }, numeric(1)),
    infit_p_adj = .ppp_maxt(obs, rep_mat),
    row.names = NULL
  )
}

run_one <- function(cond_row, rep) {
  f <- file.path(out_dir, cond_row$id, sprintf("rep_%03d.rds", rep))
  if (file.exists(f)) return(invisible(NULL))
  seed <- quick_seed(cond_row$model_type, cond_row$cond, rep)
  set.seed(seed)
  dat <- simulate_data(cond_row$model_type, cond_row$cond, cond_row$n)
  t0 <- Sys.time()
  fit <- suppressMessages(update(
    base_fits[[cond_row$key]],
    newdata = dat, recompile = FALSE, chains = CHAINS, cores = 1,
    iter = ITER, warmup = WARMUP, refresh = 0, seed = seed
  ))
  fit_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rhat_max <- max(brms::rhat(fit), na.rm = TRUE)
  n_div <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)
  info <- data.frame(id = cond_row$id, arm = cond_row$arm,
                     model_type = cond_row$model_type, cond = cond_row$cond,
                     n = cond_row$n, rep = rep)

  # Q3, with the quick run's seed for the draws
  t1 <- Sys.time()
  set.seed(seed + 1)
  q3 <- suppressMessages(q3_statistic(fit, facet_var = rater,
                                      ndraws_use = NDRAWS))
  hs <- halo_stats(q3)
  q3s <- q3_post(q3)$summary
  pair <- data.frame(
    info,
    min_p_adj_same = min(q3s$q3_p_adj[q3s$pair_type == "same rater"],
                         na.rm = TRUE),
    min_p_adj_all = min(q3s$q3_p_adj, na.rm = TRUE)
  )
  q3_sec <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  # Rater infit
  t2 <- Sys.time()
  set.seed(seed + 2)
  inf <- suppressMessages(infit_statistic(fit, facet_var = rater,
                                          by = "facet", ndraws_use = NDRAWS))
  ir <- infit_by_rater(inf)
  infit_sec <- as.numeric(difftime(Sys.time(), t2, units = "secs"))

  # Rater thresholds (thres arm)
  thresholds <- NULL
  if (cond_row$arm == "thres") {
    th <- facet_parameters(fit, facet_var = rater)$thresholds
    thresholds <- data.frame(info, rater = as.character(th$rater),
                             threshold = th$threshold,
                             location = th$location)
  }

  rs <- rater_setup(cond_row$cond, cond_row$model_type)
  truth <- function(r) {
    data.frame(halo_sd_true = unname(rs$halo_sd[r]),
               disc_true = unname(rs$disc[r]),
               stretch_true = unname(rs$stretch[r]))
  }
  per_rater <- merge(hs$per_rater, ir, by = "rater")
  per_rater <- data.frame(info, per_rater, truth(per_rater$rater))
  res <- list(
    info = data.frame(info, rhat_max = rhat_max, n_div = n_div,
                      fit_sec = fit_sec, q3_sec = q3_sec,
                      infit_sec = infit_sec),
    pooled = data.frame(info, hs$pooled),
    per_rater = per_rater,
    pair = pair,
    thresholds = thresholds
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
parts <- c("info", "pooled", "per_rater", "pair", "thresholds")
all_res <- lapply(parts, function(el) dplyr::bind_rows(lapply(reps, `[[`, el)))
names(all_res) <- parts
saveRDS(all_res, file.path("dev", "results", paste0(name, ".rds")))
saveRDS(all_res, file.path("dev", "results",
                           paste0(name, Sys.Date(), ".rds")))
cat("saved", length(reps), "replications\n")
