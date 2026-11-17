# Null calibration of the many-facet fit statistics (easyRaschBayes
# 0.4.0.9001, MFRM step 5, see dev/mfrm-design-plan.md)
#
# Questions:
#   1. Under data generated from the fitted MFRM, are the conditional
#      restscore ppp and p_adj per item x rater cell close to uniform in a
#      sparse rating design with many cells (40)? The standard-design null
#      run (null_restscore, 2026-11-06) found them close to uniform.
#   2. Conditional infit per cell, pooled by rater and pooled by item: how
#      conservative? (Standard designs: two-sided rej05 0.5-2%.)
#   3. Conditional Q3 per cell pair by pair type, with q3_p_adj.
#   4. With a halo term (1 | rater:id) in data and model, are 1-3 still
#      calibrated? This settles the halo message (plan, open questions).
#   5. Power check, not null: halo in the data but not in the model. Does
#      same-rater Q3 detect it, and what do restscore and infit do?
#   6. facet_parameters(): coverage of the 95% HDCI for rater severity
#      (centred truth), and fair-average sanity.
#
# Design: 2 models x 3 halo settings x 2 sample sizes, REPS replications.
#   model  rm   dichotomous, response ~ 1 + item + (1 | rater) + (1 | id)
#          rsm  rating scale, acat, 4 categories, same formula
#   halo   none        no halo in data or model
#          modelled    halo SD 0.6 in data, + (1 | rater:id) in the model
#          unmodelled  halo SD 0.6 in data, model without the halo term
#                      (power check; leave out of the null summaries)
#   n      300 and 600 persons, theta ~ N(0, 1)
#   5 items, 8 raters; each person is rated by 2 adjacent raters (cyclic),
#   so 10 of the 40 cells per person, connected design.
#   Items are fixed effects with default (treatment) contrasts, the most
#   common user setup; the statistics do not depend on the coding.
#
# Run from the package root on the host:
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 nohup Rscript dev/null_mfrm_calibration.R > dev/results/null_mfrm.log 2>&1 &
# Smoke test (2 reps, n = 150, short chains, 2 workers):
#   SMOKE=1 Rscript dev/null_mfrm_calibration.R
#
# Each replication is saved on its own (one list per rep with tables for
# cells, raters, items, Q3 pairs and facet parameters), so an interrupted
# run resumes where it stopped. Summarise with
# dev/null_mfrm_calibration_summary.R.
#
# Rough cost: 2,400 fits. Halo models have 600-1,200 extra parameters, and
# with a halo term every person is its own weight group in the conditional
# engine, so those replications are the slow ones. Measured in the
# container (2026-11-13, one core, simulated probabilities, no fit): one
# conditional draw for n = 600 with every person its own group takes about
# 0.2 s per statistic call, so 5 calls x NDRAWS = 1000 draws is about 17
# minutes of statistics per halo-modelled replication at n = 600, on top of
# the fit. Run the smoke test first and extrapolate from its fit_sec and
# stat_sec columns.

stopifnot(file.exists("DESCRIPTION"),
          read.dcf("DESCRIPTION")[, "Package"] == "easyRaschBayes")
smoke <- Sys.getenv("SMOKE") == "1"

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(brms)
})

N_WORKERS <- if (smoke) 2L else 10L   # max 10 cores on the M4 Max
REPS      <- if (smoke) 2L else 200L
ITER      <- if (smoke) 400L else 1500L
WARMUP    <- if (smoke) 200L else 500L
CHAINS    <- 2L
NDRAWS    <- if (smoke) 200L else 1000L   # draws used by the statistics
SIZES     <- if (smoke) c(150L) else c(300L, 600L)

out_dir <- file.path("dev", "results",
                     if (smoke) "null_mfrm_smoke" else "null_mfrm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Generating values, fixed across replications ---
items  <- sprintf("I%d", 1:5)
raters <- sprintf("R%d", 1:8)
b_item   <- setNames(seq(-0.8, 0.8, length.out = 5), items)
severity <- setNames(seq(-0.5, 0.5, length.out = 8), raters)
tau_rsm  <- c(-1.2, 0, 1.2)
HALO_SD  <- 0.6

simulate_data <- function(model_type, halo, n) {
  theta <- rnorm(n)
  d <- do.call(rbind, lapply(seq_len(n), function(p) {
    r <- raters[((p - 1L) + 0:1) %% 8L + 1L]
    expand.grid(id = p, item = items, rater = r, stringsAsFactors = FALSE)
  }))
  h <- if (halo == "none") rep(0, n * 8L) else rnorm(n * 8L, 0, HALO_SD)
  loc <- theta[d$id] + h[(d$id - 1L) * 8L + match(d$rater, raters)] -
    b_item[d$item] - severity[d$rater]
  if (model_type == "rsm") {
    lp <- cbind(0, outer(loc, tau_rsm, "-"))
    lp <- t(apply(lp, 1L, cumsum))
    p <- exp(lp - apply(lp, 1L, max))
    p <- p / rowSums(p)
    u <- runif(nrow(d))
    d$response <- 1L + rowSums(u > t(apply(p, 1L, cumsum))[, -4L,
                                                           drop = FALSE])
  } else {
    d$response <- rbinom(nrow(d), 1L, plogis(loc))
  }
  d$id <- factor(d$id)
  d$item <- factor(d$item, levels = items)
  d$rater <- factor(d$rater, levels = raters)
  d
}

halo_in_model <- function(halo) halo == "modelled"
model_key <- function(model_type, halo) {
  paste0(model_type, if (halo_in_model(halo)) "_halo" else "")
}
model_formula <- function(model_type, halo) {
  f <- response ~ 1 + item + (1 | rater) + (1 | id)
  if (halo_in_model(halo)) f <- update(f, . ~ . + (1 | rater:id))
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

conditions <- expand.grid(model_type = c("rm", "rsm"),
                          halo = c("none", "modelled", "unmodelled"),
                          n = SIZES, stringsAsFactors = FALSE)
conditions$cond <- paste0(conditions$model_type, "_", conditions$halo,
                          "_n", conditions$n)

# --- Compile each model once, then reuse with update(recompile = FALSE) ---
# chains = 0 compiles without sampling. rstan then prints "the number of
# chains is less than 1; sampling not done", which is expected, so the
# output is captured and replaced by one line per model.
base_fits <- list()
for (i in seq_len(nrow(conditions))) {
  key <- model_key(conditions$model_type[i], conditions$halo[i])
  if (!is.null(base_fits[[key]])) next
  set.seed(1)
  t0 <- Sys.time()
  invisible(utils::capture.output(suppressMessages(
    base_fits[[key]] <- brm(
      model_formula(conditions$model_type[i], conditions$halo[i]),
      data = simulate_data(conditions$model_type[i], conditions$halo[i], 40),
      family = model_family(conditions$model_type[i]),
      prior = model_prior(), chains = 0, refresh = 0
    )
  )))
  cat(format(Sys.time()), "compiled", key, "in",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}

# Mid-p ppp per group from long draws (obs, rep columns).
ppp_by <- function(d, key, obs, rep) {
  d <- as.data.frame(d)
  s <- split(d, d[[key]])
  data.frame(
    key = names(s),
    ppp = vapply(s, function(x) {
      mean(x[[obs]] > x[[rep]]) + 0.5 * mean(x[[obs]] == x[[rep]])
    }, numeric(1)),
    obs_mean = vapply(s, function(x) mean(x[[obs]]), numeric(1)),
    rep_mean = vapply(s, function(x) mean(x[[rep]]), numeric(1)),
    row.names = NULL
  )
}

# Max-T p-value for infit groups from long draws, the same procedure as
# p_adj in item_restscore_statistic() (two-sided).
infit_maxt <- function(d, key) {
  d <- as.data.frame(d)
  keys <- sort(unique(as.character(d[[key]])))
  rep_mat <- sapply(keys, function(k) {
    x <- d[d[[key]] == k, ]
    x$infit_rep[order(x$draw)]
  })
  obs <- sapply(keys, function(k) mean(d$infit[d[[key]] == k]))
  setNames(.ppp_maxt(obs, rep_mat), keys)
}

run_one <- function(cond_row, rep) {
  f <- file.path(out_dir, cond_row$cond, sprintf("rep_%03d.rds", rep))
  if (file.exists(f)) return(invisible(NULL))
  seed <- 1e5 * match(cond_row$cond, conditions$cond) + rep
  set.seed(seed)
  dat <- simulate_data(cond_row$model_type, cond_row$halo, cond_row$n)
  t0 <- Sys.time()
  fit <- suppressMessages(update(
    base_fits[[model_key(cond_row$model_type, cond_row$halo)]],
    newdata = dat, recompile = FALSE, chains = CHAINS, cores = 1,
    iter = ITER, warmup = WARMUP, refresh = 0, seed = seed
  ))
  fit_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rhat_max <- max(brms::rhat(fit), na.rm = TRUE)
  n_div <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)

  t1 <- Sys.time()
  info <- data.frame(cond = cond_row$cond, model_type = cond_row$model_type,
                     halo = cond_row$halo, n = cond_row$n, rep = rep)

  # 1. Restscore per cell
  set.seed(seed + 1)
  rs <- suppressMessages(item_restscore_statistic(
    fit, facet_var = rater, ndraws_use = NDRAWS))$result
  cells <- data.frame(info, cell = as.character(rs$cell),
                      item = as.character(rs$item),
                      rater = as.character(rs$rater),
                      gamma_obs = rs$gamma_obs, gamma_exp = rs$gamma_exp,
                      rs_ppp = rs$ppp, rs_p_adj = rs$p_adj)

  # 2. Infit per cell, by rater, by item (conditional)
  set.seed(seed + 2)
  inf_cell <- suppressMessages(infit_statistic(
    fit, facet_var = rater, ndraws_use = NDRAWS))
  pc <- ppp_by(inf_cell, "cell", "infit", "infit_rep")
  cells$infit_ppp <- pc$ppp[match(cells$cell, pc$key)]
  cells$infit_obs <- pc$obs_mean[match(cells$cell, pc$key)]
  set.seed(seed + 3)
  inf_rater <- suppressMessages(infit_statistic(
    fit, facet_var = rater, by = "facet", ndraws_use = NDRAWS))
  pr <- ppp_by(inf_rater, "rater", "infit", "infit_rep")
  adj_r <- infit_maxt(inf_rater, "rater")
  set.seed(seed + 4)
  inf_item <- suppressMessages(infit_statistic(
    fit, facet_var = rater, by = "item", ndraws_use = NDRAWS))
  pit <- ppp_by(inf_item, "item", "infit", "infit_rep")
  adj_i <- infit_maxt(inf_item, "item")
  groups <- rbind(
    data.frame(info, level = "rater", key = pr$key, infit_ppp = pr$ppp,
               infit_obs = pr$obs_mean, infit_p_adj = adj_r[pr$key]),
    data.frame(info, level = "item", key = pit$key, infit_ppp = pit$ppp,
               infit_obs = pit$obs_mean, infit_p_adj = adj_i[pit$key])
  )

  # 3. Q3 per cell pair (conditional), summarised with q3_post()
  set.seed(seed + 5)
  q3 <- suppressMessages(q3_statistic(fit, facet_var = rater,
                                      ndraws_use = NDRAWS))
  q3s <- q3_post(q3)$summary
  q3_pairs <- data.frame(info, item_pair = q3s$item_pair,
                         pair_type = q3s$pair_type, q3_obs = q3s$q3_obs,
                         q3_rep = q3s$q3_rep, q3_ppp = q3s$q3_ppp,
                         q3_p_adj = q3s$q3_p_adj)
  q3_pairs <- q3_pairs[!is.na(q3_pairs$q3_obs), ]

  # 6. Facet parameters against the centred truth
  fp <- facet_parameters(fit, facet_var = rater)$summary
  truth <- severity - mean(severity)
  facet <- data.frame(info, rater = fp$rater,
                      truth = unname(truth[fp$rater]),
                      location = fp$location, hdci_lower = fp$hdci_lower,
                      hdci_upper = fp$hdci_upper,
                      observed_avg = fp$observed_avg, fair_avg = fp$fair_avg)

  stat_sec <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  res <- list(
    info = data.frame(info, rhat_max = rhat_max, n_div = n_div,
                      fit_sec = fit_sec, stat_sec = stat_sec),
    cells = cells, groups = groups, q3 = q3_pairs, facet = facet
  )
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, f)
  invisible(NULL)
}

jobs <- expand.grid(ci = seq_len(nrow(conditions)), rep = seq_len(REPS))
cat(format(Sys.time()), "starting", nrow(jobs), "fits on", N_WORKERS,
    "workers\n")
status <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  cond <- conditions$cond[jobs$ci[j]]
  t0 <- Sys.time()
  out <- tryCatch({
    run_one(conditions[jobs$ci[j], ], jobs$rep[j])
    "ok"
  }, error = function(e) paste("error:", conditionMessage(e)))
  cat(format(Sys.time()), cond, "rep", jobs$rep[j], out, "in",
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
all_res <- lapply(c("info", "cells", "groups", "q3", "facet"), function(el) {
  dplyr::bind_rows(lapply(reps, `[[`, el))
})
names(all_res) <- c("info", "cells", "groups", "q3", "facet")
name <- if (smoke) "null_mfrm_smoke" else "null_mfrm"
saveRDS(all_res, file.path("dev", "results", paste0(name, ".rds")))
saveRDS(all_res, file.path("dev", "results",
                           paste0(name, Sys.Date(), ".rds")))
cat("saved", length(reps), "replications\n")
