# Null calibration of a max-T adjusted p-value for infit_statistic()
# (easyRaschBayes 0.3.0.9000)
#
# Question: item_restscore_statistic() has p_adj and q3_post() has q3_p_adj.
# Would the same two-sided max-T over items work for infit? The conditional
# infit ppp is known to be conservative (dev/null_restscore_calibration.R:
# two-sided rejection at .05 for 0.5% to 2% of items), because the item
# parameters are fitted to the same responses. Does the max-T inherit this?
#
# Per item, from the draws of infit_statistic() (observed infit and
# replicated infit per draw):
#   ppp         mid-p P(observed > replicated), as infit_post()
#   p_marg      two-sided per-item p from the max-T construction without the
#               maximum: the posterior mean observed infit, standardised by
#               the mean and SD of the replicated infit, against the
#               standardised replicated values of the same item. Separates
#               the effect of this construction from that of the maximum.
#   p_adj       two-sided single-step max-T, as p_adj in
#               item_restscore_statistic() (.ppp_maxt): the same standardised
#               observed value against the largest absolute standardised
#               replicated value over items in each draw.
#   p_adj_draw  the same maximum, but each draw's observed value is
#               standardised and compared with that draw's maximum, averaged
#               over draws. A ppp-style variant of p_adj.
# Each for both methods, conditional (default) and theta.
#
# Design as in dev/null_restscore_calibration.R (3 models x 2 sample sizes,
# REPS replications, 10 dichotomous or 8 PCM items), with new seeds.
# Expected run time about 3 hours on 10 workers.
#
# Run from the package root on the host:
#   nohup Rscript dev/null_infit_adj_calibration.R > dev/results/null_infit_adj.log 2>&1 &
# Smoke test (2 reps, small fits, 2 workers):
#   SMOKE=1 Rscript dev/null_infit_adj_calibration.R
#
# Each replication is saved on its own, so an interrupted run resumes where
# it stopped. Summarise with dev/null_infit_adj_calibration_summary.R.

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
NDRAWS    <- if (smoke) 100L else 1000L
SIZES     <- if (smoke) c(150L) else c(300L, 1000L)

name    <- if (smoke) "null_infit_adj_smoke" else "null_infit_adj"
out_dir <- file.path("dev", "results", name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Item parameters, as in dev/null_restscore_calibration.R ---
k_dich <- 10L
b_dich <- seq(-1.5, 1.5, length.out = k_dich)
k_pcm  <- 8L
loc_pcm <- seq(-1, 1, length.out = k_pcm)
tau_pcm <- lapply(loc_pcm, function(l) l + c(-1.2, 0, 1.2))

simulate_data <- function(model_type, n) {
  theta <- rnorm(n)
  if (model_type == "pcm") {
    X <- sapply(tau_pcm, function(tau) {
      ct <- c(0, cumsum(tau))
      num <- exp(outer(theta, 0:length(tau)) - matrix(ct, n, length(ct), byrow = TRUE))
      p <- num / rowSums(num)
      u <- runif(n)
      rowSums(u > t(apply(p, 1, cumsum))[, -ncol(p), drop = FALSE])
    })
    resp <- as.vector(X) + 1L   # acat responses are coded from 1
  } else {
    X <- sapply(b_dich, function(b) rbinom(n, 1, plogis(theta - b)))
    resp <- as.vector(X)
  }
  k <- ncol(X)
  data.frame(
    id       = factor(rep(seq_len(n), times = k)),
    item     = factor(rep(sprintf("I%02d", seq_len(k)), each = n)),
    response = resp
  )
}

model_formula <- function(model_type) {
  switch(model_type,
         dich_ri = bf(response ~ 1 + (1 | item) + (1 | id)),
         dich_fi = bf(response ~ 0 + item + (1 | id)),
         pcm     = bf(response | thres(gr = item) ~ 1 + (1 | id)))
}
model_family <- function(model_type) {
  if (model_type == "pcm") acat() else bernoulli()
}

# --- p-values from draws x items matrices of observed and replicated infit ---
infit_pvalues <- function(d, method) {
  d <- as.data.frame(d)
  items <- sort(unique(as.character(d$item)))
  to_mat <- function(col) {
    sapply(items, function(i) {
      x <- d[d$item == i, ]
      x[[col]][order(x$draw)]
    })
  }
  obs <- to_mat("infit")
  rep <- to_mat("infit_rep")
  m <- colMeans(rep, na.rm = TRUE)
  s <- apply(rep, 2L, stats::sd, na.rm = TRUE)
  z_rep <- abs(sweep(sweep(rep, 2L, m), 2L, s, "/"))
  z_obs_draw <- abs(sweep(sweep(obs, 2L, m), 2L, s, "/"))
  z_obs <- abs(colMeans(obs, na.rm = TRUE) - m) / s
  t_max <- apply(z_rep, 1L, max, na.rm = TRUE)
  data.frame(
    item = items, method = method,
    infit_obs = colMeans(obs, na.rm = TRUE), infit_rep = m,
    sd_obs = apply(obs, 2L, stats::sd, na.rm = TRUE), sd_rep = s,
    ppp = colMeans(obs > rep, na.rm = TRUE) + 0.5 * colMeans(obs == rep, na.rm = TRUE),
    p_marg = vapply(seq_along(items), function(j) mean(z_rep[, j] >= z_obs[j]), 0),
    p_adj = vapply(z_obs, function(z) mean(t_max >= z), 0),
    p_adj_draw = colMeans(t_max >= z_obs_draw, na.rm = TRUE),
    # the package helper, as a check that p_adj matches it
    p_adj_pkg = .ppp_maxt(colMeans(obs, na.rm = TRUE), rep),
    row.names = NULL
  )
}

conditions <- expand.grid(model_type = c("dich_ri", "dich_fi", "pcm"),
                          n = SIZES, stringsAsFactors = FALSE)
conditions$cond <- paste0(conditions$model_type, "_n", conditions$n)

# --- Compile each model once, then reuse with update(recompile = FALSE) ---
base_fits <- list()
for (mt in unique(conditions$model_type)) {
  set.seed(1)
  base_fits[[mt]] <- brm(model_formula(mt), data = simulate_data(mt, 100),
                         family = model_family(mt), chains = 0, refresh = 0)
}

run_one <- function(cond_row, rep) {
  f <- file.path(out_dir, cond_row$cond, sprintf("rep_%03d.rds", rep))
  if (file.exists(f)) return(invisible(NULL))
  seed <- 1e5 * match(cond_row$cond, conditions$cond) + rep + 9e6
  set.seed(seed)
  dat <- simulate_data(cond_row$model_type, cond_row$n)
  t0 <- Sys.time()
  fit <- suppressMessages(update(
    base_fits[[cond_row$model_type]], newdata = dat, recompile = FALSE,
    chains = CHAINS, cores = 1, iter = ITER, warmup = WARMUP,
    refresh = 0, seed = seed
  ))
  fit_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rhat_max <- max(brms::rhat(fit), na.rm = TRUE)
  n_div <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)

  res <- list()
  for (m in c("conditional", "theta")) {
    set.seed(seed + if (m == "conditional") 1 else 2)
    inf <- infit_statistic(fit, ndraws_use = NDRAWS, method = m)
    res[[m]] <- infit_pvalues(inf, m)
  }
  res <- do.call(rbind, res)
  res <- cbind(data.frame(cond = cond_row$cond, model_type = cond_row$model_type,
                          n = cond_row$n, rep = rep),
               res,
               data.frame(rhat_max = rhat_max, n_div = n_div, fit_sec = fit_sec))
  rownames(res) <- NULL
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, f)
  invisible(NULL)
}

jobs <- expand.grid(ci = seq_len(nrow(conditions)), rep = seq_len(REPS))
cat(format(Sys.time()), "starting", nrow(jobs), "fits on", N_WORKERS, "workers\n")
status <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  tryCatch({
    run_one(conditions[jobs$ci[j], ], jobs$rep[j])
    "ok"
  }, error = function(e) paste("error:", conditionMessage(e)))
}, mc.cores = N_WORKERS, mc.preschedule = FALSE)
status <- unlist(status)
cat(format(Sys.time()), "finished.", sum(status == "ok"), "ok,",
    sum(status != "ok"), "errors\n")
if (any(status != "ok")) print(table(status[status != "ok"]))

# --- Combine (dated and undated copies) ---
files <- list.files(out_dir, pattern = "^rep_.*\\.rds$", recursive = TRUE,
                    full.names = TRUE)
all_res <- dplyr::bind_rows(lapply(files, readRDS))
saveRDS(all_res, file.path("dev", "results", paste0(name, ".rds")))
saveRDS(all_res, file.path("dev", "results", paste0(name, Sys.Date(), ".rds")))
cat("saved", nrow(all_res), "rows\n")
