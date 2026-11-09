# Null calibration of item_restscore_statistic() and infit_statistic()
# (easyRaschBayes 0.3.0.9000)
#
# Question: under data generated from the fitted Rasch model, are the new
# total-score-conditional ppp and ppp_adj close to uniform, and how do they
# compare with the old theta-based ppp (0.3.0)? The same fits are used for
# infit, conditional (new default) against theta-based (method = "theta").
#
# Design: 3 models x 2 sample sizes, REPS replications each.
#   dich_ri  dichotomous, response ~ 1 + (1 | item) + (1 | id), 10 items
#   dich_fi  dichotomous, response ~ 0 + item + (1 | id),       10 items
#   pcm      partial credit, response | thres(gr = item) ~ 1 + (1 | id),
#            8 items with 4 categories
#   n = 300 and 1000 persons, theta ~ N(0, 1)
# dich_ri is the specification in the package examples. Its item prior
# shrinks difficulties, which moves gamma_exp slightly (about 0.004 in the
# raschdat1 fit), so it is compared with fixed items.
#
# Run from the package root on the host:
#   nohup Rscript dev/null_restscore_calibration.R > dev/results/null_restscore.log 2>&1 &
# Smoke test (2 reps, small fits, 2 workers):
#   SMOKE=1 Rscript dev/null_restscore_calibration.R
#
# Each replication is saved on its own, so an interrupted run resumes where
# it stopped. Summarise with dev/null_restscore_calibration_summary.R.

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
SIZES     <- if (smoke) c(150L) else c(300L, 1000L)

out_dir <- file.path("dev", "results", if (smoke) "null_restscore_smoke" else "null_restscore")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Item parameters, fixed across replications ---
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

# The 0.3.0 statistic: replicated data from theta draws, ppp = P(obs > rep).
old_restscore_ppp <- function(fit, draw_ids) {
  yrep <- posterior_predict(fit, draw_ids = draw_ids)
  d <- fit$data
  pid <- as.integer(d$id)
  iid <- as.integer(d$item)
  n <- max(pid)
  k <- max(iid)
  lin <- (iid - 1L) * n + pid
  obs <- matrix(NA_real_, n, k)
  obs[lin] <- d$response
  tot <- rowSums(obs)
  g_obs <- vapply(seq_len(k), function(i) gk_gamma_vec(obs[, i], tot - obs[, i]), 0)
  g_rep <- matrix(NA_real_, nrow(yrep), k)
  for (s in seq_len(nrow(yrep))) {
    w <- matrix(NA_real_, n, k)
    w[lin] <- yrep[s, ]
    ts <- rowSums(w)
    for (i in seq_len(k)) g_rep[s, i] <- gk_gamma_vec(w[, i], ts - w[, i])
  }
  colMeans(matrix(g_obs, nrow(g_rep), k, byrow = TRUE) > g_rep)
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
  seed <- 1e5 * match(cond_row$cond, conditions$cond) + rep
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

  draw_ids <- sample(ndraws(fit), NDRAWS)
  set.seed(seed + 1)
  t1 <- Sys.time()
  new <- item_restscore_statistic(fit, ndraws_use = NDRAWS)$result
  new_sec <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  ppp_old <- old_restscore_ppp(fit, draw_ids)
  set.seed(seed + 2)
  inf_c <- infit_statistic(fit, ndraws_use = NDRAWS)
  inf_t <- infit_statistic(fit, ndraws_use = NDRAWS, method = "theta")
  infit_ppp <- function(d) {
    d <- as.data.frame(d)
    sapply(split(d, d$item), function(x) {
      mean(x$infit > x$infit_rep) + 0.5 * mean(x$infit == x$infit_rep)
    })
  }
  ppp_inf_c <- infit_ppp(inf_c)
  ppp_inf_t <- infit_ppp(inf_t)

  new <- new[order(new$item), ]
  res <- data.frame(
    cond = cond_row$cond, model_type = cond_row$model_type, n = cond_row$n,
    rep = rep, item = as.character(new$item),
    gamma_obs = new$gamma_obs, gamma_exp = new$gamma_exp,
    ppp = new$ppp, ppp_adj = new$p_adj,
    ppp_old = ppp_old[match(as.character(new$item), levels(dat$item))],
    infit_ppp_cond = ppp_inf_c[as.character(new$item)],
    infit_ppp_theta = ppp_inf_t[as.character(new$item)],
    rhat_max = rhat_max, n_div = n_div, fit_sec = fit_sec, new_sec = new_sec
  )
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
name <- if (smoke) "null_restscore_smoke" else "null_restscore"
saveRDS(all_res, file.path("dev", "results", paste0(name, ".rds")))
saveRDS(all_res, file.path("dev", "results", paste0(name, Sys.Date(), ".rds")))
cat("saved", nrow(all_res), "rows\n")
