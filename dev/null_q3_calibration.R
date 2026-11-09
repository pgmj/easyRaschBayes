# Null calibration of q3_statistic() and a conditional Q3 prototype
# (easyRaschBayes 0.3.0.9000)
#
# Question: under data generated from the fitted Rasch model, how
# conservative is the Q3 ppp from q3_statistic() (residuals from posterior
# draws of the person locations), and is a Q3 conditional on the total score
# close to uniform, as the conditional restscore ppp was?
#
# Statistics per item pair, from the same fit:
#   theta   q3_statistic(): residuals X - E(X | theta draw), replicated data
#           from posterior_predict() with the same draw. ppp = mid-p
#           P(observed Q3 > replicated Q3).
#   cond    prototype: residuals X - E(X | total score), with E from the
#           ESF engine (R/conditional_score.R) for the item parameters of
#           each draw, replicated data drawn given each person's total
#           score. Persons with the minimum or maximum total have zero
#           conditional variance and drop out, as in conditional infit.
#           ppp = mid-p P(observed Q3 > replicated Q3). Also a one-sided
#           max-T adjusted ppp over pairs (cond_ppp_adj).
# Local dependence is one-sided (Q3 above the model's), so the summary
# reports one-sided rejection (ppp > .95) and two-sided rates.
#
# Design as in dev/null_restscore_calibration.R, so results can be read
# together: 3 models x 2 sample sizes, REPS replications each.
#   dich_ri  dichotomous, response ~ 1 + (1 | item) + (1 | id), 10 items
#   dich_fi  dichotomous, response ~ 0 + item + (1 | id),       10 items
#   pcm      partial credit, response | thres(gr = item) ~ 1 + (1 | id),
#            8 items with 4 categories
#   n = 300 and 1000 persons, theta ~ N(0, 1)
# The previous run took about 3 hours on 10 workers. The Q3 statistics add
# somewhat less than the restscore and infit statistics did.
#
# Run from the package root on the host:
#   nohup Rscript dev/null_q3_calibration.R > dev/results/null_q3.log 2>&1 &
# Smoke test (2 reps, small fits, 2 workers):
#   SMOKE=1 Rscript dev/null_q3_calibration.R
#
# Each replication is saved on its own, so an interrupted run resumes where
# it stopped. Summarise with dev/null_q3_calibration_summary.R.

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
NDRAWS    <- if (smoke) 100L else 1000L   # draws used by the statistics
SIZES     <- if (smoke) c(150L) else c(300L, 1000L)

name    <- if (smoke) "null_q3_smoke" else "null_q3"
out_dir <- file.path("dev", "results", name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Item parameters, fixed across replications (as in the restscore run) ---
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

# --- Conditional Q3 prototype -------------------------------------------
# One draw: Q3 (upper triangle, row-major pair order) for observed and
# replicated residuals, both relative to E(X | total score).
cond_q3_draw <- function(prob, person_idx, item_idx, x_obs) {
  k <- ncol(x_obs)
  st <- .cond_draw_setup(prob, person_idx, item_idx, x_obs)
  mo <- .cond_moments(st, k)
  x_rep <- .cond_sample(st$W, st$suffix, st$key, st$r)
  x_rep[is.na(x_obs)] <- NA
  use <- !is.na(mo$V) & mo$V > 1e-12
  r_obs <- ifelse(use, x_obs - mo$E, NA_real_)
  r_rep <- ifelse(use, x_rep - mo$E, NA_real_)
  ut <- upper.tri(diag(k))
  c_obs <- suppressWarnings(stats::cor(r_obs, use = "pairwise.complete.obs"))
  c_rep <- suppressWarnings(stats::cor(r_rep, use = "pairwise.complete.obs"))
  list(obs = t(c_obs)[t(ut)], rep = t(c_rep)[t(ut)])
}

cond_q3 <- function(fit, draw_ids) {
  dd <- .cond_prepare_data(fit, "response", "item", "id")
  per_draw <- .cond_map_draws(fit, draw_ids, dd$n_cat_max, cond_q3_draw,
                              person_idx = dd$person_idx,
                              item_idx = dd$item_idx, x_obs = dd$x_obs)
  obs <- do.call(rbind, lapply(per_draw, `[[`, "obs"))
  rep <- do.call(rbind, lapply(per_draw, `[[`, "rep"))
  k <- length(dd$unique_items)
  pg <- which(upper.tri(diag(k)), arr.ind = TRUE)
  pg <- pg[order(pg[, "row"], pg[, "col"]), , drop = FALSE]
  labels <- paste0(dd$unique_items[pg[, "row"]], " : ",
                   dd$unique_items[pg[, "col"]])
  ppp <- colMeans(obs > rep) + 0.5 * colMeans(obs == rep)
  # one-sided max-T over pairs, observed at its posterior mean
  m <- colMeans(rep)
  s <- apply(rep, 2, stats::sd)
  t_max <- apply(sweep(sweep(rep, 2, m), 2, s, "/"), 1, max)
  z_obs <- (colMeans(obs) - m) / s
  # p = P(max_T >= z): small = dependence; stored as 1 - p to share the
  # "near 1 flags a pair" direction with ppp
  ppp_adj <- 1 - vapply(z_obs, function(z) mean(t_max >= z), numeric(1))
  data.frame(item_pair = labels,
             cond_q3_obs = colMeans(obs), cond_q3_rep = colMeans(rep),
             cond_ppp = ppp, cond_ppp_adj = ppp_adj,
             cond_sd_obs = apply(obs, 2, stats::sd), cond_sd_rep = s)
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
  seed <- 1e5 * match(cond_row$cond, conditions$cond) + rep + 7e6
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

  # theta-based Q3 (the package function)
  set.seed(seed + 1)
  t1 <- Sys.time()
  q3 <- as.data.frame(q3_statistic(fit, ndraws_use = NDRAWS))
  theta_sec <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  th <- do.call(rbind, lapply(split(q3, q3$item_pair), function(d) {
    data.frame(item_pair = d$item_pair[1],
               theta_q3_obs = mean(d$q3), theta_q3_rep = mean(d$q3_rep),
               theta_ppp = mean(d$q3 > d$q3_rep) + 0.5 * mean(d$q3 == d$q3_rep),
               theta_sd_obs = stats::sd(d$q3), theta_sd_rep = stats::sd(d$q3_rep))
  }))

  # conditional prototype
  set.seed(seed + 2)
  t2 <- Sys.time()
  cq <- cond_q3(fit, sample(ndraws(fit), NDRAWS))
  cond_sec <- as.numeric(difftime(Sys.time(), t2, units = "secs"))

  res <- merge(th, cq, by = "item_pair")
  res <- cbind(data.frame(cond = cond_row$cond, model_type = cond_row$model_type,
                          n = cond_row$n, rep = rep),
               res,
               data.frame(rhat_max = rhat_max, n_div = n_div, fit_sec = fit_sec,
                          theta_sec = theta_sec, cond_sec = cond_sec))
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
