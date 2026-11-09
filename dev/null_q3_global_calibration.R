# Null calibration of max-T adjusted Q3 p-values: pair-specific vs global
# reference (easyRaschBayes 0.3.0.9000)
#
# Question: Christensen, Makransky and Horton (2017) compare each pair's Q3
# with the mean Q3 over all pairs in the same dataset (Q3* = Q3max - mean Q3)
# and use the null distribution of the maximum as a single cutoff for all
# pairs. q3_post() instead centres each pair on its own mean replicated Q3
# and standardises by its replicated SD. Under data generated from the
# fitted Rasch model, do both keep the familywise error rate at .05, and do
# both spread the false flags evenly over pairs?
#
# Four max-T variants, crossing reference and scale. In every posterior draw,
# x is the vector of Q3 over pairs (observed or replicated):
#   pair_std     (x_j - m_j) / s_j             current q3_p_adj
#   pair_raw     x_j - m_j
#   global_raw   x_j - mean over pairs of x     Christensen et al.
#   global_std   (x_j - mean over pairs of x) / s_j
# m_j is pair j's mean replicated value and s_j the SD over draws of the
# replicated values (after the global centring, for the global variants).
# The reference distribution is the maximum over pairs of the replicated
# values in each draw. The observed value of a pair is the mean over draws
# of its observed values, as in q3_post(). Adjusted p = P(max >= observed),
# small values flag local dependence.
# Each variant is computed for both q3_statistic() methods, conditional
# (the default) and theta.
#
# Design as in dev/null_q3_calibration.R (3 models x 2 sample sizes, REPS
# replications), with new seeds. Expected run time about 3 hours on 10
# workers.
#
# Run from the package root on the host:
#   nohup Rscript dev/null_q3_global_calibration.R > dev/results/null_q3_global.log 2>&1 &
# Smoke test (2 reps, small fits, 2 workers):
#   SMOKE=1 Rscript dev/null_q3_global_calibration.R
#
# Each replication is saved on its own, so an interrupted run resumes where
# it stopped. Summarise with dev/null_q3_global_calibration_summary.R.

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

name    <- if (smoke) "null_q3_global_smoke" else "null_q3_global"
out_dir <- file.path("dev", "results", name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Item parameters, as in dev/null_q3_calibration.R ---
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

# --- Max-T variants ---------------------------------------------------------
# obs, rep: draws x pairs matrices of observed and replicated Q3.
maxt_p <- function(obs, rep, reference = c("pair", "global"),
                   scale = c("std", "raw")) {
  reference <- match.arg(reference)
  scale <- match.arg(scale)
  if (reference == "global") {
    obs <- obs - rowMeans(obs, na.rm = TRUE)
    rep <- rep - rowMeans(rep, na.rm = TRUE)
    z_obs <- colMeans(obs, na.rm = TRUE)
    z_rep <- rep
  } else {
    m <- colMeans(rep, na.rm = TRUE)
    z_obs <- colMeans(obs, na.rm = TRUE) - m
    z_rep <- sweep(rep, 2L, m)
  }
  if (scale == "std") {
    s <- apply(rep, 2L, stats::sd, na.rm = TRUE)
    z_obs <- z_obs / s
    z_rep <- sweep(z_rep, 2L, s, "/")
  }
  t_max <- apply(z_rep, 1L, max, na.rm = TRUE)
  vapply(z_obs, function(z) mean(t_max >= z), numeric(1))
}

VARIANTS <- expand.grid(reference = c("pair", "global"), scale = c("std", "raw"),
                        stringsAsFactors = FALSE)

# q3_statistic() output -> draws x pairs matrices, then all variants
q3_variants <- function(q3, method) {
  q3 <- as.data.frame(q3)
  pairs <- unique(q3$item_pair)
  to_mat <- function(col) {
    sapply(pairs, function(p) {
      d <- q3[q3$item_pair == p, ]
      d[[col]][order(d$draw)]
    })
  }
  obs <- to_mat("q3")
  rep <- to_mat("q3_rep")
  out <- data.frame(item_pair = pairs, method = method,
                    q3_obs = colMeans(obs, na.rm = TRUE),
                    q3_rep = colMeans(rep, na.rm = TRUE),
                    sd_rep = apply(rep, 2L, stats::sd, na.rm = TRUE),
                    ppp = colMeans(obs > rep, na.rm = TRUE) +
                      0.5 * colMeans(obs == rep, na.rm = TRUE))
  for (v in seq_len(nrow(VARIANTS))) {
    nm <- paste0("p_", VARIANTS$reference[v], "_", VARIANTS$scale[v])
    out[[nm]] <- maxt_p(obs, rep, VARIANTS$reference[v], VARIANTS$scale[v])
  }
  rownames(out) <- NULL
  out
}

# Sanity check: pair_std reproduces q3_post()'s q3_p_adj (to rounding)
check_against_package <- function(q3, res) {
  pkg <- q3_post(q3)$summary
  m <- merge(res[, c("item_pair", "p_pair_std")],
             as.data.frame(pkg)[, c("item_pair", "q3_p_adj")], by = "item_pair")
  max(abs(m$p_pair_std - m$q3_p_adj))
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
  seed <- 1e5 * match(cond_row$cond, conditions$cond) + rep + 8e6
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
  check <- c()
  for (m in c("conditional", "theta")) {
    set.seed(seed + if (m == "conditional") 1 else 2)
    q3 <- q3_statistic(fit, ndraws_use = NDRAWS, method = m)
    res[[m]] <- q3_variants(q3, m)
    check[m] <- check_against_package(q3, res[[m]])
  }
  res <- do.call(rbind, res)
  res <- cbind(data.frame(cond = cond_row$cond, model_type = cond_row$model_type,
                          n = cond_row$n, rep = rep),
               res,
               data.frame(rhat_max = rhat_max, n_div = n_div, fit_sec = fit_sec,
                          check_cond = check[["conditional"]],
                          check_theta = check[["theta"]]))
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
