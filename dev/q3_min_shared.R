# Minimum number of shared persons for Q3 cell pairs in many-facet models
# (finding 1 of dev/code-review/code-review-2026-11-16.md). Small study, 2026-11-16.
#
# Question. q3_statistic(facet_var = ) keeps every cell pair that at least
# one person shares. With few shared persons Q3 is close to +/-1 in every
# draw. Does that (a) make q3_ppp flag such pairs too often, and (b) lower
# the power of the max-T q3_p_adj for the other pairs, by inflating the
# reference distribution? Which minimum number of shared persons, applied
# to the max-T family only, keeps the familywise rate and restores power?
#
# Shortcut. No brms fits. Each replicate uses the generating parameters as
# the "posterior": every draw has the same item, rater and threshold
# parameters, and the replicated data come from the conditional engine
# (.cond_q3_draw(), the same function q3_statistic() uses). Observed Q3 is
# then the same in every draw and only the replicated Q3 varies. This
# isolates the small-sample behaviour of Q3, but leaves out posterior
# uncertainty, which adds variation to the observed Q3 in real fits.
#
# Design (fixed across replicates). 8 raters, 4 items, rating scale with 4
# categories (acat), 32 cells. Each person is rated by 2 raters on all 4
# items. The 28 rater pairs get very different numbers of persons (2 to
# 30, `pair_n` below), so cross-rater cell pairs (same item or neither)
# share 2 to 30 persons, and same-rater pairs share 50-100.
#   cond   null       no halo
#          halo06     halo SD 0.6 for R1 only, not in the model
#          halo10     halo SD 1.0 for R1 only, not in the model
#   reps   200 null, 100 per halo condition, 400 draws each
#
# Outcomes (per minimum m, pairs with fewer than m shared persons left out
# of the max-T family, their q3_ppp kept):
#   fwer      null: any pair in the family with q3_p_adj < .05
#   power     halo: any of R1's 6 same-rater pairs with q3_p_adj < .05
#   per-pair calibration by shared-person bin (null): rate of mid-p < .05
#             and of q3_ppp > .95 (as q3_post() computes it)
#   tmax_small share of draws whose max-T maximum comes from a pair with
#             fewer than 10 shared persons (the mechanism)
#
# Run from the package root on the host (10 workers, about 5-10 min):
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 nohup Rscript dev/q3_min_shared.R \
#     > dev/results/q3_min_shared.log 2>&1 &
# Smoke test (container): SMOKE=1 Rscript dev/q3_min_shared.R
# Summary: Rscript dev/q3_min_shared_summary.R (SMOKE=1 for the smoke run)

suppressMessages(devtools::load_all(".", quiet = TRUE))
smoke <- Sys.getenv("SMOKE") == "1"
N_WORKERS <- if (smoke) 2L else 10L   # max 10 cores on the host
N_DRAWS <- if (smoke) 50L else 400L
REPS <- if (smoke) c(null = 2L, halo06 = 2L, halo10 = 2L) else
  c(null = 200L, halo06 = 100L, halo10 = 100L)
HALO_SD <- c(null = 0, halo06 = 0.6, halo10 = 1.0)
MINS <- c(0L, 3L, 5L, 8L, 10L, 15L, 20L)
name <- if (smoke) "q3_min_shared_smoke" else "q3_min_shared"
dir.create(file.path("dev", "results"), showWarnings = FALSE)

# --- Fixed design ---
n_rater <- 8L
n_item <- 4L
raters <- paste0("R", seq_len(n_rater))
rater_pairs <- t(utils::combn(n_rater, 2L))          # 28 pairs
pair_n <- c(2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10, 11, 12, 13,
            14, 15, 16, 18, 20, 22, 24, 26, 30)
set.seed(20261116)
pair_n <- sample(pair_n)                              # which raters get few
person_raters <- rater_pairs[rep(seq_along(pair_n), pair_n), , drop = FALSE]
n_person <- nrow(person_raters)                       # 300
item_loc <- c(-0.6, -0.2, 0.2, 0.6)
severity <- seq(-0.7, 0.7, length.out = n_rater)
tau <- c(-1.5, 0, 1.5)
n_cat <- length(tau) + 1L

obs_grid <- do.call(rbind, lapply(seq_len(n_person), function(p) {
  expand.grid(person = p, item = seq_len(n_item),
              rater = person_raters[p, ])
}))
obs_grid$cell <- (obs_grid$rater - 1L) * n_item + obs_grid$item
k <- n_rater * n_item
cell_rater <- rep(seq_len(n_rater), each = n_item)
cell_item <- rep(seq_len(n_item), times = n_rater)

# Pairs in the order .cond_q3_draw() returns them (upper triangle, row-major)
pg <- which(upper.tri(diag(k)), arr.ind = TRUE)
pg <- pg[order(pg[, "row"], pg[, "col"]), , drop = FALSE]
pair_type <- ifelse(cell_rater[pg[, 1]] == cell_rater[pg[, 2]], "same rater",
                    ifelse(cell_item[pg[, 1]] == cell_item[pg[, 2]],
                           "same item", "neither"))
answered <- matrix(FALSE, n_person, k)
answered[cbind(obs_grid$person, obs_grid$cell)] <- TRUE
shared_all <- crossprod(answered)
n_shared <- shared_all[cbind(pg[, 1], pg[, 2])]
keep <- n_shared > 0                                  # as q3_statistic()
target <- pair_type == "same rater" & cell_rater[pg[, 1]] == 1L

pcm_probs <- function(eta) {
  # eta: vector, returns length(eta) x n_cat matrix of acat probabilities
  num <- cbind(0, sapply(seq_along(tau), function(j) j * eta - sum(tau[1:j])))
  num <- exp(num - apply(num, 1L, max))
  num / rowSums(num)
}

# One-sided upper mid-p over draws: P(rep > obs) + 0.5 P(rep == obs)
midp_upper <- function(obs, rep) {
  mean(rep > obs, na.rm = TRUE) + 0.5 * mean(rep == obs, na.rm = TRUE)
}

one_rep <- function(cond, rep_id) {
  set.seed(1000L * match(cond, names(REPS)) + rep_id)
  theta <- rnorm(n_person)
  halo <- matrix(rnorm(n_person * n_rater, 0, HALO_SD[[cond]]), n_person)
  halo[, -1] <- 0                                     # R1 only
  eta_model <- theta[obs_grid$person] - item_loc[obs_grid$item] -
    severity[obs_grid$rater]
  eta_true <- eta_model + halo[cbind(obs_grid$person, obs_grid$rater)]
  p_true <- pcm_probs(eta_true)
  y <- apply(p_true, 1L, function(p) sample.int(n_cat, 1L, prob = p)) - 1L
  x_obs <- matrix(NA_real_, n_person, k)
  x_obs[cbind(obs_grid$person, obs_grid$cell)] <- y
  prob <- pcm_probs(eta_model)                        # model without halo

  t0 <- proc.time()[["elapsed"]]
  draws <- lapply(seq_len(N_DRAWS), function(d) {
    .cond_q3_draw(prob, obs_grid$person, obs_grid$cell, x_obs)
  })
  obs <- draws[[1]]$obs[keep]                         # same in every draw
  rep_mat <- do.call(rbind, lapply(draws, function(z) z$rep[keep]))
  sec <- proc.time()[["elapsed"]] - t0

  ns <- n_shared[keep]
  pair <- data.frame(
    cond = cond, rep = rep_id, pair = which(keep), n_shared = ns,
    pair_type = pair_type[keep], target = target[keep],
    q3_obs = obs, q3_rep = colMeans(rep_mat, na.rm = TRUE),
    midp = vapply(seq_along(obs), function(j) midp_upper(obs[j], rep_mat[, j]),
                  numeric(1)),
    # As q3_post(): P(observed > replicated), no mid-p
    q3_ppp = colMeans(sweep(rep_mat, 2L, obs, "<"), na.rm = TRUE)
  )
  fam <- lapply(MINS, function(m) {
    inc <- ns >= m
    padj <- .ppp_maxt_upper(obs[inc], rep_mat[, inc, drop = FALSE])
    # Where the max-T maximum comes from, per draw
    rm_ <- rep_mat[, inc, drop = FALSE]
    mu <- colMeans(rm_, na.rm = TRUE)
    s <- apply(rm_, 2L, stats::sd, na.rm = TRUE)
    ok <- is.finite(s) & s > 0
    z <- sweep(sweep(rm_[, ok, drop = FALSE], 2L, mu[ok]), 2L, s[ok], "/")
    z[is.na(z)] <- -Inf
    arg <- max.col(z, ties.method = "first")
    data.frame(
      cond = cond, rep = rep_id, min_shared = m, n_family = sum(inc),
      any_flag = any(padj < .05, na.rm = TRUE),
      target_flag = any(padj[target[keep][inc]] < .05, na.rm = TRUE),
      target_min_padj = min(padj[target[keep][inc]], na.rm = TRUE),
      tmax_q95 = stats::quantile(apply(z, 1L, max), .95),
      tmax_small = mean(ns[inc][ok][arg] < 10)
    )
  })
  list(pair = pair, family = do.call(rbind, fam),
       info = data.frame(cond = cond, rep = rep_id, sec = sec))
}

jobs <- do.call(rbind, lapply(names(REPS), function(cn) {
  data.frame(cond = cn, rep = seq_len(REPS[[cn]]))
}))
jobs <- jobs[order(jobs$rep, jobs$cond), ]            # interleave conditions
cat("Design:", n_person, "persons,", k, "cells,", sum(keep), "pairs with",
    "shared persons (", sum(n_shared[keep] < 10), "with fewer than 10).",
    nrow(jobs), "jobs,", N_DRAWS, "draws each.\n")
t_start <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  out <- tryCatch(one_rep(jobs$cond[j], jobs$rep[j]),
                  error = function(e) conditionMessage(e))
  if (j %% 20 == 0 || smoke) {
    cat(format(Sys.time(), "%H:%M:%S"), "job", j, "of", nrow(jobs), "\n")
  }
  out
}, mc.cores = N_WORKERS, mc.preschedule = FALSE)
failed <- !vapply(res, is.list, logical(1))
if (any(failed)) {
  cat(sum(failed), "jobs failed, first error:", res[[which(failed)[1]]], "\n")
}
res <- res[!failed]
all_res <- list(
  pair = do.call(rbind, lapply(res, `[[`, "pair")),
  family = do.call(rbind, lapply(res, `[[`, "family")),
  info = do.call(rbind, lapply(res, `[[`, "info")),
  design = list(pair_n = pair_n, rater_pairs = rater_pairs,
                n_person = n_person, mins = MINS, n_draws = N_DRAWS,
                halo_sd = HALO_SD)
)
saveRDS(all_res, file.path("dev", "results", paste0(name, ".rds")))
saveRDS(all_res, file.path("dev", "results",
                           paste0(name, "_", format(Sys.Date(), "%Y%m%d"),
                                  ".rds")))
cat("Done in", format(round(Sys.time() - t_start, 1)), "\n")
