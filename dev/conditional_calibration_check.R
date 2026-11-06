# Run from the package root: N=500 REPS=200 Rscript dev/conditional_calibration_check.R (about 3 min on 4 cores).
# Exploratory check (container, smoke scale): PAV calibration of cumulative
# events P(X_i >= k) using conditional probabilities given the total score
# versus theta-based in-sample probabilities (ML theta), both with CML item
# parameters. Bands from replicates: conditional sampling given R for the
# conditional version, independent sampling given theta-hat for the other.
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(parallel)

K <- 9; C <- 4; N <- as.integer(Sys.getenv("N", "500")); REPS <- as.integer(Sys.getenv("REPS", "100")); B <- 200
loc <- seq(-1, 1.5, length.out = K)
tau <- t(sapply(loc, function(l) l + c(-1, 0, 1)))   # K x 3

sim_pcm <- function(n, a) {
  th <- rnorm(n, -0.3, 1.3)
  x <- matrix(0L, n, K)
  for (i in 1:K) {
    eta <- sapply(0:(C - 1), function(h) a[i] * (h * th - c(0, cumsum(tau[i, ]))[h + 1]))
    p <- exp(eta - apply(eta, 1, max)); p <- p / rowSums(p)
    x[, i] <- apply(p, 1, function(pp) sample.int(C, 1, prob = pp)) - 1L
  }
  x
}

# weighted PAV on tied groups, returns fitted value per group
pav_w <- function(y, w) {
  n <- length(y); val <- y; wt <- w; len <- rep(1L, n); m <- 0L
  v <- numeric(n); ww <- numeric(n); ll <- integer(n)
  for (j in 1:n) {
    m <- m + 1L; v[m] <- val[j]; ww[m] <- wt[j]; ll[m] <- 1L
    while (m > 1L && v[m - 1L] > v[m]) {
      nw <- ww[m - 1L] + ww[m]
      v[m - 1L] <- (v[m - 1L] * ww[m - 1L] + v[m] * ww[m]) / nw
      ww[m - 1L] <- nw; ll[m - 1L] <- ll[m - 1L] + ll[m]; m <- m - 1L
    }
  }
  rep(v[1:m], ll[1:m])
}

# Calibration of one event set: p (vector), y (0/1 matrix columns = obs, rows = replicates or 1)
# grouped by distinct p. Returns observed curve, band, weights.
cal_curve <- function(p, y_obs, y_rep) {
  g <- match(p, sort(unique(p))); pu <- sort(unique(p)); w <- tabulate(g, length(pu))
  agg <- function(y) as.numeric(rowsum(as.numeric(y), g, reorder = TRUE)) / w
  obs <- pav_w(agg(y_obs), w)
  rep <- sapply(seq_len(nrow(y_rep)), function(b) pav_w(agg(y_rep[b, ]), w))
  lo <- apply(rep, 1, quantile, 0.025); hi <- apply(rep, 1, quantile, 0.975)
  list(p = pu, w = w, obs = obs, lo = lo, hi = hi, out = obs < lo | obs > hi)
}

one_rep <- function(seed, a) {
  set.seed(seed)
  x <- sim_pcm(N, a)
  colnames(x) <- paste0("I", 1:K)
  fit <- eRm::PCM(x, se = FALSE, sum0 = TRUE)
  beta <- matrix(c(rbind(0, matrix(fit$betapar, nrow = C - 1))), nrow = K, byrow = TRUE) # K x C, beta_i0 = 0
  r <- rowSums(x); keep <- r > 0 & r < K * (C - 1)
  x <- x[keep, ]; r <- r[keep]; n <- nrow(x)
  # --- conditional probabilities via the package engine
  person_idx <- rep(seq_len(n), times = K); item_idx <- rep(1:K, each = n)
  prob <- exp(beta[item_idx, ]); prob <- prob / rowSums(prob)
  st <- .cond_draw_setup(prob, person_idx, item_idx, x)
  gr <- st$groups[[1]]; ri <- match(r, gr$r_vals)
  pc <- lapply(1:K, function(i) gr$probs[[i]][ri, , drop = FALSE])   # n x C each
  yrep_c <- replicate(B, .cond_sample(st$W, st$suffix, st$key, st$r))  # n x K x B
  # --- theta-based: ML theta per raw score with CML betas
  th_r <- sapply(sort(unique(r)), function(s) uniroot(function(t) {
    e <- sapply(1:K, function(i) { q <- exp(beta[i, ] + (0:(C - 1)) * t); sum((0:(C - 1)) * q) / sum(q) })
    sum(e) - s }, c(-10, 10))$root)
  th <- th_r[match(r, sort(unique(r)))]
  pt <- lapply(1:K, function(i) { q <- exp(outer(th, 0:(C - 1)) + rep(beta[i, ], each = n)); q / rowSums(q) })
  yrep_t <- array(0L, c(n, K, B))
  for (i in 1:K) { cp <- t(apply(pt[[i]], 1, cumsum)); u <- matrix(runif(n * B), n, B)
    yrep_t[, i, ] <- (u > cp[, 1]) + (u > cp[, 2]) + (u > cp[, 3]) }
  res <- list()
  for (meth in c("cond", "theta")) {
    P <- if (meth == "cond") pc else pt; YR <- if (meth == "cond") yrep_c else yrep_t
    for (i in 1:K) for (k in 1:(C - 1)) {
      p <- round(rowSums(P[[i]][, (k + 1):C, drop = FALSE]), 10)
      cc <- cal_curve(p, matrix(as.numeric(x[, i] >= k), 1), t(YR[, i, ] >= k) * 1)
      res[[length(res) + 1]] <- data.frame(seed, meth, item = i, event = k,
        exceed = sum(cc$w * cc$out) / sum(cc$w),
        dev = sum(cc$w * (cc$obs - cc$p)) / sum(cc$w),  # mean signed deviation
        slope = coef(lm(cc$obs ~ cc$p, weights = cc$w))[2])
    }
  }
  do.call(rbind, res)
}

conds <- list(null = rep(1, K), low = c(0.5, rep(1, K - 1)), high = c(2, rep(1, K - 1)))
out <- list()
for (cn in names(conds)) {
  t0 <- Sys.time()
  rr <- mclapply(seq_len(REPS), function(s) tryCatch(one_rep(1000 * match(cn, names(conds)) + s, conds[[cn]]),
                 error = function(e) NULL), mc.cores = 4)
  out[[cn]] <- cbind(cond = cn, do.call(rbind, rr))
  message(cn, " done ", format(Sys.time() - t0), " (", sum(!sapply(rr, is.null)), " ok)")
}
res <- do.call(rbind, out)
saveRDS(res, paste0("dev/results/conditional_calibration/res_n", N, ".rds")); saveRDS(res, paste0("dev/results/conditional_calibration/res_n", N, "_", Sys.Date(), ".rds"))
agg <- aggregate(cbind(exceed, dev, slope) ~ cond + meth + I(ifelse(item == 1, "item1", "others")), data = res, FUN = mean)
names(agg)[3] <- "which"
print(agg[order(agg$cond, agg$which, agg$meth), ], digits = 3)
