# Internal machinery for checks that condition on the total score.
#
# Under a Rasch-family model (dichotomous Rasch, partial credit, rating
# scale) the total score is sufficient for the person location. Given the
# total score, the distribution of a response pattern depends only on the
# item parameters. The functions below work from the category
# probabilities that brms returns for one posterior draw:
#
# 1. `.cond_weight_groups()` removes the person location from those
#    probabilities and groups persons whose remaining item weights are the
#    same (the same answered items and the same item parameters, which
#    differ only under DIF terms or missing responses).
# 2. `.cond_esf()` computes elementary symmetric functions (ESFs) for one
#    group, as suffix and prefix tables.
# 3. `.cond_expected_probs()` gives P(X_i = x | total score) for each item.
# 4. `.cond_sample()` draws response patterns conditional on each person's
#    total score.
#
# Item weights are scaled so that each item's largest weight is 1, and each
# ESF table is rescaled by its maximum after every step. All probabilities
# are normalised over the categories of one item at a time, so these
# scalings cancel.

# Check that a brmsfit is a model in which the total score is sufficient.
# Takes plain inputs so it can be tested without a fitted model.
.check_sufficient_score_model <- function(family, link, dpar_names,
                                          disc_values, sd_names,
                                          person_name, fun_name) {
  if (!family %in% c("bernoulli", "acat") || !identical(link, "logit")) {
    stop(
      "`", fun_name, "()` conditions on the total score, which requires ",
      "a Rasch-family model: `bernoulli()` or `acat` with the logit link. ",
      "The model uses family '", family, "' with link '", link, "'.",
      call. = FALSE
    )
  }
  if ("disc" %in% dpar_names ||
      (length(disc_values) > 0L && any(abs(disc_values - 1) > 1e-8))) {
    stop(
      "`", fun_name, "()` requires equal item discrimination, but the ",
      "model estimates `disc`. The total score is then not sufficient.",
      call. = FALSE
    )
  }
  prefix <- paste0("sd_", person_name, "__")
  person_sd <- sd_names[startsWith(sd_names, prefix)]
  coefs <- substring(person_sd, nchar(prefix) + 1L)
  if (length(coefs) == 0L) {
    stop(
      "`", fun_name, "()` found no person random intercept `(1 | ",
      person_name, ")` in the model.",
      call. = FALSE
    )
  }
  if (!all(grepl("^Intercept(:|$)", coefs))) {
    stop(
      "`", fun_name, "()` requires a single person location, `(1 | ",
      person_name, ")`. The model has person-level terms '",
      paste(coefs[!grepl("^Intercept(:|$)", coefs)], collapse = "', '"),
      "', for example a multidimensional model, so the total score is ",
      "not sufficient.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Run the sufficiency check on a brmsfit.
.check_sufficient_score_brmsfit <- function(model, person_name, fun_name) {
  fam <- stats::family(model)
  vars <- brms::variables(model)
  disc_values <- numeric(0)
  if ("disc" %in% vars) {
    disc_values <- as.numeric(
      brms::as_draws_matrix(model, variable = "disc")
    )
  }
  .check_sufficient_score_model(
    family      = fam$family,
    link        = fam$link,
    dpar_names  = names(model$formula$pforms),
    disc_values = disc_values,
    sd_names    = vars[startsWith(vars, "sd_")],
    person_name = person_name,
    fun_name    = fun_name
  )
}

# Category probabilities for one draw as an n_obs x C matrix, scores 0..C-1.
.cond_prob_matrix <- function(epred, s) {
  if (length(dim(epred)) == 3L) {
    epred[s, , , drop = TRUE]
  } else {
    p <- epred[s, ]
    cbind(1 - p, p)
  }
}

# Remove the person location and group persons by their item weights.
#
# prob: n_obs x C category probabilities for one draw.
# person_idx, item_idx: integer indices (1..n_persons, 1..k) per observation.
# Returns list(key = group per person, items = list per group of the
# answered items in item order, W = list per group of n_answered x C
# matrices of scaled item weights for those items).
#
# Works on the observations only, so memory grows with the number of
# responses, not with n_persons x k. Before 0.4.0.9001 the function built
# n_persons x k x C arrays, which in many-facet designs (k = items x
# raters) were large and mostly empty. Groups, their order and the
# weights are the same as before, up to the rounding of the hash below.
.cond_weight_groups <- function(prob, person_idx, item_idx, n_persons, k) {
  n_cat <- ncol(prob)
  # Adjacent-category structure: log P(x) - log P(0) = x * theta - sum(tau).
  lo <- log(prob) - log(prob[, 1L])
  # Remove theta using the log-odds of score 1 on the person's first
  # answered item (in item order). The result depends on item parameters
  # only: x * tau_ref1 - sum_{j <= x} tau_ij.
  ord <- order(person_idx, item_idx)
  first_obs <- ord[!duplicated(person_idx[ord])]
  ref <- numeric(n_persons)
  ref[person_idx[first_obs]] <- lo[first_obs, 2L]
  lo <- lo - outer(ref[person_idx], 0:(n_cat - 1L))

  # Hash each person's weights to one number: the sum over items and
  # categories of weight x hash_w, with sentinels 1e4 for missing items,
  # -1e4 for impossible categories and 2e4 for answered items with
  # undefined weights (prob[, 1] = 0), so that those are not grouped with
  # persons who did not answer the item. Written
  # as the all-missing value plus each observation's difference from it,
  # summed in item order, so persons with the same weights get the same
  # number. Rounding absorbs floating point noise, a near-boundary case
  # only splits a group in two.
  hash_w <- matrix(sin(seq_len(k * n_cat) * 12.9898) * 43.758,
                   nrow = k, ncol = n_cat)
  lh <- lo
  lh[is.na(lh)] <- 1e4
  lh[is.infinite(lh)] <- -1e4
  # Rows handled as undefined in the weights loop below
  lh[is.na(lo[, 1L]), ] <- 2e4
  contrib <- rowSums((lh - 1e4) * hash_w[item_idx, , drop = FALSE])
  h <- rep(1e4 * sum(hash_w), n_persons)
  ps <- rowsum(contrib[ord], person_idx[ord], reorder = TRUE)
  h[as.integer(rownames(ps))] <- h[as.integer(rownames(ps))] + ps[, 1L]
  h <- round(h, 6)
  uh <- unique(h)
  key <- match(h, uh)
  rep_person <- match(seq_along(uh), key)

  # Weights of each group, from its first person's observations
  obs_by_person <- split(ord, factor(person_idx[ord], levels = seq_len(n_persons)))
  items <- vector("list", length(uh))
  W <- vector("list", length(uh))
  for (u in seq_along(uh)) {
    rows <- obs_by_person[[rep_person[u]]]
    Lu <- matrix(lo[rows, ], nrow = length(rows), ncol = n_cat)
    undefined <- is.na(Lu[, 1L])
    Lu[undefined, ] <- -Inf
    Lu[undefined, 1L] <- 0
    items[[u]] <- item_idx[rows]
    W[[u]] <- exp(Lu - apply(Lu, 1L, max))
  }
  list(key = key, items = items, W = W)
}

# Convolve two score distributions, truncated to length(out_len).
.cond_convolve <- function(a, b, out_len) {
  out <- numeric(out_len)
  nb <- length(b)
  for (j in which(a > 0)) {
    hi <- min(out_len, j + nb - 1L)
    out[j:hi] <- out[j:hi] + a[j] * b[seq_len(hi - j + 1L)]
  }
  out
}

# Suffix and prefix ESF tables for one group.
# W: k x C matrix of item weights. Returns list(suffix, prefix), each a
# (k + 1) x (M + 1) matrix with M = k * (C - 1). suffix[i, ] is the score
# distribution of items i..k, suffix[k + 1, ] that of no items. prefix[i, ]
# covers items 1..(i - 1). Rows are scaled to a maximum of 1.
.cond_esf <- function(W) {
  k <- nrow(W)
  n_cat <- ncol(W)
  M <- k * (n_cat - 1L)
  e0 <- c(1, numeric(M))
  suffix <- prefix <- matrix(0, nrow = k + 1L, ncol = M + 1L)
  suffix[k + 1L, ] <- e0
  prefix[1L, ] <- e0
  for (i in seq.int(k, 1L)) {
    s <- .cond_convolve(W[i, ], suffix[i + 1L, ], M + 1L)
    suffix[i, ] <- s / max(s)
  }
  for (i in seq_len(k)) {
    p <- .cond_convolve(W[i, ], prefix[i, ], M + 1L)
    prefix[i + 1L, ] <- p / max(p)
  }
  list(suffix = suffix, prefix = prefix)
}

# P(X_i = x | total score r) for every item, for the totals in r.
# Returns a k-length list of length(r) x C matrices.
.cond_expected_probs <- function(W, esf, r) {
  k <- nrow(W)
  n_cat <- ncol(W)
  M1 <- ncol(esf$suffix)
  # Score distribution of all items but i, G[i, t + 1], the convolution of
  # prefix[i, ] and suffix[i + 1, ]. Only t = r - x is needed, so G is
  # filled at those t, one inner product per t for all items at once.
  # With a group per person (a person x facet term) that is a few columns
  # instead of k full convolutions.
  t_need <- unique(as.vector(outer(r, seq_len(n_cat) - 1L, "-")))
  t_need <- t_need[t_need >= 0L & t_need < M1]
  pre <- esf$prefix[seq_len(k), , drop = FALSE]
  suf <- esf$suffix[seq_len(k) + 1L, , drop = FALSE]
  G <- matrix(0, nrow = k, ncol = M1)
  for (t in t_need) {
    j <- seq_len(t + 1L)
    G[, t + 1L] <- rowSums(pre[, j, drop = FALSE] *
                             suf[, t + 2L - j, drop = FALSE])
  }
  out <- vector("list", k)
  for (i in seq_len(k)) {
    g_rest <- G[i, ]
    pm <- matrix(0, nrow = length(r), ncol = n_cat)
    for (x in seq_len(n_cat)) {
      t <- r - (x - 1L)
      ok <- t >= 0L & t < M1
      pm[ok, x] <- W[i, x] * g_rest[t[ok] + 1L]
    }
    tot <- rowSums(pm)
    if (any(!is.finite(tot) | tot <= 0)) {
      stop("Conditional probabilities underflowed. The test may have ",
           "too many items or too extreme thresholds for this method.",
           call. = FALSE)
    }
    out[[i]] <- pm / tot
  }
  out
}

# Draw response patterns conditional on each person's total score.
# W: U x k x C weights, suffix: U x (k + 1) x (M + 1) suffix tables,
# key: group per person, r: total score per person.
# Returns an n_persons x k matrix of scores 0..C-1.
.cond_sample <- function(W, suffix, key, r) {
  n <- length(key)
  k <- dim(W)[2L]
  n_cat <- dim(W)[3L]
  M1 <- dim(suffix)[3L]
  out <- matrix(0L, nrow = n, ncol = k)
  remaining <- r
  for (i in seq_len(k)) {
    pm <- matrix(0, nrow = n, ncol = n_cat)
    for (x in seq_len(n_cat)) {
      t <- remaining - (x - 1L)
      ok <- t >= 0L & t < M1
      # cbind() drops a zero-length column, so an empty selection would
      # become a one-row index matrix
      if (!any(ok)) next
      pm[ok, x] <- W[cbind(key[ok], i, x)] *
        suffix[cbind(key[ok], i + 1L, t[ok] + 1L)]
    }
    tot <- rowSums(pm)
    if (any(!is.finite(tot) | tot <= 0)) {
      stop("Conditional sampling failed: a total score has zero ",
           "probability under the model weights.", call. = FALSE)
    }
    u <- stats::runif(n) * tot
    cum <- pm[, 1L]
    xi <- integer(n)
    for (x in seq_len(n_cat - 1L)) {
      xi <- xi + (u > cum)
      cum <- cum + pm[, x + 1L]
    }
    out[, i] <- xi
    remaining <- remaining - xi
  }
  out
}

# Conditional quantities shared by the statistics for one draw.
#
# prob: n_obs x C category probabilities for one draw.
# x_obs: n_persons x k observed scores (0-based, NA when missing).
#
# The ESF tables of a group cover only the items its members answered. In a
# sparse design (many-facet models, where the items are item x facet cells,
# or planned missing data) this keeps the tables at the size of one
# person's answered set instead of all k items. The answered items of group
# u sit at positions 1..n_answered of the compact arrays, in item order,
# and pos_item[u, j] gives the item column at position j (NA for padding).
# Padding positions get weights (1, 0, ..., 0), which force a score of 0.
#
# Returns list(r, key, W, suffix, pos_item, k, groups). W (U x ka x C) and
# suffix (U x (ka + 1) x (M + 1)) are the compact arrays, with ka the
# largest number of answered items in any group. Each element of groups
# describes one set of persons with the same item weights: members, the
# distinct totals r_vals and their counts, P(X_i = x | r) for those totals
# (a k-length list, NULL for unanswered items), and which items the group
# answered.
.cond_draw_setup <- function(prob, person_idx, item_idx, x_obs) {
  n_persons <- nrow(x_obs)
  k <- ncol(x_obs)
  r <- rowSums(x_obs, na.rm = TRUE)
  cw <- .cond_weight_groups(prob, person_idx, item_idx, n_persons, k)
  U <- length(cw$W)
  n_cat <- ncol(prob)
  rep_person <- match(seq_len(U), cw$key)
  answered <- !is.na(x_obs[rep_person, , drop = FALSE])
  ka <- max(rowSums(answered))
  M1 <- ka * (n_cat - 1L) + 1L

  pos_item <- matrix(NA_integer_, nrow = U, ncol = ka)
  W_c <- array(0, dim = c(U, ka, n_cat))
  W_c[, , 1L] <- 1
  suffix_all <- array(0, dim = c(U, ka + 1L, M1))
  suffix_all[, , 1L] <- 1
  groups <- vector("list", U)
  for (u in seq_len(U)) {
    items <- which(answered[u, ])
    na <- length(items)
    if (!identical(as.integer(cw$items[[u]]), as.integer(items))) {
      stop("Internal error: the item weights do not match the observed ",
           "responses.", call. = FALSE)
    }
    Wu <- cw$W[[u]]
    pos_item[u, seq_len(na)] <- items
    W_c[u, seq_len(na), ] <- Wu
    esf <- .cond_esf(Wu)
    suffix_all[u, seq_len(na + 1L), seq_len(ncol(esf$suffix))] <- esf$suffix
    members <- which(cw$key == u)
    r_tab <- table(r[members])
    r_vals <- as.integer(names(r_tab))
    probs <- vector("list", k)
    probs[items] <- .cond_expected_probs(Wu, esf, r_vals)
    groups[[u]] <- list(
      members  = members,
      r_vals   = r_vals,
      r_counts = as.numeric(r_tab),
      probs    = probs,
      answered = answered[u, ]
    )
  }
  list(r = r, key = cw$key, W = W_c, suffix = suffix_all,
       pos_item = pos_item, k = k, groups = groups)
}

# Replicated responses for one draw, conditional on each person's total.
# Samples over the compact positions of .cond_draw_setup() and maps them
# back to item columns. Returns an n_persons x k matrix, 0 for unanswered
# items.
.cond_replicate <- function(st) {
  x_pos <- .cond_sample(st$W, st$suffix, st$key, st$r)
  items <- st$pos_item[st$key, , drop = FALSE]
  ok <- !is.na(items)
  out <- matrix(0L, nrow = length(st$key), ncol = st$k)
  out[cbind(row(items)[ok], items[ok])] <- x_pos[ok]
  out
}

# Expected item x restscore table and replicated scores for one draw.
# Returns list(gamma_exp = k, gamma_rep = k).
.cond_restscore_draw <- function(prob, person_idx, item_idx, x_obs) {
  k <- ncol(x_obs)
  st <- .cond_draw_setup(prob, person_idx, item_idx, x_obs)
  r <- st$r
  n_cat <- dim(st$W)[3L]
  M1 <- dim(st$suffix)[3L]

  exp_tab <- lapply(seq_len(k), function(i) matrix(0, nrow = M1, ncol = n_cat))
  for (g in st$groups) {
    for (i in which(g$answered)) {
      for (x in seq_len(n_cat)) {
        rest <- g$r_vals - (x - 1L)
        ok <- rest >= 0L
        # cbind() drops a zero-length column, so an empty selection would
        # become a one-row index matrix
        if (!any(ok)) next
        idx <- cbind(rest[ok] + 1L, x)
        exp_tab[[i]][idx] <- exp_tab[[i]][idx] +
          g$r_counts[ok] * g$probs[[i]][ok, x]
      }
    }
  }

  x_rep <- .cond_replicate(st)
  gamma_exp <- gamma_rep <- numeric(k)
  for (i in seq_len(k)) {
    gamma_exp[i] <- gk_gamma(exp_tab[[i]])
    valid <- !is.na(x_obs[, i])
    gamma_rep[i] <- gk_gamma_vec(x_rep[valid, i], r[valid] - x_rep[valid, i])
  }
  list(gamma_exp = gamma_exp, gamma_rep = gamma_rep)
}

# Conditional expected value and variance of every response given the
# person's total score. Returns list(E, V), n_persons x k, NA when missing.
.cond_moments <- function(st, k) {
  n_persons <- length(st$r)
  n_cat <- dim(st$W)[3L]
  scores <- 0:(n_cat - 1L)
  E <- V <- matrix(NA_real_, nrow = n_persons, ncol = k)
  for (g in st$groups) {
    idx <- match(st$r[g$members], g$r_vals)
    for (i in which(g$answered)) {
      e <- as.vector(g$probs[[i]] %*% scores)
      v <- as.vector(g$probs[[i]] %*% scores^2) - e^2
      E[g$members, i] <- e[idx]
      V[g$members, i] <- pmax(v[idx], 0)
    }
  }
  list(E = E, V = V)
}

# Conditional infit (and outfit) for one draw, observed and replicated.
#
# Expectations and variances are conditional on the total score, as in
# Christensen et al. (2013) and iarm::out_infit(). Persons with the minimum
# or maximum total have zero variance on every item and drop out.
# group: optional integer vector of length k mapping items (columns) to
# groups, for many-facet models where the items are item x facet cells.
# Sums of squared residuals and variances are then pooled within groups.
# Returns list(infit, infit_rep, outfit, outfit_rep), each of length k, or
# of length max(group).
.cond_infit_draw <- function(prob, person_idx, item_idx, x_obs,
                             group = NULL) {
  k <- ncol(x_obs)
  st <- .cond_draw_setup(prob, person_idx, item_idx, x_obs)
  mo <- .cond_moments(st, k)
  x_rep <- .cond_replicate(st)
  x_rep[is.na(x_obs)] <- NA

  use <- !is.na(mo$V) & mo$V > 1e-12
  V <- ifelse(use, mo$V, NA_real_)
  res_obs <- ifelse(use, (x_obs - mo$E)^2, NA_real_)
  res_rep <- ifelse(use, (x_rep - mo$E)^2, NA_real_)
  if (is.null(group)) group <- seq_len(k)
  pool <- function(m) as.vector(rowsum(colSums(m, na.rm = TRUE), group))
  sum_v <- pool(V)
  n_use <- pool(use)
  list(
    infit      = pool(res_obs) / sum_v,
    infit_rep  = pool(res_rep) / sum_v,
    outfit     = pool(res_obs / V) / n_use,
    outfit_rep = pool(res_rep / V) / n_use
  )
}

# Data preparation shared by the conditional statistics.
# Returns the index vectors, the observed score matrix (0-based) and the
# family name.
.cond_prepare_data <- function(model, resp_var, item_name, person_name) {
  family_name <- stats::family(model)$family
  obs_response <- as.numeric(model$data[[resp_var]])
  # acat responses are coded 1..C in brms, score them from 0
  obs_score <- if (family_name == "acat") obs_response - 1 else obs_response
  items   <- model$data[[item_name]]
  persons <- model$data[[person_name]]
  unique_items   <- sort(unique(items))
  unique_persons <- sort(unique(persons))
  person_idx <- match(persons, unique_persons)
  item_idx   <- match(items, unique_items)
  x_obs <- matrix(NA_real_, nrow = length(unique_persons),
                  ncol = length(unique_items))
  x_obs[cbind(person_idx, item_idx)] <- obs_score
  list(person_idx = person_idx, item_idx = item_idx,
       unique_items = unique_items, x_obs = x_obs,
       n_cat_max = if (family_name == "acat") max(obs_response) else 2L)
}

# Apply fun(prob, ...) to every draw, fetching posterior_epred in chunks of
# draws to bound memory. Returns a list with one result per draw.
.cond_map_draws <- function(model, draw_ids, n_cat_max, fun, ...) {
  n_obs <- nrow(model$data)
  chunk_size <- max(1L, floor(2e7 / (n_obs * n_cat_max)))
  n_draws <- length(draw_ids)
  out <- vector("list", n_draws)
  chunks <- split(seq_len(n_draws), ceiling(seq_len(n_draws) / chunk_size))
  for (ch in chunks) {
    epred <- brms::posterior_epred(model, draw_ids = draw_ids[ch])
    for (j in seq_along(ch)) {
      out[[ch[j]]] <- fun(.cond_prob_matrix(epred, j), ...)
    }
  }
  out
}

# Draw subset shared by the statistics.
.select_draw_ids <- function(model, ndraws_use) {
  total_draws <- brms::ndraws(model)
  if (is.null(ndraws_use)) return(seq_len(total_draws))
  ndraws_use <- as.integer(ndraws_use)
  if (ndraws_use < 1) {
    stop("'ndraws_use' must be a positive integer.", call. = FALSE)
  }
  if (ndraws_use > total_draws) {
    warning("'ndraws_use' (", ndraws_use, ") exceeds available draws (",
            total_draws, "). Using all draws.", call. = FALSE)
    ndraws_use <- total_draws
  }
  sample(seq_len(total_draws), ndraws_use)
}

# Mid-p posterior predictive p-value: P(obs > rep) + 0.5 * P(obs == rep).
.ppp_mid <- function(obs, rep_draws) {
  vapply(seq_along(obs), function(i) {
    d <- rep_draws[, i]
    d <- d[!is.na(d)]
    if (is.na(obs[i]) || length(d) == 0L) return(NA_real_)
    mean(obs[i] > d) + 0.5 * mean(obs[i] == d)
  }, numeric(1))
}

# Two-sided single-step max-T p-values across items. Each item's deviation
# is standardised by the mean and SD of its replicated distribution. The
# reference is the maximum absolute standardised deviation over items in
# each replicated dataset.
.ppp_maxt <- function(obs, rep_draws) {
  m <- colMeans(rep_draws, na.rm = TRUE)
  s <- apply(rep_draws, 2L, stats::sd, na.rm = TRUE)
  usable <- !is.na(obs) & is.finite(s) & s > 0
  out <- rep(NA_real_, length(obs))
  if (!any(usable)) return(out)
  z_rep <- abs(sweep(sweep(rep_draws[, usable, drop = FALSE], 2L, m[usable]),
                     2L, s[usable], "/"))
  t_max <- apply(z_rep, 1L, max, na.rm = TRUE)
  z_obs <- abs(obs[usable] - m[usable]) / s[usable]
  out[usable] <- vapply(z_obs, function(z) mean(t_max >= z), numeric(1))
  out
}

# Conditional branch of infit_statistic(). Returns the same long format as
# the theta branch.
#
# group_of: optional function mapping the item labels (dd$unique_items) to
# group labels. Items in the same group are pooled, see .cond_infit_draw().
# The output then has one row per group and draw, and the first column is
# named out_name.
.infit_conditional <- function(model, resp_var, item_name, person_name,
                               ndraws_use, outfit, group_of = NULL,
                               out_name = item_name) {
  .check_unique_person_item(model$data[[person_name]], model$data[[item_name]],
                            person_name, item_name, "infit_statistic")
  .check_sufficient_score_brmsfit(model, person_name, "infit_statistic")

  draw_ids <- .select_draw_ids(model, ndraws_use)
  n_draws <- length(draw_ids)
  dd <- .cond_prepare_data(model, resp_var, item_name, person_name)
  labels <- dd$unique_items
  group <- NULL
  if (!is.null(group_of)) {
    g <- group_of(labels)
    labels <- unique(g)
    group <- match(g, labels)
  }
  k <- length(labels)

  per_draw <- .cond_map_draws(model, draw_ids, dd$n_cat_max,
                              .cond_infit_draw,
                              person_idx = dd$person_idx,
                              item_idx = dd$item_idx, x_obs = dd$x_obs,
                              group = group)
  stack <- function(name) do.call(rbind, lapply(per_draw, `[[`, name))

  result <- tibble::tibble(
    item      = rep(labels, each = n_draws),
    draw      = rep(seq_len(n_draws), times = k),
    infit     = as.vector(stack("infit")),
    infit_rep = as.vector(stack("infit_rep"))
  )
  if (outfit) {
    result$outfit     <- as.vector(stack("outfit"))
    result$outfit_rep <- as.vector(stack("outfit_rep"))
  }
  names(result)[1] <- out_name
  result <- dplyr::group_by(result, .data[[out_name]])
  dplyr::arrange(result, .data[[out_name]], .data$draw)
}

# Conditional Q3 for one draw, observed and replicated.
#
# Residuals are X - E(X | total score), with E from the item parameters of
# the draw. Replicated data keep each person's total score. Persons with the
# minimum or maximum total have zero conditional variance and drop out, as
# in conditional infit. Returns list(obs, rep): Q3 for each item pair, upper
# triangle in row-major order (1:2, 1:3, ..., 2:3, ...), or only the pairs
# where `keep` (logical, same order) is TRUE.
.cond_q3_draw <- function(prob, person_idx, item_idx, x_obs, keep = NULL) {
  k <- ncol(x_obs)
  st <- .cond_draw_setup(prob, person_idx, item_idx, x_obs)
  mo <- .cond_moments(st, k)
  x_rep <- .cond_replicate(st)
  x_rep[is.na(x_obs)] <- NA
  use <- !is.na(mo$V) & mo$V > 1e-12
  r_obs <- ifelse(use, x_obs - mo$E, NA_real_)
  r_rep <- ifelse(use, x_rep - mo$E, NA_real_)
  c_obs <- suppressWarnings(stats::cor(r_obs, use = "pairwise.complete.obs"))
  c_rep <- suppressWarnings(stats::cor(r_rep, use = "pairwise.complete.obs"))
  lower <- lower.tri(diag(k))
  obs <- t(c_obs)[lower]
  rep <- t(c_rep)[lower]
  if (!is.null(keep)) {
    obs <- obs[keep]
    rep <- rep[keep]
  }
  list(obs = obs, rep = rep)
}

# Item pairs, upper triangle in row-major order, that at least one person
# answered both. Logical vector of length k (k - 1) / 2.
.q3_shared_pairs <- function(person_idx, item_idx, n_persons, k) {
  answered <- matrix(0, n_persons, k)
  answered[cbind(person_idx, item_idx)] <- 1
  shared <- crossprod(answered)
  t(shared)[lower.tri(shared)] > 0
}

# Conditional branch of q3_statistic(). Returns the same long format as the
# theta branch. With shared_only = TRUE only item pairs that at least one
# person answered both are computed (many-facet cells).
.q3_conditional <- function(model, resp_var, item_name, person_name,
                            ndraws_use, shared_only = FALSE) {
  .check_sufficient_score_brmsfit(model, person_name, "q3_statistic")

  draw_ids <- .select_draw_ids(model, ndraws_use)
  n_draws <- length(draw_ids)
  dd <- .cond_prepare_data(model, resp_var, item_name, person_name)
  k <- length(dd$unique_items)
  keep <- if (shared_only) {
    .q3_shared_pairs(dd$person_idx, dd$item_idx, nrow(dd$x_obs), k)
  }

  per_draw <- .cond_map_draws(model, draw_ids, dd$n_cat_max,
                              .cond_q3_draw,
                              person_idx = dd$person_idx,
                              item_idx = dd$item_idx, x_obs = dd$x_obs,
                              keep = keep)
  obs <- do.call(rbind, lapply(per_draw, `[[`, "obs"))
  rep <- do.call(rbind, lapply(per_draw, `[[`, "rep"))

  pg <- which(upper.tri(diag(k)), arr.ind = TRUE)
  pg <- pg[order(pg[, "row"], pg[, "col"]), , drop = FALSE]
  if (!is.null(keep)) pg <- pg[keep, , drop = FALSE]
  n_pairs <- nrow(pg)
  item_1 <- dd$unique_items[pg[, "row"]]
  item_2 <- dd$unique_items[pg[, "col"]]

  tibble::tibble(
    draw      = rep(seq_len(n_draws), times = n_pairs),
    item_pair = rep(paste0(item_1, " : ", item_2), each = n_draws),
    item_1    = rep(item_1, each = n_draws),
    item_2    = rep(item_2, each = n_draws),
    q3        = as.vector(obs),
    q3_rep    = as.vector(rep)
  )
}

# One-sided single-step max-T p-values across item pairs, for local
# dependence (observed Q3 above the replicated). Each pair's observed mean
# is standardised by the mean and SD of its replicated distribution. The
# reference is the maximum standardised replicated value over pairs in each
# draw. Small values flag dependence.
.ppp_maxt_upper <- function(obs, rep_draws) {
  m <- colMeans(rep_draws, na.rm = TRUE)
  s <- apply(rep_draws, 2L, stats::sd, na.rm = TRUE)
  usable <- !is.na(obs) & is.finite(s) & s > 0
  out <- rep(NA_real_, length(obs))
  if (!any(usable)) return(out)
  z_rep <- sweep(sweep(rep_draws[, usable, drop = FALSE], 2L, m[usable]),
                 2L, s[usable], "/")
  t_max <- apply(z_rep, 1L, max, na.rm = TRUE)
  z_obs <- (obs[usable] - m[usable]) / s[usable]
  out[usable] <- vapply(z_obs, function(z) mean(t_max >= z), numeric(1))
  out
}
