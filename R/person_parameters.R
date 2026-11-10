#' Extract Person Parameters from a Bayesian Rasch Model
#'
#' Extracts person ability estimates from a fitted Bayesian Rasch
#' model. Returns both Bayesian EAP (expected a posteriori) estimates
#' with posterior SDs and frequentist WLE (Warm's weighted likelihood)
#' estimates with standard errors, plus a lookup table mapping ordinal
#' sum scores to both scales.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object. Supported
#'   parameterisations:
#'   \describe{
#'     \item{Polytomous ordinal (PCM)}{e.g., \code{family = acat} with
#'       \code{thres(gr = item)}, producing item-specific thresholds.}
#'     \item{Dichotomous Rasch (random items)}{e.g.,
#'       \code{response ~ 1 + (1 | item) + (1 | id)} with
#'       \code{family = bernoulli()}.}
#'     \item{Dichotomous 1PL (fixed items)}{e.g.,
#'       \code{response ~ 0 + item + (1 | id)} with
#'       \code{family = bernoulli()}.}
#'     \item{Rating scale model}{\code{family = acat} with common
#'       thresholds and item effects, e.g.,
#'       \code{response ~ 1 + (1 | item) + (1 | id)}.}
#'   }
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Default is \code{id}.
#' @param draws Logical. If \code{TRUE}, a matrix of full posterior
#'   draws (persons x draws) is included in the output. Default is
#'   \code{FALSE}.
#' @param center Logical. If \code{TRUE} (the default), person
#'   parameters and item difficulties are recentered so that the
#'   mean item difficulty is zero, matching the convention in
#'   frequentist CML Rasch estimation. If \code{FALSE}, raw brms
#'   parameterisation is used.
#' @param theta_range A numeric vector of length 2 specifying the
#'   search range for the WLE. A WLE outside it is set to the nearer
#'   bound with an \code{NA} standard error. Default is
#'   \code{c(-7, 7)}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{person_estimates}{A \code{\link[tibble]{tibble}} with one
#'     row per person containing: the person ID, \code{n_items}
#'     (number of items answered), \code{sum_score} (sum of the item
#'     scores, each scored from 0),
#'     \code{eap} (posterior mean of person random effect),
#'     \code{eap_se} (posterior SD), \code{wle} (Warm's weighted
#'     likelihood estimate), and \code{wle_se} (asymptotic SE
#'     of the WLE). Rows are ordered to match the original person
#'     order in the model data (i.e., the order of first appearance
#'     of each person ID).}
#'   \item{score_table}{A \code{\link[tibble]{tibble}} mapping each
#'     sum score observed among persons who answered every item to
#'     the number of such persons, their mean EAP and mean EAP SE,
#'     and the WLE and WLE SE for that score on the complete item
#'     set.}
#'   \item{draws_matrix}{(Only if \code{draws = TRUE}) A numeric
#'     matrix with rows = persons and columns = posterior draws.
#'     Row names are person IDs. Rows are ordered to match the
#'     original person order in the model data. Can be passed
#'     directly to \code{\link{RMUreliability}}.}
#' }
#'
#' @details
#' \strong{EAP estimates} are extracted as the posterior means of
#' the person random effects (\code{r_id[j, Intercept]}) from
#' \code{\link[brms]{as_draws_df}}, with the posterior SD serving
#' as the standard error. These are the standard Bayesian point
#' estimates and reflect shrinkage toward the population mean.
#'
#' \strong{WLE estimates} (Warm, 1989) are computed from the
#' posterior mean item parameters by solving the Warm-corrected
#' score equation
#' \deqn{r - E(\theta) + J(\theta) / (2 I(\theta)) = 0}
#' with \code{\link[stats]{uniroot}}, where \eqn{r} is the sum score,
#' \eqn{E(\theta)} the expected sum score, \eqn{I(\theta)} the test
#' information and \eqn{J(\theta) = \sum_i \sum_c (c - E_i)^3 P_{ic}}
#' the sum of third central moments. The correction reduces the
#' finite-sample bias of the MLE and keeps the estimate finite at the
#' minimum and maximum scores. The standard error is
#' \eqn{1 / \sqrt{I(\theta)}}.
#'
#' A sum score only has a WLE relative to the items that make it up.
#' Persons with missing responses therefore get the WLE for their sum
#' score on the items they answered, computed separately for each set
#' of answered items.
#'
#' When \code{center = TRUE} (the default), item difficulty
#' parameters are shifted so their mean is zero, and EAP person
#' parameters are shifted by the same constant. WLE is computed
#' from the centered item parameters. This matches the convention
#' in frequentist CML Rasch estimation.
#'
#' \strong{Row ordering:} Both \code{person_estimates} and
#' \code{draws_matrix} preserve the original person order from
#' the model data (order of first appearance of each person ID).
#' This allows direct row-binding with the source data without
#' re-matching.
#'
#' @references
#' Warm, T. A. (1989). Weighted likelihood estimation of ability
#' in item response theory. \emph{Psychometrika}, \emph{54}(3),
#' 427--450. \doi{10.1007/BF02294627}
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R
#' with brms and Stan. \emph{Journal of Statistical Software},
#' \emph{100}, 1--54. \doi{10.18637/jss.v100.i05}
#'
#' Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013).
#' \emph{Rasch Models in Health}. Iste and Wiley, pp. 63--70.
#'
#' @seealso
#' \code{\link{RMUreliability}}, \code{\link{plot_targeting}},
#' \code{\link[brms]{ranef}}, \code{\link[brms]{as_draws_df}}.
#'
#' @examples
#' \donttest{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' # --- Dichotomous Rasch Model (random items) ---
#'
#' # PHQ-9 items dichotomised (any symptom vs none), for illustration only
#' df_rm <- phq9[1:200, 1:9] %>%
#'   mutate(across(everything(), ~ as.integer(.x > 0))) %>%
#'   rownames_to_column("id") %>%
#'   pivot_longer(!id, names_to = "item", values_to = "response")
#'
#' fit_rm <- brm(
#'   response ~ 1 + (1 | item) + (1 | id),
#'   data   = df_rm,
#'   family = bernoulli(),
#'   chains = 4, cores = 2, iter = 1000 # use more iter and cores
#' )
#'
#' # Basic usage — rows match original person order
#' pp <- person_parameters(fit_rm)
#' pp$person_estimates
#' pp$score_table
#'
#' # With full posterior draws (e.g., for RMUreliability)
#' pp_draws <- person_parameters(fit_rm, draws = TRUE)
#' RMUreliability(pp_draws$draws_matrix)
#'
#' # --- Polytomous PCM (acat) ---
#'
#' df_pcm <- phq9[1:200, 1:9] %>%
#'   mutate(across(everything(), ~ .x + 1)) %>%
#'   rownames_to_column("id") %>%
#'   pivot_longer(!id, names_to = "item", values_to = "response")
#'
#' fit_pcm <- brm(
#'   response | thres(gr = item) ~ 1 + (1 | id),
#'   data   = df_pcm,
#'   family = acat,
#'   chains = 4, cores = 2, iter = 1000 # use more iter and cores
#' )
#'
#' pp_pcm <- person_parameters(fit_pcm)
#' pp_pcm$score_table
#' }
#'
#' @importFrom brms as_draws_df ranef
#' @importFrom rlang enquo as_name .data
#' @importFrom stats family sd var plogis as.formula
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr left_join group_by summarise arrange mutate n
#' @export
person_parameters <- function(
    model,
    item_var    = item,
    person_var  = id,
    draws       = FALSE,
    center      = TRUE,
    theta_range = c(-7, 7)
) {
  
  # ================================================================
  # INPUT VALIDATION
  # ================================================================
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  
  item_name   <- rlang::as_name(rlang::enquo(item_var))
  person_name <- rlang::as_name(rlang::enquo(person_var))
  
  if (!item_name %in% names(model$data)) {
    stop("Item variable '", item_name, "' not found in model data.",
         call. = FALSE)
  }
  if (!person_name %in% names(model$data)) {
    stop("Person variable '", person_name, "' not found in model data.",
         call. = FALSE)
  }
  
  resp_var <- as.character(stats::formula(model)$formula[[2]])
  if (length(resp_var) > 1) resp_var <- resp_var[2]
  
  all_draws   <- tibble::as_tibble(brms::as_draws_df(model))
  family_name <- stats::family(model)$family
  is_ordinal  <- grepl("acat|cumul|sratio|cratio",
                       family_name, ignore.case = TRUE)
  
  # ================================================================
  # ORIGINAL PERSON ORDER (first appearance in model data)
  # ================================================================
  person_col_data   <- as.character(model$data[[person_name]])
  original_order    <- unique(person_col_data)
  
  # ================================================================
  # ITEM PARAMETERS (posterior means)
  # ================================================================
  tau_draws <- .model_tau_draws(model, item_name, draws = all_draws)
  unique_items <- names(tau_draws)
  n_items <- length(unique_items)

  if (is_ordinal) {
    item_params <- list(thresholds = lapply(tau_draws, colMeans))
  } else {
    item_params <- list(difficulties = vapply(tau_draws, mean, numeric(1)))
  }

  # ================================================================
  # CENTERING
  # ================================================================
  if (center && !is_ordinal) {
    shift <- mean(item_params$difficulties)
    item_params$difficulties <- item_params$difficulties - shift
  } else if (center && is_ordinal) {
    all_thresholds <- unlist(item_params$thresholds)
    shift <- mean(all_thresholds)
    item_params$thresholds <- lapply(
      item_params$thresholds, function(x) x - shift
    )
  } else {
    shift <- 0
  }
  
  # ================================================================
  # EAP PERSON PARAMETERS
  # ================================================================
  person_re <- brms::ranef(model)[[person_name]]
  if (is.null(person_re)) {
    stop("No random effects found for '", person_name,
         "'. Check your model formula.", call. = FALSE)
  }
  
  # ranef() returns persons in its own order (typically alphabetical)
  ranef_ids   <- rownames(person_re[, , "Intercept"])
  eap_raw     <- person_re[, "Estimate", "Intercept"]
  eap_se_raw  <- person_re[, "Est.Error", "Intercept"]
  
  if (center) {
    eap <- eap_raw - shift
  } else {
    eap <- eap_raw
  }
  eap_se <- eap_se_raw
  
  # ================================================================
  # SUM SCORES AND ANSWERED ITEMS
  # ================================================================
  # Scores from 0. brms ordinal responses are categories 1..C (category 1
  # is the base of the thresholds), and an ordered factor is stored as such
  # in model$data.
  dat <- model$data
  resp <- dat[[resp_var]]
  score <- if (is.factor(resp)) as.integer(resp) else as.numeric(resp)
  if (is_ordinal) score <- score - 1

  pid <- as.character(dat[[person_name]])
  iid <- match(as.character(dat[[item_name]]), unique_items)
  persons <- unique(pid)
  pidx <- match(pid, persons)
  sum_score <- as.numeric(tapply(score, factor(pidx, levels = seq_along(persons)),
                                 sum))
  answered <- matrix(FALSE, nrow = length(persons), ncol = n_items)
  answered[cbind(pidx, iid)] <- TRUE
  n_answered <- as.integer(rowSums(answered))

  # ================================================================
  # WLE COMPUTATION, per set of answered items
  # ================================================================
  # A sum score only has a WLE relative to the items that make it up, so
  # persons are grouped by the set of items they answered and each set gets
  # its own score-to-WLE table.
  wle_table <- function(items) {
    if (is_ordinal) {
      thr <- item_params$thresholds[items]
      .compute_wle_polytomous(thr, sum(lengths(thr)), theta_range)
    } else {
      .compute_wle_dichotomous(item_params$difficulties[items], theta_range)
    }
  }
  wb <- .wle_by_pattern(answered, sum_score, wle_table)
  wle <- wb$wle
  wle_se <- wb$wle_se
  full_table <- wle_table(seq_len(n_items))

  # ================================================================
  # BUILD PERSON-LEVEL TABLE (in ranef order first, then reorder)
  # ================================================================
  m <- match(ranef_ids, persons)
  person_df <- tibble::tibble(
    person_id = ranef_ids,
    n_items   = n_answered[m],
    sum_score = sum_score[m],
    eap       = round(as.numeric(eap), 4),
    eap_se    = round(as.numeric(eap_se), 4),
    wle       = round(wle[m], 4),
    wle_se    = round(wle_se[m], 4)
  )

  # Rename the person_id column to the actual variable name
  colnames(person_df)[1] <- person_name
  
  # Reorder to match original person order in model data
  orig_idx <- match(original_order, person_df[[person_name]])
  person_df <- person_df[orig_idx, ]
  
  # ================================================================
  # BUILD SCORE TABLE (sum score -> EAP and WLE)
  # ================================================================
  # Only persons who answered every item share one sum-score scale.
  score_table <- person_df[person_df$n_items == n_items, ] |>
    dplyr::group_by(.data$sum_score) |>
    dplyr::summarise(
      n       = dplyr::n(),
      eap     = round(mean(.data$eap, na.rm = TRUE), 4),
      eap_se  = round(mean(.data$eap_se, na.rm = TRUE), 4),
      .groups = "drop"
    )
  fm <- match(score_table$sum_score, full_table$scores)
  score_table$wle    <- round(full_table$wle[fm], 4)
  score_table$wle_se <- round(full_table$wle_se[fm], 4)
  
  # ================================================================
  # FULL POSTERIOR DRAWS (optional)
  # ================================================================
  out <- list(
    person_estimates = person_df,
    score_table      = score_table
  )
  
  if (draws) {
    re_pattern <- paste0("^r_", person_name, "\\[")
    re_cols <- grep(re_pattern, names(all_draws), value = TRUE)
    
    if (length(re_cols) == 0) {
      warning("Could not find posterior draws for person random effects.",
              call. = FALSE)
    } else {
      draws_mat <- t(as.matrix(all_draws[, re_cols]))
      if (center) {
        draws_mat <- draws_mat - shift
      }
      
      # Clean row names: "r_id[person1,Intercept]" -> "person1"
      rn <- rownames(draws_mat)
      rn <- gsub(paste0("^r_", person_name, "\\["), "", rn)
      rn <- gsub(",Intercept\\]$", "", rn)
      rownames(draws_mat) <- rn
      
      # Reorder rows to match original person order in model data
      draws_idx <- match(original_order, rn)
      draws_mat <- draws_mat[draws_idx, , drop = FALSE]
      
      out$draws_matrix <- draws_mat
    }
  }
  
  out
}


# ── Internal: WLE per set of answered items ─────────────────────
# answered: persons x items logical matrix, score: sum score per person,
# table_fun(items) returns list(scores, wle, wle_se) for those items.
#' @keywords internal
.wle_by_pattern <- function(answered, score, table_fun) {
  pattern <- apply(answered, 1, function(x) paste(which(x), collapse = ","))
  wle <- wle_se <- rep(NA_real_, nrow(answered))
  for (pat in unique(pattern)) {
    who <- which(pattern == pat)
    tab <- table_fun(which(answered[who[1], ]))
    m <- match(score[who], tab$scores)
    wle[who]    <- tab$wle[m]
    wle_se[who] <- tab$wle_se[m]
  }
  list(wle = wle, wle_se = wle_se)
}


# ══════════════════════════════════════════════════════════════════
# WLE COMPUTATION HELPERS
# ══════════════════════════════════════════════════════════════════
#
# The Warm (1989) weighted-likelihood estimate solves
#   l'(theta) + J(theta) / (2 I(theta)) = 0,
# where J is the summed third central moment of the item scores. The
# bias-correction term keeps the equation solvable at the extreme
# scores, so the WLE is finite there (a sensible step beyond the next
# score) rather than diverging like the MLE. Each score's root is found
# by bracketed search (stats::uniroot) over theta_range; only a root
# that falls outside that range is clamped to the nearer boundary.
# ══════════════════════════════════════════════════════════════════

# ── Internal: WLE for dichotomous Rasch ──────────────────────────
#' @keywords internal
.compute_wle_dichotomous <- function(difficulties, theta_range, tol = 1e-6) {
  n_items <- length(difficulties)
  scores  <- 0:n_items
  lo <- theta_range[1]; hi <- theta_range[2]

  wle    <- rep(NA_real_, length(scores))
  wle_se <- rep(NA_real_, length(scores))

  # Warm score: l'(theta) + J(theta) / (2 I(theta)).
  warm_score <- function(theta, r) {
    p     <- stats::plogis(theta - difficulties)
    info  <- sum(p * (1 - p))
    third <- sum(p * (1 - p) * (1 - 2 * p))
    (r - sum(p)) + third / (2 * info)
  }
  info_at <- function(theta) {
    p <- stats::plogis(theta - difficulties)
    sum(p * (1 - p))
  }

  for (idx in seq_along(scores)) {
    r    <- scores[idx]
    g_lo <- warm_score(lo, r)
    g_hi <- warm_score(hi, r)

    if (is.finite(g_lo) && is.finite(g_hi) && g_lo * g_hi < 0) {
      theta <- stats::uniroot(warm_score, c(lo, hi), r = r, tol = tol)$root
      info  <- info_at(theta)
      wle[idx]    <- theta
      wle_se[idx] <- if (info < 1e-12) NA_real_ else 1 / sqrt(info)
    } else {
      # Root outside the search range: clamp to the nearer boundary.
      wle[idx]    <- if (r <= n_items / 2) lo else hi
      wle_se[idx] <- NA_real_
    }
  }

  list(scores = scores, wle = wle, wle_se = wle_se)
}


# ── Internal: WLE for polytomous (acat/PCM) ─────────────────────
#' @keywords internal
.compute_wle_polytomous <- function(thresholds, max_score, theta_range,
                                    tol = 1e-6) {
  scores  <- 0:max_score
  n_items <- length(thresholds)
  lo <- theta_range[1]; hi <- theta_range[2]

  wle    <- rep(NA_real_, length(scores))
  wle_se <- rep(NA_real_, length(scores))

  # Adjacent-category probabilities (log-sum-exp stable).
  cat_probs <- function(theta, tau_i) {
    cum_eta <- c(0, cumsum(theta - tau_i))
    m       <- max(cum_eta)
    exp(cum_eta - (m + log(sum(exp(cum_eta - m)))))
  }
  # Summed expected score, information and third central moment.
  moments <- function(theta) {
    E_sum <- 0; info <- 0; third <- 0
    for (i in seq_len(n_items)) {
      tau_i <- thresholds[[i]]
      cats  <- 0:length(tau_i)
      p     <- cat_probs(theta, tau_i)
      E_i   <- sum(cats * p)
      E_sum <- E_sum + E_i
      info  <- info  + sum((cats - E_i)^2 * p)
      third <- third + sum((cats - E_i)^3 * p)
    }
    list(E = E_sum, info = info, third = third)
  }
  # Warm score: l'(theta) + J(theta) / (2 I(theta)).
  warm_score <- function(theta, r) {
    m <- moments(theta)
    (r - m$E) + m$third / (2 * m$info)
  }

  for (idx in seq_along(scores)) {
    r    <- scores[idx]
    g_lo <- warm_score(lo, r)
    g_hi <- warm_score(hi, r)

    if (is.finite(g_lo) && is.finite(g_hi) && g_lo * g_hi < 0) {
      theta <- stats::uniroot(warm_score, c(lo, hi), r = r, tol = tol)$root
      info  <- moments(theta)$info
      wle[idx]    <- theta
      wle_se[idx] <- if (info < 1e-12) NA_real_ else 1 / sqrt(info)
    } else {
      # Root outside the search range: clamp to the nearer boundary.
      wle[idx]    <- if (r <= max_score / 2) lo else hi
      wle_se[idx] <- NA_real_
    }
  }

  list(scores = scores, wle = wle, wle_se = wle_se)
}