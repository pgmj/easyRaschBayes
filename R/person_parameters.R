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
#'   range for the Newton-Raphson WLE search. Default is
#'   \code{c(-7, 7)}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{person_estimates}{A \code{\link[tibble]{tibble}} with one
#'     row per person containing: the person ID, \code{sum_score},
#'     \code{eap} (posterior mean of person random effect),
#'     \code{eap_se} (posterior SD), \code{wle} (Warm's weighted
#'     likelihood estimate), and \code{wle_se} (asymptotic SE
#'     of the WLE). Rows are ordered to match the original person
#'     order in the model data (i.e., the order of first appearance
#'     of each person ID).}
#'   \item{score_table}{A \code{\link[tibble]{tibble}} mapping each
#'     observed ordinal sum score to its mean EAP, mean EAP SE,
#'     WLE, and WLE SE. Extreme scores (0 and maximum) receive
#'     WLE estimates at the boundary of \code{theta_range} with
#'     \code{NA} standard errors.}
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
#' posterior mean item parameters using Newton-Raphson iteration
#' with adaptive step damping on the Warm-corrected likelihood.
#' WLE adds a bias correction term \eqn{J(\theta) / (2 I(\theta))}
#' to the score equations, where \eqn{I(\theta)} is the test
#' information and \eqn{J(\theta) = \sum_i \sum_c (c - E_i)^3 P_{ic}}
#' is the sum of third central moments. This produces estimates with
#' reduced finite-sample bias compared to MLE, especially at
#' extreme scores (Warm, 1989).
#'
#' The Newton-Raphson algorithm uses adaptive step damping following
#' the approach in \pkg{iarm} (Mueller): the maximum allowed step
#' size shrinks by a factor of 1.05 each iteration, preventing
#' overshoot and ensuring convergence for near-extreme scores.
#'
#' For extreme scores (sum score = 0 or maximum possible), the
#' WLE is not well-defined. These cases are assigned the boundary
#' values of \code{theta_range} with \code{NA} standard errors.
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
#' df_rm <- eRm::raschdat3 %>%
#'   as.data.frame() %>%
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
#' df_pcm <- eRm::pcmdat2 %>%
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
  # DETECT MODEL PARAMETERISATION
  # ================================================================
  fe_pattern <- paste0("^b_", item_name)
  fe_cols    <- grep(fe_pattern, names(all_draws), value = TRUE)
  fe_cols    <- setdiff(fe_cols, "b_Intercept")
  
  has_fixed_items  <- !is_ordinal && length(fe_cols) > 0
  has_random_items <- !is_ordinal && any(grepl(
    paste0("^r_", item_name, "\\["), names(all_draws)
  ))
  
  if (!is_ordinal && !has_fixed_items && !has_random_items) {
    stop("Could not detect item parameters in the model. ",
         "Supported parameterisations:\n",
         "  - response ~ 0 + item + (1 | id)    [fixed items]\n",
         "  - response ~ 1 + (1 | item) + (1 | id)  [random items]\n",
         "  - response | thres(gr = item) ~ 1 + (1 | id) [PCM/acat]",
         call. = FALSE)
  }
  
  # ================================================================
  # EXTRACT ITEM PARAMETERS (posterior means)
  # ================================================================
  unique_items <- sort(unique(model$data[[item_name]]))
  n_items <- length(unique_items)
  
  if (is_ordinal) {
    item_params <- .extract_item_params_ordinal(
      all_draws, unique_items, item_name
    )
  } else if (has_fixed_items) {
    item_params <- .extract_item_params_fixed(
      all_draws, fe_cols, unique_items, item_name
    )
  } else {
    item_params <- .extract_item_params_random(
      all_draws, unique_items, item_name
    )
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
  # COMPUTE SUM SCORES
  # ================================================================
  dat <- model$data
  sum_scores <- stats::aggregate(
    stats::as.formula(paste(resp_var, "~", person_name)),
    data = dat,
    FUN  = sum,
    na.rm = TRUE
  )
  colnames(sum_scores) <- c("person_id", "sum_score")
  
  # For ordinal models, the minimum response category may be > 0
  # (brms acat uses 1-based coding). Adjust sum scores to be
  # 0-based for the WLE computation.
  min_resp <- min(dat[[resp_var]], na.rm = TRUE)
  if (min_resp > 0) {
    n_obs_per_person <- stats::aggregate(
      stats::as.formula(paste(resp_var, "~", person_name)),
      data = dat,
      FUN  = length
    )
    colnames(n_obs_per_person) <- c("person_id", "n_obs")
    sum_scores <- merge(sum_scores, n_obs_per_person, by = "person_id")
    sum_scores$sum_score <- sum_scores$sum_score -
      (min_resp * sum_scores$n_obs)
    sum_scores$n_obs <- NULL
  }
  
  # ================================================================
  # WLE COMPUTATION
  # ================================================================
  if (is_ordinal) {
    max_score <- sum(vapply(
      item_params$thresholds, length, integer(1)
    ))
    wle_results <- .compute_wle_polytomous(
      item_params$thresholds, max_score, theta_range
    )
  } else {
    max_score <- n_items
    wle_results <- .compute_wle_dichotomous(
      item_params$difficulties, theta_range
    )
  }
  
  # ================================================================
  # BUILD PERSON-LEVEL TABLE (in ranef order first, then reorder)
  # ================================================================
  person_df <- tibble::tibble(
    person_id = ranef_ids,
    sum_score = sum_scores$sum_score[
      match(ranef_ids, as.character(sum_scores$person_id))
    ],
    eap    = round(as.numeric(eap), 4),
    eap_se = round(as.numeric(eap_se), 4)
  )
  
  # Match WLE to persons via sum score
  wle_lookup <- tibble::tibble(
    sum_score = wle_results$scores,
    wle       = round(wle_results$wle, 4),
    wle_se    = round(wle_results$wle_se, 4)
  )
  
  person_df <- dplyr::left_join(person_df, wle_lookup, by = "sum_score")
  
  # Rename the person_id column to the actual variable name
  colnames(person_df)[1] <- person_name
  
  # Reorder to match original person order in model data
  orig_idx <- match(original_order, person_df[[person_name]])
  person_df <- person_df[orig_idx, ]
  
  # ================================================================
  # BUILD SCORE TABLE (sum score -> EAP and WLE)
  # ================================================================
  score_table <- person_df |>
    dplyr::group_by(.data$sum_score) |>
    dplyr::summarise(
      n       = dplyr::n(),
      eap     = round(mean(.data$eap, na.rm = TRUE), 4),
      eap_se  = round(mean(.data$eap_se, na.rm = TRUE), 4),
      wle     = round(mean(.data$wle, na.rm = TRUE), 4),
      wle_se  = round(mean(.data$wle_se, na.rm = TRUE), 4),
      .groups = "drop"
    )
  
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


# ══════════════════════════════════════════════════════════════════
# ITEM PARAMETER EXTRACTION HELPERS
# ══════════════════════════════════════════════════════════════════

# ── Internal: extract item parameters for ordinal (acat/PCM) ─────
#' @keywords internal
.extract_item_params_ordinal <- function(draws, unique_items,
                                         item_name) {
  thresholds <- list()
  
  for (item_label in unique_items) {
    thresh_pattern <- paste0(
      "^b_Intercept\\[",
      gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
      ","
    )
    thresh_cols <- grep(thresh_pattern, names(draws), value = TRUE)
    
    if (length(thresh_cols) == 0) {
      warning("No threshold parameters found for item '",
              item_label, "'.", call. = FALSE)
      next
    }
    
    thresh_nums <- as.numeric(
      gsub(".*,(\\d+)\\]$", "\\1", thresh_cols)
    )
    thresh_cols <- thresh_cols[order(thresh_nums)]
    
    thresholds[[item_label]] <- vapply(
      thresh_cols,
      function(col) mean(draws[[col]]),
      numeric(1)
    )
  }
  
  list(thresholds = thresholds)
}


# ── Internal: extract item params for fixed-effects dichotomous ──
#' @keywords internal
.extract_item_params_fixed <- function(draws, fe_cols, unique_items,
                                       item_name) {
  difficulties <- vapply(
    fe_cols,
    function(col) -mean(draws[[col]]),
    numeric(1)
  )
  names(difficulties) <- unique_items
  list(difficulties = difficulties)
}


# ── Internal: extract item params for random-effects dichotomous ─
#' @keywords internal
.extract_item_params_random <- function(draws, unique_items,
                                        item_name) {
  intercept_col <- grep("^b_Intercept$", names(draws), value = TRUE)
  intercept_mean <- if (length(intercept_col) == 1) {
    mean(draws[[intercept_col]])
  } else {
    0
  }
  
  difficulties <- vapply(unique_items, function(item_label) {
    re_col <- paste0("r_", item_name, "[", item_label, ",Intercept]")
    if (!re_col %in% names(draws)) {
      re_col_alt <- grep(
        paste0("^r_", item_name, "\\[",
               gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
               ",Intercept\\]$"),
        names(draws), value = TRUE
      )
      if (length(re_col_alt) == 1) re_col <- re_col_alt
      else return(NA_real_)
    }
    -(intercept_mean + mean(draws[[re_col]]))
  }, numeric(1))
  
  names(difficulties) <- unique_items
  list(difficulties = difficulties)
}


# ══════════════════════════════════════════════════════════════════
# WLE COMPUTATION HELPERS
# ══════════════════════════════════════════════════════════════════
#
# Newton-Raphson with adaptive step damping following the approach
# in iarm::persons_mle (Mueller): the maximum allowed step size
# starts at maxdelta_start and shrinks by a factor of
# maxdelta_shrink each iteration, preventing overshoot and ensuring
# convergence for near-extreme scores.
# ══════════════════════════════════════════════════════════════════

# ── Internal: WLE for dichotomous Rasch ──────────────────────────
#' @keywords internal
.compute_wle_dichotomous <- function(difficulties, theta_range,
                                     max_iter = 100, tol = 1e-4,
                                     maxdelta_start = 3,
                                     maxdelta_shrink = 1.05) {
  n_items  <- length(difficulties)
  scores   <- 0:n_items
  n_scores <- length(scores)
  maxval   <- max(abs(theta_range))
  
  wle    <- rep(NA_real_, n_scores)
  wle_se <- rep(NA_real_, n_scores)
  
  for (idx in seq_along(scores)) {
    r <- scores[idx]
    
    # Extreme scores: boundary values, no SE
    if (r == 0) {
      wle[idx]    <- theta_range[1]
      wle_se[idx] <- NA_real_
      next
    }
    if (r == n_items) {
      wle[idx]    <- theta_range[2]
      wle_se[idx] <- NA_real_
      next
    }
    
    # Newton-Raphson with adaptive step damping
    theta_hat <- 0
    maxdelta  <- maxdelta_start
    
    for (iter in seq_len(max_iter)) {
      theta0 <- theta_hat
      
      p_i <- stats::plogis(theta_hat - difficulties)
      q_i <- 1 - p_i
      
      # Score function and derivatives
      dll  <- r - sum(p_i)                    # first derivative
      d2ll <- -sum(p_i * q_i)                 # second derivative
      d3ll <- -sum(p_i * q_i * (q_i - p_i))  # third derivative
      
      # WLE step: NR + Warm bias correction
      delta <- -dll / d2ll - d3ll / (2 * d2ll^2)
      
      # Adaptive damping: shrink maxdelta each iteration
      maxdelta <- maxdelta / maxdelta_shrink
      if (abs(delta) > maxdelta) {
        delta <- sign(delta) * maxdelta
      }
      
      theta_hat <- theta_hat + delta
      
      # Safety clamp at boundary
      if (abs(theta_hat) > maxval) {
        theta_hat <- sign(theta_hat) * maxval
      }
      
      if (abs(theta_hat - theta0) < tol) break
    }
    
    wle[idx] <- theta_hat
    
    # SE = 1/sqrt(I) at convergence; NA if at boundary
    p_i     <- stats::plogis(theta_hat - difficulties)
    q_i     <- 1 - p_i
    I_final <- sum(p_i * q_i)
    
    if (abs(theta_hat) >= maxval || I_final < 1e-12) {
      wle_se[idx] <- NA_real_
    } else {
      wle_se[idx] <- 1 / sqrt(I_final)
    }
  }
  
  list(scores = scores, wle = wle, wle_se = wle_se)
}


# ── Internal: WLE for polytomous (acat/PCM) ─────────────────────
#' @keywords internal
.compute_wle_polytomous <- function(thresholds, max_score, theta_range,
                                    max_iter = 100, tol = 1e-4,
                                    maxdelta_start = 3,
                                    maxdelta_shrink = 1.05) {
  scores   <- 0:max_score
  n_scores <- length(scores)
  n_items  <- length(thresholds)
  maxval   <- max(abs(theta_range))
  
  wle    <- rep(NA_real_, n_scores)
  wle_se <- rep(NA_real_, n_scores)
  
  for (idx in seq_along(scores)) {
    r <- scores[idx]
    
    # Extreme scores: boundary values, no SE
    if (r == 0) {
      wle[idx]    <- theta_range[1]
      wle_se[idx] <- NA_real_
      next
    }
    if (r == max_score) {
      wle[idx]    <- theta_range[2]
      wle_se[idx] <- NA_real_
      next
    }
    
    # Newton-Raphson with adaptive step damping
    theta_hat <- 0
    maxdelta  <- maxdelta_start
    
    for (iter in seq_len(max_iter)) {
      theta0 <- theta_hat
      
      # Accumulate derivatives across items
      dll  <- r   # will subtract E_i for each item
      d2ll <- 0
      d3ll <- 0
      
      for (i in seq_len(n_items)) {
        tau_i <- thresholds[[i]]
        n_cat <- length(tau_i) + 1
        cats  <- seq(0, n_cat - 1)
        
        # Adjacent category model probabilities (log-sum-exp stable)
        cum_eta   <- c(0, cumsum(theta_hat - tau_i))
        max_eta   <- max(cum_eta)
        log_denom <- max_eta + log(sum(exp(cum_eta - max_eta)))
        p_cat     <- exp(cum_eta - log_denom)
        
        E_i <- sum(cats * p_cat)              # expected score
        V_i <- sum((cats - E_i)^2 * p_cat)   # variance = info
        K_i <- sum((cats - E_i)^3 * p_cat)   # third central moment
        
        dll  <- dll - E_i
        d2ll <- d2ll - V_i
        d3ll <- d3ll + K_i
      }
      
      # Guard against zero information
      if (abs(d2ll) < 1e-12) break
      
      # WLE step: NR + Warm bias correction
      delta <- -dll / d2ll - d3ll / (2 * d2ll^2)
      
      # Adaptive damping: shrink maxdelta each iteration
      maxdelta <- maxdelta / maxdelta_shrink
      if (abs(delta) > maxdelta) {
        delta <- sign(delta) * maxdelta
      }
      
      theta_hat <- theta_hat + delta
      
      # Safety clamp at boundary
      if (abs(theta_hat) > maxval) {
        theta_hat <- sign(theta_hat) * maxval
      }
      
      if (abs(theta_hat - theta0) < tol) break
    }
    
    wle[idx] <- theta_hat
    
    # Recompute information at convergence for SE
    I_final <- 0
    for (i in seq_len(n_items)) {
      tau_i <- thresholds[[i]]
      n_cat <- length(tau_i) + 1
      cats  <- seq(0, n_cat - 1)
      cum_eta   <- c(0, cumsum(theta_hat - tau_i))
      max_eta   <- max(cum_eta)
      log_denom <- max_eta + log(sum(exp(cum_eta - max_eta)))
      p_cat     <- exp(cum_eta - log_denom)
      E_i <- sum(cats * p_cat)
      V_i <- sum((cats - E_i)^2 * p_cat)
      I_final <- I_final + V_i
    }
    
    if (abs(theta_hat) >= maxval || I_final < 1e-12) {
      wle_se[idx] <- NA_real_
    } else {
      wle_se[idx] <- 1 / sqrt(I_final)
    }
  }
  
  list(scores = scores, wle = wle, wle_se = wle_se)
}