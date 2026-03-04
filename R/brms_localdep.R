#' Posterior Predictive Q3 Residual Correlations for Bayesian IRT Models
#'
#' Computes a Bayesian analogue of Yen's Q3 statistic (Yen, 1984) for
#' detecting local dependence between item pairs in Rasch-family models
#' fitted with \pkg{brms}. For each posterior draw, residual correlations
#' are computed for both observed and replicated data, yielding a
#' posterior predictive p-value for each item pair that is automatically
#' calibrated without requiring knowledge of the sampling distribution.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from an
#'   ordinal IRT model (e.g., \code{family = acat} for a partial credit
#'   model) or a dichotomous model (\code{family = bernoulli()}).
#' @param item_var An unquoted variable name identifying the item grouping
#'   variable in the model data (e.g., \code{item}).
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data (e.g., \code{id}).
#' @param ndraws_use Optional positive integer. If specified, a random
#'   subset of posterior draws of this size is used. If \code{NULL} (the
#'   default), all draws are used.
#'
#' @return A \code{\link[tibble]{tibble}} with the following columns:
#' \describe{
#'   \item{item_1}{First item in the pair.}
#'   \item{item_2}{Second item in the pair.}
#'   \item{q3_obs}{Posterior mean of the observed Q3 residual correlation.}
#'   \item{q3_rep}{Posterior mean of the replicated Q3 residual
#'     correlation.}
#'   \item{q3_diff}{Posterior mean of \code{q3_obs - q3_rep}. Large
#'     positive values indicate that the observed residual correlation
#'     exceeds what the model expects.}
#'   \item{ppp}{Posterior predictive p-value:
#'     \code{mean(q3_obs > q3_rep)} across draws. Values close to 1
#'     indicate local dependence (observed correlation systematically
#'     higher than replicated).}
#'   \item{q3_obs_q025, q3_obs_q975}{2.5\% and 97.5\% quantiles
#'     (95\% credible interval) of the posterior distribution of
#'     observed Q3.}
#'   \item{q3_obs_q005, q3_obs_q995}{0.5\% and 99.5\% quantiles
#'     (99\% credible interval) of the posterior distribution of
#'     observed Q3.}
#'   \item{q3_diff_q025, q3_diff_q975}{2.5\% and 97.5\% quantiles
#'     (95\% credible interval) of the posterior distribution of
#'     Q3 differences.}
#'   \item{q3_diff_q005, q3_diff_q995}{0.5\% and 99.5\% quantiles
#'     (99\% credible interval) of the posterior distribution of
#'     Q3 differences.}
#' }
#'
#' @details
#' The procedure works as follows for each posterior draw \eqn{s}:
#'
#' \enumerate{
#'   \item Compute expected values \eqn{E^{(s)}_{vi}} from the category
#'     probabilities returned by \code{\link[brms]{posterior_epred}}.
#'     For ordinal models: \eqn{E^{(s)}_{vi} = \sum_c c \cdot
#'     P^{(s)}(X_{vi} = c)}. For binary models: \eqn{E^{(s)}_{vi} =
#'     P^{(s)}(X_{vi} = 1)}.
#'   \item Compute observed residuals: \eqn{d^{(s)}_{vi} = X_{vi} -
#'     E^{(s)}_{vi}}.
#'   \item Compute replicated residuals: \eqn{d^{rep(s)}_{vi} =
#'     Y^{rep(s)}_{vi} - E^{(s)}_{vi}}, where \eqn{Y^{rep}} is drawn
#'     via \code{\link[brms]{posterior_predict}}.
#'   \item For each item pair \eqn{(i, j)}, compute Q3 as the Pearson
#'     correlation of residuals across all persons who responded to
#'     both items.
#'   \item Aggregate across draws to obtain posterior means, quantiles,
#'     and posterior predictive p-values.
#' }
#'
#' The key advantage over parametric bootstrapping is that the reference
#' distribution is obtained directly from the posterior, automatically
#' accounting for parameter uncertainty and the negative bias inherent
#' in Q3 (which depends on test length and person ability distribution).
#'
#' @references
#' Yen, W. M. (1984). Effects of local item dependence on the fit and
#' equating performance of the three-parameter logistic model.
#' \emph{Applied Psychological Measurement}, \emph{8}(2), 125--145.
#'
#' Christensen, K. B., Makransky, G. & Horton, M. (2017). Critical
#' values for Yen's Q3: Identification of local dependence in the
#' Rasch model using residual correlations.
#' \emph{Applied Psychological Measurement}, \emph{41}(3), 178--194.
#' \doi{10.1177/0146621616677520}
#'
#' Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS)
#' with Bayesian Item Response Models. \emph{Journal of Intelligence},
#' \emph{8}(1). \doi{10.3390/jintelligence8010005}
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
#' and Stan. \emph{Journal of Statistical Software}, \emph{100}, 1--54.
#' \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link{fit_statistic_pcm}} for posterior predictive fit statistics
#' with user-supplied criterion functions,
#' \code{\link{fit_statistic_rm}} for posterior predictive fit statistics
#' with user-supplied criterion functions,
#' \code{\link{infit_statistic}} for Bayesian infit/outfit,
#' \code{\link[brms]{posterior_epred}},
#' \code{\link[brms]{posterior_predict}}.
#'
#' @examples
#' \donttest{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' # --- Partial Credit Model ---
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
#'   chains = 4,
#'   cores  = 4,
#'   iter   = 2000
#' )
#'
#' # Q3 residual correlations
#' q3_results <- q3_statistic(
#'   model      = fit_pcm,
#'   item_var   = item,
#'   person_var = id,
#'   ndraws_use = 500
#' )
#'
#' # Flag item pairs with ppp > 0.95 as locally dependent
#' q3_results %>%
#'   filter(ppp > 0.95) %>%
#'   arrange(desc(q3_diff))
#'
#' # Inspect 99% credible intervals for Q3 differences
#' q3_results %>%
#'   filter(q3_diff_q005 > 0)
#'
#' # --- Dichotomous Rasch Model ---
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
#'   chains = 4,
#'   cores  = 4,
#'   iter   = 2000
#' )
#'
#' q3_rm <- q3_statistic(
#'   model      = fit_rm,
#'   item_var   = item,
#'   person_var = id,
#'   ndraws_use = 500
#' )
#'
#' q3_rm %>%
#'   filter(ppp > 0.95)
#' }
#'
#' @importFrom brms posterior_epred posterior_predict ndraws
#' @importFrom dplyr group_by summarise arrange filter desc
#' @importFrom rlang enquo !! .data as_name
#' @importFrom tibble as_tibble
#' @importFrom stats formula setNames cor quantile
#' @export
q3_statistic <- function(model, item_var = item, person_var = id,
                         ndraws_use = NULL) {
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }

  item_quo    <- rlang::enquo(item_var)
  person_quo  <- rlang::enquo(person_var)
  item_name   <- rlang::as_name(item_quo)
  person_name <- rlang::as_name(person_quo)

  # --- Extract response variable name ---
  resp_var <- as.character(formula(model)$formula[[2]])
  if (length(resp_var) > 1) {
    resp_var <- resp_var[2]
  }
  if (!resp_var %in% names(model$data)) {
    stop("Response variable '", resp_var, "' not found in model data.",
         call. = FALSE)
  }
  if (!item_name %in% names(model$data)) {
    stop("Item variable '", item_name, "' not found in model data.",
         call. = FALSE)
  }
  if (!person_name %in% names(model$data)) {
    stop("Person variable '", person_name, "' not found in model data.",
         call. = FALSE)
  }

  # --- Determine posterior draw subset ---
  draw_ids <- NULL
  if (!is.null(ndraws_use)) {
    ndraws_use <- as.integer(ndraws_use)
    if (ndraws_use < 1) {
      stop("'ndraws_use' must be a positive integer.", call. = FALSE)
    }
    total_draws <- brms::ndraws(model)
    if (ndraws_use > total_draws) {
      warning("'ndraws_use' (", ndraws_use, ") exceeds available draws (",
              total_draws, "). Using all draws.", call. = FALSE)
      ndraws_use <- total_draws
    }
    draw_ids <- sample(seq_len(total_draws), ndraws_use)
  }

  # --- Posterior predictions ---
  epred_array <- brms::posterior_epred(model, draw_ids = draw_ids)
  yrep_mat    <- brms::posterior_predict(model, draw_ids = draw_ids)

  n_draws <- dim(epred_array)[1]
  n_obs   <- dim(epred_array)[2]
  obs_response <- model$data[[resp_var]]

  # --- Compute expected values per draw x observation ---
  if (length(dim(epred_array)) == 3) {
    n_cat <- dim(epred_array)[3]
    cat_values <- seq_len(n_cat)
    E_mat <- apply(epred_array, c(1, 2), function(p) sum(cat_values * p))
  } else {
    E_mat <- epred_array
  }

  # --- Compute residuals ---
  obs_mat   <- matrix(obs_response, nrow = n_draws, ncol = n_obs, byrow = TRUE)
  resid_obs <- obs_mat - E_mat
  resid_rep <- yrep_mat - E_mat

  # --- Reshape residuals to person x item matrices per draw ---
  items   <- model$data[[item_name]]
  persons <- model$data[[person_name]]
  unique_items   <- sort(unique(items))
  unique_persons <- sort(unique(persons))
  k <- length(unique_items)
  n_persons <- length(unique_persons)

  person_idx <- match(persons, unique_persons)
  item_idx   <- match(items, unique_items)

  # Item pair indices
  pair_grid <- expand.grid(j = seq_len(k), i = seq_len(k))
  pair_grid <- pair_grid[pair_grid$i < pair_grid$j, ]
  pair_grid <- pair_grid[order(pair_grid$i, pair_grid$j), ]
  n_pairs <- nrow(pair_grid)

  q3_obs_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_pairs)
  q3_rep_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_pairs)

  for (s in seq_len(n_draws)) {
    resid_obs_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)
    resid_rep_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)

    for (obs in seq_len(n_obs)) {
      resid_obs_wide[person_idx[obs], item_idx[obs]] <- resid_obs[s, obs]
      resid_rep_wide[person_idx[obs], item_idx[obs]] <- resid_rep[s, obs]
    }

    for (p in seq_len(n_pairs)) {
      i <- pair_grid$i[p]
      j <- pair_grid$j[p]

      valid <- !is.na(resid_obs_wide[, i]) & !is.na(resid_obs_wide[, j])
      if (sum(valid) > 2) {
        q3_obs_draws[s, p] <- stats::cor(
          resid_obs_wide[valid, i], resid_obs_wide[valid, j]
        )
        q3_rep_draws[s, p] <- stats::cor(
          resid_rep_wide[valid, i], resid_rep_wide[valid, j]
        )
      }
    }
  }

  # --- Summarise across draws ---
  q3_diff_draws <- q3_obs_draws - q3_rep_draws

  result <- data.frame(
    item_1 = unique_items[pair_grid$i],
    item_2 = unique_items[pair_grid$j],
    stringsAsFactors = FALSE
  )

  result$q3_obs  <- colMeans(q3_obs_draws, na.rm = TRUE)
  result$q3_rep  <- colMeans(q3_rep_draws, na.rm = TRUE)
  result$q3_diff <- colMeans(q3_diff_draws, na.rm = TRUE)
  result$ppp     <- colMeans(q3_obs_draws > q3_rep_draws, na.rm = TRUE)

  # 95% credible intervals
  result$q3_obs_q025  <- apply(q3_obs_draws, 2, stats::quantile,
                               probs = 0.025, na.rm = TRUE)
  result$q3_obs_q975  <- apply(q3_obs_draws, 2, stats::quantile,
                               probs = 0.975, na.rm = TRUE)
  result$q3_diff_q025 <- apply(q3_diff_draws, 2, stats::quantile,
                               probs = 0.025, na.rm = TRUE)
  result$q3_diff_q975 <- apply(q3_diff_draws, 2, stats::quantile,
                               probs = 0.975, na.rm = TRUE)

  # 99% credible intervals
  result$q3_obs_q005  <- apply(q3_obs_draws, 2, stats::quantile,
                               probs = 0.005, na.rm = TRUE)
  result$q3_obs_q995  <- apply(q3_obs_draws, 2, stats::quantile,
                               probs = 0.995, na.rm = TRUE)
  result$q3_diff_q005 <- apply(q3_diff_draws, 2, stats::quantile,
                               probs = 0.005, na.rm = TRUE)
  result$q3_diff_q995 <- apply(q3_diff_draws, 2, stats::quantile,
                               probs = 0.995, na.rm = TRUE)

  result <- tibble::as_tibble(result)
  result <- dplyr::arrange(result, dplyr::desc(.data$q3_diff))

  result
}
