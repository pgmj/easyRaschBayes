#' Posterior Predictive Q3 Residual Correlations for Bayesian IRT Models
#'
#' Computes a Bayesian analogue of Yen's Q3 statistic (Yen, 1984) for
#' detecting local dependence between item pairs in Rasch-family models
#' fitted with \pkg{brms}. For each posterior draw, residual correlations
#' are computed for both observed and replicated data, yielding draw-level
#' Q3 values that can be summarized and visualized via
#' \code{\link{q3_post}}.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from an
#'   ordinal IRT model (e.g., \code{family = acat} for a partial credit
#'   model) or a dichotomous model (\code{family = bernoulli()}).
#' @param item_var An unquoted variable name identifying the item grouping
#'   variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Default is \code{id}.
#' @param ndraws_use Optional positive integer. If specified, a random
#'   subset of posterior draws of this size is used. If \code{NULL} (the
#'   default), all draws are used.
#'
#' @return A \code{\link[tibble]{tibble}} in long format with one row
#'   per draw per item pair, containing:
#' \describe{
#'   \item{draw}{Integer draw index.}
#'   \item{item_pair}{Character label of the item pair
#'     (\code{"item1 : item2"}).}
#'   \item{item_1}{First item in the pair.}
#'   \item{item_2}{Second item in the pair.}
#'   \item{q3}{Observed Q3 residual correlation for this draw.}
#'   \item{q3_rep}{Replicated Q3 residual correlation for this draw.}
#' }
#'
#' This long-format output parallels the structure of
#' \code{\link{infit_statistic}} and can be passed directly to
#' \code{\link{q3_post}} for summary tables and plots.
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
#' }
#'
#' @references
#' Yen, W. M. (1984). Effects of local item dependence on the fit and
#' equating performance of the three-parameter logistic model.
#' \emph{Applied Psychological Measurement}, \emph{8}(2), 125--145.
#' \doi{10.1177/014662168400800201}
#'
#' Christensen, K. B., Makransky, G. & Horton, M. (2017). Critical
#' values for Yen's Q3: Identification of local dependence in the
#' Rasch model using residual correlations.
#' \emph{Applied Psychological Measurement}, \emph{41}(3), 178--194.
#' \doi{10.1177/0146621616677520}
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with
#' brms and Stan. \emph{Journal of Statistical Software}, \emph{100},
#' 1--54. \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link{q3_post}} for postprocessing summaries and plots,
#' \code{\link{infit_statistic}} for Bayesian infit/outfit,
#' \code{\link[brms]{posterior_epred}},
#' \code{\link[brms]{posterior_predict}}.
#'
#' @examples
#' \dontrun{
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
#'   chains = 4, cores = 4, iter = 2000
#' )
#'
#' q3_draws <- q3_statistic(fit_pcm, ndraws_use = 500)
#'
#' # Postprocess
#' result <- q3_post(q3_draws)
#' result$summary
#' result$hdi
#' result$plot
#' }
#'
#' @importFrom brms posterior_epred posterior_predict ndraws
#' @importFrom rlang enquo as_name .data
#' @importFrom stats formula cor
#' @importFrom tibble tibble
#' @export
q3_statistic <- function(model, item_var = item, person_var = id,
                         ndraws_use = NULL) {
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }

  item_name   <- rlang::as_name(rlang::enquo(item_var))
  person_name <- rlang::as_name(rlang::enquo(person_var))

  resp_var <- as.character(formula(model)$formula[[2]])
  if (length(resp_var) > 1) resp_var <- resp_var[2]

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
  total_draws <- brms::ndraws(model)
  if (!is.null(ndraws_use)) {
    ndraws_use <- as.integer(ndraws_use)
    if (ndraws_use < 1) {
      stop("'ndraws_use' must be a positive integer.", call. = FALSE)
    }
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

  # --- Expected values ---
  if (length(dim(epred_array)) == 3) {
    n_cat <- dim(epred_array)[3]
    cat_values <- seq_len(n_cat)
    dim_orig <- dim(epred_array)
    ep_2d <- matrix(epred_array, nrow = dim_orig[1] * dim_orig[2],
                    ncol = dim_orig[3])
    E_vec <- ep_2d %*% cat_values
    E_mat <- matrix(E_vec, nrow = n_draws, ncol = n_obs)
  } else {
    E_mat <- epred_array
  }

  # --- Residuals ---
  obs_mat   <- matrix(obs_response, nrow = n_draws, ncol = n_obs,
                      byrow = TRUE)
  resid_obs <- obs_mat - E_mat
  resid_rep <- yrep_mat - E_mat

  # --- Person x item mapping ---
  items   <- model$data[[item_name]]
  persons <- model$data[[person_name]]
  unique_items   <- sort(unique(items))
  unique_persons <- sort(unique(persons))
  k <- length(unique_items)
  n_persons <- length(unique_persons)

  person_idx <- match(persons, unique_persons)
  item_idx   <- match(items, unique_items)
  lin_idx <- person_idx + (item_idx - 1L) * n_persons

  # --- Item pair indices (upper triangle) ---
  pair_grid <- expand.grid(j = seq_len(k), i = seq_len(k))
  pair_grid <- pair_grid[pair_grid$i < pair_grid$j, ]
  pair_grid <- pair_grid[order(pair_grid$i, pair_grid$j), ]
  n_pairs <- nrow(pair_grid)

  pair_labels <- paste0(unique_items[pair_grid$i], " : ",
                        unique_items[pair_grid$j])

  # --- Q3 per draw ---
  q3_obs_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_pairs)
  q3_rep_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_pairs)

  for (s in seq_len(n_draws)) {
    resid_obs_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)
    resid_rep_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)
    resid_obs_wide[lin_idx] <- resid_obs[s, ]
    resid_rep_wide[lin_idx] <- resid_rep[s, ]

    cor_obs <- stats::cor(resid_obs_wide, use = "pairwise.complete.obs")
    cor_rep <- stats::cor(resid_rep_wide, use = "pairwise.complete.obs")

    for (p in seq_len(n_pairs)) {
      q3_obs_draws[s, p] <- cor_obs[pair_grid$i[p], pair_grid$j[p]]
      q3_rep_draws[s, p] <- cor_rep[pair_grid$i[p], pair_grid$j[p]]
    }
  }

  # --- Build long-format output ---
  result_list <- vector("list", n_pairs)
  for (p in seq_len(n_pairs)) {
    result_list[[p]] <- tibble::tibble(
      draw      = seq_len(n_draws),
      item_pair = pair_labels[p],
      item_1    = unique_items[pair_grid$i[p]],
      item_2    = unique_items[pair_grid$j[p]],
      q3        = q3_obs_draws[, p],
      q3_rep    = q3_rep_draws[, p]
    )
  }

  result <- do.call(rbind, result_list)
  tibble::as_tibble(result)
}
