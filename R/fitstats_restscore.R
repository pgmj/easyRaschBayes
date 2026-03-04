#' Posterior Predictive Item-Restscore Association for Bayesian IRT Models
#'
#' Computes a Bayesian analogue of the item-restscore association test
#' (Kreiner, 2011) for Rasch-family models fitted with \pkg{brms}. For
#' each posterior draw, the Goodman-Kruskal gamma coefficient between
#' each item's score and the rest-score (total score minus that item)
#' is computed for both observed and replicated data. Posterior
#' predictive p-values indicate whether the observed association is
#' stronger than the model predicts, which signals violations of
#' local independence or unidimensionality.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from an
#'   ordinal IRT model (e.g., \code{family = acat} for a partial credit
#'   model) or a dichotomous model (\code{family = bernoulli()}).
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data (e.g., \code{item}).
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data (e.g., \code{id}).
#' @param ndraws_use Optional positive integer. If specified, a random
#'   subset of posterior draws of this size is used. If \code{NULL} (the
#'   default), all draws are used.
#'
#' @return A \code{\link[tibble]{tibble}} with the following columns:
#' \describe{
#'   \item{item}{The item identifier.}
#'   \item{gamma_obs}{Posterior mean of the observed Goodman-Kruskal
#'     gamma between this item and the rest-score.}
#'   \item{gamma_rep}{Posterior mean of the replicated gamma.}
#'   \item{gamma_diff}{Posterior mean of \code{gamma_obs - gamma_rep}.
#'     Positive values indicate the observed item-restscore association
#'     is stronger than the model expects.}
#'   \item{ppp}{Posterior predictive p-value:
#'     \code{mean(gamma_obs > gamma_rep)} across draws.
#'     Values close to 1 indicate the item discriminates more than the
#'     model predicts (too high discrimination). Values close to 0
#'     indicate the item discriminates less than expected (too low
#'     discrimination, e.g., noise or miskeyed item).}
#'   \item{gamma_obs_q025, gamma_obs_q975}{95\% credible interval for
#'     the observed gamma.}
#'   \item{gamma_obs_q005, gamma_obs_q995}{99\% credible interval for
#'     the observed gamma.}
#'   \item{gamma_diff_q025, gamma_diff_q975}{95\% credible interval for
#'     the gamma difference.}
#'   \item{gamma_diff_q005, gamma_diff_q995}{99\% credible interval for
#'     the gamma difference.}
#' }
#'
#' @details
#' The item-restscore association is a key diagnostic in Rasch
#' measurement. Under the Rasch model, each item should relate to the
#' latent trait (and hence the rest-score) only through the modelled
#' relationship. Goodman-Kruskal's gamma is a rank-based measure of
#' association for ordinal cross-tabulations that is well-suited for
#' this purpose (Kreiner, 2011).
#'
#' The procedure for each posterior draw \eqn{s} is:
#'
#' \enumerate{
#'   \item Obtain replicated responses \eqn{Y^{rep(s)}} from
#'     \code{\link[brms]{posterior_predict}}.
#'   \item For each item \eqn{i} and each person \eqn{v}, compute the
#'     rest-score: \eqn{R^{obs}_{vi} = \sum_{j \neq i} X_{vj}} for
#'     observed data and \eqn{R^{rep(s)}_{vi} = \sum_{j \neq i}
#'     Y^{rep(s)}_{vj}} for replicated data.
#'   \item Cross-tabulate item score \eqn{\times} rest-score and compute
#'     the Goodman-Kruskal gamma for both observed and replicated data.
#'   \item Compare the two gammas across draws.
#' }
#'
#' Items with \code{ppp} close to 1 have observed item-restscore
#' association that is consistently stronger than the model predicts.
#' This typically indicates that the item discriminates more than
#' assumed under the equal-discrimination Rasch model (i.e., a
#' violation of the Rasch assumption). Items with \code{ppp} close to
#' 0 discriminate less than expected.
#'
#' @references
#' Kreiner, S. (2011). A note on item-restscore association in Rasch
#' models. \emph{Applied Psychological Measurement}, \emph{35}(7),
#' 557--561.
#'
#' Goodman, L. A. & Kruskal, W. H. (1954). Measures of association for
#' cross classifications. \emph{Journal of the American Statistical
#' Association}, \emph{49}(268), 732--764.
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with
#' brms and Stan. \emph{Journal of Statistical Software}, \emph{100},
#' 1--54. \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link{fit_statistic_pcm}} for posterior predictive fit statistics,
#' \code{\link{fit_statistic_rm}} for posterior predictive fit statistics,
#' \code{\link{infit_statistic}} for Bayesian infit/outfit,
#' \code{\link{q3_statistic}} for Bayesian Q3 residual correlations,
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
#'   chains = 4,
#'   cores  = 4,
#'   iter   = 2000
#' )
#'
#' # Item-restscore association
#' irs <- item_restscore_statistic(
#'   model      = fit_pcm,
#'   ndraws_use = 500
#' )
#'
#' # Flag items with too-strong discrimination (ppp > 0.95)
#' irs$result %>% filter(ppp > 0.95)
#'
#' # Flag items with too-weak discrimination (ppp < 0.05)
#' irs$result %>% filter(ppp < 0.05)
#'
#' # --- Dichotomous Rasch Model ---
#'
#' df_rm <- eRm::rainger %>%
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
#' irs_rm <- item_restscore_statistic(
#'   model      = fit_rm,
#'   ndraws_use = 500
#' )
#'
#' irs_rm$result %>%
#'   arrange(ppp)
#' }
#'
#' @importFrom brms posterior_predict ndraws
#' @importFrom dplyr arrange desc filter
#' @importFrom rlang enquo !! .data as_name
#' @importFrom tibble as_tibble
#' @importFrom stats formula quantile
#' @export
item_restscore_statistic <- function(model, item_var = item, person_var = id,
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

  # --- Posterior predicted responses ---
  yrep_mat <- brms::posterior_predict(model, draw_ids = draw_ids)

  n_draws <- nrow(yrep_mat)
  n_obs   <- ncol(yrep_mat)
  obs_response <- model$data[[resp_var]]

  # --- Map observations to person x item structure ---
  items   <- model$data[[item_name]]
  persons <- model$data[[person_name]]
  unique_items   <- sort(unique(items))
  unique_persons <- sort(unique(persons))
  k <- length(unique_items)
  n_persons <- length(unique_persons)

  person_idx <- match(persons, unique_persons)
  item_idx   <- match(items, unique_items)

  # Pre-compute linear index for vectorized wide-matrix construction
  lin_idx <- (item_idx - 1L) * n_persons + person_idx

  # --- Build observed wide matrix (person x item) ONCE ---
  obs_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)
  obs_wide[lin_idx] <- obs_response
  obs_total <- rowSums(obs_wide, na.rm = TRUE)

  # --- Compute observed gamma ONCE (constant across draws) ---
  gamma_obs <- numeric(k)
  for (i in seq_len(k)) {
    item_score_obs <- obs_wide[, i]
    rest_score_obs <- obs_total - item_score_obs
    valid_obs <- !is.na(item_score_obs)
    if (sum(valid_obs) > 2L) {
      gamma_obs[i] <- gk_gamma_vec(item_score_obs[valid_obs],
                                   rest_score_obs[valid_obs])
    } else {
      gamma_obs[i] <- NA_real_
    }
  }

  # Broadcast observed gamma to a draws x items matrix (for summary stats)
  gamma_obs_draws <- matrix(gamma_obs, nrow = n_draws, ncol = k, byrow = TRUE)

  # --- Compute replicated gamma for each draw ---
  gamma_rep_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)

  for (s in seq_len(n_draws)) {
    # Vectorized wide-matrix construction
    rep_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)
    rep_wide[lin_idx] <- yrep_mat[s, ]
    rep_total <- rowSums(rep_wide, na.rm = TRUE)

    for (i in seq_len(k)) {
      item_score_rep <- rep_wide[, i]
      rest_score_rep <- rep_total - item_score_rep
      valid_rep <- !is.na(item_score_rep)
      if (sum(valid_rep) > 2L) {
        gamma_rep_draws[s, i] <- gk_gamma_vec(item_score_rep[valid_rep],
                                              rest_score_rep[valid_rep])
      }
    }
  }

  # --- Summarise across draws ---
  gamma_diff_draws <- gamma_obs_draws - gamma_rep_draws

  result <- data.frame(
    item = unique_items,
    stringsAsFactors = FALSE
  )
  names(result)[1] <- item_name

  result$gamma_obs  <- gamma_obs
  result$gamma_rep  <- colMeans(gamma_rep_draws, na.rm = TRUE)
  result$gamma_diff <- colMeans(gamma_diff_draws, na.rm = TRUE)
  result$ppp        <- colMeans(gamma_obs_draws > gamma_rep_draws,
                                na.rm = TRUE)

  # 95% credible intervals
  result$gamma_obs_q025  <- apply(gamma_obs_draws, 2, stats::quantile,
                                  probs = 0.025, na.rm = TRUE)
  result$gamma_obs_q975  <- apply(gamma_obs_draws, 2, stats::quantile,
                                  probs = 0.975, na.rm = TRUE)
  result$gamma_diff_q025 <- apply(gamma_diff_draws, 2, stats::quantile,
                                  probs = 0.025, na.rm = TRUE)
  result$gamma_diff_q975 <- apply(gamma_diff_draws, 2, stats::quantile,
                                  probs = 0.975, na.rm = TRUE)

  # 99% credible intervals
  result$gamma_obs_q005  <- apply(gamma_obs_draws, 2, stats::quantile,
                                  probs = 0.005, na.rm = TRUE)
  result$gamma_obs_q995  <- apply(gamma_obs_draws, 2, stats::quantile,
                                  probs = 0.995, na.rm = TRUE)
  result$gamma_diff_q005 <- apply(gamma_diff_draws, 2, stats::quantile,
                                  probs = 0.005, na.rm = TRUE)
  result$gamma_diff_q995 <- apply(gamma_diff_draws, 2, stats::quantile,
                                  probs = 0.995, na.rm = TRUE)

  result <- tibble::as_tibble(result)
  result <- dplyr::arrange(result, dplyr::desc(.data$gamma_diff))

  gamma_rep_draws_df <- gamma_rep_draws %>%
    as.data.frame() %>%
    tidyr::pivot_longer(everything(), names_to = "item", values_to = "gamma_rep")
  gamma_obs_draws_df <- gamma_obs_draws %>%
    as.data.frame() %>%
    tidyr::pivot_longer(everything(), names_to = "item", values_to = "gamma_obs")

  gamma_draws_df <- dplyr::bind_cols(
    gamma_rep_draws_df,
    gamma = gamma_obs_draws_df$gamma_obs
  )

  list(
    result = result,
    draws = tidyr::as_tibble(gamma_draws_df)
  )
}

