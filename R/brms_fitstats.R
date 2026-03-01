#' Posterior Predictive Item Fit Statistic for Bayesian IRT Models
#'
#' Computes posterior predictive item (or person) fit statistics for
#' Bayesian IRT models fitted with \pkg{brms}. For each posterior draw,
#' observed and replicated data are compared via a user-supplied criterion
#' function, grouped by item, person, or any other variable. Posterior
#' predictive p-values can then be derived from the output to assess fit.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object.
#' @param criterion A function with signature \code{function(y, p)} that
#'   computes a pointwise fit criterion. For ordinal and categorical models,
#'   \code{y} is the observed (or replicated) response category and \code{p}
#'   is the model-predicted probability of that category. For binary models,
#'   \code{y} is the binary response and \code{p} is the predicted probability
#'   of success.
#' @param group An unquoted variable name (e.g., \code{item} or \code{id})
#'   indicating the grouping variable over which the fit statistic is
#'   aggregated. Typically \code{item} for item fit or \code{id} for person
#'   fit.
#' @param ndraws_use Optional positive integer. If specified, a random subset
#'   of posterior draws of this size is used, which can speed up computation
#'   for large models. If \code{NULL} (the default), all draws are used.
#'
#' @return A \code{\link[tibble]{tibble}} with the following columns:
#' \describe{
#'   \item{\code{group}}{The grouping variable (e.g., item name or person id).}
#'   \item{draw}{Integer index of the posterior draw.}
#'   \item{crit}{The observed fit statistic (criterion applied to observed
#'     data) summed within each group and draw.}
#'   \item{crit_rep}{The replicated fit statistic (criterion applied to
#'     posterior predicted data) summed within each group and draw.}
#'   \item{crit_diff}{The difference \code{crit_rep - crit}.}
#' }
#' The output is grouped by the grouping variable. Posterior predictive
#' p-values can be obtained by computing
#' \code{mean(crit_rep > crit)} within each group.
#'
#' @details
#' The function implements the posterior predictive checking approach for
#' item fit described in Bürkner (2020). The procedure works as follows:
#'
#' \enumerate{
#'   \item Draw posterior expected category probabilities via
#'     \code{\link[brms]{posterior_epred}} and posterior predicted responses
#'     via \code{\link[brms]{posterior_predict}}.
#'   \item For ordinal or categorical models (3D array output from
#'     \code{posterior_epred}), extract the probability assigned to the
#'     observed response category and to the replicated response category
#'     for each draw and observation.
#'   \item Apply the user-supplied \code{criterion} function to compute
#'     pointwise fit values for both observed and replicated data.
#'   \item Aggregate (sum) the criterion values within each level of
#'     \code{group} and each posterior draw.
#' }
#'
#' A common choice for ordinal IRT models is the categorical
#' log-likelihood criterion \code{function(y, p) log(p)}. For binary
#' (e.g., dichotomous Rasch) models, the Bernoulli log-likelihood
#' \code{function(y, p) y * log(p) + (1 - y) * log(1 - p)} may be used
#' instead.
#'
#' @references
#' Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS)
#' with Bayesian Item Response Models. \emph{Journal of Intelligence},
#' \emph{8}(1). \doi{10.3390/jintelligence8010005}
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms and
#' Stan. \emph{Journal of Statistical Software}, \emph{100}, 1--54.
#' \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link{fit_statistic_rm}} for dichotomous Rasch models,
#' \code{\link[brms]{posterior_epred}} for expected predictions,
#' \code{\link[brms]{posterior_predict}} for posterior predictive samples,
#' \code{\link[brms]{pp_check}} for graphical posterior predictive checks.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' # --- Polytomous Rasch (Partial Credit Model) ---
#'
#' # Prepare data in long format
#' df_pcm <- eRm::pcmdat2 %>%
#'   mutate(across(everything(), ~ .x + 1)) %>%
#'   rownames_to_column("id") %>%
#'   pivot_longer(!id, names_to = "item", values_to = "response")
#'
#' # Fit a Partial Credit Model using the adjacent category family
#' fit_pcm <- brm(
#'   response | thres(gr = item) ~ 1 + (1 | id),
#'   data    = df_pcm,
#'   family  = acat,
#'   chains  = 4,
#'   cores   = 4,
#'   iter    = 2000
#' )
#'
#' # Categorical log-likelihood criterion (for polytomous models)
#' ll_categorical <- function(y, p) log(p)
#'
#' # Compute item fit statistics
#' item_fit <- fit_statistic_pcm(
#'   model      = fit_pcm,
#'   criterion  = ll_categorical,
#'   group      = item,
#'   ndraws_use = 500
#' )
#'
#' # Summarise: posterior predictive p-values per item
#' item_fit %>%
#'   group_by(item) %>%
#'   summarise(
#'     observed   = mean(crit),
#'     replicated = mean(crit_rep),
#'     ppp        = mean(crit_rep > crit)
#'   )
#'
#' # Use ggplot2 to make a histogram
#' library(ggplot2)
#' item_fit %>%
#'   ggplot(aes(crit_diff)) +
#'   geom_histogram(aes(fill = ifelse(crit_diff > 0, "above","below"))) +
#'   facet_wrap("item") +
#'   theme_bw() +
#'   theme(legend.position = "none")
#'
#' # Compute person fit statistics
#' person_fit <- fit_statistic_pcm(
#'   model      = fit_pcm,
#'   criterion  = ll_categorical,
#'   group      = id,
#'   ndraws_use = 500
#' )
#'}
#'
#' @importFrom brms posterior_epred posterior_predict ndraws
#' @importFrom dplyr mutate left_join group_by summarise arrange
#' @importFrom tidyr pivot_longer
#' @importFrom rlang enquo !! .data
#' @importFrom stats formula setNames
#' @export
fit_statistic_pcm <- function(model, criterion, group, ndraws_use = NULL) {
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  if (!is.function(criterion)) {
    stop("'criterion' must be a function with signature function(y, p).",
         call. = FALSE)
  }

  group <- rlang::enquo(group)


  # --- Extract response variable name from the model formula ---
  resp_var <- as.character(formula(model)$formula[[2]])
  # For formulas with addition terms (e.g., response | thres(gr = item)),

  # the LHS is a call and as.character() returns a vector; the actual
  # variable name is the second element.
  if (length(resp_var) > 1) {
    resp_var <- resp_var[2]
  }
  if (!resp_var %in% names(model$data)) {
    stop("Response variable '", resp_var, "' not found in model data.",
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

  # --- Posterior expected predictions and predicted responses ---
  epred_array <- brms::posterior_epred(model, draw_ids = draw_ids)
  yrep_mat    <- brms::posterior_predict(model, draw_ids = draw_ids)

  n_draws <- dim(epred_array)[1]
  n_obs   <- dim(epred_array)[2]

  # --- Compute pointwise predicted probabilities ---
  if (length(dim(epred_array)) == 3) {
    # Ordinal/categorical model: epred_array is S x N x C
    obs_response <- model$data[[resp_var]]

    # p(observed category) per draw and observation
    ppe_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
    for (n in seq_len(n_obs)) {
      ppe_mat[, n] <- epred_array[, n, obs_response[n]]
    }

    # p(replicated category) per draw and observation
    yrep_ppe_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
    for (n in seq_len(n_obs)) {
      yrep_ppe_mat[, n] <- epred_array[
        cbind(seq_len(n_draws), n, yrep_mat[, n])
      ]
    }
  } else {
    # Binary/continuous model: epred_array is S x N
    ppe_mat      <- epred_array
    yrep_ppe_mat <- epred_array
  }

  # --- Reshape matrices to long format ---
  obs_names <- seq_len(n_obs)

  ppe_long <- ppe_mat |>
    as.data.frame() |>
    stats::setNames(obs_names) |>
    dplyr::mutate(draw = dplyr::row_number()) |>
    tidyr::pivot_longer(
      -"draw", names_to = ".row", values_to = "ppe"
    ) |>
    dplyr::mutate(.row = as.integer(.data$.row))

  yrep_long <- yrep_mat |>
    as.data.frame() |>
    stats::setNames(obs_names) |>
    dplyr::mutate(draw = dplyr::row_number()) |>
    tidyr::pivot_longer(
      -"draw", names_to = ".row", values_to = "yrep"
    ) |>
    dplyr::mutate(.row = as.integer(.data$.row))

  yrep_ppe_long <- yrep_ppe_mat |>
    as.data.frame() |>
    stats::setNames(obs_names) |>
    dplyr::mutate(draw = dplyr::row_number()) |>
    tidyr::pivot_longer(
      -"draw", names_to = ".row", values_to = "yrep_ppe"
    ) |>
    dplyr::mutate(.row = as.integer(.data$.row))

  # --- Combine with model data ---
  model_data <- model$data |>
    dplyr::mutate(.row = dplyr::row_number())

  result <- ppe_long |>
    dplyr::left_join(yrep_long, by = c("draw", ".row")) |>
    dplyr::left_join(yrep_ppe_long, by = c("draw", ".row")) |>
    dplyr::left_join(model_data, by = ".row") |>
    dplyr::mutate(
      crit     = criterion(.data[[resp_var]], .data$ppe),
      crit_rep = criterion(.data$yrep, .data$yrep_ppe)
    ) |>
    dplyr::group_by(!!group, .data$draw) |>
    dplyr::summarise(
      crit      = sum(.data$crit),
      crit_rep  = sum(.data$crit_rep),
      crit_diff = .data$crit_rep - .data$crit,
      .groups   = "drop_last"
    ) |>
    dplyr::arrange(!!group, .data$draw)

  result
}


#' Posterior Predictive Item Fit Statistic for Binary Bayesian IRT Models
#'
#' Computes posterior predictive item (or person) fit statistics for
#' dichotomous Bayesian IRT models fitted with \pkg{brms}. For each
#' posterior draw, observed and replicated data are compared via a
#' user-supplied criterion function, grouped by item, person, or any other
#' variable. Posterior predictive p-values can then be derived from the
#' output to assess fit.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object with a binary
#'   response (e.g., \code{family = bernoulli()}).
#' @param criterion A function with signature \code{function(y, p)} that
#'   computes a pointwise fit criterion, where \code{y} is the binary
#'   response (0 or 1) and \code{p} is the predicted probability of
#'   success. A common choice is the Bernoulli log-likelihood:
#'   \code{function(y, p) y * log(p) + (1 - y) * log(1 - p)}.
#' @param group An unquoted variable name (e.g., \code{item} or \code{id})
#'   indicating the grouping variable over which the fit statistic is
#'   aggregated. Typically \code{item} for item fit or \code{id} for person
#'   fit.
#' @param ndraws_use Optional positive integer. If specified, a random subset
#'   of posterior draws of this size is used, which can speed up computation
#'   for large models. If \code{NULL} (the default), all draws are used.
#'
#' @return A \code{\link[tibble]{tibble}} with the following columns:
#' \describe{
#'   \item{\code{group}}{The grouping variable (e.g., item name or person id).}
#'   \item{draw}{Integer index of the posterior draw.}
#'   \item{crit}{The observed fit statistic (criterion applied to observed
#'     data) summed within each group and draw.}
#'   \item{crit_rep}{The replicated fit statistic (criterion applied to
#'     posterior predicted data) summed within each group and draw.}
#'   \item{crit_diff}{The difference \code{crit_rep - crit}.}
#' }
#' The output is grouped by the grouping variable. Posterior predictive
#' p-values can be obtained by computing
#' \code{mean(crit_rep > crit)} within each group.
#'
#' @details
#' This function is the binary-response counterpart of
#' \code{\link{fit_statistic_pcm}}, which handles polytomous (ordinal /
#' categorical) models. For dichotomous models, \code{posterior_epred()}
#' returns a 2D matrix (S x N) of success probabilities, so the criterion
#' function receives the observed binary response and the corresponding
#' probability directly.
#'
#' The procedure follows the posterior predictive checking approach
#' described in Bürkner (2020):
#'
#' \enumerate{
#'   \item Draw posterior expected success probabilities via
#'     \code{\link[brms]{posterior_epred}} and posterior predicted binary
#'     responses via \code{\link[brms]{posterior_predict}}.
#'   \item Apply the user-supplied \code{criterion} function pointwise
#'     to both observed and replicated data paired with the predicted
#'     probabilities.
#'   \item Aggregate (sum) the criterion values within each level of
#'     \code{group} and each posterior draw.
#' }
#'
#' The standard criterion for binary models is the Bernoulli log-likelihood:
#' \deqn{\ell(y, p) = y \log(p) + (1 - y) \log(1 - p).}
#'
#' @references
#' Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS)
#' with Bayesian Item Response Models. \emph{Journal of Intelligence},
#' \emph{8}(1). \doi{10.3390/jintelligence8010005}
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms and
#' Stan. \emph{Journal of Statistical Software}, \emph{100}, 1--54.
#' \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link{fit_statistic_pcm}} for polytomous (ordinal/categorical) models,
#' \code{\link[brms]{posterior_epred}} for expected predictions,
#' \code{\link[brms]{posterior_predict}} for posterior predictive samples,
#' \code{\link[brms]{pp_check}} for graphical posterior predictive checks.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' # --- Dichotomous Rasch Model ---
#'
#' # Prepare binary response data in long format
#' df_rm <- eRm::raschdat3 %>%
#'   as.data.frame() %>%
#'   rownames_to_column("id") %>%
#'   pivot_longer(!id, names_to = "item", values_to = "response")
#'
#' # Fit a dichotomous Rasch model
#' fit_rm <- brm(
#'   response ~ 1 + (1 | item) + (1 | id),
#'   data   = df_rm,
#'   family = bernoulli(),
#'   chains = 4,
#'   cores  = 4,
#'   iter   = 2000
#' )
#'
#' # Bernoulli log-likelihood criterion
#' ll_bernoulli <- function(y, p) y * log(p) + (1 - y) * log(1 - p)
#'
#' # Compute item fit statistics
#' item_fit <- fit_statistic_rm(
#'   model      = fit_rm,
#'   criterion  = ll_bernoulli,
#'   group      = item,
#'   ndraws_use = 500
#' )
#'
#' # Summarise: posterior predictive p-values per item
#' item_fit %>%
#'   group_by(item) %>%
#'   summarise(
#'     observed   = mean(crit),
#'     replicated = mean(crit_rep),
#'     ppp        = mean(crit_rep > crit)
#'   )
#'
#' # Use ggplot2 to make a histogram
#' library(ggplot2)
#' item_fit %>%
#'   ggplot(aes(crit_diff)) +
#'   geom_histogram(aes(fill = ifelse(crit_diff > 0, "above","below"))) +
#'   facet_wrap("item") +
#'   theme_bw() +
#'   theme(legend.position = "none")
#'
#' # Compute person fit statistics
#' person_fit <- fit_statistic_rm(
#'   model      = fit_rm,
#'   criterion  = ll_bernoulli,
#'   group      = id,
#'   ndraws_use = 500
#' )
#'
#' person_fit %>%
#'   group_by(id) %>%
#'   summarise(
#'     observed   = mean(crit),
#'     replicated = mean(crit_rep),
#'     ppp        = mean(crit_rep > crit)
#'   )
#'
#' # --- 1PL model with item-specific intercepts ---
#'
#' # Alternative parameterisation with fixed item effects
#' fit_1pl <- brm(
#'   response ~ 0 + item + (1 | id),
#'   data   = df_rm,
#'   family = bernoulli(),
#'   chains = 4,
#'   cores  = 4,
#'   iter   = 2000
#' )
#'
#' item_fit_1pl <- fit_statistic_rm(
#'   model      = fit_1pl,
#'   criterion  = ll_bernoulli,
#'   group      = item,
#'   ndraws_use = 500
#' )
#'
#' item_fit_1pl %>%
#'   group_by(item) %>%
#'   summarise(
#'     observed   = mean(crit),
#'     replicated = mean(crit_rep),
#'     ppp        = mean(crit_rep > crit)
#'   )
#' }
#'
#' @importFrom brms posterior_epred posterior_predict ndraws
#' @importFrom dplyr mutate left_join group_by summarise arrange row_number
#' @importFrom tidyr pivot_longer
#' @importFrom rlang enquo !! .data
#' @importFrom stats formula setNames
#' @export
fit_statistic_rm <- function(model, criterion, group,
                                 ndraws_use = NULL) {
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  if (!is.function(criterion)) {
    stop("'criterion' must be a function with signature function(y, p).",
         call. = FALSE)
  }

  group <- rlang::enquo(group)

  # --- Extract response variable name from the model formula ---
  resp_var <- as.character(formula(model)$formula[[2]])
  if (length(resp_var) > 1) {
    resp_var <- resp_var[2]
  }
  if (!resp_var %in% names(model$data)) {
    stop("Response variable '", resp_var, "' not found in model data.",
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

  # --- Posterior expected predictions and predicted responses ---
  # For binary models, posterior_epred returns an S x N matrix of P(Y = 1)
  ppe_mat  <- brms::posterior_epred(model, draw_ids = draw_ids)
  yrep_mat <- brms::posterior_predict(model, draw_ids = draw_ids)

  if (length(dim(ppe_mat)) == 3) {
    stop("Model appears to be ordinal/categorical (3D posterior_epred). ",
         "Use fit_statistic_pcm() instead.", call. = FALSE)
  }

  n_draws <- nrow(ppe_mat)
  n_obs   <- ncol(ppe_mat)

  # --- Reshape matrices to long format ---
  obs_names <- seq_len(n_obs)

  ppe_long <- ppe_mat |>
    as.data.frame() |>
    stats::setNames(obs_names) |>
    dplyr::mutate(draw = dplyr::row_number()) |>
    tidyr::pivot_longer(
      -"draw", names_to = ".row", values_to = "ppe"
    ) |>
    dplyr::mutate(.row = as.integer(.data$.row))

  yrep_long <- yrep_mat |>
    as.data.frame() |>
    stats::setNames(obs_names) |>
    dplyr::mutate(draw = dplyr::row_number()) |>
    tidyr::pivot_longer(
      -"draw", names_to = ".row", values_to = "yrep"
    ) |>
    dplyr::mutate(.row = as.integer(.data$.row))

  # --- Combine with model data ---
  model_data <- model$data |>
    dplyr::mutate(.row = dplyr::row_number())

  result <- ppe_long |>
    dplyr::left_join(yrep_long, by = c("draw", ".row")) |>
    dplyr::left_join(model_data, by = ".row") |>
    dplyr::mutate(
      crit     = criterion(.data[[resp_var]], .data$ppe),
      crit_rep = criterion(.data$yrep, .data$ppe)
    ) |>
    dplyr::group_by(!!group, .data$draw) |>
    dplyr::summarise(
      crit      = sum(.data$crit),
      crit_rep  = sum(.data$crit_rep),
      crit_diff = .data$crit_rep - .data$crit,
      .groups   = "drop_last"
    ) |>
    dplyr::arrange(!!group, .data$draw)

  result
}