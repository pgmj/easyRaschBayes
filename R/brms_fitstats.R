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
#' \donttest{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' # --- Polytomous Rasch (Partial Credit Model) ---
#'
#' # Prepare data in long format
#' df_pcm <- phq9[1:200, 1:9] %>%
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
#'   cores   = 1, # use more cores if you have
#'   iter    = 500 # use at least 2000 
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
#'   ndraws_use = 100 # use at least 500
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
#'   ndraws_use = 100 # use at least 500
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
#' \donttest{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' # --- Dichotomous Rasch Model ---
#'
#' # Prepare binary response data in long format
#' # PHQ-9 items dichotomised (any symptom vs none), for illustration only
#' df_rm <- phq9[1:200, 1:9] %>%
#'   mutate(across(everything(), ~ as.integer(.x > 0))) %>%
#'   rownames_to_column("id") %>%
#'   pivot_longer(!id, names_to = "item", values_to = "response")
#'
#' # Fit a dichotomous Rasch model
#' fit_rm <- brm(
#'   response ~ 1 + (1 | item) + (1 | id),
#'   data   = df_rm,
#'   family = bernoulli(),
#'   chains = 4,
#'   cores   = 1, # use more cores if you have
#'   iter    = 500 # use at least 2000 
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
#'   ndraws_use = 100 # use at least 500
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
#'   ndraws_use = 100 # use at least 500
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
#'   cores   = 1, # use more cores if you have
#'   iter    = 500 # use at least 2000 
#' )
#'
#' item_fit_1pl <- fit_statistic_rm(
#'   model      = fit_1pl,
#'   criterion  = ll_bernoulli,
#'   group      = item,
#'   ndraws_use = 100 # use at least 500
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


#' Posterior Predictive Infit Statistic for Bayesian IRT Models
#'
#' Computes item infit (and optionally outfit) for each posterior draw of
#' a Rasch-family model fitted with \pkg{brms}, for the observed data and
#' for replicated data. By default (\code{method = "conditional"}) the
#' expected values and variances are conditional on each person's total
#' score, as in the conditional infit of Christensen et al. (2013) and
#' Kreiner & Christensen (2011), and the replicated data keep each
#' person's observed total score. With \code{method = "theta"} they are
#' computed given each draw's person locations instead.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object. For
#'   \code{method = "conditional"} it must be a Rasch-family model in which
#'   the total score is sufficient for the person location:
#'   \code{family = acat} (partial credit or rating scale model) or
#'   \code{family = bernoulli()} (dichotomous Rasch model), both with the
#'   logit link, a person random intercept \code{(1 | id)}, and no
#'   \code{disc} term. \code{method = "theta"} accepts any ordinal or
#'   Bernoulli model.
#' @param item_var An unquoted variable name identifying the item grouping
#'   variable in the model data (e.g., \code{item}).
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data (e.g., \code{id}).
#' @param ndraws_use Optional positive integer. If specified, a random subset
#'   of posterior draws of this size is used. If \code{NULL} (the default),
#'   all draws are used.
#' @param outfit Logical. If \code{TRUE}, outfit statistics are computed
#'   alongside infit. Default is \code{FALSE} (infit only), since outfit
#'   is highly sensitive to outliers and rarely recommended for Rasch
#'   diagnostics.
#' @param method Either \code{"conditional"} (the default), with
#'   expectations conditional on the total score, or \code{"theta"}, with
#'   expectations given each draw's person locations. \code{"theta"} gives
#'   the same computation as \code{\link{infit_statistic_hpcm}} and as
#'   \code{infit_statistic()} in easyRaschBayes 0.3.0.
#'
#' @return A \code{\link[tibble]{tibble}} with the following columns:
#' \describe{
#'   \item{item}{The item identifier, named after \code{item_var}.}
#'   \item{draw}{Integer index of the posterior draw.}
#'   \item{infit}{The observed infit statistic for that item and draw.}
#'   \item{infit_rep}{The replicated infit statistic for that item and
#'     draw.}
#'   \item{outfit}{(Only if \code{outfit = TRUE}) The observed outfit
#'     statistic for that item and draw.}
#'   \item{outfit_rep}{(Only if \code{outfit = TRUE}) The replicated
#'     outfit statistic for that item and draw.}
#' }
#' The output is grouped by the item variable. Use \code{\link{infit_post}}
#' for posterior predictive p-values and a plot.
#'
#' @details
#' Both methods use the infit and outfit mean squares (Christensen et al.,
#' 2013; Müller, 2020) as discrepancy measures in a posterior predictive
#' check. For each posterior draw \eqn{s}, each observation gets an
#' expected value \eqn{E^{(s)}_{vi}} and a variance \eqn{Var^{(s)}_{vi}},
#' and
#' \deqn{Infit_i^{(s)} = \frac{\sum_v (X_{vi} - E^{(s)}_{vi})^2}
#'   {\sum_v Var^{(s)}_{vi}}, \qquad
#'   Outfit_i^{(s)} = \frac{1}{N_i} \sum_v
#'   \frac{(X_{vi} - E^{(s)}_{vi})^2}{Var^{(s)}_{vi}}.}
#' The replicated statistics use the same \eqn{E} and \eqn{Var} with
#' replicated responses in place of \eqn{X}.
#'
#' \strong{Conditional method.} Under a Rasch-family model the total score
#' \eqn{R_v} is sufficient for the person location. The function takes
#' the category probabilities from \code{\link[brms]{posterior_epred}},
#' removes the person location from them, and computes
#' \eqn{E^{(s)}_{vi} = E(X_{vi} \mid R_v)} and
#' \eqn{Var^{(s)}_{vi} = Var(X_{vi} \mid R_v)} with elementary symmetric
#' functions of that draw's item parameters. Persons with the minimum or
#' maximum total score have zero variance on every item and are left out,
#' as in \code{iarm::out_infit()}. Replicated responses are drawn
#' conditional on each person's observed total score, so the reference
#' distribution does not depend on the person locations. Given the same
#' item parameters, the observed values equal the conditional infit and
#' outfit of \code{iarm::out_infit()}. The observed infit varies across
#' draws only through the uncertainty in the item parameters. With
#' missing responses, the total score is the sum over the items a person
#' answered, scored from 0, and the conditional distribution is computed
#' for each person's set of answered items.
#'
#' \strong{Theta method.} \eqn{E} and \eqn{Var} are computed from the
#' category probabilities given the person location
#' \eqn{\theta_v^{(s)}} of each draw, and replicated responses come from
#' \code{\link[brms]{posterior_predict}}. The person locations are
#' estimated from the same data, which makes posterior predictive
#' p-values of this kind conservative. They tend to cluster around 0.5,
#' so values near 0 or 1 are rarer than p-values of the same size would
#' be (Robins et al., 2000).
#'
#' The conditional method is better calibrated than the theta method, but
#' its p-values are also conservative, because the item parameters are
#' estimated from the same data. In simulations under a fitting Rasch
#' model, two-sided p-values below .05 occurred for 0.5% to 2% of items,
#' with the lowest rates for polytomous items. Misfit can therefore go
#' undetected. The shortfall is largest for the easiest and hardest
#' dichotomous items, so misfit in an item far from the persons' locations
#' is the most likely to be missed. It is suggested to also evaluate item fit with
#' \code{\link{item_restscore_statistic}}, whose p-values were close to
#' uniform in the same simulations.
#'
#' With either method, the posterior predictive p-value is not a measure
#' of the size of misfit. For a fixed amount of misfit it moves toward 0
#' or 1 as the sample size grows.
#'
#' Replication uses random numbers, so call \code{set.seed()} first for
#' reproducible results.
#'
#' @references
#' Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013).
#' \emph{Rasch Models in Health}. Iste and Wiley, pp. 86--90.
#'
#' Kreiner, S. & Christensen, K. B. (2011). Exact evaluation of Bias in
#' Rasch model residuals. \emph{Advances in Mathematics Research}, 12,
#' 19--40.
#'
#' Müller, M. (2020). Item fit statistics for Rasch analysis: can we
#' trust them? \emph{Journal of Statistical Distributions and
#' Applications}, \emph{7}(1).
#' \doi{10.1186/s40488-020-00108-7}
#'
#' Robins, J. M., van der Vaart, A. & Ventura, V. (2000). Asymptotic
#' distribution of P values in composite null models. \emph{Journal of
#' the American Statistical Association}, \emph{95}(452), 1143--1156.
#' \doi{10.1080/01621459.2000.10474310}
#'
#' @seealso
#' \code{\link{infit_post}} for posterior predictive p-values and a plot,
#' \code{\link{item_restscore_statistic}} for the item-restscore
#' association,
#' \code{\link{fit_statistic_pcm}} and \code{\link{fit_statistic_rm}} for
#' general-purpose posterior predictive fit statistics with user-supplied
#' criterion functions,
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
#' # --- Partial Credit Model (polytomous) ---
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
#'   chains = 4,
#'   cores  = 1, # use more cores if you have
#'   iter   = 500 # use at least 2000 
#' )
#'
#' # Compute infit per item
#' item_infit <- infit_statistic(
#'   model      = fit_pcm,
#'   ndraws_use = 100 # use at least 500
#' )
#'
#' # Post-process draws
#' infit_results <- infit_post(item_infit)
#' infit_results$summary
#' infit_results$hdi
#' infit_results$plot
#'
#' # --- Dichotomous Rasch Model ---
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
#'   chains = 4,
#'   cores  = 1, # use more cores if you have
#'   iter   = 500 # use at least 2000 
#' )
#'
#' item_infit_rm <- infit_statistic(
#'   model      = fit_rm,
#'   ndraws_use = 100 # use at least 500
#' )
#' 
#' # Post-process draws
#' infit_results <- infit_post(item_infit_rm)
#' infit_results$summary
#' infit_results$hdi
#' infit_results$plot
#' 
#' }
#'
#' @importFrom brms posterior_epred posterior_predict ndraws
#' @importFrom dplyr mutate group_by summarise arrange row_number n
#' @importFrom tidyr pivot_longer
#' @importFrom rlang enquo !! .data as_name
#' @importFrom stats formula setNames
#' @export
infit_statistic <- function(model, item_var = item, person_var = id,
                            ndraws_use = NULL, outfit = FALSE,
                            method = c("conditional", "theta")) {
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  
  item_quo   <- rlang::enquo(item_var)
  person_quo <- rlang::enquo(person_var)
  item_name   <- rlang::as_name(item_quo)
  person_name <- rlang::as_name(person_quo)
  
  # --- Extract response variable name from model formula ---
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

  method <- match.arg(method)
  if (method == "conditional") {
    return(.infit_conditional(model, resp_var, item_name, person_name,
                              ndraws_use, outfit))
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
  
  # =================================================================
  # Compute E_mat [S x N] and Var_mat [S x N]  — VECTORIZED
  # =================================================================
  if (length(dim(epred_array)) == 3) {
    # Ordinal/categorical: epred_array is S x N x C
    n_cat <- dim(epred_array)[3]
    cat_values <- seq_len(n_cat)  # 1, 2, ..., C
    
    # E[X | params] = sum_c  c * P(X=c)
    # Reshape epred to (S*N) x C, multiply by cat_values, reshape back
    dim_orig <- dim(epred_array)
    ep_2d <- matrix(epred_array, nrow = dim_orig[1] * dim_orig[2],
                    ncol = dim_orig[3])
    E_vec <- ep_2d %*% cat_values          # (S*N) x 1
    E_mat <- matrix(E_vec, nrow = n_draws, ncol = n_obs)
    
    # Var[X | params] = sum_c (c - E)^2 * P(X=c)
    #                 = sum_c c^2 * P(X=c) - E^2   (computational formula)
    E2_vec <- ep_2d %*% (cat_values^2)     # (S*N) x 1
    Var_vec <- E2_vec - E_vec^2
    Var_mat <- matrix(Var_vec, nrow = n_draws, ncol = n_obs)
  } else {
    # Binary: S x N matrix with P(Y = 1)
    E_mat   <- epred_array
    Var_mat <- epred_array * (1 - epred_array)
  }
  
  # Clamp variance to avoid division by zero
  Var_mat[Var_mat < 1e-12] <- 1e-12
  
  # =================================================================
  # Squared residuals: (X - E)^2  — compute once, divide later
  # =================================================================
  obs_mat <- matrix(obs_response, nrow = n_draws, ncol = n_obs, byrow = TRUE)
  resid2_obs <- (obs_mat - E_mat)^2
  resid2_rep <- (yrep_mat - E_mat)^2
  
  # Z^2 = resid^2 / Var  (only needed for outfit)
  # Infit numerator = sum(Var * Z^2) = sum(resid^2)  — no division needed!
  # Infit denominator = sum(Var)
  
  # --- Retrieve item identifiers and pre-compute indices ---
  items <- model$data[[item_name]]
  unique_items <- unique(items)
  k <- length(unique_items)
  
  # Pre-compute column indices for each item (avoids repeated which())
  item_col_idx <- vector("list", k)
  for (idx in seq_len(k)) {
    item_col_idx[[idx]] <- which(items == unique_items[idx])
  }
  
  # =================================================================
  # Compute infit (and optionally outfit) per item per draw
  # =================================================================
  # Pre-allocate output matrices: rows = draws, cols = items
  infit_obs_mat <- matrix(NA_real_, nrow = n_draws, ncol = k)
  infit_rep_mat <- matrix(NA_real_, nrow = n_draws, ncol = k)
  if (outfit) {
    outfit_obs_mat <- matrix(NA_real_, nrow = n_draws, ncol = k)
    outfit_rep_mat <- matrix(NA_real_, nrow = n_draws, ncol = k)
  }
  
  for (idx in seq_len(k)) {
    cols <- item_col_idx[[idx]]
    
    # Infit = sum(resid^2) / sum(Var)
    #       = sum(Var * Z^2) / sum(Var)    [algebraically identical]
    # This avoids computing Z^2 entirely for infit!
    sum_var <- rowSums(Var_mat[, cols, drop = FALSE])
    sum_var[sum_var < 1e-12] <- 1e-12
    
    infit_obs_mat[, idx] <- rowSums(resid2_obs[, cols, drop = FALSE]) /
      sum_var
    infit_rep_mat[, idx] <- rowSums(resid2_rep[, cols, drop = FALSE]) /
      sum_var
    
    if (outfit) {
      # Outfit = mean(Z^2) = mean(resid^2 / Var)
      Z2_obs_i <- resid2_obs[, cols, drop = FALSE] /
        Var_mat[, cols, drop = FALSE]
      Z2_rep_i <- resid2_rep[, cols, drop = FALSE] /
        Var_mat[, cols, drop = FALSE]
      outfit_obs_mat[, idx] <- rowMeans(Z2_obs_i, na.rm = TRUE)
      outfit_rep_mat[, idx] <- rowMeans(Z2_rep_i, na.rm = TRUE)
    }
  }
  
  # =================================================================
  # Assemble output tibble
  # =================================================================
  draw_seq <- seq_len(n_draws)
  item_rep <- rep(unique_items, each = n_draws)
  draw_rep <- rep(draw_seq, times = k)
  
  if (outfit) {
    result <- tibble::tibble(
      item       = item_rep,
      draw       = draw_rep,
      infit      = as.vector(infit_obs_mat),
      infit_rep  = as.vector(infit_rep_mat),
      outfit     = as.vector(outfit_obs_mat),
      outfit_rep = as.vector(outfit_rep_mat)
    )
  } else {
    result <- tibble::tibble(
      item      = item_rep,
      draw      = draw_rep,
      infit     = as.vector(infit_obs_mat),
      infit_rep = as.vector(infit_rep_mat)
    )
  }
  
  # Rename the item column to match the user's variable name
  names(result)[names(result) == "item"] <- item_name
  
  result <- dplyr::group_by(result, .data[[item_name]])
  result <- dplyr::arrange(result, .data[[item_name]], .data$draw)
  
  result
}


#' Item-Restscore Association Conditional on the Total Score
#'
#' Computes the item-restscore association (Kreiner, 2011) for Rasch-family
#' models fitted with \pkg{brms}, with a posterior predictive check that
#' conditions on each person's total score. The observed Goodman-Kruskal
#' gamma between each item and its rest-score (total score minus that item)
#' is compared with (i) its expected value under the model, conditional on
#' the total scores, and (ii) the gammas of datasets replicated conditional
#' on the total scores. Both are computed for each posterior draw of the
#' item parameters.
#'
#' Observed gamma above the expected value means the item discriminates
#' more than the equal discrimination of the Rasch model implies, and
#' observed gamma below it means the item discriminates less.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from a
#'   Rasch-family model in which the total score is sufficient for the
#'   person location: \code{family = acat} (partial credit or rating scale
#'   model) or \code{family = bernoulli()} (dichotomous Rasch model), both
#'   with the logit link, a person random intercept \code{(1 | id)}, and no
#'   \code{disc} term. Other models stop with an error.
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data (e.g., \code{item}).
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data (e.g., \code{id}).
#' @param ndraws_use Optional positive integer. If specified, a random
#'   subset of posterior draws of this size is used. If \code{NULL} (the
#'   default), all draws are used.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{result}{A \code{\link[tibble]{tibble}} with one row per item,
#'     sorted by \code{gamma_diff} in decreasing order, and the columns:
#'     \describe{
#'       \item{item}{The item identifier, named after \code{item_var}.}
#'       \item{gamma_obs}{The observed gamma between the item and the
#'         rest-score.}
#'       \item{gamma_exp}{The expected gamma conditional on the total
#'         scores, as in Kreiner (2011), averaged over posterior draws.}
#'       \item{gamma_diff}{\code{gamma_obs - gamma_exp}.}
#'       \item{ppp}{Mid-p posterior predictive p-value, the proportion
#'         of draws where \code{gamma_obs} exceeds the replicated gamma,
#'         counting ties as one half. Values near 1 indicate stronger
#'         association than the model implies, values near 0 weaker.}
#'       \item{p_adj}{Two-sided posterior predictive p-value adjusted
#'         for the number of items by a single-step max-T procedure. See
#'         Details.}
#'       \item{gamma_rep}{Mean of the replicated gamma across draws.}
#'       \item{gamma_diff_q025, gamma_diff_q975}{95% credible interval
#'         for \code{gamma_obs} minus the expected gamma. It reflects
#'         uncertainty in the item parameters only. See Details.}
#'       \item{gamma_diff_q005, gamma_diff_q995}{99% credible interval
#'         for the same difference.}
#'       \item{gamma_rep_q025, gamma_rep_q975}{95% posterior predictive
#'         interval for the gamma of a replicated dataset.}
#'       \item{gamma_rep_q005, gamma_rep_q995}{99% posterior predictive
#'         interval for the same.}
#'     }
#'   }
#'   \item{draws}{A \code{\link[tibble]{tibble}} in long format with one
#'     row per item and draw, and the columns \code{item} (named after
#'     \code{item_var}), \code{draw}, \code{gamma_rep} (replicated gamma)
#'     and \code{gamma_exp} (expected gamma).}
#' }
#'
#' @details
#' Under a Rasch-family model the total score is sufficient for the
#' person location. Given the total score, the distribution of a
#' person's responses depends on the item parameters only. For each
#' posterior draw \eqn{s} the function:
#'
#' \enumerate{
#'   \item Takes the category probabilities from
#'     \code{\link[brms]{posterior_epred}} and removes the person
#'     location from them, which leaves item weights that depend on
#'     that draw's item parameters only.
#'   \item Computes, with elementary symmetric functions, the
#'     probability of each item score given each person's total score.
#'     Summed over persons, these give the expected item \eqn{\times}
#'     rest-score table and its gamma, \eqn{\gamma^{exp(s)}_i}.
#'   \item Draws a replicated dataset in which every person keeps their
#'     observed total score, and computes the replicated gamma
#'     \eqn{\gamma^{rep(s)}_i}.
#' }
#'
#' The observed gamma \eqn{\gamma^{obs}_i} is computed once from the data.
#' Persons with the minimum or maximum total score carry no information
#' about fit and their responses are reproduced as observed in every
#' replicated dataset.
#'
#' Because replication conditions on the total scores, the reference
#' distribution does not depend on the person locations. The data are
#' still used to estimate the item parameters, but these are few relative
#' to the number of persons. This removes the main reason why posterior
#' predictive p-values based on person parameters drawn from the same
#' data are conservative (Robins et al., 2000).
#'
#' The \code{gamma_diff} interval describes how far the observed sample
#' departs from the model expectation, given the uncertainty in the item
#' parameters. It contains no sampling variation, so it is not a test, and
#' it is not a credible interval for misfit in the population. Use
#' \code{ppp} or \code{p_adj} to judge whether a departure is larger
#' than chance would produce. As with any test, they will flag smaller
#' departures as the sample size grows.
#'
#' \code{p_adj} standardises each item's observed gamma by the mean and
#' standard deviation of its replicated gammas. Its reference is the
#' largest absolute standardised deviation across items within each
#' replicated dataset, so it accounts for testing all items at once.
#'
#' With missing responses, the total and rest-scores are sums over the
#' items a person answered, scored from 0. The conditional distribution
#' is computed for each person's set of answered items.
#'
#' Replication uses random numbers, so call \code{set.seed()} first for
#' reproducible results.
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
#' Robins, J. M., van der Vaart, A. & Ventura, V. (2000). Asymptotic
#' distribution of P values in composite null models. \emph{Journal of
#' the American Statistical Association}, \emph{95}(452), 1143--1156.
#' \doi{10.1080/01621459.2000.10474310}
#'
#' @seealso
#' \code{\link{item_restscore_post}} for a summary table and plot,
#' \code{\link{infit_statistic}} for posterior predictive infit/outfit,
#' \code{\link{q3_statistic}} for Bayesian Q3 residual correlations,
#' \code{\link[brms]{posterior_epred}}.
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
#' df_pcm <- phq9[1:200, 1:9] %>%
#'   mutate(across(everything(), ~ .x + 1)) %>%
#'   rownames_to_column("id") %>%
#'   pivot_longer(!id, names_to = "item", values_to = "response")
#'
#' fit_pcm <- brm(
#'   response | thres(gr = item) ~ 1 + (1 | id),
#'   data   = df_pcm,
#'   family = acat,
#'   chains = 4,
#'   cores  = 1, # use more cores if you have
#'   iter   = 500 # use at least 2000 
#' )
#'
#' # Item-restscore association
#' irs <- item_restscore_statistic(
#'   model      = fit_pcm,
#'   ndraws_use = 100 # use at least 500
#' )
#' 
#' # Post-process draws
#' irs_results <- item_restscore_post(irs)
#' irs_results$summary
#' irs_results$plot
#'
#' # --- Dichotomous Rasch Model ---
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
#'   chains = 4,
#'   cores  = 1, # use more cores if you have
#'   iter   = 500 # use at least 2000 
#' )
#'
#' irs_rm <- item_restscore_statistic(
#'   model      = fit_rm,
#'   ndraws_use = 100 # use at least 500
#' )
#'
#' # Post-process draws
#' irs_results <- item_restscore_post(irs_rm)
#' irs_results$summary
#' irs_results$plot
#' }
#'
#' @importFrom brms posterior_epred ndraws variables as_draws_matrix
#' @importFrom rlang enquo !! .data as_name
#' @importFrom tibble as_tibble tibble
#' @importFrom stats formula quantile runif sd
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

  .check_unique_person_item(model$data[[person_name]], model$data[[item_name]],
                            person_name, item_name, "item_restscore_statistic")
  .check_sufficient_score_brmsfit(model, person_name,
                                  "item_restscore_statistic")

  draw_ids <- .select_draw_ids(model, ndraws_use)
  n_draws <- length(draw_ids)

  # --- Map observations to person x item structure ---
  dd <- .cond_prepare_data(model, resp_var, item_name, person_name)
  unique_items <- dd$unique_items
  x_obs <- dd$x_obs
  k <- length(unique_items)
  r_obs <- rowSums(x_obs, na.rm = TRUE)

  # --- Observed gamma, computed once ---
  gamma_obs <- numeric(k)
  for (i in seq_len(k)) {
    valid <- !is.na(x_obs[, i])
    gamma_obs[i] <- if (sum(valid) > 2L) {
      gk_gamma_vec(x_obs[valid, i], r_obs[valid] - x_obs[valid, i])
    } else {
      NA_real_
    }
  }

  # --- Expected and replicated gamma per draw ---
  per_draw <- .cond_map_draws(model, draw_ids, dd$n_cat_max,
                              .cond_restscore_draw,
                              person_idx = dd$person_idx,
                              item_idx = dd$item_idx, x_obs = x_obs)
  gamma_exp_draws <- do.call(rbind, lapply(per_draw, `[[`, "gamma_exp"))
  gamma_rep_draws <- do.call(rbind, lapply(per_draw, `[[`, "gamma_rep"))

  # --- Summarise across draws ---
  diff_draws <- matrix(gamma_obs, nrow = n_draws, ncol = k, byrow = TRUE) -
    gamma_exp_draws
  q <- function(m, p) apply(m, 2L, stats::quantile, probs = p, na.rm = TRUE)

  result <- tibble::tibble(item = unique_items)
  names(result)[1] <- item_name
  result$gamma_obs       <- gamma_obs
  result$gamma_exp       <- colMeans(gamma_exp_draws, na.rm = TRUE)
  result$gamma_diff      <- gamma_obs - result$gamma_exp
  result$ppp             <- .ppp_mid(gamma_obs, gamma_rep_draws)
  result$p_adj         <- .ppp_maxt(gamma_obs, gamma_rep_draws)
  result$gamma_rep       <- colMeans(gamma_rep_draws, na.rm = TRUE)
  result$gamma_diff_q025 <- q(diff_draws, 0.025)
  result$gamma_diff_q975 <- q(diff_draws, 0.975)
  result$gamma_diff_q005 <- q(diff_draws, 0.005)
  result$gamma_diff_q995 <- q(diff_draws, 0.995)
  result$gamma_rep_q025  <- q(gamma_rep_draws, 0.025)
  result$gamma_rep_q975  <- q(gamma_rep_draws, 0.975)
  result$gamma_rep_q005  <- q(gamma_rep_draws, 0.005)
  result$gamma_rep_q995  <- q(gamma_rep_draws, 0.995)
  result <- result[order(-result$gamma_diff), ]

  draws_df <- tibble::tibble(
    item      = rep(unique_items, each = n_draws),
    draw      = rep(seq_len(n_draws), times = k),
    gamma_rep = as.vector(gamma_rep_draws),
    gamma_exp = as.vector(gamma_exp_draws)
  )
  names(draws_df)[1] <- item_name

  list(result = result, draws = draws_df)
}
