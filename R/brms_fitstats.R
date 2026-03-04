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
#' \donttest{
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


#' Posterior Predictive Infit Statistic for Bayesian IRT Models
#'
#' Computes a Bayesian analogue of the conditional item infit statistic
#' (as described in Christensen, Kreiner & Mesbah, 2013) for Rasch-family
#' models fitted with \pkg{brms}. For each posterior draw, expected values
#' and variances are derived from the category probabilities returned by
#' \code{\link[brms]{posterior_epred}}, and variance-weighted standardised
#' residuals are computed for both observed and replicated data. The result
#' can be summarised into posterior predictive p-values to assess item fit.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from an ordinal
#'   IRT model (e.g., \code{family = acat} for a partial credit model or
#'   \code{family = bernoulli()} for a dichotomous Rasch model).
#' @param item_var An unquoted variable name identifying the item grouping
#'   variable in the model data (e.g., \code{item}).
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data (e.g., \code{id}).
#' @param ndraws_use Optional positive integer. If specified, a random subset
#'   of posterior draws of this size is used. If \code{NULL} (the default),
#'   all draws are used.
#'
#' @return A \code{\link[tibble]{tibble}} with the following columns:
#' \describe{
#'   \item{item}{The item identifier.}
#'   \item{draw}{Integer index of the posterior draw.}
#'   \item{infit}{The observed infit statistic for that item and draw.}
#'   \item{infit_rep}{The replicated infit statistic (based on posterior
#'     predicted data) for that item and draw.}
#'   \item{outfit}{The observed outfit statistic for that item and draw.}
#'   \item{outfit_rep}{The replicated outfit statistic for that item and
#'     draw.}
#' }
#' The output is grouped by the item variable. Posterior predictive
#' p-values can be obtained by computing, e.g.,
#' \code{mean(infit_rep > infit)} within each item.
#'
#' @details
#' The procedure adapts the conditional infit/outfit statistics
#' (Christensen et al., 2013; Kreiner & Christensen, 2011; Müller, 2020) to the
#' Bayesian framework:
#'
#' \enumerate{
#'   \item For each posterior draw \eqn{s}, category probabilities
#'     \eqn{P^{(s)}(X_{vi} = c)} are obtained from
#'     \code{\link[brms]{posterior_epred}}.
#'   \item The conditional expected value and variance for each
#'     observation are computed as:
#'     \deqn{E^{(s)}_{vi} = \sum_c c \cdot P^{(s)}(X_{vi} = c)}
#'     \deqn{Var^{(s)}_{vi} = \sum_c (c - E^{(s)}_{vi})^2 \cdot
#'       P^{(s)}(X_{vi} = c)}
#'   \item Standardised squared residuals are:
#'     \deqn{Z^{2(s)}_{vi} = (X_{vi} - E^{(s)}_{vi})^2 / Var^{(s)}_{vi}}
#'   \item \strong{Outfit} is the unweighted mean of \eqn{Z^{2}_{vi}}
#'     across persons within each item.
#'   \item \strong{Infit} is the variance-weighted mean:
#'     \deqn{Infit^{(s)}_i = \frac{\sum_v Var^{(s)}_{vi} \cdot
#'       Z^{2(s)}_{vi}}{\sum_v Var^{(s)}_{vi}}}
#'   \item The same computations are repeated for replicated data
#'     \eqn{Y^{rep}} drawn via \code{\link[brms]{posterior_predict}}.
#' }
#'
#' Under perfect fit, both infit and outfit have an expected value of 1.
#' Values substantially above 1 indicate underfit (too much noise),
#' values below 1 indicate overfit (too little variation, e.g.,
#' redundancy). Posterior predictive p-values near 0 or 1 indicate
#' misfit.
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
#' Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013).
#' \emph{Rasch Models in Health}. Iste and Wiley, pp. 86--90.
#'
#' Kreiner, S. & Christensen, K. B. (2011). Exact evaluation of Bias in
#' Rasch model residuals. \emph{Advances in Mathematics Research}, 12,
#' 19--40.
#'
#' Müller, M. (2020). Item fit statistics for Rasch analysis: can we trust them?
#' \emph{Journal of Statistical Distributions and Applications}, \emph{7}(1).
#' \doi{10.1186/s40488-020-00108-7}
#'
#' @seealso
#' \code{\link{fit_statistic_pcm}} for a general-purpose posterior predictive
#' fit statistic with user-supplied criterion functions,
#' \code{\link{fit_statistic_rm}} for a general-purpose posterior predictive
#' fit statistic with user-supplied criterion functions,
#' \code{\link[brms]{posterior_epred}},
#' \code{\link[brms]{posterior_predict}},
#' \code{\link[brms]{pp_check}}.
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
#' # Compute infit per item
#' item_infit <- infit_statistic(
#'   model      = fit_pcm,
#'   item_var   = item,
#'   person_var = id,
#'   ndraws_use = 500
#' )
#'
#' # Summarise across draws
#' item_infit %>%
#'   group_by(item) %>%
#'   summarise(
#'     infit_obs = mean(infit),
#'     infit_rep = mean(infit_rep),
#'     infit_ppp = mean(infit_rep > infit)
#'   )
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
#' item_infit_rm <- infit_statistic(
#'   model      = fit_rm,
#'   item_var   = item,
#'   person_var = id,
#'   ndraws_use = 500
#' )
#'
#' item_infit_rm %>%
#'   group_by(item) %>%
#'   summarise(
#'     infit_obs = mean(infit),
#'     infit_rep = mean(infit_rep),
#'     infit_ppp = mean(infit_rep > infit)
#'   )
#' }
#'
#' @importFrom brms posterior_epred posterior_predict ndraws
#' @importFrom dplyr mutate group_by summarise arrange row_number n
#' @importFrom tidyr pivot_longer
#' @importFrom rlang enquo !! .data as_name
#' @importFrom stats formula setNames
#' @export
infit_statistic <- function(model, item_var = item, person_var = id,
                            ndraws_use = NULL) {
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

  # --- Compute E, Var, and Z^2 per draw and observation ---
  if (length(dim(epred_array)) == 3) {
    # Ordinal/categorical: S x N x C array of category probabilities
    n_cat <- dim(epred_array)[3]
    cat_values <- seq_len(n_cat)

    # E[X | params] = sum_c  c * P(X = c)     [S x N]
    E_mat <- apply(epred_array, c(1, 2), function(p) sum(cat_values * p))

    # Var[X | params] = sum_c (c - E)^2 * P(X = c)   [S x N]
    Var_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
    for (s in seq_len(n_draws)) {
      for (n in seq_len(n_obs)) {
        Var_mat[s, n] <- sum((cat_values - E_mat[s, n])^2 *
                               epred_array[s, n, ])
      }
    }
  } else {
    # Binary: S x N matrix with P(Y = 1)
    # E = p, Var = p * (1 - p) for Bernoulli
    E_mat   <- epred_array
    Var_mat <- epred_array * (1 - epred_array)
  }

  # Clamp variance to avoid division by zero
  Var_mat[Var_mat < 1e-12] <- 1e-12

  # Standardised squared residuals for observed data: Z^2_vi
  obs_mat <- matrix(obs_response, nrow = n_draws, ncol = n_obs, byrow = TRUE)
  Z2_obs <- (obs_mat - E_mat)^2 / Var_mat

  # Standardised squared residuals for replicated data
  Z2_rep <- (yrep_mat - E_mat)^2 / Var_mat

  # --- Retrieve item and person identifiers per observation ---
  items   <- model$data[[item_name]]
  persons <- model$data[[person_name]]
  unique_items <- unique(items)

  # --- Compute infit and outfit per item per draw ---
  result_list <- vector("list", length(unique_items))

  for (idx in seq_along(unique_items)) {
    item_label <- unique_items[idx]
    obs_idx <- which(items == item_label)

    # Extract submatrices for this item [S x n_i]
    Z2_obs_i <- Z2_obs[, obs_idx, drop = FALSE]
    Z2_rep_i <- Z2_rep[, obs_idx, drop = FALSE]
    Var_i    <- Var_mat[, obs_idx, drop = FALSE]

    # Outfit = unweighted mean of Z^2 across persons
    outfit_obs <- rowMeans(Z2_obs_i, na.rm = TRUE)
    outfit_rep <- rowMeans(Z2_rep_i, na.rm = TRUE)

    # Infit = variance-weighted mean: sum(Var * Z^2) / sum(Var)
    sum_var <- rowSums(Var_i, na.rm = TRUE)
    sum_var[sum_var < 1e-12] <- 1e-12
    infit_obs <- rowSums(Var_i * Z2_obs_i, na.rm = TRUE) / sum_var
    infit_rep <- rowSums(Var_i * Z2_rep_i, na.rm = TRUE) / sum_var

    result_list[[idx]] <- data.frame(
      item       = item_label,
      draw       = seq_len(n_draws),
      infit      = infit_obs,
      infit_rep  = infit_rep,
      outfit     = outfit_obs,
      outfit_rep = outfit_rep,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, result_list)

  # Rename the item column to match the user's variable name
  names(result)[names(result) == "item"] <- item_name

  result <- tibble::as_tibble(result)
  result <- dplyr::group_by(result, .data[[item_name]])
  result <- dplyr::arrange(result, .data[[item_name]], .data$draw)

  result
}


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
#' Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices
#' (SPM-LS) with Bayesian Item Response Models. \emph{Journal of
#' Intelligence}, \emph{8}(1). \doi{10.3390/jintelligence8010005}
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
#' # Item-restscore association
#' irs <- item_restscore_statistic(
#'   model      = fit_pcm,
#'   item_var   = item,
#'   person_var = id,
#'   ndraws_use = 500
#' )
#'
#' # Flag items with too-strong discrimination (ppp > 0.95)
#' irs %>% filter(ppp > 0.95)
#'
#' # Flag items with too-weak discrimination (ppp < 0.05)
#' irs %>% filter(ppp < 0.05)
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
#'   item_var   = item,
#'   person_var = id,
#'   ndraws_use = 500
#' )
#'
#' irs_rm %>%
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

  # --- Build observed wide matrix (person x item) ---
  obs_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)
  for (obs in seq_len(n_obs)) {
    obs_wide[person_idx[obs], item_idx[obs]] <- obs_response[obs]
  }

  # Observed rest-scores: for item i, rest = rowSum - item_i
  obs_total <- rowSums(obs_wide, na.rm = TRUE)

  # --- Compute gamma for each item and draw ---
  gamma_obs_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  gamma_rep_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)

  for (s in seq_len(n_draws)) {
    # Build replicated wide matrix for draw s
    rep_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)
    for (obs in seq_len(n_obs)) {
      rep_wide[person_idx[obs], item_idx[obs]] <- yrep_mat[s, obs]
    }
    rep_total <- rowSums(rep_wide, na.rm = TRUE)

    for (i in seq_len(k)) {
      # Observed: item score vs rest-score
      item_score_obs <- obs_wide[, i]
      rest_score_obs <- obs_total - item_score_obs

      valid_obs <- !is.na(item_score_obs)
      if (sum(valid_obs) > 2) {
        tab_obs <- table(item_score_obs[valid_obs], rest_score_obs[valid_obs])
        gamma_obs_draws[s, i] <- gk_gamma(tab_obs)
      }

      # Replicated: item score vs rest-score
      item_score_rep <- rep_wide[, i]
      rest_score_rep <- rep_total - item_score_rep

      valid_rep <- !is.na(item_score_rep)
      if (sum(valid_rep) > 2) {
        tab_rep <- table(item_score_rep[valid_rep], rest_score_rep[valid_rep])
        gamma_rep_draws[s, i] <- gk_gamma(tab_rep)
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

  result$gamma_obs  <- colMeans(gamma_obs_draws, na.rm = TRUE)
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

  result
}


#' Goodman-Kruskal Gamma from a Contingency Table
#'
#' Computes the Goodman-Kruskal gamma coefficient for an ordinal
#' contingency table. This is an internal helper function.
#'
#' @param tab A matrix or table representing a contingency table with
#'   rows and columns in ordinal order.
#'
#' @return A scalar: the gamma coefficient in \eqn{[-1, 1]}.
#'
#' @details
#' Gamma is defined as \eqn{(C - D) / (C + D)}, where \eqn{C} is the
#' number of concordant pairs and \eqn{D} is the number of discordant
#' pairs. For a contingency table, concordant pairs arise when both
#' variables increase together, and discordant pairs when one increases
#' while the other decreases.
#'
#' @references
#' Goodman, L. A. & Kruskal, W. H. (1954). Measures of association for
#' cross classifications. \emph{Journal of the American Statistical
#' Association}, \emph{49}(268), 732--764.
#'
#' @keywords internal
gk_gamma <- function(tab) {
  tab <- as.matrix(tab)
  nr <- nrow(tab)
  nc <- ncol(tab)

  if (nr < 2 || nc < 2) {
    return(NA_real_)
  }

  C <- 0  # concordant pairs
  D <- 0  # discordant pairs

  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (tab[i, j] == 0) next
      # Concordant: cells below-right
      if (i < nr && j < nc) {
        C <- C + tab[i, j] * sum(tab[(i + 1):nr, (j + 1):nc, drop = FALSE])
      }
      # Discordant: cells below-left
      if (i < nr && j > 1) {
        D <- D + tab[i, j] * sum(tab[(i + 1):nr, 1:(j - 1), drop = FALSE])
      }
    }
  }

  denom <- C + D
  if (denom == 0) {
    return(NA_real_)
  }

  (C - D) / denom
}
