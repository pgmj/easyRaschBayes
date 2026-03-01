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
#' @param outfit Logical. If \code{TRUE}, outfit statistics are computed
#'   alongside infit. Default is \code{FALSE} (infit only), since outfit
#'   is highly sensitive to outliers and rarely recommended for Rasch
#'   diagnostics.
#'
#' @return A \code{\link[tibble]{tibble}} with the following columns:
#' \describe{
#'   \item{item}{The item identifier.}
#'   \item{draw}{Integer index of the posterior draw.}
#'   \item{infit}{The observed infit statistic for that item and draw.}
#'   \item{infit_rep}{The replicated infit statistic (based on posterior
#'     predicted data) for that item and draw.}
#'   \item{outfit}{(Only if \code{outfit = TRUE}) The observed outfit
#'     statistic for that item and draw.}
#'   \item{outfit_rep}{(Only if \code{outfit = TRUE}) The replicated
#'     outfit statistic for that item and draw.}
#' }
#' The output is grouped by the item variable. Posterior predictive
#' p-values can be obtained by computing, e.g.,
#' \code{mean(infit_rep > infit)} within each item.
#'
#' @details
#' The procedure adapts the conditional infit/outfit statistics
#' (Christensen et al., 2013; Kreiner & Christensen, 2011; Müller, 2020)
#' to the Bayesian framework:
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
#'   \item The infit statistic for item \eqn{i} is the variance-weighted
#'     mean of \eqn{Z^2} across persons:
#'     \deqn{Infit_i^{(s)} = \frac{\sum_v Var_{vi}^{(s)} Z^{2(s)}_{vi}}
#'       {\sum_v Var_{vi}^{(s)}}}
#'   \item If requested, the outfit is the unweighted mean of \eqn{Z^2}.
#' }
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
#' @seealso
#' \code{\link{fit_statistic_pcm}} for a general-purpose posterior predictive
#' fit statistic with user-supplied criterion functions,
#' \code{\link{fit_statistic_rm}} for a general-purpose posterior predictive
#' fit statistic with user-supplied criterion functions,
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
#' item_infit_pcm <- infit_statistic(
#'   model      = fit_pcm,
#'   item_var   = item,
#'   person_var = id,
#'   ndraws_use = 500
#' )
#'
#' item_infit_pcm %>%
#'   group_by(item) %>%
#'   summarise(
#'     infit_obs = mean(infit),
#'     infit_rep = mean(infit_rep),
#'     infit_ppp = mean(infit_rep > infit)
#'   )
#'
#' # Including outfit
#' item_fit_full <- infit_statistic(
#'   model      = fit_pcm,
#'   item_var   = item,
#'   person_var = id,
#'   ndraws_use = 500,
#'   outfit     = TRUE
#' )
#' }
#'
#' @importFrom brms posterior_epred posterior_predict ndraws
#' @importFrom dplyr mutate group_by summarise arrange row_number n
#' @importFrom tidyr pivot_longer
#' @importFrom rlang enquo !! .data as_name
#' @importFrom stats formula setNames
#' @export
infit_statistic <- function(model, item_var = item, person_var = id,
                            ndraws_use = NULL, outfit = FALSE) {
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