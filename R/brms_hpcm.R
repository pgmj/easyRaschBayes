# ============================================================================
# Hurdle Partial Credit Model (hPCM)
# ============================================================================
#
# Custom `brms` family for a hurdle Rasch partial credit model:
#   * a Bernoulli hurdle on `1[Y == 0]` (P(Y = 0) = hu)
#   * an acat-logit partial credit model on the positive categories
#     1..K-1 (conditional on Y > 0)
#
# The custom family code (`hurdle_acat`, `hurdle_acat_stanvars`, and the
# `log_lik`, `posterior_predict`, `posterior_epred` helpers) was
# contributed by Kristoffer Magnusson and is included here with his
# permission. The companion posterior-predictive fit-diagnostic
# functions (`infit_statistic_hpcm`, `q3_statistic_hpcm`) were added by
# the easyRaschBayes maintainers.
#
# Motivation (psychometric):
#   For mental-health / addiction scales, the lowest response category
#   ("not at all") often reflects a different process from the higher
#   categories: presence/absence of the symptom vs. severity of an
#   already-present symptom. This is the structure described in
#   Magnus & Garnier-Villarreal (2022). The hurdle PCM is the
#   PCM-equivalent of their MZI-GRM: a separate "susceptibility"
#   (hurdle) trait governs Y = 0 vs. Y > 0, and a separate "severity"
#   trait (correlated with the first) governs the conditional PCM on
#   the positive categories.
#
# NOTE on brms compatibility:
#   The custom family relies on brms's native ordinal infrastructure
#   (`specials = "ordinal"` + `thres(gr = item)`). On CRAN brms 2.23.x
#   this combination emits invalid Stan code for custom ordinal
#   families with grouped thresholds. A patched branch is available:
#       devtools::install_github(
#         "rpsychologist/brms@fix-custom-ordinal-grouped-thres"
#       )
#
# Sign convention for the hurdle person trait:
#   In this brms parameterisation, `hu ~ ... + (1 | id)` adds the
#   person random effect to the logit of `hu = P(Y = 0)`, so higher
#   `r_id[, "hu"]` values mean *more* zeros (lower susceptibility).
#   This is the opposite sign of a Stan-style parameterisation where
#   `mu_hu = b_gate[item] - theta_gate[person]` would have higher
#   `theta_gate` mean *fewer* zeros (higher susceptibility). The two
#   parameterisations imply identical likelihoods; only the sign of
#   the recovered correlation `cor(r_id[, "hu"], r_id[, "mu"])` is
#   flipped relative to a "susceptibility x severity" labelling. This
#   is purely a labelling choice and does not affect any inferences
#   about person ranking, gate probabilities, or severity probabilities.
# ============================================================================


# ============================================================================
# 1. Custom brms family definition
# ============================================================================

#' Hurdle Partial Credit Model Custom brms Family
#'
#' A custom \pkg{brms} family implementing a hurdle Rasch partial credit
#' model (hPCM). A Bernoulli logit gate models \eqn{P(Y = 0)}; conditional
#' on \eqn{Y > 0}, an acat-logit partial credit process governs
#' transitions among the positive categories \eqn{1, \ldots, K - 1}.
#'
#' The two submodels can be interpreted as a \emph{hurdle} (presence /
#' absence of a symptom or behaviour, sometimes called "susceptibility")
#' and a \emph{partial credit} severity model (frequency / intensity
#' given presence), in the spirit of Magnus and Garnier-Villarreal
#' (2022). Person random effects on the two submodels can be modelled
#' as correlated (via brms's `(1 |g| id)` syntax), allowing the data to
#' inform whether susceptibility and severity are distinct latent
#' constructs or essentially the same trait.
#'
#' @section Usage:
#' \preformatted{
#' library(brms)
#' fit <- brm(
#'   bf(
#'     response | thres(gr = item) ~ 1 + (1 |g| id),
#'     hu ~ 0 + factor(item) + (1 |g| id)
#'   ),
#'   data     = dat,
#'   family   = hurdle_acat(),
#'   stanvars = hurdle_acat_stanvars(),
#'   ...
#' )
#' }
#'
#' @section brms compatibility note:
#' The custom family relies on brms's native ordinal infrastructure
#' (`specials = "ordinal"` + `thres(gr = item)`). On CRAN brms 2.23.x
#' this combination emits invalid Stan code for custom ordinal families
#' with grouped thresholds. A patched branch is available:
#' \preformatted{
#' devtools::install_github(
#'   "rpsychologist/brms@fix-custom-ordinal-grouped-thres"
#' )
#' }
#'
#' @section Sign convention:
#' The brms formula `hu ~ ... + (1 | id)` adds the person random effect
#' to the logit of `hu = P(Y = 0)`, so higher values of the person
#' random effect on `hu` mean \emph{more} zeros (lower susceptibility).
#' This is the opposite sign of a Stan-style parameterisation in which
#' higher `theta_gate` means fewer zeros (higher susceptibility). The
#' two parameterisations imply identical likelihoods; only the sign of
#' the recovered correlation between the two person random effects is
#' flipped relative to a "susceptibility x severity" labelling. No
#' inferences about person ranking, gate probabilities, or severity
#' probabilities are affected.
#'
#' @return A \code{customfamily} object suitable for the \code{family}
#'   argument of \code{\link[brms]{brm}}. The companion stanvars (the
#'   Stan code for the custom lpmf) must be passed via the
#'   \code{stanvars} argument; see \code{\link{hurdle_acat_stanvars}}.
#'
#' @author Kristoffer Magnusson
#'
#' @references
#' Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
#' zero-inflated graded response model for ordinal symptom data.
#' \emph{Psychological Methods}, \emph{27}(2), 261-279.
#' \doi{10.1037/met0000395}
#'
#' @seealso
#' \code{\link{hurdle_acat_stanvars}} for the Stan function block,
#' \code{\link{infit_statistic_hpcm}} for submodel-specific item infit,
#' \code{\link{q3_statistic_hpcm}} for submodel-specific Q3 residual
#' correlations.
#'
#' @importFrom brms custom_family
#' @export
hurdle_acat <- function() {
  brms::custom_family(
    "hurdle_acat",
    dpars     = c("mu", "hu"),
    links     = c("identity", "logit"),
    type      = "int",
    specials  = c("ordinal", "extra_cat"),
    threshold = "flexible"
  )
}


#' Stan Function Block for `hurdle_acat`
#'
#' Returns a \code{\link[brms]{stanvar}} containing the Stan
#' implementation of the hurdle-acat log-PMF used by
#' \code{\link{hurdle_acat}}. Pass it as the \code{stanvars} argument
#' to \code{\link[brms]{brm}}.
#'
#' @return A \code{stanvars} object adding the
#'   \code{hurdle_acat_merged_lpmf} (and helper \code{acat_logit_lpmf})
#'   to the \code{functions} block of the generated Stan model.
#'
#' @author Kristoffer Magnusson
#'
#' @seealso \code{\link{hurdle_acat}}
#'
#' @importFrom brms stanvar
#' @keywords internal
#' @export
hurdle_acat_stanvars <- function() {
  fun_block <- "
  /* acat-logit log-PMF for a single response.
   * Args:
   *   y: response category
   *   mu: latent mean parameter
   *   disc: discrimination parameter
   *   thres: ordinal thresholds (already sliced to this observation's item)
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real acat_logit_lpmf(int y, real mu, real disc, vector thres) {
    int nthres = num_elements(thres);
    vector[nthres + 1] p =
      append_row(0, cumulative_sum(disc * (mu - thres)));
    return p[y] - log_sum_exp(p);
  }

  /* Hurdle acat-logit log-PMF for a single response with merged thresholds.
   * Args:
   *   y: response category (0 = gate not crossed; 1..K-1 = severity)
   *   mu: severity latent mean
   *   hu: hurdle (zero) probability
   *   thres: merged ordinal thresholds (concatenation across items)
   *   j: start/end index for this observation's item within 'thres'
   * Returns:
   *   a scalar to be added to the log posterior
   */
  real hurdle_acat_merged_lpmf(int y, real mu, real hu,
                                vector thres, array[] int j) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    }
    return bernoulli_lpmf(0 | hu)
           + acat_logit_lpmf(y | mu, 1.0, thres[j[1]:j[2]]);
  }
  "
  brms::stanvar(scode = fun_block, block = "functions")
}


# ============================================================================
# 2. Post-processing helpers for the custom family
# ============================================================================
#
# These functions are exported because brms looks them up by name (it
# matches `log_lik_<family>`, `posterior_predict_<family>`,
# `posterior_epred_<family>`) when the user calls the corresponding
# top-level brms methods on a model fitted with `hurdle_acat()`. They
# are not meant to be called directly by end-users.

# ---------------------------------------------------------------------------
# Internal acat-logit kernel
# ---------------------------------------------------------------------------
#
# The post-processing helpers below previously relied on unexported brms
# internals (`brms:::dacat`, `brms:::pordinal`, `brms:::first_greater`,
# `brms:::subset_thres`). Using `:::` triggers an R CMD check NOTE and is
# fragile: those functions can change or disappear between brms releases.
#
# For the *logit* link -- the only link these helpers ever use -- the
# adjacent-category (acat) model has a simple closed form, so we
# reimplement the needed pieces locally. The implementations below were
# verified to reproduce the corresponding brms internals to < 1e-10 over
# hundreds of random inputs (category probabilities, CDF, and the
# inverse-CDF index). The only remaining brms coupling is via the public
# `brms::get_dpar()` and the `prep$thres` list elements, both already
# used elsewhere in these helpers.

#' Adjacent-category (logit) probabilities for the hurdle PCM severity part
#'
#' Local replacement for the relevant behaviour of \code{brms:::dacat()}
#' with \code{link = "logit"}. For category scores \eqn{k = 0, \ldots, K-1}
#' (with \eqn{K - 1} thresholds), the adjacent-category logit model gives
#' \deqn{P(Y = k) \propto \exp\left(\sum_{j \le k} disc \cdot
#'   (\eta - \tau_j)\right), \quad P(Y = 0) \propto 1.}
#'
#' @param eta Numeric vector (length S) of latent values per draw.
#' @param thres Numeric matrix (S x nthres) of thresholds per draw.
#' @param disc Scalar discrimination (default 1).
#' @return An S x (nthres + 1) matrix of category probabilities, columns
#'   ordered as categories 1..K (1-based, matching \code{brms:::dacat()}).
#' @keywords internal
#' @noRd
.hpcm_acat_probs <- function(eta, thres, disc = 1) {
  thres <- as.matrix(thres)
  # disc * (eta - thres): eta (length S) recycles down the rows so that
  # element [s, j] = disc * (eta[s] - thres[s, j]).
  x  <- disc * (eta - thres)
  ex <- exp(x)
  # Cumulative product across thresholds (columns), then prepend the
  # category-0 reference column of 1s.
  cum <- ex
  if (ncol(ex) > 1L) {
    for (k in 2:ncol(ex)) cum[, k] <- cum[, k - 1L] * ex[, k]
  }
  out <- cbind(1, cum)
  out / rowSums(out)
}

#' Select category probabilities (local replacement for brms:::dacat)
#'
#' @param x Integer vector of (1-based) category indices to return.
#' @inheritParams .hpcm_acat_probs
#' @return An S x length(x) matrix, matching \code{brms:::dacat(..., }
#'   \code{link = "logit")[, x, drop = FALSE]}.
#' @keywords internal
#' @noRd
.hpcm_dacat <- function(x, eta, thres, disc = 1) {
  probs <- .hpcm_acat_probs(eta, thres, disc = disc)
  probs[, x, drop = FALSE]
}

#' Cumulative acat probabilities (local replacement for brms:::pordinal)
#'
#' @param q Integer vector of (1-based) category indices whose CDF values
#'   to return.
#' @inheritParams .hpcm_acat_probs
#' @return An S x length(q) matrix of cumulative category probabilities
#'   \eqn{P(Y \le q)}, matching \code{brms:::pordinal(..., }
#'   \code{family = "acat", link = "logit")}.
#' @keywords internal
#' @noRd
.hpcm_pordinal <- function(q, eta, thres, disc = 1) {
  probs <- .hpcm_acat_probs(eta, thres, disc = disc)
  cdf <- probs
  if (ncol(probs) > 1L) {
    for (k in 2:ncol(probs)) cdf[, k] <- cdf[, k - 1L] + probs[, k]
  }
  cdf[, q, drop = FALSE]
}

#' First column index exceeding a target (local replacement for
#' brms:::first_greater)
#'
#' For each row of \code{A}, returns the index of the first column whose
#' value is \eqn{\ge} \code{target}; if none, returns \code{ncol(A)}.
#' Used for inverse-CDF sampling. Matches \code{brms:::first_greater()}
#' in value (returned as integer rather than double).
#'
#' @param A Numeric matrix (S x K), typically a per-draw CDF.
#' @param target Numeric vector (length S) of thresholds, e.g. uniform
#'   draws.
#' @return Integer vector (length S) of selected column indices.
#' @keywords internal
#' @noRd
.hpcm_first_greater <- function(A, target) {
  A <- as.matrix(A)
  n <- nrow(A)
  K <- ncol(A)
  out   <- rep.int(K, n)
  found <- logical(n)
  for (i in seq_len(K)) {
    hit <- !found & (target <= A[, i])
    out[hit] <- i
    found <- found | hit
  }
  out
}

#' Per-observation thresholds (local replacement for brms:::subset_thres)
#'
#' Identical to \code{brms:::subset_thres()} but accesses the
#' \code{prep$thres} list elements directly (no \code{:::}).
#'
#' @param prep A brms prep object.
#' @param i Observation index.
#' @return An S x nthres_i matrix of thresholds for observation \code{i}.
#' @keywords internal
#' @noRd
.hpcm_subset_thres <- function(prep, i) {
  thres  <- prep$thres$thres
  Jthres <- prep$thres$Jthres
  if (!is.null(Jthres)) {
    thres <- thres[, Jthres[i, 1]:Jthres[i, 2], drop = FALSE]
  }
  thres
}

#' Log-Likelihood Method for `hurdle_acat`
#'
#' Not called directly. brms invokes this when \code{\link[brms]{log_lik}}
#' is called on a model fitted with \code{\link{hurdle_acat}}.
#'
#' @param i Observation index (supplied by brms).
#' @param prep A brms prep object (supplied by brms).
#'
#' @return A numeric vector of length \code{prep$ndraws}.
#'
#' @author Kristoffer Magnusson
#'
#' @keywords internal
#' @export
log_lik_hurdle_acat <- function(i, prep) {
  hu <- brms::get_dpar(prep, "hu", i = i)
  y <- prep$data$Y[i]
  if (y == 0L) {
    return(log(hu))
  }
  mu <- brms::get_dpar(prep, "mu", i = i)
  thres <- .hpcm_subset_thres(prep, i)
  sev_p <- .hpcm_dacat(
    x = y, eta = mu, thres = thres, disc = 1
  )
  log1p(-hu) + log(as.vector(sev_p))
}


#' Posterior-Predict Method for `hurdle_acat`
#'
#' Not called directly. brms invokes this when
#' \code{\link[brms]{posterior_predict}} is called on a model fitted
#' with \code{\link{hurdle_acat}}.
#'
#' @param i Observation index (supplied by brms).
#' @param prep A brms prep object (supplied by brms).
#' @param ... Currently unused.
#'
#' @return An integer vector of length \code{prep$ndraws} containing
#'   replicated category labels in \eqn{0, 1, \ldots, K - 1}.
#'
#' @author Kristoffer Magnusson
#'
#' @keywords internal
#' @export
posterior_predict_hurdle_acat <- function(i, prep, ...) {
  hu <- brms::get_dpar(prep, "hu", i = i)
  mu <- brms::get_dpar(prep, "mu", i = i)
  thres <- .hpcm_subset_thres(prep, i)
  p <- .hpcm_pordinal(
    seq_len(ncol(thres) + 1L),
    eta = mu, disc = 1, thres = thres
  )
  draws <- .hpcm_first_greater(p, target = stats::runif(prep$ndraws))
  draws[stats::runif(prep$ndraws) < hu] <- 0L
  draws
}


#' Posterior-Epred Method for `hurdle_acat`
#'
#' Not called directly. brms invokes this when
#' \code{\link[brms]{posterior_epred}} is called on a model fitted with
#' \code{\link{hurdle_acat}}. Returns an \eqn{S \times N \times K}
#' array of category probabilities, where slot 1 is \eqn{hu = P(Y = 0)}
#' and slots \eqn{2..K} are \eqn{(1 - hu) \cdot P_{sev}(k)} for
#' \eqn{k = 1, \ldots, K - 1}.
#'
#' @param prep A brms prep object (supplied by brms).
#'
#' @return A numeric \eqn{S \times N \times K_{total}} array of
#'   category probabilities.
#'
#' @author Kristoffer Magnusson
#'
#' @keywords internal
#' @export
posterior_epred_hurdle_acat <- function(prep) {
  hu_mat <- brms::get_dpar(prep, "hu")
  mu_mat <- brms::get_dpar(prep, "mu")
  nobs <- ncol(mu_mat)
  K_total <- max(prep$thres$nthres) + 2L
  out <- array(0, dim = c(nrow(mu_mat), nobs, K_total))
  for (n in seq_len(nobs)) {
    # Items with fewer thresholds than the maximum have fewer categories.
    # Their probabilities fill the first columns and the rest stay 0, as
    # in brms's posterior_epred for ordinal families.
    sev <- .hpcm_acat_probs(
      eta = mu_mat[, n],
      thres = .hpcm_subset_thres(prep, n),
      disc = 1
    )
    out[, n, 1L] <- hu_mat[, n]
    out[, n, 1L + seq_len(ncol(sev))] <- (1 - hu_mat[, n]) * sev
  }
  out
}


# Response scores 0..K-1 as brms codes them for hurdle_acat: numeric
# responses as they are, ordered factors as their integer code minus 1.
#' @keywords internal
#' @noRd
.hpcm_response_scores <- function(x) {
  if (is.factor(x)) as.integer(x) - 1L else as.numeric(x)
}

# ============================================================================
# 3. Posterior-predictive infit (per-submodel)
# ============================================================================

#' Posterior Predictive Infit Statistic for the Hurdle Partial Credit Model
#'
#' Computes posterior predictive item infit statistics separately for
#' the two submodels of a hurdle partial credit model fitted with \pkg{brms}
#' using the \code{\link{hurdle_acat}} custom family: (i) the
#' \emph{hurdle} submodel for \eqn{P(Y > 0)} (Bernoulli) and (ii) the
#' \emph{partial credit} severity submodel for
#' \eqn{P(Y = k \mid Y > 0)} on the positive categories. For each
#' posterior draw, expected values and variances are derived from the
#' submodel-specific category probabilities, and variance-weighted
#' standardised residuals are computed for both observed and
#' replicated data. Expectations are given each draw's person
#' parameters, as in \code{infit_statistic(method = "theta")}, since the
#' hurdle model has no sufficient total score. Use that method for a
#' Rasch-family model when results should be comparable with this
#' function.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object using the
#'   \code{\link{hurdle_acat}} custom family (i.e.,
#'   \code{posterior_epred} returns an \code{S x N x K_total} array
#'   whose first category is the hurdle / zero probability).
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data (e.g., \code{item}).
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data (e.g., \code{id}).
#' @param ndraws_use Optional positive integer. If specified, a random
#'   subset of posterior draws of this size is used. If \code{NULL}
#'   (the default), all draws are used.
#' @param outfit Logical. If \code{TRUE}, outfit statistics are
#'   computed alongside infit. Default is \code{FALSE}.
#'
#' @return A list with two elements, each a \code{\link[tibble]{tibble}}
#'   in the same format as the output of \code{\link{infit_statistic}}
#'   (and directly compatible with \code{\link{infit_post}}):
#' \describe{
#'   \item{\code{hurdle}}{Item infit for the Bernoulli hurdle submodel
#'     on \eqn{1[Y > 0]}, evaluated on \emph{all} observations.}
#'   \item{\code{pcm}}{Item infit for the partial credit severity
#'     submodel on \eqn{P(Y = k \mid Y > 0)}, evaluated only on the
#'     observations with \eqn{Y_{obs} > 0}.}
#' }
#'
#' @details
#' The hurdle PCM splits the generative process into:
#' \enumerate{
#'   \item A Bernoulli hurdle with \eqn{hu = P(Y = 0)}.
#'   \item A partial credit / acat-logit severity process over the
#'     positive categories \eqn{1, \ldots, K - 1}, applied only when
#'     the hurdle is crossed.
#' }
#'
#' \code{posterior_epred} for the \code{hurdle_acat} family returns an
#' \code{S x N x K_total} array whose first category is \eqn{hu} and
#' whose remaining categories are \eqn{(1 - hu) \cdot P_{sev}(k)}. The
#' two submodel infits are computed as follows:
#'
#' \strong{Hurdle submodel.} All observations contribute. The Bernoulli
#' moments are \eqn{E_{hurdle} = 1 - hu} and
#' \eqn{Var_{hurdle} = hu \cdot (1 - hu)}. Observed and replicated
#' scores are \eqn{1[Y_{obs} > 0]} and \eqn{1[Y_{rep} > 0]} with
#' \eqn{Y_{rep}} obtained from \code{\link[brms]{posterior_predict}}.
#'
#' \strong{Partial credit submodel.} Only observations with
#' \eqn{Y_{obs} > 0} contribute. Conditional probabilities are
#' \deqn{P(Y = k \mid Y > 0) = epred[, , k+1] / (1 - hu),
#'   \quad k = 1, \ldots, K - 1.}
#' The expected value and variance under these probabilities use
#' category scores
#' \eqn{1, \ldots, K - 1}. Replicated severity values are drawn
#' independently for each (draw, observation) from these conditional
#' probabilities via inverse-CDF sampling, so the partial credit PPC
#' is not contaminated by hurdle misfit.
#'
#' Within each submodel the infit per item is
#' \deqn{Infit_i^{(s)} = \sum_v (X_{vi} - E_{vi}^{(s)})^2 /
#'   \sum_v Var_{vi}^{(s)},}
#' with the sum restricted to the rows the submodel applies to
#' (all rows for the hurdle; rows with \eqn{Y_{obs} > 0} for partial
#' credit).
#'
#' @references
#' Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013).
#' \emph{Rasch Models in Health}. Iste and Wiley, pp. 86--90.
#'
#' Kreiner, S. & Christensen, K. B. (2011). Exact evaluation of bias in
#' Rasch model residuals. \emph{Advances in Mathematics Research},
#' \emph{12}, 19--40.
#'
#' Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
#' zero-inflated graded response model for ordinal symptom data.
#' \emph{Psychological Methods}, \emph{27}(2), 261-279.
#' \doi{10.1037/met0000395}
#'
#' @seealso
#' \code{\link{infit_statistic}} for the single-submodel version,
#' \code{\link{infit_post}} for summarising and plotting the draws,
#' \code{\link{q3_statistic_hpcm}} for hPCM Q3 residual correlations.
#'
#' @importFrom brms posterior_epred posterior_predict ndraws
#' @importFrom dplyr arrange group_by
#' @importFrom rlang enquo !! .data as_name
#' @importFrom stats formula runif
#' @importFrom tibble tibble
#' @export
infit_statistic_hpcm <- function(model,
                                 item_var = item,
                                 person_var = id,
                                 ndraws_use = NULL,
                                 outfit = FALSE) {
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }

  item_quo    <- rlang::enquo(item_var)
  person_quo  <- rlang::enquo(person_var)
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

  if (length(dim(epred_array)) != 3) {
    stop("Expected a 3D posterior_epred array (S x N x K_total) for a ",
         "hurdle_acat model.", call. = FALSE)
  }

  n_draws <- dim(epred_array)[1]
  n_obs   <- dim(epred_array)[2]
  K_total <- dim(epred_array)[3]
  K_sev   <- K_total - 1L  # number of positive (severity) categories

  obs_response <- .hpcm_response_scores(model$data[[resp_var]])

  # =================================================================
  # HURDLE submodel: Bernoulli on 1[Y > 0]
  # =================================================================
  hu_mat       <- matrix(epred_array[, , 1L], nrow = n_draws)  # S x N: P(Y = 0)
  E_hurdle     <- 1 - hu_mat                       # S x N: P(Y > 0)
  Var_hurdle   <- hu_mat * (1 - hu_mat)            # S x N
  Var_hurdle[Var_hurdle < 1e-12] <- 1e-12

  hurdle_obs_vec <- as.integer(obs_response > 0)
  hurdle_obs_mat <- matrix(hurdle_obs_vec,
                           nrow = n_draws, ncol = n_obs, byrow = TRUE)
  hurdle_rep_mat <- (yrep_mat > 0) * 1L

  resid2_hurdle_obs <- (hurdle_obs_mat - E_hurdle)^2
  resid2_hurdle_rep <- (hurdle_rep_mat - E_hurdle)^2

  # =================================================================
  # PCM submodel: PCM on Y | Y > 0
  # =================================================================
  one_minus_hu <- 1 - hu_mat
  one_minus_hu[one_minus_hu < 1e-12] <- 1e-12

  sev_p_array <- epred_array[, , 2:K_total, drop = FALSE]   # S x N x K_sev
  divisor_arr <- array(rep(one_minus_hu, K_sev),
                       dim = c(n_draws, n_obs, K_sev))
  cond_p_array <- sev_p_array / divisor_arr
  cond_p_array[!is.finite(cond_p_array)] <- 0

  cat_values_sev <- seq_len(K_sev)  # severity scores 1, 2, ..., K-1

  cond_p_2d  <- matrix(cond_p_array,
                       nrow = n_draws * n_obs, ncol = K_sev)
  E_pcm_vec  <- cond_p_2d %*% cat_values_sev
  E_pcm_mat  <- matrix(E_pcm_vec,  nrow = n_draws, ncol = n_obs)
  E2_pcm_vec <- cond_p_2d %*% (cat_values_sev^2)
  Var_pcm_mat <- matrix(E2_pcm_vec, nrow = n_draws, ncol = n_obs) -
    E_pcm_mat^2
  Var_pcm_mat[Var_pcm_mat < 1e-12] <- 1e-12

  # Replicated severity values: inverse-CDF draws from conditional
  # probabilities for the rows that contribute (Y_obs > 0).
  yrep_pcm_mat <- matrix(NA_integer_, nrow = n_draws, ncol = n_obs)
  pcm_rows <- which(hurdle_obs_vec == 1L)
  if (length(pcm_rows) > 0L) {
    for (n in pcm_rows) {
      probs_sn <- matrix(cond_p_array[, n, ], nrow = n_draws)  # S x K_sev
      # Column-wise cumulative sum, valid for any number of draws and
      # categories
      cum_p    <- probs_sn
      if (ncol(cum_p) > 1L) {
        for (j in 2:ncol(cum_p)) cum_p[, j] <- cum_p[, j - 1L] + cum_p[, j]
      }
      u        <- stats::runif(n_draws)
      yrep_pcm_mat[, n] <- max.col(cum_p >= u, ties.method = "first")
    }
  }

  obs_pcm_mat <- matrix(obs_response,
                        nrow = n_draws, ncol = n_obs, byrow = TRUE)
  mask_pcm <- hurdle_obs_mat == 1L

  resid2_pcm_obs <- (obs_pcm_mat - E_pcm_mat)^2
  resid2_pcm_obs[!mask_pcm] <- 0
  resid2_pcm_rep <- (yrep_pcm_mat - E_pcm_mat)^2
  resid2_pcm_rep[!mask_pcm] <- 0
  Var_pcm_masked <- Var_pcm_mat
  Var_pcm_masked[!mask_pcm] <- 0

  # =================================================================
  # Aggregate by item per draw
  # =================================================================
  items <- model$data[[item_name]]
  unique_items <- unique(items)
  k_items <- length(unique_items)

  item_col_idx <- vector("list", k_items)
  for (idx in seq_len(k_items)) {
    item_col_idx[[idx]] <- which(items == unique_items[idx])
  }

  infit_hurdle_obs <- matrix(NA_real_, nrow = n_draws, ncol = k_items)
  infit_hurdle_rep <- matrix(NA_real_, nrow = n_draws, ncol = k_items)
  infit_pcm_obs    <- matrix(NA_real_, nrow = n_draws, ncol = k_items)
  infit_pcm_rep    <- matrix(NA_real_, nrow = n_draws, ncol = k_items)
  if (outfit) {
    outfit_hurdle_obs <- matrix(NA_real_, nrow = n_draws, ncol = k_items)
    outfit_hurdle_rep <- matrix(NA_real_, nrow = n_draws, ncol = k_items)
    outfit_pcm_obs    <- matrix(NA_real_, nrow = n_draws, ncol = k_items)
    outfit_pcm_rep    <- matrix(NA_real_, nrow = n_draws, ncol = k_items)
  }

  for (idx in seq_len(k_items)) {
    cols <- item_col_idx[[idx]]

    # --- Hurdle: all observations contribute ---
    sum_var_h <- rowSums(Var_hurdle[, cols, drop = FALSE])
    sum_var_h[sum_var_h < 1e-12] <- 1e-12
    infit_hurdle_obs[, idx] <- rowSums(resid2_hurdle_obs[, cols, drop = FALSE]) /
      sum_var_h
    infit_hurdle_rep[, idx] <- rowSums(resid2_hurdle_rep[, cols, drop = FALSE]) /
      sum_var_h

    # --- PCM: only Y_obs > 0 rows contribute ---
    sum_var_p <- rowSums(Var_pcm_masked[, cols, drop = FALSE])
    has_pcm <- sum_var_p >= 1e-12
    sum_var_p_safe <- sum_var_p
    sum_var_p_safe[!has_pcm] <- NA_real_
    infit_pcm_obs[, idx] <- rowSums(resid2_pcm_obs[, cols, drop = FALSE]) /
      sum_var_p_safe
    infit_pcm_rep[, idx] <- rowSums(resid2_pcm_rep[, cols, drop = FALSE]) /
      sum_var_p_safe

    if (outfit) {
      Z2_h_obs <- resid2_hurdle_obs[, cols, drop = FALSE] /
        Var_hurdle[, cols, drop = FALSE]
      Z2_h_rep <- resid2_hurdle_rep[, cols, drop = FALSE] /
        Var_hurdle[, cols, drop = FALSE]
      outfit_hurdle_obs[, idx] <- rowMeans(Z2_h_obs, na.rm = TRUE)
      outfit_hurdle_rep[, idx] <- rowMeans(Z2_h_rep, na.rm = TRUE)

      Z2_p_obs <- resid2_pcm_obs[, cols, drop = FALSE] /
        Var_pcm_mat[, cols, drop = FALSE]
      Z2_p_rep <- resid2_pcm_rep[, cols, drop = FALSE] /
        Var_pcm_mat[, cols, drop = FALSE]
      pcm_mask_cols <- mask_pcm[, cols, drop = FALSE]
      Z2_p_obs[!pcm_mask_cols] <- NA_real_
      Z2_p_rep[!pcm_mask_cols] <- NA_real_
      outfit_pcm_obs[, idx] <- rowMeans(Z2_p_obs, na.rm = TRUE)
      outfit_pcm_rep[, idx] <- rowMeans(Z2_p_rep, na.rm = TRUE)
    }
  }

  # =================================================================
  # Assemble output: list(hurdle = tbl, pcm = tbl)
  # =================================================================
  draw_seq <- seq_len(n_draws)
  item_rep <- rep(unique_items, each = n_draws)
  draw_rep <- rep(draw_seq, times = k_items)

  build_tbl <- function(infit_obs, infit_rep,
                        out_obs = NULL, out_rep = NULL) {
    tbl <- tibble::tibble(
      item      = item_rep,
      draw      = draw_rep,
      infit     = as.vector(infit_obs),
      infit_rep = as.vector(infit_rep)
    )
    if (!is.null(out_obs)) {
      tbl$outfit     <- as.vector(out_obs)
      tbl$outfit_rep <- as.vector(out_rep)
    }
    names(tbl)[names(tbl) == "item"] <- item_name
    tbl <- dplyr::arrange(tbl, .data[[item_name]], .data$draw)
    dplyr::group_by(tbl, .data[[item_name]])
  }

  if (outfit) {
    list(
      hurdle = build_tbl(infit_hurdle_obs, infit_hurdle_rep,
                         outfit_hurdle_obs, outfit_hurdle_rep),
      pcm    = build_tbl(infit_pcm_obs,    infit_pcm_rep,
                         outfit_pcm_obs,    outfit_pcm_rep)
    )
  } else {
    list(
      hurdle = build_tbl(infit_hurdle_obs, infit_hurdle_rep),
      pcm    = build_tbl(infit_pcm_obs,    infit_pcm_rep)
    )
  }
}


# ============================================================================
# 4. Posterior-predictive Q3 (per-submodel)
# ============================================================================

#' Posterior Predictive Q3 Residual Correlations for the Hurdle PCM
#'
#' Computes Yen's Q3 residual correlations separately for the two
#' submodels of a hurdle partial credit model fitted with \pkg{brms}
#' using the \code{\link{hurdle_acat}} custom family: (i) the
#' \emph{hurdle} submodel for \eqn{P(Y > 0)} (Bernoulli) and (ii) the
#' \emph{partial credit} severity submodel for
#' \eqn{P(Y = k \mid Y > 0)} on the positive categories. For each
#' posterior draw, residuals are computed from each submodel's
#' expected values and Pearson correlations between item-pair
#' residuals are returned for both observed and replicated data.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object using the
#'   \code{\link{hurdle_acat}} custom family.
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Default is \code{id}.
#' @param ndraws_use Optional positive integer. If specified, a random
#'   subset of posterior draws of this size is used. If \code{NULL}
#'   (the default), all draws are used.
#'
#' @return A list with two elements, each a \code{\link[tibble]{tibble}}
#'   in the same long format as \code{\link{q3_statistic}} (and
#'   directly compatible with \code{\link{q3_post}}):
#' \describe{
#'   \item{\code{hurdle}}{Q3 correlations of hurdle residuals
#'     \eqn{1[Y > 0] - (1 - hu)}, using \emph{all} observations.}
#'   \item{\code{pcm}}{Q3 correlations of partial credit residuals
#'     \eqn{Y - E_{pcm}}, using only observations with
#'     \eqn{Y_{obs} > 0}. Rows with \eqn{Y_{obs} = 0} are set to
#'     \code{NA} and excluded pairwise.}
#' }
#'
#' @details
#' The hurdle PCM has a Bernoulli hurdle \eqn{hu = P(Y = 0)} and a
#' partial credit severity process on the positive categories
#' \eqn{1, \ldots, K - 1}. \code{posterior_epred} for
#' \code{\link{hurdle_acat}} returns an \code{S x N x K_total} array
#' whose first category is \eqn{hu} and whose remaining categories are
#' \eqn{(1 - hu) \cdot P_{sev}(k)}.
#'
#' \strong{Hurdle residuals.} For each (draw, observation):
#' \deqn{d_{hurdle} = 1[Y > 0] - (1 - hu).}
#' All observations contribute. Replicated residuals use
#' \eqn{1[Y_{rep} > 0]} from \code{\link[brms]{posterior_predict}}.
#'
#' \strong{Partial credit residuals.} For each (draw, observation)
#' with \eqn{Y_{obs} > 0}:
#' \deqn{d_{pcm} = Y - \sum_{k=1}^{K-1} k \cdot P(Y = k \mid Y > 0),}
#' with conditional probabilities computed as
#' \eqn{epred[, , k+1] / (1 - hu)}. Replicated partial credit values
#' are drawn independently per (draw, observation) from these
#' conditional probabilities via inverse-CDF sampling, so the partial
#' credit PPC is not contaminated by hurdle misfit. Rows with
#' \eqn{Y_{obs} = 0} are set to \code{NA} in the wide residual matrix
#' and excluded via \code{use = "pairwise.complete.obs"} when
#' correlations are computed.
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
#'
#' Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
#' zero-inflated graded response model for ordinal symptom data.
#' \emph{Psychological Methods}, \emph{27}(2), 261-279.
#' \doi{10.1037/met0000395}
#'
#' @seealso
#' \code{\link{q3_statistic}} for the single-submodel version,
#' \code{\link{q3_post}} for summarising and plotting the draws,
#' \code{\link{infit_statistic_hpcm}} for hPCM infit.
#'
#' @importFrom brms posterior_epred posterior_predict ndraws
#' @importFrom rlang enquo as_name .data
#' @importFrom stats formula cor runif
#' @importFrom tibble tibble as_tibble
#' @export
q3_statistic_hpcm <- function(model,
                              item_var = item,
                              person_var = id,
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

  .check_unique_person_item(model$data[[person_name]], model$data[[item_name]],
                            person_name, item_name, "q3_statistic_hpcm")

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

  if (length(dim(epred_array)) != 3) {
    stop("Expected a 3D posterior_epred array (S x N x K_total) for a ",
         "hurdle_acat model.", call. = FALSE)
  }

  n_draws <- dim(epred_array)[1]
  n_obs   <- dim(epred_array)[2]
  K_total <- dim(epred_array)[3]
  K_sev   <- K_total - 1L

  obs_response <- .hpcm_response_scores(model$data[[resp_var]])

  # =================================================================
  # HURDLE submodel residuals (Bernoulli on 1[Y > 0])
  # =================================================================
  hu_mat   <- matrix(epred_array[, , 1L], nrow = n_draws)  # S x N: P(Y = 0)
  E_hurdle <- 1 - hu_mat                     # S x N: P(Y > 0)

  hurdle_obs_vec <- as.integer(obs_response > 0)
  hurdle_obs_mat <- matrix(hurdle_obs_vec,
                           nrow = n_draws, ncol = n_obs, byrow = TRUE)
  hurdle_rep_mat <- (yrep_mat > 0) * 1L

  resid_hurdle_obs <- hurdle_obs_mat - E_hurdle
  resid_hurdle_rep <- hurdle_rep_mat - E_hurdle

  # =================================================================
  # PCM submodel residuals (Y | Y > 0)
  # =================================================================
  one_minus_hu <- 1 - hu_mat
  one_minus_hu[one_minus_hu < 1e-12] <- 1e-12

  sev_p_array <- epred_array[, , 2:K_total, drop = FALSE]
  divisor_arr <- array(rep(one_minus_hu, K_sev),
                       dim = c(n_draws, n_obs, K_sev))
  cond_p_array <- sev_p_array / divisor_arr
  cond_p_array[!is.finite(cond_p_array)] <- 0

  cat_values_sev <- seq_len(K_sev)
  cond_p_2d <- matrix(cond_p_array,
                      nrow = n_draws * n_obs, ncol = K_sev)
  E_pcm_vec <- cond_p_2d %*% cat_values_sev
  E_pcm_mat <- matrix(E_pcm_vec, nrow = n_draws, ncol = n_obs)

  yrep_pcm_mat <- matrix(NA_integer_, nrow = n_draws, ncol = n_obs)
  pcm_rows <- which(hurdle_obs_vec == 1L)
  if (length(pcm_rows) > 0L) {
    for (n in pcm_rows) {
      probs_sn <- matrix(cond_p_array[, n, ], nrow = n_draws)
      # Column-wise cumulative sum, valid for any number of draws and
      # categories
      cum_p    <- probs_sn
      if (ncol(cum_p) > 1L) {
        for (j in 2:ncol(cum_p)) cum_p[, j] <- cum_p[, j - 1L] + cum_p[, j]
      }
      u        <- stats::runif(n_draws)
      yrep_pcm_mat[, n] <- max.col(cum_p >= u, ties.method = "first")
    }
  }

  obs_pcm_mat <- matrix(obs_response,
                        nrow = n_draws, ncol = n_obs, byrow = TRUE)
  mask_pcm <- hurdle_obs_mat == 1L

  resid_pcm_obs <- obs_pcm_mat - E_pcm_mat
  resid_pcm_rep <- yrep_pcm_mat - E_pcm_mat
  resid_pcm_obs[!mask_pcm] <- NA_real_
  resid_pcm_rep[!mask_pcm] <- NA_real_

  # =================================================================
  # Person x item mapping (shared by both submodels)
  # =================================================================
  items   <- model$data[[item_name]]
  persons <- model$data[[person_name]]
  unique_items   <- sort(unique(items))
  unique_persons <- sort(unique(persons))
  k <- length(unique_items)
  n_persons <- length(unique_persons)

  person_idx <- match(persons, unique_persons)
  item_idx   <- match(items, unique_items)
  lin_idx <- person_idx + (item_idx - 1L) * n_persons

  pair_grid <- expand.grid(j = seq_len(k), i = seq_len(k))
  pair_grid <- pair_grid[pair_grid$i < pair_grid$j, ]
  pair_grid <- pair_grid[order(pair_grid$i, pair_grid$j), ]
  n_pairs <- nrow(pair_grid)

  pair_labels <- paste0(unique_items[pair_grid$i], " : ",
                        unique_items[pair_grid$j])
  pair_i_idx <- pair_grid$i
  pair_j_idx <- pair_grid$j

  # =================================================================
  # Compute Q3 per draw for each submodel
  # =================================================================
  compute_q3 <- function(resid_obs_mat, resid_rep_mat) {
    q3_obs <- matrix(NA_real_, nrow = n_draws, ncol = n_pairs)
    q3_rep <- matrix(NA_real_, nrow = n_draws, ncol = n_pairs)

    for (s in seq_len(n_draws)) {
      wide_obs <- matrix(NA_real_, nrow = n_persons, ncol = k)
      wide_rep <- matrix(NA_real_, nrow = n_persons, ncol = k)
      wide_obs[lin_idx] <- resid_obs_mat[s, ]
      wide_rep[lin_idx] <- resid_rep_mat[s, ]

      cor_obs <- suppressWarnings(
        stats::cor(wide_obs, use = "pairwise.complete.obs")
      )
      cor_rep <- suppressWarnings(
        stats::cor(wide_rep, use = "pairwise.complete.obs")
      )

      for (p in seq_len(n_pairs)) {
        q3_obs[s, p] <- cor_obs[pair_i_idx[p], pair_j_idx[p]]
        q3_rep[s, p] <- cor_rep[pair_i_idx[p], pair_j_idx[p]]
      }
    }

    result_list <- vector("list", n_pairs)
    for (p in seq_len(n_pairs)) {
      result_list[[p]] <- tibble::tibble(
        draw      = seq_len(n_draws),
        item_pair = pair_labels[p],
        item_1    = unique_items[pair_i_idx[p]],
        item_2    = unique_items[pair_j_idx[p]],
        q3        = q3_obs[, p],
        q3_rep    = q3_rep[, p]
      )
    }
    tibble::as_tibble(do.call(rbind, result_list))
  }

  list(
    hurdle = compute_q3(resid_hurdle_obs, resid_hurdle_rep),
    pcm    = compute_q3(resid_pcm_obs,    resid_pcm_rep)
  )
}


# ============================================================================
# 5. Parameter mapping conventions for the brms hurdle_acat model
# ============================================================================
#
# Given the recommended formula
#
#     bf(
#       response | thres(gr = item) ~ 1 + (1 |g| id),
#       hu ~ 0 + factor(item) + (1 |g| id)
#     ),
#     family = hurdle_acat()
#
# brms produces these posterior variables:
#
#   PCM submodel (severity):
#     b_Intercept[<item>, <k>]       per-item PCM thresholds tau_{ik}
#     r_id[<v>, Intercept]           person random effect on severity
#     sd_id__Intercept               SD of person severity
#
#   Hurdle submodel (presence/absence):
#     b_hu_factoritem<item>          per-item hurdle intercepts on logit(hu);
#                                    hu = P(Y == 0), so higher value = more
#                                    zeros = harder to cross. Used directly
#                                    as the hurdle item difficulty.
#     r_id[<v>, hu_Intercept]        person random effect on logit(hu)
#     sd_id__hu_Intercept            SD of person hurdle effect
#
#   Cross-submodel:
#     cor_id__Intercept__hu_Intercept  correlation between the two person
#                                      random effects
#
# To turn the brms hurdle person random effect into a Rasch-style person
# trait (higher = more presence / susceptibility), we negate the brms
# random effect:
#
#     theta_hurdle[v] = -r_id[v, hu_Intercept]
#
# Under this sign convention, the hurdle submodel reads as
#     P(Y > 0) = plogis(theta_hurdle[v] - delta_hurdle[i]),
# with delta_hurdle[i] = b_hu_factoritem<i>. The trait correlation
# reported by the functions below is the brms correlation with its sign
# flipped, so it is read as cor(theta_hurdle, theta_pcm).
# ============================================================================


# ============================================================================
# 6. Item parameters (per-submodel)
# ============================================================================

#' Extract Item Parameters from a Hurdle Partial Credit Model
#'
#' Extracts item parameters from a fitted hurdle Rasch partial credit
#' model (\code{family = \link{hurdle_acat}()}), returning a per-submodel
#' breakdown: hurdle item difficulties (Bernoulli logit on \eqn{P(Y > 0)})
#' and partial credit thresholds. Each submodel's output mirrors the
#' shape of \code{\link{item_parameters}} so existing plotting and
#' downstream functions can be applied directly to \code{res$hurdle} or
#' \code{res$pcm}.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object using the
#'   \code{\link{hurdle_acat}} custom family. The recommended formula
#'   is \preformatted{
#' bf(
#'   response | thres(gr = item) ~ 1 + (1 |g| id),
#'   hu ~ 0 + factor(item) + (1 |g| id)
#' )}
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Default is \code{id}.
#' @param draws Logical. If \code{TRUE}, a draws matrix of full
#'   posterior draws is included in each submodel's output. Default
#'   is \code{FALSE}.
#' @param center Logical. If \code{TRUE} (the default), item parameters
#'   are shifted within each submodel so their mean is zero, matching
#'   the convention in \code{\link{item_parameters}}. Person parameters
#'   reported by \code{\link{person_parameters_hpcm}} use the same
#'   shifts.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the highest density
#'   continuous interval (HDCI) reported in the summary. Default is
#'   0.95.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{hurdle}}{A list with the same structure as the output
#'     of \code{\link{item_parameters}} applied to a dichotomous Rasch
#'     model: \code{locations}, \code{locations_wide}, \code{summary},
#'     \code{item_information}, \code{person_sd}, optionally
#'     \code{draws_matrix}. \code{location} is the brms posterior mean
#'     of \code{b_hu_factoritem<i>}; higher values mean more zeros
#'     (a harder hurdle to cross).}
#'   \item{\code{pcm}}{A list with the same structure as
#'     \code{\link{item_parameters}} applied to a PCM: \code{locations}
#'     (long), \code{locations_wide} (t1, t2, ..., location),
#'     \code{summary}, \code{item_information}, \code{threshold_order},
#'     \code{person_sd}, optionally \code{draws_matrix}. \code{location}
#'     columns are posterior means of \code{b_Intercept[<item>, <k>]}.}
#'   \item{\code{correlation}}{A \code{\link[tibble]{tibble}} with
#'     \code{mean}, \code{sd}, \code{hdci_lower}, \code{hdci_upper}
#'     summarising the posterior of
#'     \eqn{\rho(\theta_{hurdle}, \theta_{pcm})}. This is the brms
#'     correlation \code{cor_id__Intercept__hu_Intercept} with its
#'     sign flipped to match the "higher = more presence" convention
#'     used for the hurdle person trait (see Details).}
#' }
#'
#' @details
#' \strong{Hurdle submodel.} The hurdle linear predictor is
#' \eqn{logit(hu) = \delta_{hurdle, i} + \tilde{r}_v} with
#' \eqn{hu = P(Y = 0)}; higher \eqn{\delta_{hurdle, i}} means more
#' zeros, so this is reported directly as the item "location" (harder
#' = higher value, standard Rasch convention for the Bernoulli on
#' \eqn{P(Y > 0)}). Hurdle item information is the Bernoulli variance
#' \eqn{p(1-p)} evaluated at \eqn{\theta = \delta_{hurdle, i}}, which
#' is exactly \eqn{1/4}.
#'
#' \strong{Partial credit submodel.} Thresholds \eqn{\tau_{ik}} are
#' the per-item PCM threshold parameters. Item information uses the
#' standard PCM formula; threshold ordering diagnostics are the same
#' as for \code{\link{item_parameters}}.
#'
#' \strong{Sign of the trait correlation.} The brms model reports
#' \code{cor_id__Intercept__hu_Intercept}, which is the correlation
#' between the brms random effects \code{r_id[, Intercept]} and
#' \code{r_id[, hu_Intercept]}. The latter has the opposite sign of
#' the conventional "susceptibility" person trait \eqn{\theta_{hurdle}}
#' (because higher values of the brms random effect mean more zeros,
#' i.e., lower susceptibility). The correlation reported here is
#' therefore \eqn{-\text{cor}_{brms}}, so that positive values mean
#' higher susceptibility goes with higher severity. The marginal SDs
#' are unchanged by the sign flip.
#'
#' \strong{Centering.} When \code{center = TRUE}, the hurdle item
#' difficulties are shifted by their mean and the PCM thresholds are
#' shifted by the grand mean of all PCM thresholds. The same shifts
#' are applied to the corresponding person traits in
#' \code{\link{person_parameters_hpcm}}, preserving the underlying
#' likelihood.
#'
#' @references
#' Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
#' zero-inflated graded response model for ordinal symptom data.
#' \emph{Psychological Methods}, \emph{27}(2), 261-279.
#' \doi{10.1037/met0000395}
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with
#' brms and Stan. \emph{Journal of Statistical Software}, \emph{100},
#' 1--54. \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link{item_parameters}} for the single-submodel version,
#' \code{\link{person_parameters_hpcm}} for the person-side counterpart,
#' \code{\link{hurdle_acat}} for the custom brms family.
#'
#' @importFrom brms as_draws_df
#' @importFrom dplyr arrange
#' @importFrom ggdist mean_hdci
#' @importFrom rlang enquo as_name .data
#' @importFrom stats sd plogis optimize
#' @importFrom tibble tibble as_tibble
#' @export
item_parameters_hpcm <- function(model,
                                 item_var   = item,
                                 person_var = id,
                                 draws      = FALSE,
                                 center     = TRUE,
                                 prob       = 0.95) {

  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  if (!is.numeric(prob) || prob <= 0 || prob >= 1) {
    stop("'prob' must be a numeric value between 0 and 1.",
         call. = FALSE)
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

  all_draws    <- tibble::as_tibble(brms::as_draws_df(model))
  unique_items <- sort(unique(model$data[[item_name]]))

  # ---------------------------------------------------------------
  # PCM submodel: reuse the existing ordinal extractor
  # ---------------------------------------------------------------
  pcm_param_draws <- .extract_item_param_draws(
    all_draws, unique_items, item_name,
    is_ordinal = TRUE, has_fixed_items = FALSE,
    item_levels = .item_levels(model, item_name)
  )

  pcm_shift <- if (center) {
    mean(vapply(unlist(pcm_param_draws, recursive = FALSE),
                mean, numeric(1)))
  } else {
    0
  }
  pcm_param_draws <- lapply(pcm_param_draws, function(item_list) {
    lapply(item_list, function(d) d - pcm_shift)
  })

  pcm_out <- .build_param_tables(
    param_draws = pcm_param_draws,
    submodel    = "pcm",
    prob        = prob,
    person_sd_col = paste0("sd_", person_name, "__Intercept"),
    all_draws   = all_draws
  )

  # ---------------------------------------------------------------
  # Hurdle submodel: extract b_hu_factoritem<X> (or b_hu_<X>)
  # ---------------------------------------------------------------
  hurdle_draws_per_item <- .extract_hurdle_item_draws(
    all_draws, unique_items, item_name
  )

  hu_shift <- if (center) {
    mean(vapply(hurdle_draws_per_item, mean, numeric(1)))
  } else {
    0
  }
  hurdle_param_draws <- lapply(hurdle_draws_per_item,
                               function(d) list(d - hu_shift))
  names(hurdle_param_draws) <- unique_items

  hurdle_out <- .build_param_tables(
    param_draws = hurdle_param_draws,
    submodel    = "hurdle",
    prob        = prob,
    person_sd_col = paste0("sd_", person_name, "__hu_Intercept"),
    all_draws   = all_draws
  )

  # ---------------------------------------------------------------
  # Optional: draws matrices
  # ---------------------------------------------------------------
  if (draws) {
    pcm_out$draws_matrix    <- .pack_draws(pcm_param_draws,
                                           include_threshold_index = TRUE)
    hurdle_out$draws_matrix <- .pack_draws(hurdle_param_draws,
                                           include_threshold_index = FALSE)
  }

  # ---------------------------------------------------------------
  # Cross-submodel correlation (sign-flipped to match theta_hurdle
  # convention)
  # ---------------------------------------------------------------
  cor_col <- grep(
    paste0("^cor_", person_name,
           "__(hu_)?Intercept__(hu_)?Intercept$"),
    names(all_draws), value = TRUE
  )
  if (length(cor_col) >= 1) {
    cor_vals  <- -as.numeric(all_draws[[cor_col[1]]])  # sign flip
    cor_hdci  <- ggdist::mean_hdci(cor_vals, .width = prob)
    correlation <- tibble::tibble(
      mean       = round(mean(cor_vals), 4),
      sd         = round(stats::sd(cor_vals), 4),
      hdci_lower = round(cor_hdci$ymin, 4),
      hdci_upper = round(cor_hdci$ymax, 4)
    )
  } else {
    warning("Could not find brms correlation parameter ",
            "'cor_", person_name, "__Intercept__hu_Intercept'. ",
            "Submodel-trait correlation will not be reported.",
            call. = FALSE)
    correlation <- NULL
  }

  list(hurdle = hurdle_out, pcm = pcm_out, correlation = correlation)
}


# ============================================================================
# 7. Person parameters (per-submodel)
# ============================================================================

#' Extract Person Parameters from a Hurdle Partial Credit Model
#'
#' Extracts person trait estimates from a fitted hurdle Rasch partial
#' credit model (\code{family = \link{hurdle_acat}()}), returning a
#' per-submodel breakdown: a presence trait \eqn{\theta_{hurdle}}
#' (susceptibility) on the Bernoulli gate and a severity trait
#' \eqn{\theta_{pcm}} on the partial credit submodel. Each submodel's
#' output mirrors the shape of \code{\link{person_parameters}} so
#' \code{res$hurdle} and \code{res$pcm} can be passed to existing
#' downstream functions (\code{\link{RMUreliability}},
#' \code{\link{plot_targeting}}, etc.).
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object using the
#'   \code{\link{hurdle_acat}} custom family.
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Default is \code{id}.
#' @param draws Logical. If \code{TRUE}, a matrix of full posterior
#'   draws (persons x draws) is included in each submodel's output.
#'   Default is \code{FALSE}.
#' @param center Logical. If \code{TRUE} (the default), item parameters
#'   and person parameters are shifted within each submodel so that
#'   the mean item parameter is zero. Person traits are shifted by the
#'   same constant as the items in their submodel.
#' @param theta_range A numeric vector of length 2 giving the search
#'   range for the WLE on the hurdle submodel. Default is
#'   \code{c(-7, 7)}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{hurdle}}{A list with the same structure as
#'     \code{\link{person_parameters}} for a dichotomous Rasch model:
#'     \code{person_estimates} (one row per person; \code{sum_score}
#'     is the number of items with \eqn{Y > 0} and \code{n_items} the
#'     number of items answered; both EAP and Warm's WLE are provided,
#'     and with missing responses the WLE uses the items the person
#'     answered), \code{score_table} (persons who answered every
#'     item), and optionally
#'     \code{draws_matrix}. The hurdle EAP is the sign-flipped brms
#'     random effect on \code{hu}, so higher values mean greater
#'     presence / susceptibility.}
#'   \item{\code{pcm}}{A list with \code{person_estimates} (columns
#'     \code{id}, \code{sum_score}, \code{n_active}, \code{eap},
#'     \code{eap_se}), \code{score_table}, and optionally
#'     \code{draws_matrix}. \code{sum_score} is the sum of \eqn{Y}
#'     across items with \eqn{Y_{vi} > 0} for person \eqn{v};
#'     \code{n_active} is the count of such items. WLE columns are
#'     omitted intentionally (see Details).}
#'   \item{\code{correlation}}{A \code{\link[tibble]{tibble}}
#'     summarising the posterior of
#'     \eqn{\rho(\theta_{hurdle}, \theta_{pcm})} with a 95% interval,
#'     the same as the \code{correlation} element returned by
#'     \code{\link{item_parameters_hpcm}} with \code{prob = 0.95}.}
#' }
#'
#' @details
#' \strong{Hurdle person trait.} Defined as
#' \eqn{\theta_{hurdle, v} = -r_{id}(v, \texttt{hu\_Intercept})},
#' i.e., the brms random effect on \code{hu} with its sign flipped.
#' Under this convention, the hurdle submodel reads
#' \eqn{P(Y > 0) = \mathrm{plogis}(\theta_{hurdle, v} -
#' \delta_{hurdle, i})}, so higher \eqn{\theta_{hurdle}} means greater
#' presence / susceptibility. WLE is computed by treating the gate
#' as a dichotomous Rasch model on \eqn{1[Y > 0]} with item
#' difficulties \eqn{\delta_{hurdle, i}} from
#' \code{\link{item_parameters_hpcm}}.
#'
#' \strong{Partial credit person trait.} Defined as
#' \eqn{\theta_{pcm, v} = r_{id}(v, \texttt{Intercept})}. Higher values
#' mean greater severity given presence.
#'
#' \strong{Why no WLE on the partial credit submodel.} Conditional on
#' the hurdle, the severity submodel is a PCM, but the set of items
#' contributing to person \eqn{v}'s severity score depends on
#' \eqn{v}'s own pattern of \eqn{1[Y_{vi} > 0]}: only items with
#' \eqn{Y_{vi} > 0} provide PCM information about \eqn{\theta_{pcm, v}}.
#' The sum of partial credit scores is therefore not a sufficient
#' statistic across the sample, and a standard score-table WLE is not
#' well defined. The Bayesian EAP, which integrates over the
#' (correlated) multivariate prior on
#' \eqn{(\theta_{hurdle}, \theta_{pcm})}, remains well defined for all
#' persons (including those with all-zero patterns) and is therefore
#' the recommended point estimate; see Magnus and Garnier-Villarreal
#' (2022, p. 272ff) for discussion in the closely-related MZI-GRM.
#'
#' \strong{Centering.} When \code{center = TRUE}, the same shifts that
#' \code{\link{item_parameters_hpcm}} applies to the item parameters
#' are applied here to the corresponding person traits. The likelihood
#' is unchanged.
#'
#' \strong{Row ordering.} \code{person_estimates} and
#' \code{draws_matrix} preserve the order of first appearance of each
#' person ID in the model data, matching \code{\link{person_parameters}}.
#'
#' @references
#' Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
#' zero-inflated graded response model for ordinal symptom data.
#' \emph{Psychological Methods}, \emph{27}(2), 261-279.
#' \doi{10.1037/met0000395}
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in
#' item response theory. \emph{Psychometrika}, \emph{54}(3), 427--450.
#' \doi{10.1007/BF02294627}
#'
#' @seealso
#' \code{\link{person_parameters}} for the single-submodel version,
#' \code{\link{item_parameters_hpcm}} for the item-side counterpart,
#' \code{\link{hurdle_acat}} for the custom brms family.
#'
#' @importFrom brms as_draws_df ranef
#' @importFrom dplyr group_by summarise n arrange
#' @importFrom ggdist mean_hdci
#' @importFrom rlang enquo as_name .data
#' @importFrom stats sd as.formula aggregate
#' @importFrom tibble tibble as_tibble
#' @export
person_parameters_hpcm <- function(model,
                                   item_var    = item,
                                   person_var  = id,
                                   draws       = FALSE,
                                   center      = TRUE,
                                   theta_range = c(-7, 7)) {

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
  if (!resp_var %in% names(model$data)) {
    stop("Response variable '", resp_var, "' not found in model data.",
         call. = FALSE)
  }

  all_draws    <- tibble::as_tibble(brms::as_draws_df(model))
  unique_items <- sort(unique(model$data[[item_name]]))
  dat          <- model$data
  obs_response <- .hpcm_response_scores(dat[[resp_var]])

  # ---------------------------------------------------------------
  # Item parameters (needed for centering shifts and hurdle WLE)
  # ---------------------------------------------------------------
  pcm_param_draws <- .extract_item_param_draws(
    all_draws, unique_items, item_name,
    is_ordinal = TRUE, has_fixed_items = FALSE,
    item_levels = .item_levels(model, item_name)
  )
  pcm_thresh_means_uncentered <- lapply(pcm_param_draws, function(item_list) {
    vapply(item_list, mean, numeric(1))
  })

  hurdle_item_draws <- .extract_hurdle_item_draws(
    all_draws, unique_items, item_name
  )
  hurdle_diff_uncentered <- vapply(hurdle_item_draws, mean, numeric(1))

  pcm_shift <- if (center) {
    mean(unlist(pcm_thresh_means_uncentered))
  } else 0
  hu_shift <- if (center) {
    mean(hurdle_diff_uncentered)
  } else 0

  pcm_thresholds <- lapply(pcm_thresh_means_uncentered,
                           function(x) x - pcm_shift)
  hurdle_diff    <- hurdle_diff_uncentered - hu_shift

  # ---------------------------------------------------------------
  # Original person order (first appearance in model data)
  # ---------------------------------------------------------------
  person_col_data <- as.character(dat[[person_name]])
  original_order  <- unique(person_col_data)

  # ---------------------------------------------------------------
  # EAP for both submodels via ranef()
  # ---------------------------------------------------------------
  person_re <- brms::ranef(model)[[person_name]]
  if (is.null(person_re)) {
    stop("No random effects found for '", person_name,
         "'. Check your model formula.", call. = FALSE)
  }
  re_params <- dimnames(person_re)[[3]]
  pcm_param <- intersect(re_params, c("Intercept"))
  hu_param  <- grep("hu_?[Ii]ntercept|hu_Intercept",
                    re_params, value = TRUE)
  if (length(pcm_param) == 0 || length(hu_param) == 0) {
    stop("Could not locate both PCM and hurdle random effects in ",
         "ranef(model). Found: ",
         paste(re_params, collapse = ", "), call. = FALSE)
  }
  ranef_ids <- rownames(person_re[, , pcm_param[1]])

  eap_pcm    <- as.numeric(person_re[, "Estimate",  pcm_param[1]]) - pcm_shift
  eap_pcm_se <- as.numeric(person_re[, "Est.Error", pcm_param[1]])
  # Sign-flip hurdle random effect, then apply the same shift the
  # hurdle item difficulties received (since theta_h - delta_h is the
  # quantity that must be invariant).
  eap_hu     <- -as.numeric(person_re[, "Estimate", hu_param[1]]) - hu_shift
  eap_hu_se  <- as.numeric(person_re[, "Est.Error", hu_param[1]])

  # ---------------------------------------------------------------
  # Sum scores (per submodel)
  # ---------------------------------------------------------------
  hurdle_obs_long <- data.frame(
    person_id = person_col_data,
    is_pos    = as.integer(obs_response > 0),
    pcm_score = ifelse(obs_response > 0,
                       as.integer(obs_response), 0L),
    stringsAsFactors = FALSE
  )

  hurdle_sum <- stats::aggregate(
    is_pos ~ person_id, data = hurdle_obs_long, FUN = sum, na.rm = TRUE
  )
  colnames(hurdle_sum) <- c("person_id", "sum_score")

  pcm_sum <- stats::aggregate(
    cbind(pcm_score, is_pos) ~ person_id,
    data = hurdle_obs_long, FUN = sum, na.rm = TRUE
  )
  colnames(pcm_sum) <- c("person_id", "sum_score", "n_active")

  # ---------------------------------------------------------------
  # Hurdle WLE (dichotomous Rasch on 1[Y > 0])
  # ---------------------------------------------------------------
  # A sum score only has a WLE relative to the items that make it up, so
  # each set of answered items gets its own score-to-WLE table.
  n_it <- length(unique_items)
  persons_u <- unique(person_col_data)
  answered <- matrix(FALSE, nrow = length(persons_u), ncol = n_it)
  answered[cbind(match(person_col_data, persons_u),
                 match(as.character(dat[[item_name]]),
                       as.character(unique_items)))] <- TRUE
  hu_score <- hurdle_sum$sum_score[match(persons_u,
                                         as.character(hurdle_sum$person_id))]
  wb <- .wle_by_pattern(answered, hu_score, function(items) {
    .compute_wle_dichotomous(difficulties = hurdle_diff[items],
                             theta_range  = theta_range)
  })
  hu_wle <- wb$wle
  hu_wle_se <- wb$wle_se
  hurdle_wle_lookup <- tibble::tibble(
    person_id = persons_u,
    n_items   = as.integer(rowSums(answered)),
    wle       = round(hu_wle, 4),
    wle_se    = round(hu_wle_se, 4)
  )

  # ---------------------------------------------------------------
  # Build person tables
  # ---------------------------------------------------------------
  hurdle_df <- tibble::tibble(
    person_id = ranef_ids,
    sum_score = hurdle_sum$sum_score[match(ranef_ids,
                                           as.character(hurdle_sum$person_id))],
    eap       = round(eap_hu,    4),
    eap_se    = round(eap_hu_se, 4)
  )
  hurdle_df <- dplyr::left_join(hurdle_df, hurdle_wle_lookup,
                                by = "person_id")
  hurdle_df <- hurdle_df[, c("person_id", "n_items", "sum_score", "eap",
                             "eap_se", "wle", "wle_se")]
  colnames(hurdle_df)[1] <- person_name
  hurdle_df <- hurdle_df[match(original_order,
                               hurdle_df[[person_name]]), ]

  pcm_df <- tibble::tibble(
    person_id = ranef_ids,
    sum_score = pcm_sum$sum_score[match(ranef_ids,
                                        as.character(pcm_sum$person_id))],
    n_active  = pcm_sum$n_active[match(ranef_ids,
                                       as.character(pcm_sum$person_id))],
    eap       = round(eap_pcm,    4),
    eap_se    = round(eap_pcm_se, 4)
  )
  colnames(pcm_df)[1] <- person_name
  pcm_df <- pcm_df[match(original_order, pcm_df[[person_name]]), ]

  # ---------------------------------------------------------------
  # Score tables
  # ---------------------------------------------------------------
  # Only persons who answered every item share one sum-score scale.
  hurdle_score_table <- hurdle_df[hurdle_df$n_items == n_it, ] |>
    dplyr::group_by(.data$sum_score) |>
    dplyr::summarise(
      n      = dplyr::n(),
      eap    = round(mean(.data$eap,    na.rm = TRUE), 4),
      eap_se = round(mean(.data$eap_se, na.rm = TRUE), 4),
      wle    = round(mean(.data$wle,    na.rm = TRUE), 4),
      wle_se = round(mean(.data$wle_se, na.rm = TRUE), 4),
      .groups = "drop"
    )

  pcm_score_table <- pcm_df |>
    dplyr::group_by(.data$sum_score) |>
    dplyr::summarise(
      n        = dplyr::n(),
      n_active = round(mean(.data$n_active, na.rm = TRUE), 2),
      eap      = round(mean(.data$eap,      na.rm = TRUE), 4),
      eap_se   = round(mean(.data$eap_se,   na.rm = TRUE), 4),
      .groups = "drop"
    )

  hurdle_out <- list(
    person_estimates = hurdle_df,
    score_table      = hurdle_score_table
  )
  pcm_out <- list(
    person_estimates = pcm_df,
    score_table      = pcm_score_table
  )

  # ---------------------------------------------------------------
  # Optional: full draws matrices
  # ---------------------------------------------------------------
  if (draws) {
    hurdle_out$draws_matrix <- .pack_person_draws(
      all_draws, person_name, original_order,
      flip_sign = TRUE, shift = hu_shift, dpar = "hu"
    )
    pcm_out$draws_matrix <- .pack_person_draws(
      all_draws, person_name, original_order,
      flip_sign = FALSE, shift = pcm_shift, dpar = NULL
    )
  }

  # ---------------------------------------------------------------
  # Correlation (sign-flipped to match theta_hurdle convention)
  # ---------------------------------------------------------------
  cor_col <- grep(
    paste0("^cor_", person_name,
           "__(hu_)?Intercept__(hu_)?Intercept$"),
    names(all_draws), value = TRUE
  )
  if (length(cor_col) >= 1) {
    cor_vals <- -as.numeric(all_draws[[cor_col[1]]])
    cor_hdci <- ggdist::mean_hdci(cor_vals, .width = 0.95)
    correlation <- tibble::tibble(
      mean       = round(mean(cor_vals), 4),
      sd         = round(stats::sd(cor_vals), 4),
      hdci_lower = round(cor_hdci$ymin, 4),
      hdci_upper = round(cor_hdci$ymax, 4)
    )
  } else {
    warning("Could not find brms correlation parameter ",
            "'cor_", person_name, "__Intercept__hu_Intercept'. ",
            "Submodel-trait correlation will not be reported.",
            call. = FALSE)
    correlation <- NULL
  }

  list(hurdle = hurdle_out, pcm = pcm_out, correlation = correlation)
}


# ============================================================================
# 8. Internal helpers for hpcm parameter extraction
# ============================================================================

# ── Extract hurdle item draws keyed by item label ────────────────
#' @keywords internal
.extract_hurdle_item_draws <- function(draws, unique_items, item_name) {
  # Hurdle item effects from `hu ~ 0 + factor(item)` or `hu ~ 0 + item`.
  # brms orders them by factor level and may alter labels in the names, so
  # they are matched to items by position (see .check_label_map()).
  levs <- levels(factor(unique_items))
  nm <- names(draws)
  cols <- character(0)
  for (prefix in c(paste0("b_hu_factor", item_name), paste0("b_hu_", item_name))) {
    pat <- paste0("^", prefix, "(.+)$")
    cols <- grep(pat, nm, value = TRUE)
    if (length(cols)) break
  }
  if (length(cols) != length(levs)) {
    stop("Could not find one hurdle parameter per item (found ",
         length(cols), " for ", length(levs), " items). Is the hurdle ",
         "submodel specified as `hu ~ 0 + factor(", item_name, ") + ...` ?",
         call. = FALSE)
  }
  .check_label_map(sub(pat, "\\1", cols), levs, "hurdle item parameters")
  out <- lapply(as.character(unique_items), function(it) {
    as.numeric(draws[[cols[match(it, levs)]]])
  })
  names(out) <- unique_items
  out
}


# ── Build item-level tables shared by both submodels ─────────────
#
# `param_draws` is a named list: one entry per item, each itself a
# list of length 1 (dichotomous / hurdle) or K-1 (PCM) of numeric
# draw vectors.
#' @keywords internal
.build_param_tables <- function(param_draws, submodel, prob,
                                person_sd_col, all_draws) {

  # ---- locations (long) ----
  loc_rows <- list()
  for (item_label in names(param_draws)) {
    thresh_list <- param_draws[[item_label]]
    for (k in seq_along(thresh_list)) {
      loc_rows[[length(loc_rows) + 1]] <- data.frame(
        item      = item_label,
        threshold = as.integer(k),
        location  = round(mean(thresh_list[[k]]), 4),
        stringsAsFactors = FALSE
      )
    }
  }
  locations_long <- tibble::as_tibble(do.call(rbind, loc_rows))

  # ---- locations (wide) ----
  max_thresh <- max(vapply(param_draws, length, integer(1)))
  wide_rows  <- list()
  for (item_label in names(param_draws)) {
    thresh_list  <- param_draws[[item_label]]
    thresh_means <- vapply(thresh_list, mean, numeric(1))
    row <- list(item = item_label)
    if (max_thresh == 1) {
      row[["location"]] <- round(thresh_means[1], 4)
    } else {
      for (k in seq_len(max_thresh)) {
        col_name <- paste0("t", k)
        if (k <= length(thresh_means)) {
          row[[col_name]] <- round(thresh_means[k], 4)
        } else {
          row[[col_name]] <- NA_real_
        }
      }
      row[["location"]] <- round(mean(thresh_means), 4)
    }
    wide_rows[[length(wide_rows) + 1]] <- as.data.frame(
      row, stringsAsFactors = FALSE
    )
  }
  locations_wide <- tibble::as_tibble(do.call(rbind, wide_rows))
  locations_wide <- dplyr::arrange(locations_wide, .data$location)

  # ---- summary (SE, HDCI, n_eff) ----
  sum_rows <- list()
  for (item_label in names(param_draws)) {
    thresh_list <- param_draws[[item_label]]
    for (k in seq_along(thresh_list)) {
      vals  <- thresh_list[[k]]
      hdci  <- ggdist::mean_hdci(vals, .width = prob)
      n_eff <- .compute_ess(vals, all_draws$.chain)
      sum_rows[[length(sum_rows) + 1]] <- data.frame(
        item       = item_label,
        threshold  = as.integer(k),
        location   = round(mean(vals), 4),
        se         = round(stats::sd(vals), 4),
        hdci_lower = round(hdci$ymin, 4),
        hdci_upper = round(hdci$ymax, 4),
        n_eff      = round(n_eff),
        stringsAsFactors = FALSE
      )
    }
  }
  summary_tbl <- tibble::as_tibble(do.call(rbind, sum_rows))

  # ---- item information ----
  info_rows <- list()
  for (item_label in names(param_draws)) {
    thresh_list  <- param_draws[[item_label]]
    thresh_means <- vapply(thresh_list, mean, numeric(1))
    item_loc     <- mean(thresh_means)
    if (length(thresh_means) == 1) {
      p_at_loc    <- stats::plogis(item_loc - thresh_means[1])
      info_at_loc <- p_at_loc * (1 - p_at_loc)
      max_info    <- 0.25
    } else {
      info_at_loc <- .pcm_item_info(item_loc, thresh_means)
      opt <- stats::optimize(
        function(th) .pcm_item_info(th, thresh_means),
        interval = c(min(thresh_means) - 3, max(thresh_means) + 3),
        maximum  = TRUE
      )
      max_info <- opt$objective
    }
    info_rows[[length(info_rows) + 1]] <- data.frame(
      item             = item_label,
      location         = round(item_loc, 4),
      info_at_location = round(info_at_loc, 4),
      max_info         = round(max_info, 4),
      stringsAsFactors = FALSE
    )
  }
  item_info <- tibble::as_tibble(do.call(rbind, info_rows))

  # ---- threshold ordering (PCM only) ----
  threshold_order <- NULL
  if (submodel == "pcm") {
    ord_rows <- list()
    for (item_label in names(param_draws)) {
      thresh_list <- param_draws[[item_label]]
      n_thresh    <- length(thresh_list)
      if (n_thresh < 2) {
        ord_rows[[length(ord_rows) + 1]] <- data.frame(
          item            = item_label,
          n_thresholds    = n_thresh,
          ordered         = TRUE,
          prob_disordered = 0,
          stringsAsFactors = FALSE
        )
        next
      }
      thresh_means <- vapply(thresh_list, mean, numeric(1))
      is_ordered   <- all(diff(thresh_means) > 0)
      n_draws_total <- length(thresh_list[[1]])
      disordered_count <- 0
      for (s in seq_len(n_draws_total)) {
        draw_thresholds <- vapply(thresh_list,
                                  function(x) x[s], numeric(1))
        if (any(diff(draw_thresholds) <= 0)) {
          disordered_count <- disordered_count + 1
        }
      }
      ord_rows[[length(ord_rows) + 1]] <- data.frame(
        item            = item_label,
        n_thresholds    = n_thresh,
        ordered         = is_ordered,
        prob_disordered = round(disordered_count / n_draws_total, 4),
        stringsAsFactors = FALSE
      )
    }
    threshold_order <- tibble::as_tibble(do.call(rbind, ord_rows))
  }

  # ---- person SD ----
  person_sd <- NULL
  if (person_sd_col %in% names(all_draws)) {
    sd_vals <- as.numeric(all_draws[[person_sd_col]])
    sd_hdci <- ggdist::mean_hdci(sd_vals, .width = prob)
    person_sd <- tibble::tibble(
      mean       = round(mean(sd_vals), 4),
      sd         = round(stats::sd(sd_vals), 4),
      hdci_lower = round(sd_hdci$ymin, 4),
      hdci_upper = round(sd_hdci$ymax, 4)
    )
  } else {
    warning("Could not find person-level SD parameter '",
            person_sd_col, "' in the posterior draws.",
            call. = FALSE)
  }

  out <- list(
    locations        = locations_long,
    locations_wide   = locations_wide,
    summary          = summary_tbl,
    item_information = item_info,
    person_sd        = person_sd
  )
  if (submodel == "pcm") out$threshold_order <- threshold_order
  out
}


# ── Pack item draws into a matrix (thresholds x draws) ──────────
#' @keywords internal
.pack_draws <- function(param_draws, include_threshold_index) {
  draw_list <- list()
  row_names <- character(0)
  for (item_label in names(param_draws)) {
    thresh_list <- param_draws[[item_label]]
    for (k in seq_along(thresh_list)) {
      draw_list[[length(draw_list) + 1]] <- thresh_list[[k]]
      row_names <- c(
        row_names,
        if (include_threshold_index && length(thresh_list) > 1) {
          paste0(item_label, "[", k, "]")
        } else {
          item_label
        }
      )
    }
  }
  draws_mat <- do.call(rbind, draw_list)
  rownames(draws_mat) <- row_names
  draws_mat
}


# ── Pack person draws (persons x draws), with optional sign flip ─
#
# `dpar = "hu"` looks for r_<person>[<level>,hu_Intercept] columns;
# `dpar = NULL` looks for r_<person>[<level>,Intercept] columns.
#' @keywords internal
.pack_person_draws <- function(all_draws, person_name, original_order,
                               flip_sign, shift, dpar) {
  if (is.null(dpar)) {
    re_pattern <- paste0("^r_", person_name, "\\[.+,Intercept\\]$")
  } else {
    re_pattern <- paste0("^r_", person_name, "\\[.+,",
                         dpar, "_Intercept\\]$")
  }
  re_cols <- grep(re_pattern, names(all_draws), value = TRUE)
  if (length(re_cols) == 0) {
    # brms sometimes splits correlated random effects into r_<id>
    # vs r_<id>__hu columns (for uncorrelated formulations).
    alt_pattern <- if (is.null(dpar)) {
      paste0("^r_", person_name, "\\[.+,Intercept\\]$")
    } else {
      paste0("^r_", person_name, "__hu\\[.+,Intercept\\]$")
    }
    re_cols <- grep(alt_pattern, names(all_draws), value = TRUE)
    if (length(re_cols) == 0) {
      warning("Could not find posterior draws for person random ",
              "effects on submodel '",
              if (is.null(dpar)) "pcm" else dpar, "'.",
              call. = FALSE)
      return(NULL)
    }
  }

  draws_mat <- t(as.matrix(all_draws[, re_cols, drop = FALSE]))
  if (flip_sign) draws_mat <- -draws_mat
  if (shift != 0) draws_mat <- draws_mat - shift

  rn <- rownames(draws_mat)
  rn <- sub(paste0("^r_", person_name, "(__hu)?\\["), "", rn)
  rn <- sub(",.*\\]$", "", rn)
  rownames(draws_mat) <- rn

  idx <- match(original_order, rn)
  draws_mat <- draws_mat[idx, , drop = FALSE]
  draws_mat
}


# ============================================================================
# 9. Per-submodel RMU reliability
# ============================================================================

#' Relative Measurement Uncertainty Reliability for the Hurdle PCM
#'
#' Computes posterior reliability via Relative Measurement Uncertainty
#' (Bignardi, Kievit & Bürkner, 2026) separately for each of the two
#' person traits in a hurdle partial credit model fitted with the
#' \code{\link{hurdle_acat}} custom family: the presence /
#' susceptibility trait \eqn{\theta_{hurdle}} and the severity trait
#' \eqn{\theta_{pcm}}. Internally extracts per-submodel person draws
#' via \code{\link{person_parameters_hpcm}} and applies
#' \code{\link{RMUreliability}} to each.
#'
#' @param x Either:
#'   \itemize{
#'     \item a fitted \code{\link[brms]{brmsfit}} object using the
#'       \code{\link{hurdle_acat}} family, or
#'     \item the result of
#'       \code{\link{person_parameters_hpcm}(model, draws = TRUE)},
#'       i.e., a list with \code{$hurdle$draws_matrix} and
#'       \code{$pcm$draws_matrix}.
#'   }
#'   Passing a pre-computed \code{person_parameters_hpcm()} result
#'   avoids re-extracting the posterior draws from \code{brms}.
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Used only when \code{x} is a
#'   \code{brmsfit}; ignored otherwise. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Used only when \code{x} is a
#'   \code{brmsfit}; ignored otherwise. Default is \code{id}.
#' @param level Numeric in \eqn{(0, 1)}. Credibility level for the
#'   highest density continuous interval on each reliability. Default
#'   is \code{0.95}.
#' @param verbose Logical. If \code{TRUE}, print the number of subjects
#'   and posterior draws used for each submodel. Default is
#'   \code{FALSE}.
#' @param center Logical. If \code{TRUE} (the default), use centered
#'   person draws (mean item parameter = 0 in each submodel). The
#'   centering shift is a constant per draw and does not affect the
#'   correlation across draw halves, so reliability is invariant to
#'   this choice; the argument is exposed only for consistency with
#'   \code{\link{person_parameters_hpcm}}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{hurdle}}{A data frame with one row containing
#'     \code{rmu_estimate}, \code{hdci_lowerbound}, and
#'     \code{hdci_upperbound} for reliability of
#'     \eqn{\theta_{hurdle}}, as returned by
#'     \code{\link{RMUreliability}}.}
#'   \item{\code{pcm}}{Same format, for reliability of
#'     \eqn{\theta_{pcm}}.}
#' }
#'
#' @details
#' \strong{Two reliabilities, two traits.} The hurdle PCM gives every
#' person a posterior on \eqn{\theta_{hurdle}} (presence /
#' susceptibility) and on \eqn{\theta_{pcm}} (severity given
#' presence). Each posterior has its own measurement uncertainty, and
#' the RMU framework therefore yields one reliability per trait. The
#' two values address different questions: \eqn{\theta_{hurdle}}
#' reliability is about how well the items rank persons by overall
#' endorsement (a function of how many zeros appear and how
#' discriminating the items are at the gate), while
#' \eqn{\theta_{pcm}} reliability is about how well the conditional
#' severity is recovered (a function of how many positive responses
#' each person provides and how spread out the PCM thresholds are).
#' Reporting both is necessary, as emphasised in Magnus and
#' Garnier-Villarreal (2022, p. 277).
#'
#' \strong{What about all-zero responders?} For persons with
#' \eqn{Y = 0} on every item, no items provide information about
#' \eqn{\theta_{pcm}} directly. Their posterior on
#' \eqn{\theta_{pcm}} is shaped by the multivariate-normal prior
#' linking it to \eqn{\theta_{hurdle}} (and is therefore wide but
#' well-defined). These persons still contribute to the RMU
#' calculation; their contribution is appropriately downweighted
#' through their large posterior variance.
#'
#' \strong{Computational note.} When \code{x} is a \code{brmsfit},
#' this function calls \code{\link{person_parameters_hpcm}(..., draws
#' = TRUE)} internally, which extracts the full posterior draws of
#' the person random effects. For large models you may prefer to call
#' \code{\link{person_parameters_hpcm}} once and pass its result to
#' \code{RMUreliability_hpcm}, \code{\link{plot_targeting}}, etc.
#'
#' @references
#' Bignardi, G., Kievit, R., & Bürkner, P.-C. (2026). A general
#' method for estimating reliability using Bayesian Measurement
#' Uncertainty. \emph{PsyArXiv}. \url{https://osf.io/preprints/psyarxiv/h54k8}
#'
#' Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
#' zero-inflated graded response model for ordinal symptom data.
#' \emph{Psychological Methods}, \emph{27}(2), 261-279.
#' \doi{10.1037/met0000395}
#'
#' @seealso
#' \code{\link{RMUreliability}} for the single-trait version,
#' \code{\link{person_parameters_hpcm}} for the underlying person
#' draws extraction.
#'
#' @examples
#' \dontrun{
#' # Fit a hurdle PCM with the hurdle_acat() family
#' fit <- brms::brm(
#'   brms::bf(
#'     response | brms::thres(gr = item) ~ 1 + (1 |g| id),
#'     hu ~ 0 + factor(item) + (1 |g| id)
#'   ),
#'   data     = dat,
#'   family   = hurdle_acat(),
#'   stanvars = hurdle_acat_stanvars()
#' )
#'
#' # Option 1: pass the brmsfit directly
#' rel <- RMUreliability_hpcm(fit)
#' rel$hurdle
#' rel$pcm
#'
#' # Option 2: reuse a cached person_parameters_hpcm() result
#' pp  <- person_parameters_hpcm(fit, draws = TRUE)
#' rel <- RMUreliability_hpcm(pp)
#' }
#'
#' @importFrom rlang enquo inject !!
#' @export
RMUreliability_hpcm <- function(x,
                                item_var   = item,
                                person_var = id,
                                level      = 0.95,
                                verbose    = FALSE,
                                center     = TRUE) {

  if (!is.numeric(level) || length(level) != 1 ||
      level <= 0 || level >= 1) {
    stop("'level' must be a single numeric value in (0, 1).",
         call. = FALSE)
  }

  # ---------------------------------------------------------------
  # Resolve `x` to a person_parameters_hpcm-style list
  # ---------------------------------------------------------------
  if (inherits(x, "brmsfit")) {
    item_quo   <- rlang::enquo(item_var)
    person_quo <- rlang::enquo(person_var)
    pp <- rlang::inject(person_parameters_hpcm(
      model       = x,
      item_var    = !!item_quo,
      person_var  = !!person_quo,
      draws       = TRUE,
      center      = center
    ))
  } else if (is.list(x) &&
             !is.null(x$hurdle) && !is.null(x$pcm) &&
             !is.null(x$hurdle$draws_matrix) &&
             !is.null(x$pcm$draws_matrix)) {
    pp <- x
  } else {
    stop("'x' must be either a brmsfit object (fitted with the ",
         "hurdle_acat() family) or the list returned by ",
         "person_parameters_hpcm(model, draws = TRUE).",
         call. = FALSE)
  }

  if (is.null(pp$hurdle$draws_matrix) ||
      is.null(pp$pcm$draws_matrix)) {
    stop("Draws matrices for the hurdle and pcm submodels could not ",
         "be located. Re-run person_parameters_hpcm() with ",
         "`draws = TRUE` and verify the brms parameter names.",
         call. = FALSE)
  }

  # ---------------------------------------------------------------
  # Per-submodel reliability via RMUreliability()
  # ---------------------------------------------------------------
  list(
    hurdle = RMUreliability(pp$hurdle$draws_matrix,
                            verbose = verbose, level = level),
    pcm    = RMUreliability(pp$pcm$draws_matrix,
                            verbose = verbose, level = level)
  )
}


# ============================================================================
# 10. Item characteristic curves with class intervals (per-submodel)
# ============================================================================

#' Item Characteristic Curves with Class Intervals for Hurdle PCM
#'
#' Plots Item Characteristic Curves (ICCs) with class-interval overlays
#' separately for each submodel of a hurdle partial credit model fitted
#' with the \code{\link{hurdle_acat}} custom family. Returns two
#' \pkg{ggplot2} plots — one for the Bernoulli hurdle on
#' \eqn{P(Y > 0)}, one for the conditional partial credit
#' \eqn{E(Y \mid Y > 0)} — that can be inspected separately or
#' combined with \pkg{patchwork}. Conceptually the Bayesian /
#' hurdle-PCM analogue of \code{\link{plot_icc}}, which only handles
#' single-submodel ordinal / dichotomous IRT.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object using the
#'   \code{\link{hurdle_acat}} custom family.
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Default is \code{id}.
#' @param items An optional character vector of item names to plot.
#'   If \code{NULL} (the default), all items are plotted.
#' @param class_intervals Integer >= 2. The number of class intervals
#'   into which persons are binned within each submodel. Default is 5.
#' @param theta_range A numeric vector of length 2 specifying the
#'   range of theta for the expected curves. Default is
#'   \code{c(-4, 4)}.
#' @param n_points Integer. Number of evenly spaced theta values for
#'   computing the expected curves. Default is 200.
#' @param center Logical. If \code{TRUE} (the default), each submodel
#'   is recentered so that the mean item parameter is zero, matching
#'   \code{\link{item_parameters_hpcm}} and
#'   \code{\link{person_parameters_hpcm}}.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the credible interval
#'   ribbon around each expected curve. Default is 0.95.
#' @param ncol Integer. Number of columns in the faceted layout. If
#'   \code{NULL}, chosen automatically.
#' @param line_size Numeric. Line width for the expected curves.
#'   Default is 0.8.
#' @param ribbon_alpha Numeric in \eqn{[0, 1]}. Transparency of the
#'   credible interval ribbons. Default is 0.3.
#' @param point_size Numeric. Size of observed score points. Default
#'   is 2.5.
#' @param n_intervals Deprecated. Use \code{class_intervals}.
#' @param min_n Integer. Minimum number of observations required in
#'   a class interval (per item, per submodel) for the observed mean
#'   to be plotted. Default is 5.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{hurdle}}{A \code{\link[ggplot2]{ggplot}} object
#'     showing per-item Bernoulli ICCs on the \eqn{P(Y > 0)} scale,
#'     with class-interval observed empirical hurdle-crossing rates
#'     overlaid.}
#'   \item{\code{pcm}}{A \code{\link[ggplot2]{ggplot}} object showing
#'     per-item conditional partial credit ICCs on the
#'     \eqn{E(Y \mid Y > 0)} scale, with class-interval observed mean
#'     severities overlaid (restricted to persons with \eqn{Y > 0} on
#'     that item).}
#' }
#'
#' @details
#' \strong{Hurdle submodel.} For each item the expected curve is
#' \eqn{P(Y_{vi} > 0 \mid \theta_{hurdle}) = \mathrm{plogis}
#' (\theta_{hurdle} - \delta_{hurdle, i})}, computed across posterior
#' draws of the hurdle item difficulty. Persons are binned into class
#' intervals by their hurdle sum score (count of items with
#' \eqn{Y > 0}); each interval is positioned on the
#' \eqn{\theta_{hurdle}} axis by inverting the posterior-mean
#' expected total \eqn{\sum_i P(Y_i > 0 \mid \theta_{hurdle})}, exactly
#' analogous to the procedure in \code{\link{plot_icc}}. The y-axis
#' ranges from 0 to 1.
#'
#' \strong{Partial credit submodel.} For each item the expected curve
#' is the conditional expected severity
#' \eqn{E(Y_{vi} \mid Y_{vi} > 0, \theta_{pcm}) =
#' \sum_{c = 1}^{K - 1} c \cdot P(Y_{vi} = c \mid Y_{vi} > 0)},
#' with conditional PCM category probabilities computed from the
#' posterior draws of \eqn{\tau_{ik}}. The y-axis ranges from 1 to
#' \eqn{K - 1}.
#'
#' Persons are binned into class intervals by their \strong{EAP}
#' \eqn{\theta_{pcm}} (taken from \code{\link[brms]{ranef}}), not by
#' a sum score. This is intentional: under the hurdle PCM the sum of
#' \eqn{Y} over items with \eqn{Y_{vi} > 0} is not a sufficient
#' statistic for \eqn{\theta_{pcm}} (because the number of contributing
#' items varies across persons; see
#' \code{\link{person_parameters_hpcm}}). EAP-based binning sidesteps
#' the inverse-CDF construction and avoids privileging any particular
#' summary score. Within each interval, the observed mean response for
#' an item is computed \emph{only} over persons with \eqn{Y > 0} on
#' that item, the cells the PCM submodel applies to.
#'
#' EAPs are shrunk toward the mean and are estimated from the same
#' responses that are plotted, so a class interval placed at its mean
#' EAP would put the observed points off the curve under good fit.
#' Instead, for each item and interval the model-expected mean severity
#' for the same responses is computed from \code{posterior_epred()},
#' averaged over each person's posterior (up to 200 draws), and the
#' point is placed at the \eqn{\theta_{pcm}} where the item's expected
#' curve equals that mean. This is the analogue of placing intervals by
#' inverting the expected total score in \code{\link{plot_icc}}. The
#' same interval can therefore sit at slightly different
#' \eqn{\theta_{pcm}} values for different items.
#'
#' \strong{Uncertainty.} Each expected curve carries a credible
#' interval ribbon at width \code{prob}. Observed points have
#' \eqn{\pm 1.96\,SE} error bars where SE is the within-bin sampling
#' standard error of the mean (binomial for the hurdle, ordinal for
#' the partial credit submodel).
#'
#' \strong{Why two separate plots and not patchwork-combined.} The two
#' submodels have qualitatively different y-axes (probability vs.
#' expected severity) and different x-axis interpretations (presence
#' trait vs. severity trait). Returning them as a list preserves the
#' option to inspect each in isolation; combine them yourself with
#' \pkg{patchwork} when a side-by-side or stacked layout is wanted:
#' \preformatted{
#' plots <- plot_icc_hpcm(fit)
#' patchwork::wrap_plots(plots$hurdle, plots$pcm, ncol = 1)
#' }
#'
#' @references
#' Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
#' zero-inflated graded response model for ordinal symptom data.
#' \emph{Psychological Methods}, \emph{27}(2), 261-279.
#' \doi{10.1037/met0000395}
#'
#' @seealso
#' \code{\link{plot_icc}} for the single-submodel ICC plot,
#' \code{\link{item_parameters_hpcm}},
#' \code{\link{person_parameters_hpcm}}.
#'
#' @importFrom brms as_draws_df ranef
#' @importFrom rlang enquo as_name .data
#' @importFrom stats family formula quantile plogis sd aggregate
#' @importFrom tibble as_tibble
#' @export
plot_icc_hpcm <- function(model,
                           item_var     = item,
                           person_var   = id,
                           items        = NULL,
                           class_intervals = 5,
                           theta_range  = c(-4, 4),
                           n_points     = 200,
                           center       = TRUE,
                           prob         = 0.95,
                           ncol         = NULL,
                           line_size    = 0.8,
                           ribbon_alpha = 0.3,
                           point_size   = 2.5,
                           min_n        = 5,
                           n_intervals  = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  if (!is.null(n_intervals)) {
    warning("`n_intervals` is deprecated. Use `class_intervals`.",
            call. = FALSE)
    class_intervals <- n_intervals
  }
  if (!is.numeric(class_intervals) || length(class_intervals) != 1L ||
      !is.finite(class_intervals) || class_intervals < 2) {
    stop("`class_intervals` must be a single integer >= 2.", call. = FALSE)
  }

  # Some brms versions report a custom family's family() as the literal
  # name (e.g., "hurdle_acat"); others report it as "custom". Accept both.
  family_name <- stats::family(model)$family
  if (!family_name %in% c("hurdle_acat", "custom")) {
    stop("plot_icc_hpcm() requires a model fitted with the ",
         "hurdle_acat() custom family. Got: '", family_name, "'.",
         call. = FALSE)
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

  model_data        <- model$data
  obs_response_raw  <- model_data[[resp_var]]
  # Responses 0..K-1 with 0 the hurdle category, coded as brms does
  obs_response <- .hpcm_response_scores(obs_response_raw)
  unique_items <- sort(unique(model_data[[item_name]]))

  if (!is.null(items)) {
    invalid <- setdiff(items, unique_items)
    if (length(invalid) > 0) {
      stop("Items not found in model data: ",
           paste(invalid, collapse = ", "), call. = FALSE)
    }
    plot_items <- unique_items[unique_items %in% items]
  } else {
    plot_items <- unique_items
  }

  # ---------------------------------------------------------------
  # Posterior draws and centering shifts
  # ---------------------------------------------------------------
  draws <- tibble::as_tibble(brms::as_draws_df(model))

  hu_draws_per_item <- .extract_hurdle_item_draws(
    draws, unique_items, item_name
  )
  hu_diff_means <- vapply(hu_draws_per_item, mean, numeric(1))
  hu_shift      <- if (center) mean(hu_diff_means) else 0
  hu_diff_centered <- hu_diff_means - hu_shift

  # PCM thresholds per item, matched to items by brms level order
  pcm_tau <- .item_tau_draws(draws, levels(factor(model_data[[item_name]])),
                             item_name, ordinal = TRUE)
  pcm_thresh_means_uncentered <- lapply(as.character(unique_items),
                                        function(it) colMeans(pcm_tau[[it]]))
  names(pcm_thresh_means_uncentered) <- unique_items
  pcm_shift <- if (center) {
    mean(unlist(pcm_thresh_means_uncentered))
  } else 0

  # Maximum severity score K - 1 (max number of PCM thresholds across
  # items + 1). Used as the upper clamp for the conditional-severity
  # error bars; the lower clamp is 1 because Y > 0 by construction in
  # the PCM submodel.
  pcm_max_cat <- max(vapply(
    pcm_thresh_means_uncentered, length, integer(1)
  )) + 1L

  lower_prob <- (1 - prob) / 2
  upper_prob <- 1 - lower_prob

  theta_grid      <- seq(theta_range[1], theta_range[2],
                         length.out = n_points)
  theta_grid_fine <- seq(theta_range[1], theta_range[2],
                         length.out = 1000)

  n_draws_total <- nrow(draws)
  max_draws     <- min(n_draws_total, 500)
  draw_sample   <- sample(seq_len(n_draws_total), max_draws)

  # ---------------------------------------------------------------
  # Per-person summary statistics
  # ---------------------------------------------------------------
  person_col <- as.character(model_data[[person_name]])
  unique_persons <- unique(person_col)

  hurdle_sum_per <- vapply(unique_persons, function(pid) {
    sum(obs_response[person_col == pid] > 0)
  }, integer(1))
  hurdle_sum_df <- data.frame(
    person_id  = unique_persons,
    sum_hurdle = hurdle_sum_per,
    stringsAsFactors = FALSE
  )
  # Hurdle class intervals compare counts of Y > 0 across persons, which
  # requires the same items for everyone: complete responders only.
  n_obs_per <- tabulate(match(person_col, unique_persons),
                        nbins = length(unique_persons))
  complete_hu <- n_obs_per == length(unique_items)
  if (!all(complete_hu)) {
    message("plot_icc_hpcm() uses the ", sum(complete_hu), " persons who ",
            "answered every item for the hurdle panel. ", sum(!complete_hu),
            " with missing responses are left out.")
    hurdle_sum_df <- hurdle_sum_df[complete_hu, , drop = FALSE]
  }

  pcm_n_active <- vapply(unique_persons, function(pid) {
    sum(obs_response[person_col == pid] > 0)
  }, integer(1))

  # Person EAP on the PCM trait (for class-interval binning)
  person_re <- brms::ranef(model)[[person_name]]
  if (is.null(person_re)) {
    stop("No random effects found for '", person_name, "'.",
         call. = FALSE)
  }
  re_params <- dimnames(person_re)[[3]]
  if (!"Intercept" %in% re_params) {
    stop("Could not find PCM random effect 'Intercept' in ranef().",
         call. = FALSE)
  }
  pcm_eap_raw <- as.numeric(person_re[, "Estimate", "Intercept"])
  pcm_eap <- setNames(pcm_eap_raw - pcm_shift,
                      rownames(person_re[, , "Intercept"]))

  # ---------------------------------------------------------------
  # Expected curves (per item) and total expected hurdle count
  # ---------------------------------------------------------------
  hu_expected_list  <- list()
  pcm_expected_list <- list()

  for (item_label in plot_items) {

    # --- Hurdle expected curve ---
    delta_draws <- hu_draws_per_item[[item_label]][draw_sample] -
      hu_shift
    e_draws_hu <- matrix(NA_real_, nrow = max_draws, ncol = n_points)
    for (t in seq_along(theta_grid)) {
      e_draws_hu[, t] <- stats::plogis(theta_grid[t] - delta_draws)
    }
    hu_expected_list[[item_label]] <- data.frame(
      item    = item_label,
      theta   = theta_grid,
      e_score = colMeans(e_draws_hu),
      e_lower = apply(e_draws_hu, 2, stats::quantile,
                      probs = lower_prob),
      e_upper = apply(e_draws_hu, 2, stats::quantile,
                      probs = upper_prob),
      stringsAsFactors = FALSE
    )

    # --- PCM expected curve (conditional on Y > 0) ---
    thresh_sub <- pcm_tau[[as.character(item_label)]][draw_sample, ,
                                                      drop = FALSE] - pcm_shift
    n_thresh   <- ncol(thresh_sub)
    cats_pos   <- seq_len(n_thresh + 1L)   # severity scores 1..K-1

    e_draws_pcm <- matrix(NA_real_, nrow = max_draws, ncol = n_points)
    for (t in seq_along(theta_grid)) {
      theta <- theta_grid[t]
      for (s in seq_len(max_draws)) {
        p_cat <- .compute_cat_probs(
          theta, thresh_sub[s, ],
          is_acat   = TRUE,  is_cumul = FALSE,
          is_sratio = FALSE, is_cratio = FALSE
        )
        e_draws_pcm[s, t] <- sum(cats_pos * p_cat)
      }
    }
    pcm_expected_list[[item_label]] <- data.frame(
      item    = item_label,
      theta   = theta_grid,
      e_score = colMeans(e_draws_pcm),
      e_lower = apply(e_draws_pcm, 2, stats::quantile,
                      probs = lower_prob),
      e_upper = apply(e_draws_pcm, 2, stats::quantile,
                      probs = upper_prob),
      stringsAsFactors = FALSE
    )
  }
  hu_expected_data  <- do.call(rbind, hu_expected_list)
  pcm_expected_data <- do.call(rbind, pcm_expected_list)

  # Expected total hurdle count across theta_grid_fine
  e_total_hu <- vapply(theta_grid_fine, function(t) {
    sum(stats::plogis(t - hu_diff_centered))
  }, numeric(1))

  # ---------------------------------------------------------------
  # Hurdle class intervals (sum-score-based, with inverse mapping)
  # ---------------------------------------------------------------
  max_hu_score <- length(unique_items)
  hu_ci <- hurdle_sum_df[
    hurdle_sum_df$sum_hurdle > 0 &
      hurdle_sum_df$sum_hurdle < max_hu_score, ,
    drop = FALSE
  ]

  hu_observed_data <- .build_hpcm_observed(
    person_ci      = hu_ci,
    sum_col        = "sum_hurdle",
    bin_to_theta   = function(s) {
      idx <- which.min(abs(e_total_hu - s))
      theta_grid_fine[idx]
    },
    class_intervals = class_intervals,
    plot_items     = plot_items,
    model_data     = model_data,
    item_name      = item_name,
    person_name    = person_name,
    obs_response   = obs_response,
    response_fn    = function(y) as.numeric(y > 0),
    min_n          = min_n,
    restrict_to_y_pos = FALSE
  )

  # ---------------------------------------------------------------
  # PCM class intervals (EAP theta_pcm based)
  # ---------------------------------------------------------------
  pcm_ci <- data.frame(
    person_id = unique_persons,
    n_active  = pcm_n_active,
    eap_theta = pcm_eap[unique_persons],
    stringsAsFactors = FALSE
  )
  pcm_ci <- pcm_ci[pcm_ci$n_active > 0 & !is.na(pcm_ci$eap_theta), ,
                   drop = FALSE]

  # Model-expected severity for the same observations, averaged over
  # each person's posterior: E(Y | Y > 0) per draw from posterior_epred,
  # then the posterior mean. Comparing the observed bin means with these
  # expected bin means does not depend on where the bin is placed on the
  # x-axis, so EAP shrinkage does not distort the comparison.
  epred_draws <- draw_sample[seq_len(min(length(draw_sample), 200L))]
  epred <- brms::posterior_epred(model, draw_ids = epred_draws)
  n_cat_total <- dim(epred)[3L]
  p_pos <- 1 - epred[, , 1L, drop = FALSE][, , 1L]
  e_sev <- 0
  for (k in seq_len(n_cat_total - 1L)) {
    e_sev <- e_sev + k * epred[, , k + 1L]
  }
  pcm_expected_obs <- colMeans(matrix(e_sev / p_pos, nrow = length(epred_draws)))
  rm(epred)

  pcm_observed_data <- .build_hpcm_observed(
    person_ci      = pcm_ci,
    sum_col        = "eap_theta",
    bin_to_theta   = identity,  # bins already in theta_pcm scale
    class_intervals = class_intervals,
    plot_items     = plot_items,
    model_data     = model_data,
    item_name      = item_name,
    person_name    = person_name,
    obs_response   = obs_response,
    response_fn    = function(y) as.numeric(y),
    min_n          = min_n,
    restrict_to_y_pos = TRUE,
    expected_obs   = pcm_expected_obs
  )

  # Place each severity point where the item's expected curve equals the
  # model-expected mean for the same responses, as plot_icc() places
  # class intervals by inverting the expected total score. The expected
  # mean then lies on the curve, and the observed mean is compared with
  # the curve. Bins on shrunk EAPs would otherwise sit at the wrong x.
  if (!is.null(pcm_observed_data) && nrow(pcm_observed_data) > 0) {
    for (item_label in unique(pcm_observed_data$item)) {
      rows  <- which(pcm_observed_data$item == item_label)
      curve <- pcm_expected_data[pcm_expected_data$item == item_label, ,
                                 drop = FALSE]
      pcm_observed_data$theta[rows] <- .hpcm_invert_curve(
        pcm_observed_data$expected[rows], curve$theta, curve$e_score
      )
    }
  }

  # ---------------------------------------------------------------
  # Build ggplot objects
  # ---------------------------------------------------------------
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(length(plot_items)))
  }

  build_plot <- function(expected_data, observed_data,
                         y_label,
                         y_lim_min   = NA_real_,
                         y_lim_max   = NA_real_,
                         error_floor = -Inf,
                         error_ceil  =  Inf) {
    if (is.null(expected_data) || nrow(expected_data) == 0) {
      return(ggplot2::ggplot() +
               ggplot2::labs(title = "(no data)"))
    }
    expected_data$item <- factor(expected_data$item, levels = plot_items)

    # Precompute error-bar bounds and drop anything outside the
    # legitimate response range. The bounds are stored as data columns
    # (rather than evaluated inside aes() from the build_plot scope)
    # so ggplot's quosure-based aes evaluation never has to look up
    # `error_floor` / `error_ceil` lexically, which sidesteps a class
    # of subtle scope bugs that can otherwise leave a stray cap drawn
    # at y = 0 in sparse class intervals.
    if (!is.null(observed_data) && nrow(observed_data) > 0) {
      observed_data$item <- factor(observed_data$item,
                                   levels = plot_items)
      observed_data <- observed_data[
        is.finite(observed_data$response) &
          observed_data$response >= error_floor &
          observed_data$response <= error_ceil, ,
        drop = FALSE
      ]
      observed_data$ymin <- pmax(
        observed_data$response - 1.96 * observed_data$se, error_floor
      )
      observed_data$ymax <- pmin(
        observed_data$response + 1.96 * observed_data$se, error_ceil
      )
    }

    p <- ggplot2::ggplot()
    if (ribbon_alpha > 0) {
      p <- p + ggplot2::geom_ribbon(
        data = expected_data,
        mapping = ggplot2::aes(
          x = .data$theta,
          ymin = .data$e_lower,
          ymax = .data$e_upper
        ),
        fill = "grey70", alpha = ribbon_alpha
      )
    }
    p <- p + ggplot2::geom_line(
      data = expected_data,
      mapping = ggplot2::aes(x = .data$theta, y = .data$e_score),
      colour = "grey40", linewidth = line_size
    )

    if (!is.null(observed_data) && nrow(observed_data) > 0) {
      p <- p +
        ggplot2::geom_errorbar(
          data = observed_data,
          mapping = ggplot2::aes(
            x    = .data$theta,
            ymin = .data$ymin,
            ymax = .data$ymax
          ),
          colour = "#1F78B4", width = 0.05, linewidth = 0.4
        ) +
        ggplot2::geom_point(
          data = observed_data,
          mapping = ggplot2::aes(x = .data$theta, y = .data$response),
          colour = "#1F78B4", size = point_size
        )
    }

    p +
      ggplot2::facet_wrap(~ item, ncol = ncol, scales = "free_y") +
      ggplot2::scale_y_continuous(limits = c(y_lim_min, y_lim_max)) +
      ggplot2::labs(x = expression(theta), y = y_label) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        strip.text = ggplot2::element_text(face = "bold")
      )
  }

  hurdle_plot <- build_plot(
    expected_data = hu_expected_data,
    observed_data = hu_observed_data,
    y_label       = expression(P(Y > 0)),
    y_lim_min     = 0,
    y_lim_max     = 1,
    error_floor   = 0,
    error_ceil    = 1
  )

  pcm_plot <- build_plot(
    expected_data = pcm_expected_data,
    observed_data = pcm_observed_data,
    y_label       = expression(E(Y * "|" * Y > 0)),
    # Let the y-axis auto-scale (with ggplot's default 5% expansion)
    # rather than hard-clipping at the response floor of 1. A hard
    # `limits = c(1, NA)` glues observed points at the floor — common
    # in the leftmost class interval where every active responder
    # gave Y = 1 — to the bottom edge of the panel, which reads
    # visually as "y = 0". The defensive observed-data filter inside
    # build_plot (response >= error_floor) plus the error-bar clamps
    # (pmax/pmin to [error_floor, error_ceil]) already prevent any
    # spurious below-floor values, so removing the hard scale clip is
    # safe.
    y_lim_min     = NA_real_,
    y_lim_max     = NA_real_,
    error_floor   = 1,
    error_ceil    = pcm_max_cat
  ) +
    ggplot2::labs(caption = .erb_caption(paste(
      "Points: observed mean severity in each class interval of the",
      "severity EAP, among responses above 0, with 95% intervals. Each",
      "point is placed where the curve equals the model-expected mean",
      "for the same responses."))) +
    .erb_caption_theme()

  list(hurdle = hurdle_plot, pcm = pcm_plot)
}


# Invert a monotone expected-score curve: the theta at which the curve
# equals each value in `target`. Values outside the curve's range are
# placed at the nearest end of the theta grid.
#' @keywords internal
.hpcm_invert_curve <- function(target, theta, e_score) {
  ord <- order(theta)
  theta <- theta[ord]
  e_score <- cummax(e_score[ord])  # guard against numerical wiggles
  out <- rep(NA_real_, length(target))
  ok <- is.finite(target)
  out[ok] <- stats::approx(e_score, theta, xout = target[ok],
                           rule = 2, ties = mean)$y
  out
}


# ── Build observed-points data frame for one submodel ────────────
#
# person_ci: data.frame with `person_id`, plus a column named by
#   `sum_col` used to define the class intervals.
# bin_to_theta: function mapping the bin mean (in sum_col units) to
#   a theta-axis position.
# response_fn: function applied to the raw response before computing
#   per-bin means (e.g., 1[y > 0] for hurdle; identity for PCM).
# restrict_to_y_pos: if TRUE, only person x item cells with Y > 0
#   contribute (PCM submodel only).
#' @keywords internal
.build_hpcm_observed <- function(person_ci, sum_col, bin_to_theta,
                                 class_intervals, plot_items, model_data,
                                 item_name, person_name, obs_response,
                                 response_fn, min_n,
                                 restrict_to_y_pos,
                                 expected_obs = NULL) {
  if (nrow(person_ci) == 0) return(NULL)

  values <- person_ci[[sum_col]]
  breaks <- stats::quantile(values,
                            probs = seq(0, 1, length.out = class_intervals + 1))
  breaks <- unique(breaks)

  if (length(breaks) >= 3) {
    # Standard path: enough distinct quantile breaks to define
    # 2+ class intervals. Used by the PCM submodel (EAP-theta is
    # near-continuous) and by hurdle submodels with many items.
    person_ci$interval <- cut(values, breaks = breaks,
                              include.lowest = TRUE, labels = FALSE)
    bin_means <- stats::aggregate(
      person_ci[[sum_col]], by = list(interval = person_ci$interval),
      FUN = mean
    )
    colnames(bin_means) <- c("interval", "bin_mean")
  } else {
    # Fallback: `values` has too few distinct values for quantile-
    # based binning. Typical for the hurdle submodel of a short
    # instrument: with K items the non-extreme hurdle sum scores are
    # {1, ..., K-1}, so for small K (e.g., K = 3) we can have only
    # 1-2 unique scores and quantile breaks collapse. Use each
    # unique value as its own class interval.
    unique_vals <- sort(unique(values))
    if (length(unique_vals) < 2L) return(NULL)
    person_ci$interval <- match(values, unique_vals)
    bin_means <- data.frame(
      interval = seq_along(unique_vals),
      bin_mean = unique_vals,
      stringsAsFactors = FALSE
    )
  }

  bin_means$theta <- vapply(bin_means$bin_mean,
                            function(s) bin_to_theta(s),
                            numeric(1))

  person_ci <- merge(person_ci,
                     bin_means[, c("interval", "theta")],
                     by = "interval")

  obs_list <- list()
  for (item_label in plot_items) {
    item_rows  <- which(model_data[[item_name]] == item_label)
    if (length(item_rows) == 0) next
    item_persons <- as.character(model_data[[person_name]][item_rows])
    item_y       <- as.numeric(obs_response[item_rows])
    item_response <- response_fn(item_y)
    item_expected <- if (is.null(expected_obs)) {
      rep(NA_real_, length(item_rows))
    } else {
      expected_obs[item_rows]
    }

    if (restrict_to_y_pos) {
      keep <- item_y > 0
      if (!any(keep)) next
      item_persons  <- item_persons[keep]
      item_response <- item_response[keep]
      item_expected <- item_expected[keep]
    }

    item_df <- data.frame(
      person_id = item_persons,
      response  = item_response,
      expected  = item_expected,
      stringsAsFactors = FALSE
    )
    item_df <- merge(item_df, person_ci, by = "person_id")
    if (nrow(item_df) == 0) next

    agg <- do.call(rbind, lapply(
      split(item_df, list(item_df$interval, item_df$theta),
            drop = TRUE),
      function(d) {
        n_obs <- sum(!is.na(d$response))
        data.frame(
          interval = d$interval[1],
          theta    = d$theta[1],
          response = mean(d$response, na.rm = TRUE),
          expected = mean(d$expected),
          n        = n_obs,
          se       = if (n_obs > 1L) {
            stats::sd(d$response, na.rm = TRUE) / sqrt(n_obs)
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      }
    ))
    agg <- agg[agg$n >= min_n, , drop = FALSE]
    if (nrow(agg) == 0) next
    agg$item <- item_label
    obs_list[[item_label]] <- agg
  }

  if (length(obs_list) == 0) return(NULL)
  do.call(rbind, obs_list)
}


# ============================================================================
# 11. Person-item targeting plot (per-submodel)
# ============================================================================

#' Person-Item Targeting Plot for the Hurdle Partial Credit Model
#'
#' Builds two person-item targeting plots, one for each part of a hurdle
#' partial credit model fitted with the \code{\link{hurdle_acat}} custom
#' family, with the same design as \code{\link{plot_targeting}}: a person
#' histogram, an inverted histogram of the item locations, and a bottom
#' panel with response-category bands or a dot-and-whisker plot of the
#' thresholds. The presence part \eqn{(\theta_{hurdle}, \delta_{hurdle})}
#' and the severity part \eqn{(\theta_{pcm}, \tau_{ik})} are different
#' latent variables, so each gets its own plot and axis.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object using the
#'   \code{\link{hurdle_acat}} custom family.
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Default is \code{id}.
#' @param robust Logical. If \code{TRUE}, the histogram panels mark the
#'   median and MAD instead of the mean and SD. Default is \code{FALSE}.
#' @param center Logical. If \code{TRUE} (the default), each part is
#'   shifted so that its mean item location is zero, with the person
#'   locations shifted by the same constant, as in
#'   \code{\link{item_parameters_hpcm}} and
#'   \code{\link{person_parameters_hpcm}}.
#' @param sort_items Character. \code{"data"} (the default) keeps the
#'   data order, first item at the top. \code{"location"} sorts by mean
#'   item location, easiest item at the top.
#' @param bins Integer. Number of histogram bins. Default as in
#'   \code{\link{plot_targeting}}.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the credible intervals of
#'   the item locations and of the intervals in the returned tables.
#'   Default is 0.95. \code{NULL} hides the intervals in the plots.
#' @param palette An optional character vector of colours, one per
#'   response category. Overrides the viridis colours.
#' @param person_fill,threshold_fill Fill colours for the person and item
#'   location histograms.
#' @param height_ratios Numeric vector of length 3 giving the relative
#'   heights of the three panels. Default is \code{c(3, 2, 5)}.
#' @param panel Character. \code{"categories"} (the default) draws each
#'   item as a bar partitioned into response-category bands.
#'   \code{"thresholds"} draws a dot-and-whisker plot.
#' @param xlim Numeric vector of length 2. Initial limits of the x-axis,
#'   widened to whole numbers if needed. Default is \code{c(-4, 4)}.
#' @param category_labels Optional list with elements \code{hurdle} (two
#'   labels) and/or \code{pcm} (one label per severity category, starting
#'   at category 1). Defaults are \code{c("0", "1+")} for the hurdle part
#'   and the category scores for the severity part.
#' @param viridis_option,viridis_begin,viridis_end Viridis palette for the
#'   category colours, as in \code{\link{plot_targeting}}.
#' @param row_gap Numeric. Vertical spacing between item rows in the
#'   \code{"categories"} panel. Default \code{NULL} chooses it from the
#'   number of collapsed categories and disordered thresholds.
#' @param output Character. \code{"patchwork"} (the default) returns the
#'   two combined plots. \code{"list"} returns, for each part, the three
#'   panels and the posterior summaries of the category structure.
#'
#' @return With \code{output = "patchwork"}, a list with elements
#'   \code{hurdle} and \code{pcm}, each a \pkg{patchwork} object. With
#'   \code{output = "list"}, the same two elements, each a list with
#'   \code{p1}, \code{p2}, \code{p3}, \code{categories} and
#'   \code{disordered} as described in \code{\link{plot_targeting}}. In the
#'   severity part, categories are numbered from 1.
#'
#' @details
#' \strong{Hurdle part.} The person location is the posterior mean of the
#' brms random effect on \code{hu} with its sign flipped, so that higher
#' values mean a symptom is more likely to be present. The item location
#' is the hurdle item effect, where higher values mean more zeros. Then
#' \eqn{P(Y_{vi} > 0) = \mathrm{plogis}(\theta_{hurdle, v} -
#' \delta_{hurdle, i})}, a dichotomous Rasch model, and each item has two
#' bands: 0 and any response above 0.
#'
#' \strong{Severity part.} The person location is the posterior mean of
#' the random intercept, and the thresholds are those of the partial
#' credit model for the responses above 0. The bands show the severity
#' categories 1 and up. Persons who answered 0 on every item gave no
#' severity responses, so their severity location comes from the prior
#' and the correlation with the hurdle location only. They are left out
#' of the severity person histogram, and the caption gives their number.
#'
#' Category bands, red marks for categories that are never the most
#' likely response and for disordered thresholds, and their posterior
#' probabilities are computed as in \code{\link{plot_targeting}}.
#'
#' \strong{Combining the two plots.} Each element is a \pkg{patchwork}
#' object:
#' \preformatted{
#' plots <- plot_targeting_hpcm(fit)
#' patchwork::wrap_plots(plots$hurdle, plots$pcm, ncol = 2)
#' }
#'
#' @references
#' Wright, B. D. & Stone, M. H. (1979). \emph{Best Test Design}.
#' MESA Press.
#'
#' Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
#' zero-inflated graded response model for ordinal symptom data.
#' \emph{Psychological Methods}, \emph{27}(2), 261-279.
#' \doi{10.1037/met0000395}
#'
#' @seealso
#' \code{\link{plot_targeting}} for the standard Rasch models,
#' \code{\link{item_parameters_hpcm}},
#' \code{\link{person_parameters_hpcm}}.
#'
#' @importFrom brms as_draws_df ranef
#' @importFrom rlang enquo as_name .data
#' @importFrom stats family
#' @export
plot_targeting_hpcm <- function(model,
                                item_var        = item,
                                person_var      = id,
                                robust          = FALSE,
                                center          = TRUE,
                                sort_items      = c("data", "location"),
                                bins,
                                prob            = 0.95,
                                palette         = NULL,
                                person_fill     = "#0072B2",
                                threshold_fill  = "#D55E00",
                                height_ratios   = c(3, 2, 5),
                                panel           = c("categories", "thresholds"),
                                xlim            = c(-4, 4),
                                category_labels = NULL,
                                viridis_option  = "G",
                                viridis_begin   = 0.9,
                                viridis_end     = 0.2,
                                row_gap         = NULL,
                                output          = c("patchwork", "list")) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.", call. = FALSE)
  }
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }

  # Some brms versions report a custom family's family() as the literal
  # name (e.g., "hurdle_acat"); others report it as "custom". Accept both.
  family_name <- stats::family(model)$family
  if (!family_name %in% c("hurdle_acat", "custom")) {
    stop("plot_targeting_hpcm() requires a model fitted with the ",
         "hurdle_acat() custom family. Got: '", family_name, "'.",
         call. = FALSE)
  }

  sort_items  <- match.arg(sort_items)
  panel       <- match.arg(panel)
  output      <- match.arg(output)
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
  if (!is.null(category_labels) &&
      (!is.list(category_labels) ||
       !all(names(category_labels) %in% c("hurdle", "pcm")))) {
    stop("`category_labels` must be a list with elements `hurdle` and/or ",
         "`pcm`.", call. = FALSE)
  }
  show_ci <- !is.null(prob)
  ci_prob <- if (show_ci) prob else 0.95
  if (missing(bins)) bins <- NULL
  data_order <- as.character(unique(model$data[[item_name]]))

  # ---------------------------------------------------------------
  # Item parameter draws per part, centred
  # ---------------------------------------------------------------
  draws <- brms::as_draws_df(model)
  item_levels <- .item_levels(model, item_name)
  tau_pcm <- .item_tau_draws(draws, item_levels, item_name, TRUE)
  hu_draws <- .extract_hurdle_item_draws(draws, item_levels, item_name)
  tau_hu <- lapply(hu_draws, function(v) matrix(v, ncol = 1L))

  pcm_shift <- if (center) mean(unlist(lapply(tau_pcm, colMeans))) else 0
  hu_shift  <- if (center) mean(vapply(hu_draws, mean, numeric(1))) else 0
  tau_pcm <- lapply(tau_pcm, function(m) m - pcm_shift)
  tau_hu  <- lapply(tau_hu,  function(m) m - hu_shift)

  # ---------------------------------------------------------------
  # Person locations (EAP)
  # ---------------------------------------------------------------
  person_re <- brms::ranef(model)[[person_name]]
  if (is.null(person_re)) {
    stop("No random effects found for '", person_name, "'.",
         call. = FALSE)
  }
  re_params <- dimnames(person_re)[[3]]
  hu_param  <- grep("hu_?[Ii]ntercept", re_params, value = TRUE)
  if (!("Intercept" %in% re_params) || length(hu_param) == 0L) {
    stop("Could not locate both PCM and hurdle random effects in ",
         "ranef(model). Found: ",
         paste(re_params, collapse = ", "), call. = FALSE)
  }
  pcm_theta <- as.numeric(person_re[, "Estimate", "Intercept"]) - pcm_shift
  hu_theta  <- -as.numeric(person_re[, "Estimate", hu_param[1]]) - hu_shift

  # Persons who answered 0 on every item gave no severity responses, so
  # their severity EAP comes from the prior and the correlation with the
  # hurdle location only. Leave them out of the severity histogram.
  resp_var <- as.character(stats::formula(model)$formula[[2]])
  if (length(resp_var) > 1) resp_var <- resp_var[2]
  y_scores <- .hpcm_response_scores(model$data[[resp_var]])
  n_active <- tapply(y_scores > 0, as.character(model$data[[person_name]]),
                     sum, na.rm = TRUE)
  active <- n_active[rownames(person_re)] > 0
  n_inactive <- sum(!active)
  pcm_theta <- pcm_theta[active]

  # ---------------------------------------------------------------
  # Panels
  # ---------------------------------------------------------------
  common <- list(
    person_label = "EAP", data_order = data_order, center = center,
    robust = robust, sort_items = sort_items, bins = bins,
    show_ci = show_ci, ci_prob = ci_prob, palette = palette,
    person_fill = person_fill, threshold_fill = threshold_fill,
    panel = panel, xlim = xlim, viridis_option = viridis_option,
    viridis_begin = viridis_begin, viridis_end = viridis_end,
    row_gap = row_gap
  )
  hu_labels <- if (!is.null(category_labels$hurdle)) {
    category_labels$hurdle
  } else {
    c("0", "1+")
  }
  hurdle <- do.call(.targeting_build, c(common, list(
    tau = tau_hu, person_theta = hu_theta[is.finite(hu_theta)],
    n_persons = length(hu_theta), category_labels = hu_labels,
    cat_offset = 0L,
    extra_caption = paste(
      "Hurdle part: the bands show a response of 0 and any response",
      "above 0 (1+).")
  )))
  pcm_extra <- "Severity part: categories above 0."
  if (n_inactive > 0) {
    pcm_extra <- paste0(
      pcm_extra, " The person histogram leaves out ", n_inactive,
      " persons who answered 0 on every item, since they gave no severity ",
      "responses.")
  }
  pcm <- do.call(.targeting_build, c(common, list(
    tau = tau_pcm, person_theta = pcm_theta[is.finite(pcm_theta)],
    n_persons = sum(active), category_labels = category_labels$pcm,
    cat_offset = 1L, extra_caption = pcm_extra
  )))

  if (output == "list") {
    return(list(hurdle = hurdle, pcm = pcm))
  }
  stack <- function(b) {
    b$p1 / b$p2 / b$p3 + patchwork::plot_layout(heights = height_ratios)
  }
  list(hurdle = stack(hurdle), pcm = stack(pcm))
}






