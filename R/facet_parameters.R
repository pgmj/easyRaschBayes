#' Facet Parameters from a Many-Facet Rasch Model
#'
#' Extracts the measures of a third facet, for example raters, from a
#' many-facet Rasch model fitted with \pkg{brms}, with posterior summaries,
#' the number of ratings, the observed average rating and a fair average
#' rating for each element.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from a
#'   many-facet model with the facet as a random intercept, a fixed effect,
#'   or a threshold group. Examples:
#'   \describe{
#'     \item{Rating scale, random raters}{\code{rating ~ 1 + item +
#'       (1 | rater) + (1 | id)}, \code{family = acat}.}
#'     \item{Partial credit by item}{\code{rating | thres(gr = item) ~ 1 +
#'       (1 | rater) + (1 | id)}.}
#'     \item{Rater-specific rating scales}{\code{rating | thres(gr = rater)
#'       ~ 1 + item + (1 | id)}.}
#'     \item{Dichotomous}{\code{y ~ 1 + item + (1 | rater) + (1 | id)},
#'       \code{family = bernoulli()}.}
#'   }
#' @param facet_var An unquoted variable name identifying the facet, for
#'   example \code{rater}.
#' @param item_var An unquoted variable name identifying the items.
#'   Default is \code{item}.
#' @param person_var An unquoted variable name identifying the persons.
#'   Default is \code{id}.
#' @param draws Logical. If \code{TRUE}, the centred posterior draws of
#'   the facet locations are included as a matrix (elements x draws).
#'   Default is \code{FALSE}.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the highest density
#'   continuous interval (HDCI). Default is 0.95.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{summary}{A \code{\link[tibble]{tibble}} with one row per facet
#'     element and the columns: the facet variable, \code{n} (number of
#'     ratings), \code{observed_avg} (mean observed rating),
#'     \code{fair_avg}, \code{fair_avg_lower} and \code{fair_avg_upper}
#'     (posterior mean and HDCI of the fair average, see Details),
#'     \code{location} (posterior mean location, positive = lower
#'     ratings), \code{se} (posterior SD), \code{hdci_lower},
#'     \code{hdci_upper} and \code{n_eff} (bulk effective sample size).}
#'   \item{thresholds}{(Only with \code{thres(gr = facet)}) A
#'     \code{\link[tibble]{tibble}} with one row per element and
#'     threshold: the facet variable, \code{threshold}, \code{location},
#'     \code{se}, \code{hdci_lower} and \code{hdci_upper}, on the same
#'     scale as \code{summary}. \code{NULL} otherwise.}
#'   \item{facet_sd}{(Only for a random facet) A
#'     \code{\link[tibble]{tibble}} with the posterior mean, SD and HDCI of
#'     the facet standard deviation. \code{NULL} otherwise.}
#'   \item{draws_matrix}{(Only if \code{draws = TRUE}) A numeric matrix of
#'     centred location draws with one row per element and one column per
#'     posterior draw.}
#' }
#'
#' @details
#' \strong{Locations.} The model is read as
#' \eqn{\theta_v - \delta_i - \lambda_j} on the logit scale, where
#' \eqn{\delta_i} is the item location and \eqn{\lambda_j} the location
#' of facet element \eqn{j}. A positive location means lower ratings, so
#' for raters it is severity, on the same side as item difficulty in
#' \code{\link{item_parameters}}. The locations are centred within each
#' posterior draw so that they average 0. The result is the same whether
#' the facet is a random effect or a fixed effect with any contrasts,
#' including \code{\link{contr_sum_zero}}. With facet-specific thresholds
#' (\code{thres(gr = rater)}) the location is the mean of the element's
#' thresholds, and the \code{thresholds} table shows each threshold.
#' brms sets the number of thresholds per element from the categories that
#' element used, so set it in the formula, as in
#' \code{thres(3, gr = rater)} for four categories. The function stops if
#' the elements have different numbers of thresholds. Widely spaced thresholds for a rater mean the middle categories cover a
#' wide range of locations and are used often (central tendency). Closely
#' spaced thresholds mean the extreme categories are used often.
#'
#' \strong{Fair average.} The observed average rating of a rater depends on
#' which persons and items the rater happened to rate. The fair average is
#' the expected rating the element would give a person at the mean
#' location of the persons in the sample, averaged over the items, computed
#' in each posterior draw. It corresponds to the "Fair-M average" of the
#' Facets program. Person x facet (halo) effects are set to 0. Ratings are
#' on the scale brms uses: 0 and 1 for \code{bernoulli()}, categories
#' 1 to \eqn{C} for ordinal families.
#'
#' \strong{Reliability.} With \code{draws = TRUE}, the draws matrix can be
#' passed to \code{\link{RMUreliability}}. For raters, a high value means
#' that the raters differ reliably in severity, which is a problem for the
#' ratings, not a sign of their quality. With few elements, for example 3
#' to 5 raters, the correlation across elements that RMU is based on is
#' very uncertain.
#'
#' @references
#' Linacre, J. M. (1989). \emph{Many-facet Rasch measurement}. MESA Press.
#'
#' Eckes, T. (2015). \emph{Introduction to many-facet Rasch measurement:
#' Analyzing and evaluating rater-mediated assessments} (2nd ed.). Peter
#' Lang.
#'
#' @seealso
#' \code{\link{item_parameters}}, \code{\link{contr_sum_zero}},
#' \code{\link{infit_statistic}} (with \code{facet_var}),
#' \code{\link{RMUreliability}}.
#'
#' @examples
#' \donttest{
#' library(brms)
#'
#' # Simulated ratings: 60 persons, 3 items, 4 raters, two raters per person
#' set.seed(1)
#' d <- do.call(rbind, lapply(1:60, function(p) {
#'   expand.grid(id = p, item = c("i1", "i2", "i3"),
#'               rater = paste0("R", c((p - 1) %% 4 + 1, p %% 4 + 1)),
#'               stringsAsFactors = FALSE)
#' }))
#' severity <- c(R1 = -0.5, R2 = 0, R3 = 0.2, R4 = 0.6)
#' theta <- rnorm(60)
#' d$y <- rbinom(nrow(d), 1, plogis(theta[d$id] - severity[d$rater]))
#'
#' fit <- brm(y ~ 1 + item + (1 | rater) + (1 | id), data = d,
#'            family = bernoulli(), chains = 2, iter = 1000)
#'
#' fp <- facet_parameters(fit, facet_var = rater, draws = TRUE)
#' fp$summary
#' RMUreliability(fp$draws_matrix)
#' }
#'
#' @importFrom rlang enquo as_name
#' @export
facet_parameters <- function(model, facet_var, item_var = item,
                             person_var = id, draws = FALSE, prob = 0.95) {
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  if (!is.numeric(prob) || length(prob) != 1L || prob <= 0 || prob >= 1) {
    stop("'prob' must be a single number between 0 and 1.", call. = FALSE)
  }
  facet_name  <- rlang::as_name(rlang::enquo(facet_var))
  item_name   <- rlang::as_name(rlang::enquo(item_var))
  person_name <- rlang::as_name(rlang::enquo(person_var))

  all_draws <- as.data.frame(brms::as_draws_df(model))
  .facet_parameters_core(model, all_draws, facet_name, item_name,
                         person_name, draws, prob)
}

# Everything after the draws are extracted, separated so that tests can
# supply draws for a model shell.
.facet_parameters_core <- function(model, all_draws, facet_name, item_name,
                                   person_name, draws, prob) {
  fr <- .resolve_facet(model, person_name, item_name, facet_name,
                       "facet_parameters")
  if (!any(unlist(fr$roles))) {
    stop("`facet_parameters()`: '", facet_name, "' is not in the model as ",
         "a fixed effect, a random intercept or a threshold group.",
         call. = FALSE)
  }

  fam <- stats::family(model)
  # Exact names: grepl() would also match the hurdle_acat custom family,
  # whose zero-hurdle the expected rating below would ignore.
  ordinal <- fam$family %in% c("acat", "cumulative", "sratio", "cratio")
  if (!ordinal && fam$family != "bernoulli") {
    stop("`facet_parameters()` supports the ordinal families acat(), ",
         "cumulative(), sratio() and cratio(), and bernoulli(). The model ",
         "uses '", fam$family, "' (hurdle and other custom families are ",
         "not supported).", call. = FALSE)
  }
  linkinv <- if (ordinal) {
    .ordinal_linkinv(model, "facet_parameters")
  } else {
    if (!identical(fam$link, "logit")) {
      stop("`facet_parameters()` supports bernoulli() with the logit link.",
           call. = FALSE)
    }
    stats::plogis
  }

  pr <- .facet_draws(all_draws, model, facet_name, item_name, person_name,
                     fr, ordinal)
  facet_lev <- pr$levels
  n_el <- length(facet_lev)

  # --- Fair average per element and draw ---
  fair <- .fair_average_draws(pr, fam$family, linkinv)

  # --- Observed ratings ---
  resp_var <- as.character(stats::formula(model)$formula[[2]])
  if (length(resp_var) > 1) resp_var <- resp_var[2]
  y <- as.numeric(model$data[[resp_var]])
  f <- factor(as.character(model$data[[facet_name]]), levels = facet_lev)
  n_obs <- as.integer(table(f))
  obs_avg <- as.numeric(tapply(y, f, mean))

  # --- Summaries ---
  chain <- all_draws$.chain
  summ <- lapply(seq_len(n_el), function(j) {
    loc <- pr$location[, j]
    h <- ggdist::mean_hdci(loc, .width = prob)
    fa <- ggdist::mean_hdci(fair[, j], .width = prob)
    tibble::tibble(
      element        = facet_lev[j],
      n              = n_obs[j],
      observed_avg   = round(obs_avg[j], 3),
      fair_avg       = round(fa$y, 3),
      fair_avg_lower = round(fa$ymin, 3),
      fair_avg_upper = round(fa$ymax, 3),
      location       = round(h$y, 4),
      se             = round(stats::sd(loc), 4),
      hdci_lower     = round(h$ymin, 4),
      hdci_upper     = round(h$ymax, 4),
      n_eff          = round(.compute_ess(loc, chain))
    )
  })
  summary_tbl <- dplyr::bind_rows(summ)
  names(summary_tbl)[1] <- facet_name

  thr_tbl <- NULL
  if (!is.null(pr$thresholds)) {
    rows <- list()
    for (j in seq_len(n_el)) {
      m <- pr$thresholds[[j]]
      for (k in seq_len(ncol(m))) {
        h <- ggdist::mean_hdci(m[, k], .width = prob)
        rows[[length(rows) + 1L]] <- tibble::tibble(
          element = facet_lev[j], threshold = k,
          location = round(h$y, 4), se = round(stats::sd(m[, k]), 4),
          hdci_lower = round(h$ymin, 4), hdci_upper = round(h$ymax, 4)
        )
      }
    }
    thr_tbl <- dplyr::bind_rows(rows)
    names(thr_tbl)[1] <- facet_name
  }

  facet_sd <- NULL
  sd_col <- paste0("sd_", facet_name, "__Intercept")
  if (sd_col %in% names(all_draws)) {
    v <- all_draws[[sd_col]]
    h <- ggdist::mean_hdci(v, .width = prob)
    facet_sd <- tibble::tibble(mean = round(h$y, 4),
                               sd = round(stats::sd(v), 4),
                               hdci_lower = round(h$ymin, 4),
                               hdci_upper = round(h$ymax, 4))
  }

  out <- list(summary = summary_tbl, thresholds = thr_tbl,
              facet_sd = facet_sd)
  if (isTRUE(draws)) {
    dm <- t(pr$location)
    rownames(dm) <- facet_lev
    out$draws_matrix <- dm
  }
  out
}

# Facet locations and the threshold sets needed for the fair average.
#
# Returns list(levels, location = draws x elements (centred, positive =
# lower ratings), thresholds = list per element of draws x K centred
# thresholds or NULL, tau = list over items of lists over elements of
# draws x K thresholds on the logit scale (uncentred), theta = draws-long
# vector of mean person location).
.facet_draws <- function(all_draws, model, facet_name, item_name,
                         person_name, fr, ordinal) {
  facet_lev <- .item_levels(model, facet_name)
  item_lev <- .item_levels(model, item_name)
  n_draws <- nrow(all_draws)
  n_el <- length(facet_lev)
  vars <- names(model$data)

  thresholds <- NULL
  if (isTRUE(fr$roles$thres)) {
    # Facet-specific thresholds; item effects shift them
    tau_f <- .item_tau_draws(all_draws, facet_lev, facet_name, ordinal = TRUE,
                             other_vars = vars)
    # brms sets the number of thresholds per group from the categories
    # used, so a rater who never used the top category gets fewer. Means
    # of different numbers of thresholds are not comparable locations.
    n_thr <- vapply(tau_f, ncol, integer(1))
    if (length(unique(n_thr)) > 1L) {
      stop("`facet_parameters()`: the levels of '", facet_name, "' have ",
           "different numbers of thresholds (",
           paste0(names(tau_f), ": ", n_thr, collapse = ", "), "), ",
           "since brms counts the categories each level used. Their ",
           "locations cannot be compared. Refit with the number of ",
           "thresholds set, thres(", max(n_thr), ", gr = ", facet_name,
           "), so that every level has all of them.", call. = FALSE)
    }
    item_eff <- .item_effect_draws(all_draws, item_lev, item_name, vars)
    if (is.null(item_eff)) {
      # Zero item effects only when the model has none
      item_roles <- .formula_roles(model, item_name)
      if (item_roles$fixed || item_roles$random) {
        stop("`facet_parameters()`: could not find the parameters of '",
             item_name, "' in the model, although it is in the formula. ",
             "Please report this with your model formula.", call. = FALSE)
      }
      item_eff <- matrix(0, n_draws, length(item_lev))
    }
    raw_loc <- vapply(tau_f, rowMeans, numeric(n_draws))
    shift <- rowMeans(raw_loc)
    thresholds <- lapply(tau_f, function(m) unname(m - shift))
    tau <- lapply(seq_along(item_lev), function(i) {
      lapply(tau_f, function(m) m - item_eff[, i])
    })
  } else {
    eff <- .item_effect_draws(all_draws, facet_lev, facet_name, vars)
    if (is.null(eff)) {
      stop("`facet_parameters()`: no parameters found for '", facet_name,
           "'.", call. = FALSE)
    }
    raw_loc <- -eff
    tau_i <- .item_tau_draws(all_draws, item_lev, item_name, ordinal,
                             other_vars = vars)
    tau <- lapply(tau_i, function(t) {
      lapply(seq_len(n_el), function(j) t - eff[, j])
    })
  }
  location <- raw_loc - rowMeans(raw_loc)
  colnames(location) <- as.character(facet_lev)

  p_cols <- grep(paste0("^r_", person_name, "\\[.+,Intercept\\]$"),
                 names(all_draws), value = TRUE)
  theta <- if (length(p_cols)) {
    rowMeans(as.matrix(all_draws[, p_cols, drop = FALSE]))
  } else {
    rep(0, n_draws)
  }

  list(levels = as.character(facet_lev), location = location,
       thresholds = thresholds, tau = tau, theta = theta)
}

# Fair average (draws x elements): expected rating at the mean person
# location, averaged over items. Ratings 0/1 for bernoulli, 1..C for
# ordinal families.
.fair_average_draws <- function(pr, family_name, linkinv) {
  n_draws <- length(pr$theta)
  n_el <- length(pr$levels)
  n_items <- length(pr$tau)
  out <- matrix(0, n_draws, n_el)
  for (i in seq_len(n_items)) {
    for (j in seq_len(n_el)) {
      t_ij <- pr$tau[[i]][[j]]
      out[, j] <- out[, j] + .expected_rating(pr$theta, t_ij, family_name,
                                              linkinv)
    }
  }
  out / n_items
}

# Expected rating for each draw given location theta (vector) and
# thresholds tau (draws x K).
.expected_rating <- function(theta, tau, family_name, linkinv) {
  tau <- as.matrix(tau)
  if (family_name == "bernoulli") {
    return(linkinv(theta - tau[, 1L]))
  }
  K <- ncol(tau)
  if (grepl("acat", family_name)) {
    # log P(x) - log P(1) = sum_{k < x} (theta - tau_k)
    lp <- matrix(0, nrow(tau), K + 1L)
    for (k in seq_len(K)) lp[, k + 1L] <- lp[, k] + theta - tau[, k]
    lp <- lp - apply(lp, 1L, max)
    p <- exp(lp)
    p <- p / rowSums(p)
    return(as.vector(p %*% seq_len(K + 1L)))
  }
  # Same probabilities as .compute_cat_probs(), for all draws at once
  theta <- rep_len(theta, nrow(tau))
  if (grepl("cumul", family_name)) {
    cum <- linkinv(tau - theta)
    p <- cbind(cum, 1) - cbind(0, cum)
    p <- pmax(p, 0)
    p <- p / rowSums(p)
  } else {
    sratio <- grepl("sratio", family_name)
    if (!sratio && !grepl("cratio", family_name)) {
      stop("Unsupported family '", family_name, "'.", call. = FALSE)
    }
    p <- matrix(0, nrow(tau), K + 1L)
    remaining <- rep(1, nrow(tau))
    for (k in seq_len(K)) {
      # Probabilities of stopping at category k and of going on, given
      # that k is reached
      if (sratio) {
        stop_k <- linkinv(tau[, k] - theta)
        go_k <- 1 - stop_k
      } else {
        go_k <- linkinv(theta - tau[, k])
        stop_k <- 1 - go_k
      }
      p[, k] <- remaining * stop_k
      remaining <- remaining * go_k
    }
    p[, K + 1L] <- remaining
  }
  as.vector(p %*% seq_len(K + 1L))
}
