#' Summarize and Plot Posterior Predictive Infit Statistics
#'
#' Postprocesses the output of \code{\link{infit_statistic}} to produce summary tables
#' of posterior predictive infit statistics and a combined slab +
#' interval plot comparing observed infit values to the posterior
#' predictive distribution.
#'
#' @param infit_draws A data frame (or tibble) as returned by
#'   \code{\link{infit_statistic}} containing at minimum the columns
#'   \code{item}, \code{infit} (observed infit per draw), and
#'   \code{infit_rep} (replicated infit per draw).
#' @param ci Numeric in \eqn{(0, 1)}. Width of the credible
#'   interval used for the posterior predictive HDI and the slab
#'   display. Default is 0.84.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{summary}{A \code{\link[tibble]{tibble}} with one row per
#'     item containing: \code{item}, \code{infit_obs} (posterior
#'     mean of observed infit), \code{infit_rep} (posterior mean of
#'     replicated infit), and \code{infit_ppp} (posterior predictive
#'     p-value: proportion of draws where the replicated infit
#'     exceeds the observed infit). Values near 0.5 indicate good
#'     fit; values near 0 suggest underfit; values near 1 suggest
#'     overfit.}
#'   \item{hdi}{A \code{\link[tibble]{tibble}} with one row per
#'     item containing: \code{item}, \code{underfit} (posterior
#'     probability that the observed infit exceeds the upper HDI
#'     bound of the replicated distribution), and \code{overfit}
#'     (posterior probability that the observed infit falls below
#'     the lower HDI bound).}
#'   \item{plot}{A \code{\link[ggplot2]{ggplot}} object showing
#'     the posterior predictive distribution of replicated infit
#'     (grey filled slab) overlaid with the observed infit
#'     distribution (coloured slab + interval), with a dashed
#'     reference line at 1 (perfect fit).}
#' }
#'
#' @details
#' Two complementary summary tables are provided:
#'
#' \describe{
#'   \item{\code{summary}}{Reports the posterior mean observed and
#'     replicated infit values alongside the posterior predictive
#'     p-value (ppp). The ppp is the proportion of draws where the
#'     replicated infit exceeds the observed infit. Under good fit,
#'     the ppp should be near 0.5. A ppp near 0 indicates the
#'     observed infit is consistently larger than expected
#'     (underfit); a ppp near 1 indicates it is consistently
#'     smaller (overfit).}
#'   \item{\code{hdi}}{Reports the probability that the observed
#'     infit falls above (underfit) or below (overfit) the HDI of
#'     the replicated distribution. This provides a more
#'     distributional assessment than the ppp alone.}
#' }
#'
#' The plot uses two layers from the \pkg{ggdist} package:
#' \describe{
#'   \item{\code{stat_slab}}{Displays the posterior predictive
#'     (replicated) infit distribution as a filled density slab
#'     per item, shaded by credible interval level.}
#'   \item{\code{stat_slabinterval}}{Displays the observed infit
#'     distribution per item as a semi-transparent slab with
#'     point and interval summaries.}
#' }
#'
#' Under good model fit, the observed infit distribution should
#' overlap substantially with the replicated distribution. Items
#' where the observed distribution sits systematically above the
#' replicated HDI indicate underfit (more variation than expected);
#' items below indicate overfit (less variation than expected).
#'
#' @seealso
#' \code{\link{infit_statistic}},
#' \code{\link{item_restscore_post}}.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(ggplot2)
#'
#' # Assuming fit_pcm is a fitted brmsfit object
#' infit_draws <- infit_statistic(fit_pcm)
#'
#' result <- infit_post(infit_draws)
#' result$summary
#' result$hdi
#' result$plot
#' }
#'
#' @importFrom dplyr group_by summarise left_join mutate
#' @importFrom ggdist stat_slab stat_slabinterval median_hdi
#' @importFrom ggplot2 ggplot aes after_stat geom_vline labs
#'   theme_bw theme scale_fill_manual scale_fill_brewer
#' @importFrom forcats fct_rev
#' @importFrom rlang .data
#' @export
infit_post <- function(infit_draws, ci = 0.84) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required.", call. = FALSE)
  }
  if (!all(c("item", "infit", "infit_rep") %in% names(infit_draws))) {
    stop("'infit_draws' must contain columns 'item', 'infit', ",
         "and 'infit_rep'.", call. = FALSE)
  }
  
  # --- Summary table: means and ppp ---
  summary_tbl <- infit_draws |>
    dplyr::group_by(.data$item) |>
    dplyr::summarise(
      infit_obs = round(mean(.data$infit), 3),
      infit_rep = round(mean(.data$infit_rep), 3),
      infit_ppp = round(mean(.data$infit_rep > .data$infit), 3),
      .groups = "drop"
    )
  
  # --- HDI table: over/underfit probabilities ---
  infit_rep_hdi <- infit_draws |>
    dplyr::group_by(.data$item) |>
    dplyr::summarise(
      ggdist::median_hdi(.data$infit_rep, .width = ci),
      .groups = "drop"
    )
  
  hdi_tbl <- infit_draws |>
    dplyr::left_join(
      infit_rep_hdi[, c("item", "ymin", "ymax")],
      by = "item"
    ) |>
    dplyr::group_by(.data$item) |>
    dplyr::summarise(
      underfit = round(mean(.data$infit > .data$ymax), 3),
      overfit  = round(mean(.data$infit < .data$ymin), 3),
      .groups = "drop"
    )
  
  # --- Plot ---
  p <- ggplot2::ggplot(
    infit_draws,
    ggplot2::aes(x = .data$infit, y = forcats::fct_rev(.data$item))
  ) +
    ggdist::stat_slab(
      ggplot2::aes(x = .data$infit_rep,
                   fill = ggplot2::after_stat(level)),
      .width = c(ci, 0.95)
    ) +
    ggdist::stat_slabinterval(
      ggplot2::aes(slab_fill = ggplot2::after_stat(level)),
      .width = c(ci, 0.95),
      alpha = 0.8
    ) +
    ggplot2::scale_fill_manual(
      aesthetics = "fill",
      values = c("darkgrey", "lightgrey")
    ) +
    ggplot2::scale_fill_brewer(
      aesthetics = "slab_fill",
      na.value = "lightblue1"
    ) +
    ggplot2::geom_vline(
      xintercept = 1, linetype = 2
    ) +
    ggplot2::labs(
      x = "Observed (blue) and expected conditional item infit",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  list(summary = summary_tbl, hdi = hdi_tbl, plot = p)
}


#' Summarize and Plot Posterior Predictive Item-Restscore Associations
#'
#' Postprocesses the output of \code{\link{item_restscore_statistic}}
#' to produce a summary table and a slab plot comparing observed
#' item-restscore gamma associations to the posterior predictive
#' distribution.
#'
#' @param item_restscore A list as returned by
#'   \code{\link{item_restscore_statistic}}, containing at minimum:
#'   \describe{
#'     \item{result}{A data frame with columns including \code{item}
#'       and summary statistics (first 5 columns are used).}
#'     \item{draws}{A data frame with columns \code{item},
#'       \code{gamma} (observed gamma per draw), and
#'       \code{gamma_rep} (replicated gamma per draw).}
#'   }
#'
#' @return A list with two elements:
#' \describe{
#'   \item{summary}{A \code{\link[tibble]{tibble}} with the first 5
#'     columns of the result table, rounded to 3 decimal places and
#'     sorted by item.}
#'   \item{plot}{A \code{\link[ggplot2]{ggplot}} object showing
#'     the posterior predictive distribution of replicated gamma
#'     (grey filled slab) with the observed gamma values overlaid
#'     as orange diamond points.}
#' }
#'
#' @details
#' The item-restscore gamma association measures the strength of
#' the relationship between each item's responses and the rest
#' score (total score excluding that item). Under good fit, the
#' observed gamma should fall within the posterior predictive
#' distribution.
#'
#' The plot displays:
#' \describe{
#'   \item{Grey slab}{The posterior predictive distribution of
#'     replicated gamma values, shaded by 84\% and 95\% credible
#'     interval levels.}
#'   \item{Orange diamonds}{The observed gamma values per draw,
#'     plotted as points on top of the slab.}
#' }
#'
#' Items where the observed gamma (orange) falls consistently
#' outside the replicated distribution (grey) indicate poor fit
#' in terms of item discrimination.
#'
#' @seealso
#' \code{\link{item_restscore_statistic}},
#' \code{\link{infit_post}}.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(ggplot2)
#'
#' # Assuming fit_pcm is a fitted brmsfit object
#' irs <- item_restscore_statistic(fit_pcm)
#'
#' result <- item_restscore_post(irs)
#' result$summary
#' result$plot
#' }
#'
#' @importFrom dplyr mutate across where arrange
#' @importFrom ggdist stat_slab
#' @importFrom ggplot2 ggplot aes after_stat geom_point labs
#'   theme_bw theme scale_fill_manual scale_fill_brewer
#' @importFrom forcats fct_rev
#' @export
item_restscore_post <- function(item_restscore) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required.", call. = FALSE)
  }
  if (!is.list(item_restscore) ||
      !all(c("result", "draws") %in% names(item_restscore))) {
    stop("'item_restscore' must be a list with elements 'result' ",
         "and 'draws'.", call. = FALSE)
  }
  if (!all(c("item", "gamma", "gamma_rep") %in%
           names(item_restscore$draws))) {
    stop("'item_restscore$draws' must contain columns 'item', ",
         "'gamma', and 'gamma_rep'.", call. = FALSE)
  }
  
  # --- Summary table ---
  summary_tbl <- item_restscore$result[, 1:5] |>
    dplyr::mutate(
      dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))
    ) |>
    dplyr::arrange(.data$item)
  
  # --- Plot ---
  draws_df <- item_restscore$draws
  
  p <- ggplot2::ggplot(
    draws_df,
    ggplot2::aes(x = .data$gamma, y = forcats::fct_rev(.data$item))
  ) +
    ggdist::stat_slab(
      ggplot2::aes(x = .data$gamma_rep,
                   fill = ggplot2::after_stat(level)),
      .width = c(0.84, 0.95)
    ) +
    ggplot2::geom_point(
      size = 6, shape = 18, colour = "sienna1",
      position = ggplot2::position_nudge(y = 0.08)
    ) +
    ggplot2::scale_fill_manual(
      aesthetics = "fill",
      values = c("grey", "lightgrey")
    ) +
    ggplot2::scale_fill_brewer(
      aesthetics = "slab_fill",
      na.value = "lightblue1"
    ) +
    ggplot2::labs(
      x = "Observed (orange) and expected item-restscore gamma",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  list(summary = summary_tbl, plot = p)
}

#' Summarize and Plot Posterior Predictive Q3 Residual Correlations
#'
#' Postprocesses the output of \code{\link{q3_statistic}} to produce
#' summary tables of posterior predictive Q3 statistics and a combined
#' slab + interval plot comparing observed Q3 residual correlations
#' to the posterior predictive distribution for each item pair.
#'
#' @param q3_draws A data frame (or tibble) as returned by
#'   \code{\link{q3_statistic}} containing at minimum the columns
#'   \code{item_pair}, \code{q3} (observed Q3 per draw), and
#'   \code{q3_rep} (replicated Q3 per draw).
#' @param ci Numeric in \eqn{(0, 1)}. Width of the credible
#'   interval used for the posterior predictive HDI and the slab
#'   display. Default is 0.84.
#' @param n_pairs Integer. Maximum number of item pairs to display
#'   in the plot, selected by largest absolute \code{q3_diff}. If
#'   \code{NULL} (the default), all pairs are shown. Useful when
#'   the number of item pairs is large.
#' @param sort_by Character. How to order item pairs on the y-axis.
#'   \code{"q3_diff"} (the default) sorts by the posterior mean
#'   difference between observed and replicated Q3 (largest at top).
#'   \code{"q3_obs"} sorts by the posterior mean observed Q3.
#'   \code{"ppp"} sorts by the posterior predictive p-value.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{summary}{A \code{\link[tibble]{tibble}} with one row per
#'     item pair containing: \code{item_pair}, \code{item_1},
#'     \code{item_2}, \code{q3_obs} (posterior mean observed Q3),
#'     \code{q3_rep} (posterior mean replicated Q3),
#'     \code{q3_diff} (posterior mean difference), and \code{q3_ppp}
#'     (posterior predictive p-value: proportion of draws where the
#'     observed Q3 exceeds the replicated Q3).}
#'   \item{hdi}{A \code{\link[tibble]{tibble}} with one row per
#'     item pair containing: \code{item_pair}, \code{item_1},
#'     \code{item_2}, \code{ld} (local dependence probability:
#'     proportion of draws where observed Q3 exceeds upper HDI
#'     bound of replicated distribution), and \code{lr} (local
#'     repulsion probability: proportion of draws where observed
#'     Q3 falls below lower HDI bound).}
#'   \item{plot}{A \code{\link[ggplot2]{ggplot}} object showing
#'     the posterior predictive distribution of replicated Q3
#'     (grey filled slab) overlaid with the observed Q3
#'     distribution (coloured slab + interval), with a dashed
#'     reference line at 0.}
#' }
#'
#' @details
#' Two complementary summary tables are provided:
#'
#' \describe{
#'   \item{\code{summary}}{Reports posterior mean observed and
#'     replicated Q3 values alongside the posterior predictive
#'     p-value (ppp). The ppp is the proportion of draws where
#'     the observed Q3 exceeds the replicated Q3. Under good fit
#'     (no local dependence), the ppp should be near 0.5. A ppp
#'     close to 1 indicates the observed correlation is
#'     systematically higher than expected (local dependence);
#'     a ppp close to 0 indicates it is systematically lower
#'     (local repulsion, e.g., speed-accuracy tradeoffs).}
#'   \item{\code{hdi}}{Reports the probability that the observed
#'     Q3 falls above (local dependence) or below (local repulsion)
#'     the HDI of the replicated distribution. This provides a
#'     more distributional assessment than the ppp alone.}
#' }
#'
#' The plot uses two layers from the \pkg{ggdist} package:
#' \describe{
#'   \item{\code{stat_slab}}{Displays the posterior predictive
#'     (replicated) Q3 distribution as a filled density slab per
#'     item pair, shaded by credible interval level.}
#'   \item{\code{stat_slabinterval}}{Displays the observed Q3
#'     distribution per item pair as a semi-transparent slab with
#'     point and interval summaries.}
#' }
#'
#' @seealso
#' \code{\link{q3_statistic}},
#' \code{\link{infit_post}},
#' \code{\link{item_restscore_post}}.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(ggplot2)
#'
#' # Assuming fit_pcm is a fitted brmsfit object
#' q3_draws <- q3_statistic(fit_pcm, ndraws_use = 500)
#'
#' result <- q3_post(q3_draws)
#' result$summary
#' result$hdi
#' result$plot
#'
#' # Show only top 10 pairs by Q3 difference
#' result_top <- q3_post(q3_draws, n_pairs = 10)
#' result_top$plot
#' }
#'
#' @importFrom dplyr group_by summarise left_join mutate arrange
#'   desc slice_head filter
#' @importFrom ggdist stat_slab stat_slabinterval median_hdi
#' @importFrom ggplot2 ggplot aes after_stat geom_vline labs
#'   theme_bw theme scale_fill_manual scale_fill_brewer
#' @importFrom forcats fct_rev fct_reorder
#' @importFrom rlang .data
#' @importFrom tibble tibble
#' @export
q3_post <- function(q3_draws,
                    ci = 0.84,
                    n_pairs = NULL,
                    sort_by = c("q3_diff", "q3_obs", "ppp")) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("ggdist", quietly = TRUE)) {
    stop("Package 'ggdist' is required.", call. = FALSE)
  }
  if (!all(c("item_pair", "q3", "q3_rep") %in% names(q3_draws))) {
    stop("'q3_draws' must contain columns 'item_pair', 'q3', ",
         "and 'q3_rep'.", call. = FALSE)
  }
  
  sort_by <- match.arg(sort_by)
  
  # --- Summary table: means and ppp ---
  summary_tbl <- q3_draws |>
    dplyr::group_by(.data$item_pair) |>
    dplyr::summarise(
      item_1  = .data$item_1[1],
      item_2  = .data$item_2[1],
      q3_obs  = round(mean(.data$q3, na.rm = TRUE), 3),
      q3_rep  = round(mean(.data$q3_rep, na.rm = TRUE), 3),
      q3_diff = round(mean(.data$q3 - .data$q3_rep, na.rm = TRUE), 3),
      q3_ppp  = round(mean(.data$q3 > .data$q3_rep, na.rm = TRUE), 3),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$q3_diff))
  
  # --- HDI table: LD / LR probabilities ---
  q3_rep_hdi <- q3_draws |>
    dplyr::group_by(.data$item_pair) |>
    dplyr::summarise(
      ggdist::median_hdi(.data$q3_rep, .width = ci),
      .groups = "drop"
    )
  
  hdi_tbl <- q3_draws |>
    dplyr::left_join(
      q3_rep_hdi[, c("item_pair", "ymin", "ymax")],
      by = "item_pair"
    ) |>
    dplyr::group_by(.data$item_pair) |>
    dplyr::summarise(
      item_1 = .data$item_1[1],
      item_2 = .data$item_2[1],
      ld     = round(mean(.data$q3 > .data$ymax, na.rm = TRUE), 3),
      lr     = round(mean(.data$q3 < .data$ymin, na.rm = TRUE), 3),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$ld))
  
  # --- Determine which pairs to plot ---
  if (sort_by == "q3_diff") {
    sort_var <- abs(summary_tbl$q3_diff)
  } else if (sort_by == "q3_obs") {
    sort_var <- abs(summary_tbl$q3_obs)
  } else {
    sort_var <- summary_tbl$q3_ppp
  }
  pair_order <- summary_tbl$item_pair[order(sort_var, decreasing = TRUE)]
  
  if (!is.null(n_pairs)) {
    n_pairs <- min(as.integer(n_pairs), length(pair_order))
    pairs_to_plot <- pair_order[seq_len(n_pairs)]
  } else {
    pairs_to_plot <- pair_order
  }
  
  plot_data <- q3_draws |>
    dplyr::filter(.data$item_pair %in% pairs_to_plot) |>
    dplyr::mutate(
      item_pair = factor(.data$item_pair, levels = rev(pairs_to_plot))
    )
  
  # --- Plot ---
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$q3, y = .data$item_pair)
  ) +
    ggdist::stat_slab(
      ggplot2::aes(x = .data$q3_rep,
                   fill = ggplot2::after_stat(level)),
      .width = c(ci, 0.95)
    ) +
    ggdist::stat_slabinterval(
      ggplot2::aes(slab_fill = ggplot2::after_stat(level)),
      .width = c(ci, 0.95),
      alpha = 0.8
    ) +
    ggplot2::scale_fill_manual(
      aesthetics = "fill",
      values = c("darkgrey", "lightgrey")
    ) +
    ggplot2::scale_fill_brewer(
      aesthetics = "slab_fill",
      na.value = "lightblue1"
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = 2
    ) +
    ggplot2::labs(
      x = "Observed and expected Q3 residual correlation",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  list(summary = summary_tbl, hdi = hdi_tbl, plot = p)
}

#' Extract Informative Priors from a Fitted Bayesian IRT Model
#'
#' Takes a fitted \code{\link[brms]{brmsfit}} object and constructs a
#' \code{\link[brms]{brmsprior}} object in which each item parameter
#' receives a \code{normal(mean, sd)} prior derived from its posterior
#' distribution. The person-level random effect SD prior is also
#' updated. The returned prior can be passed directly to
#' \code{\link[stats]{update}} (or \code{\link[brms]{brm}}) to refit
#' the model with empirical Bayes / informative priors — useful for
#' anchoring scales, warm-starting a model on new data, or
#' regularising estimation with small samples.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object. Supported
#'   parameterisations:
#'   \describe{
#'     \item{Polytomous ordinal}{e.g., \code{family = acat} with
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
#' @param mult Numeric multiplier applied to each posterior SD before
#'   it is used as the prior SD. Values > 1 widen the priors
#'   (less informative); values < 1 tighten them. Default is 1
#'   (use posterior SD directly).
#'
#' @return A \code{\link[brms]{brmsprior}} object that can be supplied
#'   to the \code{prior} argument of \code{\link[brms]{brm}} or
#'   \code{\link[stats]{update}}.
#'
#' @details
#' The function extracts all posterior draws via
#' \code{\link[brms]{as_draws_df}}, computes the mean and SD of each
#' parameter's marginal posterior, and constructs
#' \code{normal(mean, sd * mult)} priors.
#'
#' \strong{Polytomous ordinal models} with grouped thresholds
#' (\code{thres(gr = item)}): each threshold receives its own prior
#' via \code{brms::set_prior("normal(...)", class = "Intercept",
#' group = item, coef = threshold_index)}.
#'
#' \strong{Dichotomous Rasch models} parameterised as
#' \code{response ~ 1 + (1 | item) + (1 | id)}: priors are set on
#' the global intercept (\code{class = "Intercept"}), the item-level
#' SD (\code{class = "sd", group = item_var}), and the person-level
#' SD.
#'
#' \strong{Dichotomous 1PL models} parameterised as
#' \code{response ~ 0 + item + (1 | id)}: each item-specific fixed
#' effect (e.g., \code{b_itemI1}) receives its own
#' \code{normal(mean, sd)} prior via
#' \code{brms::set_prior(..., class = "b", coef = "itemI1")}.
#'
#' In all cases the person-level SD receives a
#' \code{normal(mean, sd * mult)} prior (brms applies the lower
#' bound of zero automatically for SD parameters).
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
#'   chains = 4, cores = 2, iter = 1000 # use more iter (and cores if you have)
#' )
#'
#' # Extract posterior-informed priors
#' new_priors <- posterior_to_prior(fit_pcm)
#' new_priors
#'
#' # Narrow the prior's sd by a factor of 0.5
#' wide_priors <- posterior_to_prior(fit_pcm, mult = 0.5)
#'
#' # --- Dichotomous 1PL (fixed item effects) ---
#'
#' df_rm <- eRm::raschdat3 %>%
#'   as.data.frame() %>%
#'   rownames_to_column("id") %>%
#'   pivot_longer(!id, names_to = "item", values_to = "response")
#'
#' fit_1pl <- brm(
#'   response ~ 0 + item + (1 | id),
#'   data   = df_rm,
#'   family = bernoulli(),
#'   chains = 4, cores = 2, iter = 1000 # use more iter (and cores if you have)
#' )
#'
#' priors_1pl <- posterior_to_prior(fit_1pl)
#' priors_1pl
#' }
#'
#' @importFrom brms as_draws_df set_prior
#' @importFrom rlang enquo as_name
#' @importFrom stats sd family
#' @export
posterior_to_prior <- function(model,
                               item_var = item,
                               person_var = id,
                               mult = 1) {
  
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  if (!is.numeric(mult) || length(mult) != 1 || mult <= 0) {
    stop("'mult' must be a single positive number.", call. = FALSE)
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
  
  draws       <- brms::as_draws_df(model)
  family_name <- stats::family(model)$family
  is_ordinal  <- grepl("acat|cumul|sratio|cratio",
                       family_name, ignore.case = TRUE)
  
  prior_list <- list()
  
  # ==================================================================
  # ORDINAL MODELS (grouped thresholds: thres(gr = item))
  # ==================================================================
  if (is_ordinal) {
    has_grouped <- any(grepl("^b_Intercept\\[.+,\\d+\\]$", names(draws)))
    
    if (!has_grouped) {
      stop("Ordinal model does not appear to use grouped thresholds ",
           "(thres(gr = ...)). This function currently requires ",
           "item-specific thresholds.", call. = FALSE)
    }
    
    unique_items <- sort(unique(model$data[[item_name]]))
    
    for (item_label in unique_items) {
      thresh_pattern <- paste0(
        "^b_Intercept\\[",
        gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
        ","
      )
      thresh_cols <- grep(thresh_pattern, names(draws), value = TRUE)
      
      if (length(thresh_cols) == 0) {
        warning("No threshold parameters found for item '", item_label,
                "'. Skipping.", call. = FALSE)
        next
      }
      
      # Sort numerically
      thresh_nums <- as.numeric(
        gsub(".*,(\\d+)\\]$", "\\1", thresh_cols)
      )
      thresh_cols <- thresh_cols[order(thresh_nums)]
      
      for (idx in seq_along(thresh_cols)) {
        vals   <- draws[[thresh_cols[idx]]]
        p_mean <- mean(vals)
        p_sd   <- stats::sd(vals) * mult
        
        prior_list[[length(prior_list) + 1]] <- brms::set_prior(
          paste0("normal(", round(p_mean, 4), ", ", round(p_sd, 4), ")"),
          class = "Intercept",
          coef  = as.character(idx),
          group = item_label
        )
      }
    }
    
    # ==================================================================
    # DICHOTOMOUS MODELS
    # ==================================================================
  } else {
    
    # --- Detect parameterisation ---
    # Fixed item effects: columns like b_itemI1, b_itemI2, ...
    # Random item effects: b_Intercept + r_item[I1,Intercept], ...
    fe_pattern <- paste0("^b_", item_name)
    fe_cols    <- grep(fe_pattern, names(draws), value = TRUE)
    # Exclude b_Intercept from the match
    fe_cols    <- setdiff(fe_cols, "b_Intercept")
    
    has_fixed_items  <- length(fe_cols) > 0
    has_random_items <- any(grepl(
      paste0("^r_", item_name, "\\["), names(draws)
    ))
    
    # ── Fixed item effects: response ~ 0 + item + (1 | id) ──
    if (has_fixed_items) {
      
      for (col_name in fe_cols) {
        vals   <- draws[[col_name]]
        p_mean <- mean(vals)
        p_sd   <- stats::sd(vals) * mult
        
        # Column name is e.g. "b_itemI1" -> coef is "itemI1"
        coef_name <- sub("^b_", "", col_name)
        
        prior_list[[length(prior_list) + 1]] <- brms::set_prior(
          paste0("normal(", round(p_mean, 4), ", ", round(p_sd, 4), ")"),
          class = "b",
          coef  = coef_name
        )
      }
      
      # ── Random item effects: response ~ 1 + (1 | item) + (1 | id) ──
    } else if (has_random_items) {
      
      # Global intercept
      intercept_col <- grep("^b_Intercept$", names(draws), value = TRUE)
      if (length(intercept_col) == 1) {
        vals   <- draws[[intercept_col]]
        p_mean <- mean(vals)
        p_sd   <- stats::sd(vals) * mult
        
        prior_list[[length(prior_list) + 1]] <- brms::set_prior(
          paste0("normal(", round(p_mean, 4), ", ", round(p_sd, 4), ")"),
          class = "Intercept"
        )
      }
      
      # Item-level SD
      item_sd_col <- grep(
        paste0("^sd_", item_name, "__Intercept$"),
        names(draws), value = TRUE
      )
      if (length(item_sd_col) == 1) {
        vals   <- draws[[item_sd_col]]
        p_mean <- mean(vals)
        p_sd   <- stats::sd(vals) * mult
        
        prior_list[[length(prior_list) + 1]] <- brms::set_prior(
          paste0("normal(", round(p_mean, 4), ", ", round(p_sd, 4), ")"),
          class = "sd",
          group = item_name
        )
      }
      
    } else {
      stop("Could not detect item parameters in the model. ",
           "Supported dichotomous parameterisations:\n",
           "  - response ~ 0 + item + (1 | id)    [fixed items]\n",
           "  - response ~ 1 + (1 | item) + (1 | id)  [random items]",
           call. = FALSE)
    }
  }
  
  # ==================================================================
  # PERSON-LEVEL SD (common to all parameterisations)
  # ==================================================================
  person_sd_col <- grep(
    paste0("^sd_", person_name, "__Intercept$"),
    names(draws), value = TRUE
  )
  if (length(person_sd_col) == 1) {
    vals   <- draws[[person_sd_col]]
    p_mean <- mean(vals)
    p_sd   <- stats::sd(vals) * mult
    
    prior_list[[length(prior_list) + 1]] <- brms::set_prior(
      paste0("normal(", round(p_mean, 4), ", ", round(p_sd, 4), ")"),
      class = "sd",
      group = person_name
    )
  }
  
  # ==================================================================
  # COMBINE INTO A SINGLE brmsprior OBJECT
  # ==================================================================
  if (length(prior_list) == 0) {
    stop("Could not extract any parameters from the model.",
         call. = FALSE)
  }
  
  combined_prior <- prior_list[[1]]
  if (length(prior_list) > 1) {
    for (i in 2:length(prior_list)) {
      combined_prior <- combined_prior + prior_list[[i]]
    }
  }
  
  combined_prior
}