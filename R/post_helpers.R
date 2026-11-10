#' Summarize and Plot Posterior Predictive Infit Statistics
#'
#' Postprocesses the output of \code{\link{infit_statistic}} to produce summary tables
#' of posterior predictive infit statistics and a combined slab +
#' interval plot comparing observed infit values to the posterior
#' predictive distribution.
#'
#' @param infit_draws A data frame (or tibble) as returned by
#'   \code{\link{infit_statistic}} or one element of the output of
#'   \code{\link{infit_statistic_hpcm}}, containing the columns
#'   \code{infit} (observed infit per draw) and \code{infit_rep}
#'   (replicated infit per draw). The item identifier is taken from a
#'   column named \code{item} if there is one, and otherwise from the
#'   first column.
#' @param ci Numeric in \eqn{(0, 1)}. Width of the highest density
#'   continuous interval of the replicated (posterior predictive)
#'   infit distribution, used in the \code{hdi} table and the slab
#'   display. Default is 0.84.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{summary}{A \code{\link[tibble]{tibble}} with one row per
#'     item containing: \code{item}, \code{infit_obs} (posterior
#'     mean of observed infit), \code{infit_rep} (mean of
#'     replicated infit across draws), and \code{infit_ppp} (posterior
#'     predictive p-value: proportion of draws where the observed infit
#'     exceeds the replicated infit, counting ties as one half). Values
#'     near 0.5 indicate good fit, values near 1 suggest underfit and
#'     values near 0 suggest overfit.}
#'   \item{hdi}{A \code{\link[tibble]{tibble}} with one row per
#'     item containing: \code{item}, \code{underfit} (proportion
#'     of draws where the observed infit exceeds the upper HDI
#'     bound of the replicated distribution), and \code{overfit}
#'     (proportion of draws where the observed infit falls below
#'     the lower HDI bound). See Details for how to read these.}
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
#'     observed infit exceeds the replicated infit, with ties counted
#'     as one half. Under good fit, the ppp should be near 0.5. A ppp
#'     near 1 indicates the observed infit is consistently larger
#'     than expected (underfit), and a ppp near 0 that it is
#'     consistently smaller (overfit). The direction is the same as
#'     for the ppp of \code{\link{q3_post}} and
#'     \code{\link{item_restscore_statistic}}.}
#'   \item{\code{hdi}}{Reports the proportion of draws where the
#'     observed infit falls above (underfit) or below (overfit) the
#'     HDI of the replicated distribution. The HDI is computed over
#'     all draws, so each draw's observed value is compared to the
#'     same fixed interval. The observed infit varies across draws
#'     only through parameter uncertainty, while the replicated
#'     infit also varies through the sampling of new responses. The
#'     observed distribution is therefore narrower, and these
#'     proportions move quickly from near 0 to near 1 as the observed
#'     values pass an interval bound. They are exceedance
#'     proportions, not posterior probabilities of misfit.}
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
#' The infit ppp is conservative. Values near 0 or 1 are rarer than for a
#' uniform p-value, and misfit can go undetected (see
#' \code{\link{infit_statistic}}).
#' It is suggested to also evaluate item fit with
#' \code{\link{item_restscore_statistic}} and
#' \code{\link{item_restscore_post}}.
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
#' @importFrom ggdist stat_slab stat_slabinterval median_hdci
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
  if (!all(c("infit", "infit_rep") %in% names(infit_draws))) {
    stop("'infit_draws' must contain columns 'infit' and 'infit_rep'.",
         call. = FALSE)
  }
  item_col <- if ("item" %in% names(infit_draws)) "item" else names(infit_draws)[1]
  infit_draws <- dplyr::ungroup(infit_draws)
  infit_draws$item <- infit_draws[[item_col]]
  
  # --- Summary table: means and ppp ---
  summary_tbl <- infit_draws |>
    dplyr::group_by(.data$item) |>
    dplyr::summarise(
      # The ppp must be computed before infit_rep is overwritten by its
      # mean, since later expressions in summarise() see the new column.
      infit_ppp = round(mean(.data$infit > .data$infit_rep) +
                          0.5 * mean(.data$infit == .data$infit_rep), 3),
      infit_obs = round(mean(.data$infit), 3),
      infit_rep = round(mean(.data$infit_rep), 3),
      .groups = "drop"
    ) |>
    dplyr::select("item", "infit_obs", "infit_rep", "infit_ppp")
  
  # --- HDI table: over/underfit probabilities ---
  infit_rep_hdi <- infit_draws |>
    dplyr::group_by(.data$item) |>
    dplyr::summarise(
      # Use the continuous-interval variant (`_hdci`) rather than
      # `_hdi`. The plain `median_hdi()` returns a *set* of disjoint
      # HDIs when the posterior is multimodal at the requested width,
      # i.e., one row per disjoint segment. Inside dplyr::summarise()
      # that violates the one-row-per-group contract and errors out
      # with "must be size 1, not 2". The continuous version always
      # returns a single [lower, upper] interval, which is the right
      # quantity to use as a cutoff for the LD / LR classification.
      ggdist::median_hdci(.data$infit_rep, .width = ci),
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
      x = "Observed (blue) and replicated item infit",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  names(summary_tbl)[names(summary_tbl) == "item"] <- item_col
  names(hdi_tbl)[names(hdi_tbl) == "item"] <- item_col
  list(summary = summary_tbl, hdi = hdi_tbl, plot = p)
}


#' Summarize and Plot Item-Restscore Associations
#'
#' Postprocesses the output of \code{\link{item_restscore_statistic}}
#' to produce a summary table and a slab plot comparing observed
#' item-restscore gamma associations to the distribution of gamma in
#' datasets replicated conditional on the total scores.
#'
#' @param item_restscore A list as returned by
#'   \code{\link{item_restscore_statistic}}, with elements
#'   \code{result} and \code{draws}. The first column of both is the
#'   item identifier.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{summary}{A \code{\link[tibble]{tibble}} with one row per
#'     item and the columns item, \code{gamma_obs}, \code{gamma_exp},
#'     \code{gamma_diff}, \code{ppp} and \code{p_adj}, rounded to 3
#'     decimals and sorted by item.}
#'   \item{plot}{A \code{\link[ggplot2]{ggplot}} object showing the
#'     distribution of replicated gamma (grey slab) with the observed
#'     gamma overlaid as an orange diamond.}
#' }
#'
#' @details
#' The item-restscore gamma association measures the strength of
#' the relationship between each item's responses and the rest
#' score (total score excluding that item). Under good fit, the
#' observed gamma should fall within the distribution of replicated
#' gamma.
#'
#' The plot displays:
#' \describe{
#'   \item{Grey slab}{The distribution of replicated gamma values,
#'     shaded at the 84% and 95% interval levels.}
#'   \item{Orange diamond}{The observed gamma.}
#' }
#'
#' Items where the observed gamma (orange) falls outside the
#' replicated distribution (grey) indicate poor fit in terms of
#' item discrimination. See \code{\link{item_restscore_statistic}} for
#' the definitions of the summary columns.
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
#' @importFrom ggdist stat_slab
#' @importFrom ggplot2 ggplot aes after_stat geom_point labs
#'   theme_bw theme scale_fill_manual
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
  summary_cols <- c("gamma_obs", "gamma_exp", "gamma_diff", "ppp", "p_adj")
  if (!all(summary_cols %in% names(item_restscore$result)) ||
      !"gamma_rep" %in% names(item_restscore$draws)) {
    stop("'item_restscore' must be the output of ",
         "item_restscore_statistic() from easyRaschBayes 0.3.0.9000 or ",
         "later.", call. = FALSE)
  }

  item_col <- names(item_restscore$result)[1]

  # --- Summary table ---
  summary_tbl <- item_restscore$result[, c(item_col, summary_cols)]
  summary_tbl[summary_cols] <- lapply(summary_tbl[summary_cols], round, 3)
  summary_tbl <- summary_tbl[order(summary_tbl[[item_col]]), ]

  # --- Plot ---
  draws_df <- data.frame(
    item      = as.character(item_restscore$draws[[1]]),
    gamma_rep = item_restscore$draws$gamma_rep
  )
  obs_df <- data.frame(
    item      = as.character(item_restscore$result[[item_col]]),
    gamma_obs = item_restscore$result$gamma_obs
  )
  item_levels <- sort(unique(obs_df$item))
  draws_df$item <- factor(draws_df$item, levels = rev(item_levels))
  obs_df$item   <- factor(obs_df$item, levels = rev(item_levels))

  p <- ggplot2::ggplot(
    draws_df,
    ggplot2::aes(x = .data$gamma_rep, y = .data$item)
  ) +
    ggdist::stat_slab(
      ggplot2::aes(fill = ggplot2::after_stat(level)),
      .width = c(0.84, 0.95)
    ) +
    ggplot2::geom_point(
      data = obs_df,
      ggplot2::aes(x = .data$gamma_obs, y = .data$item),
      size = 6, shape = 18, colour = "sienna1",
      position = ggplot2::position_nudge(y = 0.08)
    ) +
    ggplot2::scale_fill_manual(
      aesthetics = "fill",
      values = c("grey", "lightgrey")
    ) +
    ggplot2::labs(
      x = "Observed (orange) and replicated item-restscore gamma",
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
#' @param ci Numeric in \eqn{(0, 1)}. Width of the highest density
#'   continuous interval of the replicated (posterior predictive) Q3
#'   distribution, used in the \code{hdi} table and the slab
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
#'     \code{q3_diff} (posterior mean difference), \code{q3_ppp}
#'     (posterior predictive p-value: proportion of draws where the
#'     observed Q3 exceeds the replicated Q3), and \code{q3_p_adj}
#'     (one-sided p-value for local dependence, adjusted for the number
#'     of item pairs; small values flag dependence).}
#'   \item{hdi}{A \code{\link[tibble]{tibble}} with one row per
#'     item pair containing: \code{item_pair}, \code{item_1},
#'     \code{item_2}, \code{ld} (proportion of draws where observed
#'     Q3 exceeds the upper HDI bound of the replicated
#'     distribution), and \code{lr} (proportion of draws where
#'     observed Q3 falls below the lower HDI bound). See Details for
#'     how to read these.}
#'   \item{plot}{A \code{\link[ggplot2]{ggplot}} object showing
#'     the posterior predictive distribution of replicated Q3
#'     (grey filled slab) overlaid with the observed Q3
#'     distribution (coloured slab + interval), with a dashed
#'     reference line at the mean replicated Q3 over all pairs (about
#'     \eqn{-1/(k-1)} for the conditional method, 0 for the theta
#'     method).}
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
#'     (local repulsion, e.g., speed-accuracy tradeoffs).
#'     \code{q3_p_adj} is a single-step max-T adjustment over item
#'     pairs, as \code{p_adj} in \code{\link{item_restscore_statistic}},
#'     but one-sided: each pair's mean observed Q3 is standardised by the
#'     mean and SD of its replicated Q3, and compared with the largest
#'     standardised replicated value over all pairs in each draw. Unlike
#'     \code{q3_ppp}, small values flag local dependence.}
#'   \item{\code{hdi}}{Reports the proportion of draws where the
#'     observed Q3 falls above (local dependence) or below (local
#'     repulsion) the HDI of the replicated distribution. The HDI is
#'     computed over all draws, so each draw's observed value is
#'     compared to the same fixed interval. The observed Q3 varies
#'     across draws only through parameter uncertainty, while the
#'     replicated Q3 also varies through the sampling of new
#'     responses. The observed distribution is therefore narrower,
#'     and these proportions move quickly from near 0 to near 1 as
#'     the observed values pass an interval bound. They are
#'     exceedance proportions, not posterior probabilities of local
#'     dependence.}
#' }
#'
#' With \code{q3_statistic(method = "conditional")}, the default,
#' \code{q3_ppp} was close to uniform in simulations under a fitting Rasch
#' model, and \code{q3_p_adj} kept the proportion of datasets with any
#' pair flagged at .05 near or below .05. With \code{method = "theta"},
#' \code{q3_ppp} is conservative, and so is \code{q3_p_adj} (2\% to 3\% of
#' datasets with any pair flagged at .05). See \code{\link{q3_statistic}}.
#'
#' \strong{Reading \code{q3_p_adj}.} For the pair with the smallest value,
#' \code{q3_p_adj} is a test of local independence over all pairs, and it
#' keeps the familywise error rate at .05 when no pair is dependent. The
#' values for the other pairs do not have this guarantee when some pair is
#' dependent, since all pairs are computed from the same residuals
#' (Ferreira, 2024, Example 2). The suggested use is to act on the most
#' deviant pair, for example by removing an item, combining the two items
#' into a superitem or modelling the dependence, then refit the model and
#' recompute Q3 (Christensen et al., 2017). The adjustment covers the Q3
#' pairs of one call, not other fit statistics computed on the same data.
#'
#' Each pair is compared with its own mean replicated Q3 rather than with
#' the mean over all pairs. Given the total score, expected Q3 differs
#' between pairs, and is highest for pairs of the easiest and hardest
#' items. In simulations, a reference at the mean over all pairs made such
#' a pair the most deviant one in 11\% to 13\% of datasets without local
#' dependence, against about 2\% if every pair were equally likely.
#'
#' @references
#' Christensen, K. B., Makransky, G. & Horton, M. (2017). Critical
#' values for Yen's Q3: Identification of local dependence in the
#' Rasch model using residual correlations.
#' \emph{Applied Psychological Measurement}, \emph{41}(3), 178--194.
#' \doi{10.1177/0146621616677520}
#'
#' Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
#' number of hypotheses simultaneously: An account focusing on the
#' control of the probability of at least one incorrect rejection and of
#' the false discovery rate. \emph{Journal of Statistical Theory and
#' Practice}, \emph{19}, 6. \doi{10.1007/s42519-024-00412-4}
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
#' @importFrom ggdist stat_slab stat_slabinterval median_hdci
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
  if (!all(c("draw", "item_pair", "q3", "q3_rep") %in% names(q3_draws))) {
    stop("'q3_draws' must contain columns 'draw', 'item_pair', 'q3', ",
         "and 'q3_rep'.", call. = FALSE)
  }
  n_per_pair <- table(q3_draws$item_pair)
  if (length(unique(n_per_pair)) != 1L) {
    stop("Every item pair in 'q3_draws' must have the same number of ",
         "draws.", call. = FALSE)
  }
  
  sort_by <- match.arg(sort_by)
  
  # --- Summary table: means and ppp ---
  summary_tbl <- q3_draws |>
    dplyr::group_by(.data$item_pair) |>
    dplyr::summarise(
      item_1  = .data$item_1[1],
      item_2  = .data$item_2[1],
      # q3_diff and q3_ppp must be computed before q3_rep is overwritten
      # by its mean, since later expressions in summarise() see the new
      # column.
      q3_diff = round(mean(.data$q3 - .data$q3_rep, na.rm = TRUE), 3),
      q3_ppp  = round(mean(.data$q3 > .data$q3_rep, na.rm = TRUE), 3),
      q3_obs  = round(mean(.data$q3, na.rm = TRUE), 3),
      q3_rep  = round(mean(.data$q3_rep, na.rm = TRUE), 3),
      .groups = "drop"
    )

  # --- Max-T adjusted one-sided p-value over pairs ---
  wide <- q3_draws |>
    dplyr::ungroup() |>
    dplyr::select("draw", "item_pair", "q3", "q3_rep")
  pairs <- summary_tbl$item_pair
  rep_mat <- vapply(pairs, function(p) {
    d <- wide[wide$item_pair == p, ]
    d$q3_rep[order(d$draw)]
  }, numeric(sum(wide$item_pair == pairs[1])))
  rep_mat <- matrix(rep_mat, ncol = length(pairs))
  obs_mean <- vapply(pairs, function(p) {
    mean(wide$q3[wide$item_pair == p], na.rm = TRUE)
  }, numeric(1))
  summary_tbl$q3_p_adj <- round(.ppp_maxt_upper(obs_mean, rep_mat), 3)

  summary_tbl <- summary_tbl |>
    dplyr::select("item_pair", "item_1", "item_2", "q3_obs", "q3_rep",
                  "q3_diff", "q3_ppp", "q3_p_adj") |>
    dplyr::arrange(dplyr::desc(.data$q3_diff))
  ref_line <- mean(q3_draws$q3_rep, na.rm = TRUE)
  
  # --- HDI table: LD / LR probabilities ---
  q3_rep_hdi <- q3_draws |>
    dplyr::group_by(.data$item_pair) |>
    dplyr::summarise(
      # Use the continuous-interval variant (`_hdci`) rather than
      # `_hdi`. The plain `median_hdi()` returns a *set* of disjoint
      # HDIs when the q3_rep posterior is multimodal at the requested
      # width (e.g., bimodal item-pair residual correlations), i.e.,
      # one row per disjoint segment. Inside dplyr::summarise() that
      # violates the one-row-per-group contract and errors out with
      # "must be size 1, not 2". The continuous version always returns
      # a single [lower, upper] interval, which is the right quantity
      # to use as a cutoff for the LD / LR classification below.
      ggdist::median_hdci(.data$q3_rep, .width = ci),
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
      xintercept = ref_line, linetype = 2
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
#' @param target_link Character string specifying the link function of
#'   the model the priors will be used with. One of \code{"logit"},
#'   \code{"probit"}, or \code{"source"} (the default). When
#'   \code{"source"}, the link function of the fitted model is used
#'   and no transformation is applied. When different from the source
#'   model's link, all location and scale parameters are rescaled
#'   using the approximation \eqn{\beta_{\text{probit}} \approx
#'   \beta_{\text{logit}} / 1.7}. This is useful when transferring
#'   priors from a logit-fitted model to a probit model or vice versa.
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
#' \strong{Independent priors.} Each prior is built from one
#' parameter's marginal posterior, so the posterior correlations
#' between parameters are not carried over. Thresholds are usually
#' positively correlated through the shared location of the scale, so
#' the independent priors are more informative about overall location
#' than the posterior was. For example, the SD of the mean threshold
#' can be noticeably smaller under the priors than under the posterior.
#' Use \code{mult > 1} to widen them if this matters.
#'
#' \strong{Link function transformation:} When \code{target_link}
#' differs from the source model's link function, all parameters
#' (means and SDs) are rescaled by a factor of approximately 1.7.
#' This uses the well-known approximation that
#' \eqn{\Phi(x) \approx \text{logistic}(1.7 \, x)}, so
#' logit-scale parameters can be converted to probit-scale by
#' dividing by 1.7, and vice versa. The approximation is excellent
#' for parameters in the range \eqn{|\beta| < 3} and adequate
#' beyond that range.
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
#'   chains = 4, cores = 2, iter = 1000 # use more iter (and cores if you have)
#' )
#'
#' # Extract posterior-informed priors (same link)
#' new_priors <- posterior_to_prior(fit_pcm)
#' new_priors
#'
#' # Narrow the prior's sd by a factor of 0.5
#' wide_priors <- posterior_to_prior(fit_pcm, mult = 0.5)
#'
#' # Extract priors for use with a probit model
#' probit_priors <- posterior_to_prior(fit_pcm, target_link = "probit")
#'
#' # --- Dichotomous 1PL (fixed item effects) ---
#'
#' # PHQ-9 items dichotomised (any symptom vs none), for illustration only
#' df_rm <- phq9[1:200, 1:9] %>%
#'   mutate(across(everything(), ~ as.integer(.x > 0))) %>%
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
#'
#' # Transfer logit priors to a probit refit
#' priors_probit <- posterior_to_prior(fit_1pl, target_link = "probit")
#' }
#'
#' @importFrom brms as_draws_df set_prior
#' @importFrom rlang enquo as_name
#' @importFrom stats sd family
#' @importFrom tibble as_tibble
#' @export
posterior_to_prior <- function(model,
                               item_var = item,
                               person_var = id,
                               mult = 1,
                               target_link = c("source", "logit", "probit")) {
  
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  if (!is.numeric(mult) || length(mult) != 1 || mult <= 0) {
    stop("'mult' must be a single positive number.", call. = FALSE)
  }
  
  target_link <- match.arg(target_link)
  
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
  
  draws       <- tibble::as_tibble(brms::as_draws_df(model))
  family_name <- stats::family(model)$family
  is_ordinal  <- grepl("acat|cumul|sratio|cratio",
                       family_name, ignore.case = TRUE)
  
  # ================================================================
  # DETERMINE LINK TRANSFORMATION FACTOR
  # ================================================================
  source_link <- stats::family(model)$link
  if (is.null(source_link) || source_link == "") {
    source_link <- "logit"  # brms default for bernoulli/acat
  }
  
  if (target_link == "source") {
    link_factor <- 1
  } else if (source_link == target_link) {
    link_factor <- 1
  } else if (source_link == "logit" && target_link == "probit") {
    link_factor <- 1 / 1.7
  } else if (source_link == "probit" && target_link == "logit") {
    link_factor <- 1.7
  } else {
    warning("Link transformation from '", source_link, "' to '",
            target_link, "' is not supported. No rescaling applied.",
            call. = FALSE)
    link_factor <- 1
  }
  
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
    
    # Thresholds per item, matched to items by brms level order
    tau <- .item_tau_draws(draws, .item_levels(model, item_name), item_name,
                           ordinal = TRUE)
    unique_items <- names(tau)

    for (item_label in unique_items) {
      thresh_mat <- tau[[item_label]]
      for (idx in seq_len(ncol(thresh_mat))) {
        vals   <- thresh_mat[, idx]
        p_mean <- mean(vals) * link_factor
        p_sd   <- stats::sd(vals) * mult * abs(link_factor)
        
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
    fe_pattern <- paste0("^b_", item_name)
    fe_cols    <- grep(fe_pattern, names(draws), value = TRUE)
    fe_cols    <- setdiff(fe_cols, "b_Intercept")
    
    has_fixed_items  <- length(fe_cols) > 0
    has_random_items <- any(grepl(
      paste0("^r_", item_name, "\\["), names(draws)
    ))
    
    # ── Fixed item effects: response ~ 0 + item + (1 | id) ──
    if (has_fixed_items) {
      
      for (col_name in fe_cols) {
        vals   <- draws[[col_name]]
        p_mean <- mean(vals) * link_factor
        p_sd   <- stats::sd(vals) * mult * abs(link_factor)
        
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
        p_mean <- mean(vals) * link_factor
        p_sd   <- stats::sd(vals) * mult * abs(link_factor)
        
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
        p_mean <- mean(vals) * abs(link_factor)
        p_sd   <- stats::sd(vals) * mult * abs(link_factor)
        
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
    p_mean <- mean(vals) * abs(link_factor)
    p_sd   <- stats::sd(vals) * mult * abs(link_factor)
    
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
  
  if (link_factor != 1) {
    message("Priors rescaled from '", source_link, "' to '",
            target_link, "' link (factor = ",
            round(link_factor, 4), ").")
  }
  
  combined_prior
}