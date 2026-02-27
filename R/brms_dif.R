#' Differential Item Functioning (DIF) Analysis for Bayesian IRT Models
#'
#' Tests for differential item functioning (DIF) in Bayesian
#' Rasch-family models fitted with \pkg{brms} by comparing item
#' parameters across subgroups defined by an exogenous variable.
#' The function fits a DIF model that includes group-by-item
#' interactions and summarizes the posterior distribution of the
#' DIF effects.
#'
#' For polytomous models, two types of DIF can be tested:
#' \describe{
#'   \item{Uniform DIF (\code{dif_type = "uniform"}, default)}{
#'     A single location shift per item across groups, modelled as
#'     a \code{group:item} fixed-effect interaction. This tests
#'     whether the average item difficulty differs between groups.}
#'   \item{Non-uniform / threshold-level DIF
#'     (\code{dif_type = "non-uniform"})}{Each item receives
#'     group-specific thresholds via
#'     \code{thres(gr = interaction(item, group))}. DIF effects
#'     are computed as the difference in each threshold between
#'     groups, revealing whether DIF affects specific response
#'     categories.}
#' }
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from
#'   the baseline (no-DIF) model.
#' @param group_var An unquoted variable name identifying the
#'   grouping variable for DIF testing (e.g., \code{gender}).
#'   Must be a factor or character variable with exactly 2 levels
#'   in the current implementation.
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the
#'   person grouping variable. Default is \code{id}.
#' @param data An optional data frame containing all variables
#'   needed for the DIF model, including the group variable. If
#'   \code{NULL} (the default), the function attempts to use
#'   \code{model$data}. Since the baseline model formula typically
#'   does not include the group variable, brms will have dropped
#'   it from the stored model data. In that case, you must supply
#'   the original data frame here.
#' @param dif_type Character. For polytomous ordinal models only.
#'   \code{"uniform"} (the default) tests for a uniform location
#'   shift per item via a \code{group:item} fixed-effect
#'   interaction. \code{"non-uniform"} fits group-specific
#'   thresholds per item and computes per-threshold DIF effects
#'   as the difference between groups. Ignored for dichotomous
#'   models.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the credible
#'   intervals. Default is 0.95.
#' @param rope Numeric. Half-width of the Region of Practical
#'   Equivalence (ROPE) around zero for DIF effects, on the logit
#'   scale. Default is 0.5, corresponding to a practically
#'   negligible DIF effect. Set to 0 to skip ROPE analysis.
#' @param refit Logical. If \code{TRUE} (the default), the DIF
#'   model is fitted automatically by updating the baseline model
#'   via \code{\link[stats]{update}}, which reuses the compiled
#'   Stan code for faster sampling. If \code{FALSE}, only the DIF
#'   model formula is returned (useful for manual fitting with
#'   custom settings).
#' @param ... Additional arguments passed to
#'   \code{\link[brms]{update.brmsfit}} when refitting the DIF
#'   model (e.g., \code{cores}, \code{control}).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{summary}{A \code{\link[tibble]{tibble}} with one row per
#'     item (for uniform DIF) or per item × threshold (for
#'     non-uniform DIF) containing: \code{item}, optionally
#'     \code{threshold}, \code{dif_estimate} (posterior mean),
#'     \code{dif_lower}, \code{dif_upper} (credible interval),
#'     \code{dif_sd} (posterior SD), \code{pd} (probability of
#'     direction), \code{rope_percentage} (proportion inside ROPE),
#'     and \code{flag} (classification).}
#'   \item{dif_draws}{A matrix of posterior draws for the DIF
#'     effects (draws × effects), for further analysis.}
#'   \item{dif_model}{The fitted DIF \code{brmsfit} object (if
#'     \code{refit = TRUE}), or \code{NULL}.}
#'   \item{dif_formula}{The \code{brmsformula} used for the DIF
#'     model.}
#'   \item{baseline_model}{The original baseline model.}
#'   \item{plot}{A \code{\link[ggplot2]{ggplot}} forest plot of
#'     DIF effects with credible intervals and ROPE.}
#' }
#'
#' @details
#' The function constructs a DIF model by adding a group-by-item
#' interaction to the baseline model:
#'
#' \itemize{
#'   \item \strong{Dichotomous models}
#'     (\code{family = bernoulli()}): The baseline
#'     \code{response ~ 1 + (1 | item) + (1 | id)} becomes
#'     \code{response ~ 1 + group + (1 + group | item) + (1 | id)},
#'     where the group slope varying by item captures item-specific
#'     DIF.
#'   \item \strong{Polytomous uniform DIF}
#'     (\code{dif_type = "uniform"}): The baseline
#'     \code{response | thres(gr = item) ~ 1 + (1 | id)} becomes
#'     \code{response | thres(gr = item) ~ 1 + group:item + (1 | id)}.
#'   \item \strong{Polytomous non-uniform DIF}
#'     (\code{dif_type = "non-uniform"}): The baseline becomes
#'     \code{response | thres(gr = item_group) ~ 1 + (1 | id)},
#'     where \code{item_group = interaction(item, group)}. Each
#'     item × group combination gets its own thresholds. DIF
#'     effects are the differences between group-specific
#'     thresholds for each item, computed draw-by-draw from the
#'     posterior.
#' }
#'
#' DIF effects are summarized using:
#' \describe{
#'   \item{Probability of Direction (pd)}{The proportion of the
#'     posterior on the dominant side of zero. Values > 0.975
#'     indicate strong directional evidence.}
#'   \item{ROPE}{The Region of Practical Equivalence (Kruschke,
#'     2018). If > 95\% of the posterior falls inside the ROPE,
#'     the DIF effect is practically negligible. If > 95\% falls
#'     outside, the effect is practically significant.}
#'   \item{Credible Interval}{If the CI excludes zero, there is
#'     evidence of DIF at the specified credibility level.}
#' }
#'
#' @references
#' Kruschke, J. K. (2018). Rejecting or accepting parameter values
#' in Bayesian estimation. \emph{Advances in Methods and Practices
#' in Psychological Science}, \emph{1}(2), 270--280.
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with
#' brms and Stan. \emph{Journal of Statistical Software}, \emph{100},
#' 1--54. \doi{10.18637/jss.v100.i05}
#'
#' Holland, P. W. & Wainer, H. (1993). \emph{Differential Item
#' Functioning}. Lawrence Erlbaum Associates.
#'
#' @seealso
#' \code{\link{infit_statistic}} for item fit,
#' \code{\link{q3_statistic}} for local dependence,
#' \code{\link[brms]{brm}},
#' \code{\link[brms]{hypothesis}}.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#'
#' # --- Dichotomous Rasch with DIF testing ---
#'
#' set.seed(123)
#' df <- expand.grid(id = 1:200, item = paste0("I", 1:10)) %>%
#'   mutate(
#'     gender = rep(sample(c("M", "F"), 200, TRUE), each = 10),
#'     theta  = rep(rnorm(200), each = 10),
#'     delta  = rep(seq(-2, 2, length.out = 10), 200),
#'     dif    = ifelse(item == "I3" & gender == "F", 1.0,
#'              ifelse(item == "I7" & gender == "F", -0.8, 0)),
#'     p      = plogis(theta - delta - dif),
#'     response = rbinom(n(), 1, p)
#'   )
#'
#' fit_base <- brm(
#'   response ~ 1 + (1 | item) + (1 | id),
#'   data   = df,
#'   family = bernoulli(),
#'   chains = 4, cores = 4, iter = 2000
#' )
#'
#' dif_result <- dif_statistic(
#'   model     = fit_base,
#'   group_var = gender,
#'   data      = df
#' )
#'
#' dif_result$summary
#' dif_result$plot
#'
#' # --- Partial Credit Model: uniform DIF ---
#'
#' df_pcm <- eRm::pcmdat2 %>%
#'   mutate(across(everything(), ~ .x + 1)) %>%
#'   rownames_to_column("id") %>%
#'   mutate(gender = sample(c("M", "F"), n(), TRUE)) %>%
#'   pivot_longer(!c(id, gender),
#'                names_to = "item", values_to = "response")
#'
#' fit_pcm <- brm(
#'   response | thres(gr = item) ~ 1 + (1 | id),
#'   data   = df_pcm,
#'   family = acat,
#'   chains = 4, cores = 4, iter = 2000
#' )
#'
#' # Uniform DIF (default): one shift per item
#' dif_uni <- dif_statistic(fit_pcm, group_var = gender, data = df_pcm)
#' dif_uni$plot
#'
#' # Non-uniform DIF: threshold-level effects
#' dif_nu <- dif_statistic(fit_pcm, group_var = gender, data = df_pcm,
#'                          dif_type = "non-uniform")
#' dif_nu$summary
#' dif_nu$plot
#' }
#'
#' @importFrom brms brm bf as_draws_df ndraws
#' @importFrom rlang enquo as_name .data
#' @importFrom stats formula quantile sd update family as.formula
#' @importFrom tibble tibble
#' @export
dif_statistic <- function(
    model,
    group_var,
    item_var = item,
    person_var = id,
    data = NULL,
    dif_type = c("uniform", "non-uniform"),
    prob = 0.95,
    rope = 0.5,
    refit = TRUE,
    ...
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }

  dif_type    <- match.arg(dif_type)
  item_name   <- rlang::as_name(rlang::enquo(item_var))
  person_name <- rlang::as_name(rlang::enquo(person_var))
  group_name  <- rlang::as_name(rlang::enquo(group_var))

  # --- Resolve data source ---
  if (is.null(data)) {
    dif_data <- model$data
  } else {
    dif_data <- as.data.frame(data)
  }

  if (!item_name %in% names(dif_data)) {
    stop("Item variable '", item_name, "' not found in data.",
         call. = FALSE)
  }
  if (!person_name %in% names(dif_data)) {
    stop("Person variable '", person_name, "' not found in data.",
         call. = FALSE)
  }
  if (!group_name %in% names(dif_data)) {
    stop("Group variable '", group_name, "' not found in data. ",
         "If '", group_name, "' was not part of the baseline model ",
         "formula, brms will have dropped it from model$data. ",
         "Please supply the original data frame via the 'data' ",
         "argument.", call. = FALSE)
  }

  # Ensure group is factor
  if (!is.factor(dif_data[[group_name]])) {
    dif_data[[group_name]] <- as.factor(dif_data[[group_name]])
  }
  group_vals   <- dif_data[[group_name]]
  group_levels <- levels(group_vals)
  if (length(group_levels) != 2) {
    stop("'", group_name, "' must have exactly 2 levels. Found: ",
         length(group_levels), ".", call. = FALSE)
  }

  lower_prob <- (1 - prob) / 2
  upper_prob <- 1 - lower_prob

  family_name  <- stats::family(model)$family
  is_ordinal   <- grepl("acat|cumul|sratio|cratio",
                        family_name, ignore.case = TRUE)
  unique_items <- unique(dif_data[[item_name]])

  # ================================================================
  # CONSTRUCT DIF FORMULA
  # ================================================================
  base_formula <- formula(model)$formula
  resp_var <- as.character(base_formula[[2]])
  if (length(resp_var) > 1) resp_var <- resp_var[2]

  if (is_ordinal) {
    resp_text <- deparse(base_formula[[2]], width.cutoff = 500)
    rhs_text  <- deparse(base_formula[[3]], width.cutoff = 500)

    if (dif_type == "uniform") {
      dif_formula_text <- paste0(
        resp_text, " ~ ", rhs_text, " + ",
        group_name, ":", item_name
      )
      dif_formula <- brms::bf(stats::as.formula(dif_formula_text))
    } else {
      # Non-uniform: create item_group interaction for thres(gr = ...)
      ig_var <- ".item_group_dif"
      dif_data[[ig_var]] <- interaction(
        dif_data[[item_name]], dif_data[[group_name]],
        drop = TRUE, sep = "___"
      )

      dif_formula_text <- paste0(
        resp_var, " | thres(gr = ", ig_var, ") ~ 1 + (1 | ",
        person_name, ")"
      )
      dif_formula <- brms::bf(stats::as.formula(dif_formula_text))
    }
  } else {
    # Dichotomous: always random slope
    dif_formula <- brms::bf(stats::as.formula(paste0(
      deparse(base_formula[[2]], width.cutoff = 500),
      " ~ 1 + ", group_name,
      " + (1 + ", group_name, " | ", item_name, ")",
      " + (1 | ", person_name, ")"
    )))
  }

  # ================================================================
  # FIT DIF MODEL
  # ================================================================
  dif_model <- NULL

  if (refit) {
    dif_model <- stats::update(
      model,
      formula. = dif_formula,
      newdata  = dif_data,
      ...
    )
  }

  if (is.null(dif_model) && !refit) {
    message("DIF formula constructed but model not fitted ",
            "(refit = FALSE).")
    message("Formula: ",
            deparse(dif_formula$formula, width.cutoff = 500))
    return(list(
      summary        = NULL,
      dif_draws      = NULL,
      dif_model      = NULL,
      dif_formula    = dif_formula,
      baseline_model = model,
      plot           = NULL
    ))
  }

  # ================================================================
  # EXTRACT DIF EFFECTS
  # ================================================================
  draws <- brms::as_draws_df(dif_model)
  group_level <- group_levels[2]

  if (is_ordinal && dif_type == "uniform") {
    dif_info <- .extract_dif_uniform(
      draws, unique_items, item_name, group_name, group_level
    )
  } else if (is_ordinal && dif_type == "non-uniform") {
    dif_info <- .extract_dif_threshold(
      draws, unique_items, group_levels
    )
  } else {
    dif_info <- .extract_dif_dichotomous(
      draws, unique_items, item_name, group_name, group_level
    )
  }

  dif_draws_mat <- dif_info$draws_mat
  effect_labels <- dif_info$labels

  # ================================================================
  # SUMMARIZE DIF EFFECTS
  # ================================================================
  n_effects <- ncol(dif_draws_mat)

  dif_mean  <- colMeans(dif_draws_mat, na.rm = TRUE)
  dif_sd    <- apply(dif_draws_mat, 2, stats::sd, na.rm = TRUE)
  dif_lower <- apply(dif_draws_mat, 2, stats::quantile,
                     probs = lower_prob, na.rm = TRUE)
  dif_upper <- apply(dif_draws_mat, 2, stats::quantile,
                     probs = upper_prob, na.rm = TRUE)

  pd <- apply(dif_draws_mat, 2, function(x) {
    p_pos <- mean(x > 0, na.rm = TRUE)
    max(p_pos, 1 - p_pos)
  })

  if (rope > 0) {
    rope_pct <- apply(dif_draws_mat, 2, function(x) {
      mean(abs(x) <= rope, na.rm = TRUE)
    })
  } else {
    rope_pct <- rep(NA_real_, n_effects)
  }

  flag <- character(n_effects)
  for (i in seq_len(n_effects)) {
    ci_excludes_zero <- (dif_lower[i] > 0) | (dif_upper[i] < 0)
    if (rope > 0) {
      mostly_outside_rope <- rope_pct[i] < 0.05
      mostly_inside_rope  <- rope_pct[i] > 0.95
    } else {
      mostly_outside_rope <- FALSE
      mostly_inside_rope  <- FALSE
    }

    if (ci_excludes_zero && mostly_outside_rope) {
      flag[i] <- "DIF detected"
    } else if (ci_excludes_zero) {
      flag[i] <- "DIF likely"
    } else if (mostly_inside_rope) {
      flag[i] <- "Negligible"
    } else {
      flag[i] <- "Inconclusive"
    }
  }

  summary_df <- tibble::tibble(
    item            = effect_labels$item,
    dif_estimate    = dif_mean,
    dif_sd          = dif_sd,
    dif_lower       = dif_lower,
    dif_upper       = dif_upper,
    pd              = pd,
    rope_percentage = rope_pct,
    flag            = flag
  )

  if ("threshold" %in% names(effect_labels)) {
    summary_df$threshold <- effect_labels$threshold
    summary_df <- summary_df[, c("item", "threshold",
                                 setdiff(names(summary_df),
                                         c("item", "threshold")))]
  }

  # ================================================================
  # FOREST PLOT
  # ================================================================
  has_threshold <- "threshold" %in% names(summary_df)

  if (has_threshold) {
    summary_df$label <- paste0(summary_df$item, " [\u03C4",
                               summary_df$threshold, "]")
    summary_df$label <- factor(
      summary_df$label,
      levels = rev(unique(summary_df$label))
    )
    y_var <- "label"

    # Map thresholds to shapes
    unique_thresh <- sort(unique(summary_df$threshold))
    shape_vals <- c(16, 17, 15, 18, 8, 3, 4, 1, 2, 0)
    shape_vals <- shape_vals[seq_along(unique_thresh)]
    names(shape_vals) <- unique_thresh
    summary_df$threshold_f <- factor(summary_df$threshold,
                                     levels = unique_thresh)
  } else {
    summary_df$item_f <- factor(
      summary_df$item,
      levels = summary_df$item[order(summary_df$dif_estimate)]
    )
    y_var <- "item_f"
  }

  flag_colors <- c(
    "DIF detected"  = "#D55E00",
    "DIF likely"    = "#E69F00",
    "Inconclusive"  = "#56B4E9",
    "Negligible"    = "#009E73"
  )

  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x      = .data$dif_estimate,
      y      = .data[[y_var]],
      colour = .data$flag
    )
  )

  # Alternating item background bands
  if (has_threshold) {
    band_items <- unique(summary_df$item)
    label_levels <- levels(summary_df$label)
    for (b in seq_along(band_items)) {
      if (b %% 2 == 0) {
        rows <- summary_df$item == band_items[b]
        positions <- match(summary_df$label[rows], label_levels)
        ymin_b <- min(positions) - 0.5
        ymax_b <- max(positions) + 0.5
        p <- p +
          ggplot2::annotate(
            "rect",
            xmin = -Inf, xmax = Inf,
            ymin = ymin_b, ymax = ymax_b,
            fill = "grey95", alpha = 0.6
          )
      }
    }
  } else {
    item_levels <- levels(summary_df$item_f)
    for (b in seq_along(item_levels)) {
      if (b %% 2 == 0) {
        p <- p +
          ggplot2::annotate(
            "rect",
            xmin = -Inf, xmax = Inf,
            ymin = b - 0.5, ymax = b + 0.5,
            fill = "grey95", alpha = 0.6
          )
      }
    }
  }

  if (rope > 0) {
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin = -rope, xmax = rope,
        ymin = -Inf, ymax = Inf,
        fill = "grey85", alpha = 0.4
      )
  }

  p <- p +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", colour = "grey50"
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$dif_lower, xmax = .data$dif_upper),
      linewidth = 0.5
    )

  if (has_threshold) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(shape = .data$threshold_f),
        size = 2.5
      ) +
      ggplot2::scale_shape_manual(
        values = shape_vals,
        name   = "Threshold"
      )
  } else {
    p <- p +
      ggplot2::geom_point(size = 2.5)
  }

  p <- p +
    ggplot2::scale_colour_manual(
      values = flag_colors,
      name   = "DIF status"
    ) +
    ggplot2::labs(
      x = paste0("DIF effect (", group_levels[2],
                 " vs. ", group_levels[1], ")"),
      y = NULL,
      subtitle = if (rope > 0) {
        paste0(
          if (has_threshold) "Threshold-level DIF | " else "",
          "ROPE = [", -rope, ", ", rope,
          "], shaded region = negligible DIF"
        )
      } else {
        if (has_threshold) "Threshold-level DIF" else NULL
      }
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background  = ggplot2::element_rect(fill = "white",
                                               colour = NA),
      legend.position  = "bottom"
    )

  # Clean up temporary plot columns
  summary_df$label       <- NULL
  summary_df$item_f      <- NULL
  summary_df$threshold_f <- NULL

  list(
    summary        = summary_df,
    dif_draws      = dif_draws_mat,
    dif_model      = dif_model,
    dif_formula    = dif_formula,
    baseline_model = model,
    plot           = p
  )
}


# ── Internal: extract uniform DIF (fixed interaction) ────────────
#' @keywords internal
.extract_dif_uniform <- function(draws, unique_items, item_name,
                                 group_name, group_level) {
  n_draws <- nrow(draws)
  k <- length(unique_items)
  dif_draws_mat <- matrix(0, nrow = n_draws, ncol = k)
  colnames(dif_draws_mat) <- unique_items

  for (i in seq_along(unique_items)) {
    item_label <- unique_items[i]
    item_esc <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label)

    patterns <- c(
      paste0("^b_", group_name, group_level, ":",
             item_name, item_label, "$"),
      paste0("^b_", item_name, item_label, ":",
             group_name, group_level, "$"),
      paste0("^b_", group_name, group_level, ":",
             item_name, item_esc, "$"),
      paste0("^b_", item_name, item_esc, ":",
             group_name, group_level, "$")
    )

    col_found <- NULL
    for (pat in patterns) {
      candidates <- grep(pat, names(draws), value = TRUE)
      if (length(candidates) == 1) {
        col_found <- candidates
        break
      }
    }

    if (!is.null(col_found)) {
      dif_draws_mat[, i] <- draws[[col_found]]
    }
  }

  list(
    draws_mat = dif_draws_mat,
    labels    = data.frame(item = unique_items,
                           stringsAsFactors = FALSE)
  )
}


# ── Internal: extract threshold-level DIF ────────────────────────
# For non-uniform DIF with thres(gr = interaction(item, group)):
# Each item × group gets its own thresholds. DIF = difference
# between group-specific thresholds for each item, per threshold.
#
# brms names these: b_Intercept[ItemLabel___GroupLevel,ThreshNum]
#' @keywords internal
.extract_dif_threshold <- function(draws, unique_items, group_levels) {

  n_draws <- nrow(draws)
  sep <- "___"

  # Find all threshold columns
  thresh_cols <- grep("^b_Intercept\\[", names(draws), value = TRUE)

  if (length(thresh_cols) == 0) {
    stop("Could not find threshold parameters in DIF model.",
         call. = FALSE)
  }

  # Parse: b_Intercept[ItemLabel___GroupLevel,ThreshNum]
  parsed <- data.frame(
    col   = thresh_cols,
    inner = gsub("^b_Intercept\\[(.+)\\]$", "\\1", thresh_cols),
    stringsAsFactors = FALSE
  )
  parsed$ig_label <- sub(",\\d+$", "", parsed$inner)
  parsed$thresh   <- as.integer(sub("^.+,(\\d+)$", "\\1", parsed$inner))

  # Split ig_label into item and group
  parsed$item  <- sub(paste0(sep, ".*$"), "", parsed$ig_label)
  parsed$group <- sub(paste0("^.*", sep), "", parsed$ig_label)

  # Separate by group
  items_g1 <- parsed[parsed$group == group_levels[1], ]
  items_g2 <- parsed[parsed$group == group_levels[2], ]

  # Build DIF draws: for each item × threshold, group2 minus group1
  result_items  <- character(0)
  result_thresh <- character(0)
  result_draws  <- list()

  for (item_label in unique_items) {
    cols_g1 <- items_g1[items_g1$item == item_label, ]
    cols_g2 <- items_g2[items_g2$item == item_label, ]

    if (nrow(cols_g1) == 0 || nrow(cols_g2) == 0) next

    cols_g1 <- cols_g1[order(cols_g1$thresh), ]
    cols_g2 <- cols_g2[order(cols_g2$thresh), ]

    n_thresh <- min(nrow(cols_g1), nrow(cols_g2))

    for (t in seq_len(n_thresh)) {
      dif_draw <- draws[[cols_g2$col[t]]] - draws[[cols_g1$col[t]]]
      result_items  <- c(result_items, item_label)
      result_thresh <- c(result_thresh, as.character(t))
      result_draws[[length(result_draws) + 1]] <- dif_draw
    }
  }

  if (length(result_draws) == 0) {
    stop("Could not compute threshold-level DIF effects. ",
         "Check that thres(gr = ...) produced item x group ",
         "threshold parameters.", call. = FALSE)
  }

  dif_draws_mat <- do.call(cbind, result_draws)
  colnames(dif_draws_mat) <- paste0(result_items, "_t", result_thresh)

  list(
    draws_mat = dif_draws_mat,
    labels    = data.frame(
      item      = result_items,
      threshold = result_thresh,
      stringsAsFactors = FALSE
    )
  )
}


# ── Internal: extract dichotomous DIF (random slopes) ────────────
#' @keywords internal
.extract_dif_dichotomous <- function(draws, unique_items, item_name,
                                     group_name, group_level) {
  n_draws <- nrow(draws)
  k <- length(unique_items)
  dif_draws_mat <- matrix(0, nrow = n_draws, ncol = k)
  colnames(dif_draws_mat) <- unique_items

  for (i in seq_along(unique_items)) {
    item_label <- unique_items[i]
    item_esc <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label)

    re_pattern <- paste0(
      "^r_", item_name, "\\[", item_esc,
      ",", group_name, group_level, "\\]$"
    )
    re_col <- grep(re_pattern, names(draws), value = TRUE)

    if (length(re_col) == 1) {
      dif_draws_mat[, i] <- draws[[re_col]]
    } else {
      re_col2 <- grep(
        paste0("^r_", item_name, "\\[", item_esc, ",.*",
               group_name, ".*\\]$"),
        names(draws), value = TRUE
      )
      if (length(re_col2) == 1) {
        dif_draws_mat[, i] <- draws[[re_col2]]
      }
    }
  }

  list(
    draws_mat = dif_draws_mat,
    labels    = data.frame(item = unique_items,
                           stringsAsFactors = FALSE)
  )
}
