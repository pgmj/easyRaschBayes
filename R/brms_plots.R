#' Item Category Probability Function Curves for Polytomous IRT Models
#'
#' Plots item category probability functions (ICPFs) for polytomous
#' Bayesian IRT models fitted with \pkg{brms}. For each item, the
#' probability of endorsing each response category is plotted as a
#' function of the latent variable (theta), with separate colored
#' curves per category. All items are displayed in a combined faceted
#' plot, similar to the trace plots produced by
#' \code{\link[mirt]{itemplot}} in the \pkg{mirt} package.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from a
#'   polytomous IRT model (e.g., \code{family = acat} for a partial
#'   credit model or \code{family = cumulative} for a graded response
#'   model).
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data (e.g., \code{item}).
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data (e.g., \code{id}).
#' @param items An optional character vector of item names to plot.
#'   If \code{NULL} (the default), all items in the model are plotted.
#' @param theta_range A numeric vector of length 2 specifying the range
#'   of the latent variable (theta) for the x-axis. Default is
#'   \code{c(-4, 4)}.
#' @param n_points Integer. Number of evenly spaced theta values at
#'   which to evaluate the category probabilities. Default is 100.
#' @param ncol Integer. Number of columns in the faceted plot layout.
#'   If \code{NULL} (the default), an appropriate number is chosen
#'   automatically.
#' @param line_size Numeric. Line width for the probability curves.
#'   Default is 0.8.
#' @param ribbon_alpha Numeric in \eqn{[0, 1]}. Transparency of the
#'   credible interval ribbons. Default is 0.15. Set to 0 to hide
#'   ribbons.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the credible interval
#'   for the ribbons. Default is 0.95.
#' @param category_labels An optional character vector of labels for
#'   the response categories. If \code{NULL} (the default), categories
#'   are labelled as integers starting from 1.
#' @param palette An optional character vector of colors, one per
#'   response category. If \code{NULL} (the default), the
#'   \code{viridis} discrete scale from \pkg{ggplot2} is used.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The plot can be
#'   further customised using standard \pkg{ggplot2} functions.
#'
#' @details
#' The function computes category probabilities directly from the
#' posterior draws of the item threshold parameters. For the brms
#' \code{acat} (adjacent category / partial credit) family with
#' logit link, the density is:
#' \deqn{P(Y = y | \eta) = \frac{\exp\bigl(\sum_{k=1}^{y}(\eta -
#'   \tau_k)\bigr)}{\sum_{k=0}^{K} \exp\bigl(\sum_{j=1}^{k}(\eta -
#'   \tau_j)\bigr)}}
#' where \eqn{\eta} is the linear predictor (i.e., theta for a Rasch
#' model with no additional fixed effects) and \eqn{\tau_k} are the
#' item thresholds. Analogous formulas are used for the
#' \code{cumulative}, \code{sratio}, and \code{cratio} families.
#'
#' Posterior uncertainty in the thresholds propagates into credible
#' interval ribbons around the category probability curves — a
#' Bayesian advantage over point-estimate-based plots from packages
#' like \pkg{mirt} or \pkg{eRm}.
#'
#' @references
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with
#' brms and Stan. \emph{Journal of Statistical Software}, \emph{100},
#' 1--54. \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link[brms]{posterior_epred}},
#' \code{\link[brms]{conditional_effects}},
#' \code{\link[mirt]{itemplot}}.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(ggplot2)
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
#' # Plot all items
#' plot_ipf(fit_pcm, item_var = item, person_var = id)
#'
#' # Plot a subset of items
#' plot_ipf(fit_pcm, item_var = item, person_var = id,
#'          items = c("I1", "I2", "I3"))
#'
#' # Customise appearance
#' plot_ipf(fit_pcm, item_var = item, person_var = id,
#'          theta_range = c(-6, 6), ncol = 3, prob = 0.90) +
#'   theme_minimal() +
#'   labs(title = "Item Category Probability Functions")
#' }
#'
#' @importFrom brms as_draws_df
#' @importFrom rlang enquo as_name
#' @importFrom stats formula quantile plogis family
#' @importFrom tidyr as_tibble
#' @export
plot_ipf <- function(
    model,
    item_var = item,
    person_var = id,
    items = NULL,
    theta_range = c(-4, 4),
    n_points = 100,
    ncol = NULL,
    line_size = 0.8,
    ribbon_alpha = 0.15,
    prob = 0.95,
    category_labels = NULL,
    palette = NULL
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.", call. = FALSE)
  }
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }

  item_name   <- rlang::as_name(rlang::enquo(item_var))
  person_name <- rlang::as_name(rlang::enquo(person_var))

  # --- Extract response variable name ---
  resp_var <- as.character(formula(model)$formula[[2]])
  if (length(resp_var) > 1) {
    resp_var <- resp_var[2]
  }

  # --- Validate inputs ---
  if (!item_name %in% names(model$data)) {
    stop("Item variable '", item_name, "' not found in model data.",
         call. = FALSE)
  }
  if (!person_name %in% names(model$data)) {
    stop("Person variable '", person_name, "' not found in model data.",
         call. = FALSE)
  }

  unique_items <- sort(unique(model$data[[item_name]]))

  if (!is.null(items)) {
    invalid <- setdiff(items, unique_items)
    if (length(invalid) > 0) {
      stop("Items not found in model data: ",
           paste(invalid, collapse = ", "), call. = FALSE)
    }
    unique_items <- unique_items[unique_items %in% items]
  }

  # --- Create theta grid ---
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = n_points)

  # --- Extract posterior draws ---
  draws <- tidyr::as_tibble(brms::as_draws_df(model))
  n_draws_total <- nrow(draws)
  family_name <- stats::family(model)$family
  is_acat   <- grepl("acat", family_name, ignore.case = TRUE)
  is_cumul  <- grepl("cumul", family_name, ignore.case = TRUE)
  is_sratio <- grepl("sratio", family_name, ignore.case = TRUE)
  is_cratio <- grepl("cratio", family_name, ignore.case = TRUE)

  if (!any(is_acat, is_cumul, is_sratio, is_cratio)) {
    stop("This function is designed for ordinal models (acat, cumulative, ",
         "sratio, cratio). Your model uses family '", family_name, "'.",
         call. = FALSE)
  }

  # Use a subsample of draws for speed
  max_draws <- min(n_draws_total, 500)
  draw_sample <- sample(seq_len(n_draws_total), max_draws)

  lower_prob <- (1 - prob) / 2
  upper_prob <- 1 - lower_prob

  # --- Compute category probabilities per item ---
  plot_data_list <- list()

  for (item_label in unique_items) {
    # Find threshold columns for this item
    # Grouped thresholds: b_Intercept[item_label,1], etc.
    thresh_pattern <- paste0(
      "^b_Intercept\\[",
      gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
      ","
    )
    thresh_cols <- grep(thresh_pattern, names(draws), value = TRUE)

    if (length(thresh_cols) == 0) {
      # Ungrouped thresholds: b_Intercept[1], b_Intercept[2], ...
      thresh_cols <- grep("^b_Intercept\\[\\d+\\]$", names(draws),
                          value = TRUE)
    }

    if (length(thresh_cols) == 0) {
      warning("Could not find threshold parameters for item '",
              item_label, "'. Skipping.", call. = FALSE)
      next
    }

    # Sort thresholds numerically
    thresh_nums <- as.numeric(
      gsub(".*,(\\d+)\\]$|.*\\[(\\d+)\\]$", "\\1\\2", thresh_cols)
    )
    thresh_cols <- thresh_cols[order(thresh_nums)]

    n_thresholds <- length(thresh_cols)
    n_cat <- n_thresholds + 1
    thresh_sub <- as.matrix(draws[draw_sample, thresh_cols, drop = FALSE])

    # Preallocate
    prob_mean  <- matrix(NA_real_, nrow = n_points, ncol = n_cat)
    prob_lower <- matrix(NA_real_, nrow = n_points, ncol = n_cat)
    prob_upper <- matrix(NA_real_, nrow = n_points, ncol = n_cat)

    for (t in seq_along(theta_grid)) {
      theta <- theta_grid[t]
      prob_draws <- matrix(NA_real_, nrow = max_draws, ncol = n_cat)

      if (is_acat) {
        # Adjacent category model (brms formula, logit link):
        #   P(Y = y) = exp(sum_{k=1}^{y} (eta - tau_k)) /
        #              sum_{k=0}^{K} exp(sum_{j=1}^{k} (eta - tau_j))
        # where eta = theta (the person parameter)
        for (s in seq_len(max_draws)) {
          tau_s <- thresh_sub[s, ]
          # For y = 0 (first category): cumsum = 0
          # For y = c: cumsum of (theta - tau_1), ..., (theta - tau_c)
          cum_eta <- c(0, cumsum(theta - tau_s))
          # Log-sum-exp for numerical stability
          max_eta <- max(cum_eta)
          log_denom <- max_eta + log(sum(exp(cum_eta - max_eta)))
          prob_draws[s, ] <- exp(cum_eta - log_denom)
        }
      } else if (is_cumul) {
        # Cumulative model:
        #   P(Y <= c) = logistic(tau_c - eta)
        for (s in seq_len(max_draws)) {
          tau_s <- thresh_sub[s, ]
          cum_probs <- stats::plogis(tau_s - theta)
          probs <- diff(c(0, cum_probs, 1))
          probs <- pmax(probs, 0)
          prob_draws[s, ] <- probs / sum(probs)
        }
      } else if (is_sratio) {
        # Stopping ratio model:
        #   f(y) = g(tau_{y+1} - eta) * prod_{k=1}^{y} (1 - g(tau_k - eta))
        for (s in seq_len(max_draws)) {
          tau_s <- thresh_sub[s, ]
          probs <- numeric(n_cat)
          remaining <- 1
          for (cc in seq_len(n_thresholds)) {
            p_stop <- stats::plogis(tau_s[cc] - theta)
            probs[cc] <- remaining * p_stop
            remaining <- remaining * (1 - p_stop)
          }
          probs[n_cat] <- remaining
          prob_draws[s, ] <- probs
        }
      } else if (is_cratio) {
        # Continuation ratio model:
        #   f(y) = (1 - g(eta - tau_{y+1})) * prod_{k=1}^{y} g(eta - tau_k)
        for (s in seq_len(max_draws)) {
          tau_s <- thresh_sub[s, ]
          probs <- numeric(n_cat)
          remaining <- 1
          for (cc in seq_len(n_thresholds)) {
            p_cont <- stats::plogis(theta - tau_s[cc])
            probs[cc] <- remaining * (1 - p_cont)
            remaining <- remaining * p_cont
          }
          probs[n_cat] <- remaining
          prob_draws[s, ] <- probs
        }
      }

      prob_mean[t, ]  <- colMeans(prob_draws)
      prob_lower[t, ] <- apply(prob_draws, 2, stats::quantile,
                               probs = lower_prob)
      prob_upper[t, ] <- apply(prob_draws, 2, stats::quantile,
                               probs = upper_prob)
    }

    # Category labels
    if (is.null(category_labels)) {
      cat_labs <- as.character(seq_len(n_cat))
    } else {
      if (length(category_labels) != n_cat) {
        stop("Length of 'category_labels' (", length(category_labels),
             ") does not match the number of categories (", n_cat,
             ") for item '", item_label, "'.", call. = FALSE)
      }
      cat_labs <- category_labels
    }

    for (cc in seq_len(n_cat)) {
      plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
        item     = item_label,
        theta    = theta_grid,
        category = cat_labs[cc],
        prob     = prob_mean[, cc],
        lower    = prob_lower[, cc],
        upper    = prob_upper[, cc],
        stringsAsFactors = FALSE
      )
    }
  }

  plot_data <- do.call(rbind, plot_data_list)
  plot_data$item <- factor(plot_data$item, levels = unique_items)
  if (is.null(category_labels)) {
    plot_data$category <- factor(
      plot_data$category,
      levels = as.character(seq_len(max(as.integer(plot_data$category))))
    )
  } else {
    plot_data$category <- factor(plot_data$category, levels = category_labels)
  }

  # --- Layout ---
  n_items <- length(unique_items)
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_items))
  }

  # --- Build plot ---
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x      = .data$theta,
      y      = .data$prob,
      colour = .data$category,
      fill   = .data$category
    )
  )

  if (ribbon_alpha > 0) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      alpha = ribbon_alpha,
      linewidth = 0
    )
  }

  p <- p +
    ggplot2::geom_line(linewidth = line_size) +
    ggplot2::facet_wrap(~ item, ncol = ncol) +
    ggplot2::labs(
      x      = expression(theta),
      y      = "Probability",
      colour = "Category",
      fill   = "Category"
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text       = ggplot2::element_text(face = "bold")
    )

  if (!is.null(palette)) {
    p <- p +
      ggplot2::scale_colour_manual(values = palette) +
      ggplot2::scale_fill_manual(values = palette)
  } else {
    p <- p +
      ggplot2::scale_colour_viridis_d(end = 0.9) +
      ggplot2::scale_fill_viridis_d(end = 0.9)
  }

  p
}


#' Person-Item Map (Targeting Plot) for Bayesian IRT Models
#'
#' Plots a person-item map (also known as a Wright map or targeting
#' plot) for Bayesian IRT models fitted with \pkg{brms}. The plot
#' consists of three vertically stacked panels sharing the same
#' latent variable (theta / logit) x-axis:
#'
#' \enumerate{
#'   \item \strong{Top}: A histogram of person ability estimates,
#'     with a reference line for the mean (or median) and shading
#'     for ±1 SD (or ±1 MAD).
#'   \item \strong{Middle}: An inverted histogram of item threshold
#'     locations, with a reference line for the mean (or median) and
#'     shading for ±1 SD (or ±1 MAD), mirroring the top panel to
#'     visualise the overlap between person abilities and item
#'     difficulties.
#'   \item \strong{Bottom}: A dot-and-whisker plot of item thresholds
#'     by item, with credible intervals and color-coded response
#'     categories.
#' }
#'
#' Together, the top and middle panels form a half-moon (or
#' back-to-back histogram) display that makes it easy to assess
#' whether the test is well-targeted to the sample.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from an
#'   ordinal IRT model (e.g., \code{family = acat}) or a dichotomous
#'   model (\code{family = bernoulli()}).
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Default is \code{id}.
#' @param robust Logical. If \code{FALSE} (the default), the
#'   histogram annotations use mean ± SD. If \code{TRUE}, median ±
#'   MAD is used instead.
#' @param center Logical. If \code{TRUE} (the default), the scale is
#'   recentered so that the grand mean of all item threshold
#'   locations is zero, following the convention in frequentist
#'   Rasch analysis. Person estimates are shifted by the same
#'   constant. If \code{FALSE}, the raw brms parameterisation is
#'   used.
#' @param sort_items Character. How to order items on the y-axis of
#'   the bottom panel. \code{"data"} (the default) preserves the
#'   order in which items first appear in the model data, with the
#'   first item at the top. \code{"location"} sorts items by their
#'   mean threshold location (easiest at top, hardest at bottom).
#' @param bins Integer. Number of bins for both histograms. Default
#'   is 30.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the credible
#'   intervals for the item threshold whiskers. Default is 0.95.
#' @param palette An optional character vector of colors for the
#'   response categories. If \code{NULL} (the default), the
#'   \code{viridis} discrete scale is used.
#' @param person_fill Fill color for the person histogram. Default is
#'   \code{"#0072B2"} (blue).
#' @param threshold_fill Fill color for the threshold histogram.
#'   Default is \code{"#D55E00"} (vermillion).
#' @param height_ratios Numeric vector of length 3 specifying the
#'   relative heights of the top (person), middle (threshold), and
#'   bottom (dot-whisker) panels. Default is \code{c(3, 2, 5)}.
#'
#' @return A \code{patchwork} object (combined \code{ggplot}).
#'
#' @details
#' \strong{Person estimates} are obtained as the posterior means of
#' the person random effects from the fitted model via
#' \code{\link[brms]{ranef}}.
#'
#' \strong{Item thresholds} are extracted from the posterior draws.
#' For models with grouped thresholds (\code{thres(gr = item)}),
#' each item has its own set of threshold parameters. For models
#' with a single set of thresholds (e.g., dichotomous Rasch with
#' \code{(1 | item)}), the item random effects are subtracted from
#' the global thresholds to obtain item-specific locations.
#'
#' When \code{center = TRUE} (the default), the grand mean of all
#' item threshold posterior means is computed and subtracted from
#' every threshold estimate, its credible interval bounds, and every
#' person estimate. This is a uniform translation of the entire
#' scale that preserves all relative distances and matches the
#' zero-centered item difficulty convention used in frequentist CML
#' estimation.
#'
#' @references
#' Wright, B. D. & Stone, M. H. (1979). \emph{Best Test Design}.
#' MESA Press.
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with
#' brms and Stan. \emph{Journal of Statistical Software}, \emph{100},
#' 1--54. \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link{plot_ipf}} for item category probability curves,
#' \code{\link[brms]{ranef}},
#' \code{\link[brms]{as_draws_df}}.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(ggplot2)
#' library(patchwork)
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
#' # Default: centered, mean ± SD, items in data order
#' plot_targeting(fit_pcm)
#'
#' # Uncentered (raw brms parameterisation)
#' plot_targeting(fit_pcm, center = FALSE)
#'
#' # Robust: median ± MAD, items sorted by location
#' plot_targeting(fit_pcm, robust = TRUE, sort_items = "location")
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
#' plot_targeting(fit_rm, sort_items = "location")
#' }
#'
#' @importFrom brms ranef as_draws_df
#' @importFrom rlang enquo as_name .data
#' @importFrom stats formula quantile median mad sd aggregate family
#' @export
plot_targeting <- function(
    model,
    item_var = item,
    person_var = id,
    robust = FALSE,
    center = TRUE,
    sort_items = c("data", "location"),
    bins = 30,
    prob = 0.95,
    palette = NULL,
    person_fill = "#0072B2",
    threshold_fill = "#D55E00",
    height_ratios = c(3, 2, 5)
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.", call. = FALSE)
  }
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }

  sort_items  <- match.arg(sort_items)
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

  lower_prob <- (1 - prob) / 2
  upper_prob <- 1 - lower_prob

  # ================================================================
  # PERSON ESTIMATES
  # ================================================================
  person_re <- brms::ranef(model)[[person_name]]
  if (is.null(person_re)) {
    stop("No random effects found for '", person_name,
         "'. Check your model formula.", call. = FALSE)
  }
  person_theta <- person_re[, "Estimate", "Intercept"]

  # ================================================================
  # ITEM THRESHOLDS
  # ================================================================
  draws <- brms::as_draws_df(model)
  unique_items_data_order <- unique(model$data[[item_name]])
  unique_items <- sort(unique(model$data[[item_name]]))
  family_name <- stats::family(model)$family
  is_ordinal <- grepl("acat|cumul|sratio|cratio",
                      family_name, ignore.case = TRUE)

  threshold_data <- .extract_threshold_data(
    draws, model, unique_items, item_name, person_name,
    is_ordinal, lower_prob, upper_prob
  )

  # ================================================================
  # CENTERING (optional)
  # ================================================================
  # The entire scale (persons AND thresholds) is shifted by the same

  # constant so that the grand mean of item thresholds = 0.
  # Both move in the same direction to preserve relative positions.
  if (center) {
    shift <- mean(threshold_data$estimate)

    threshold_data$estimate <- threshold_data$estimate - shift
    threshold_data$lower    <- threshold_data$lower - shift
    threshold_data$upper    <- threshold_data$upper - shift

    person_theta <- person_theta - shift
  }

  # ================================================================
  # SUMMARY STATISTICS (computed after centering)
  # ================================================================
  if (robust) {
    p_center <- stats::median(person_theta)
    p_spread <- stats::mad(person_theta)
    center_label <- "Median"
    spread_label <- "MAD"
  } else {
    p_center <- mean(person_theta)
    p_spread <- stats::sd(person_theta)
    center_label <- "Mean"
    spread_label <- "SD"
  }

  thresh_locations <- threshold_data$estimate
  if (robust) {
    t_center <- stats::median(thresh_locations)
    t_spread <- stats::mad(thresh_locations)
  } else {
    t_center <- mean(thresh_locations)
    t_spread <- stats::sd(thresh_locations)
  }

  # ================================================================
  # SHARED X-AXIS LIMITS
  # ================================================================
  all_values <- c(
    person_theta,
    threshold_data$estimate,
    threshold_data$lower,
    threshold_data$upper
  )
  x_range <- range(all_values, na.rm = TRUE)
  x_pad   <- diff(x_range) * 0.05
  x_lim   <- c(x_range[1] - x_pad, x_range[2] + x_pad)

  if (center) {
    x_label <- expression("Centered latent variable" ~ (theta))
  } else {
    x_label <- expression("Latent variable" ~ (theta))
  }

  # ================================================================
  # TOP PANEL: Person histogram
  # ================================================================
  person_df <- data.frame(theta = person_theta)

  p_top <- ggplot2::ggplot(person_df, ggplot2::aes(x = .data$theta)) +
    ggplot2::geom_histogram(
      bins = bins, fill = person_fill, colour = "white", alpha = 0.85
    ) +
    ggplot2::annotate(
      "rect",
      xmin = p_center - p_spread, xmax = p_center + p_spread,
      ymin = -Inf, ymax = Inf,
      fill = person_fill, alpha = 0.12
    ) +
    ggplot2::geom_vline(
      xintercept = p_center, linewidth = 0.8,
      linetype = "dashed", colour = "grey20"
    ) +
    ggplot2::annotate(
      "text",
      x = p_center, y = Inf, vjust = -0.5,
      label = paste0(center_label, " = ", round(p_center, 2),
                     ", ", spread_label, " = ", round(p_spread, 2)),
      size = 3.2, colour = "grey20"
    ) +
    ggplot2::scale_y_continuous(
      breaks = function(lim) {
        seq(0, floor(lim[2]), by = max(1, round(lim[2] / 6)))
      }
    ) +
    ggplot2::coord_cartesian(xlim = x_lim, clip = "off") +
    ggplot2::labs(x = NULL, y = "Persons") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(5, 5, 0, 5)
    )

  # ================================================================
  # MIDDLE PANEL: Inverted threshold histogram
  # ================================================================
  thresh_hist_df <- data.frame(location = thresh_locations)

  p_mid <- ggplot2::ggplot(
    thresh_hist_df, ggplot2::aes(x = .data$location)
  ) +
    ggplot2::geom_histogram(
      bins = bins, fill = threshold_fill, colour = "white", alpha = 0.85
    ) +
    ggplot2::annotate(
      "rect",
      xmin = t_center - t_spread, xmax = t_center + t_spread,
      ymin = -Inf, ymax = Inf,
      fill = threshold_fill, alpha = 0.12
    ) +
    ggplot2::geom_vline(
      xintercept = t_center, linewidth = 0.8,
      linetype = "dashed", colour = "grey20"
    ) +
    ggplot2::annotate(
      "text",
      x = t_center, y = -Inf, vjust = 1.5,
      label = paste0(center_label, " = ", round(t_center, 2),
                     ", ", spread_label, " = ", round(t_spread, 2)),
      size = 3.2, colour = "grey20"
    ) +
    ggplot2::scale_y_reverse(
      breaks = function(lim) {
        max_val <- abs(floor(lim[1]))
        seq(0, max_val, by = max(1, round(max_val / 4)))
      }
    ) +
    ggplot2::coord_cartesian(xlim = x_lim, clip = "off") +
    ggplot2::labs(x = NULL, y = "Thresholds") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(0, 5, 0, 5)
    )

  # ================================================================
  # BOTTOM PANEL: Dot-and-whisker by item
  # ================================================================
  if (sort_items == "location") {
    item_means <- stats::aggregate(
      estimate ~ item, data = threshold_data, FUN = mean
    )
    # Reverse: easiest (lowest location) at top, hardest at bottom
    item_order <- item_means$item[order(item_means$estimate)]
  } else {
    # Reverse data order: first item in data at top of y-axis
    # ggplot places the first factor level at the bottom, so we
    # reverse so that the last level (= first item) is at the top
    item_order <- rev(unique_items_data_order)
  }
  threshold_data$item <- factor(threshold_data$item, levels = item_order)

  p_bot <- ggplot2::ggplot(
    threshold_data,
    ggplot2::aes(
      x      = .data$estimate,
      y      = .data$item,
      colour = .data$category,
      xmin   = .data$lower,
      xmax   = .data$upper
    )
  ) +
    ggplot2::geom_errorbarh(
      width = 0.25, linewidth = 0.5,
      position = ggplot2::position_dodge(width = 0.4)
    ) +
    ggplot2::geom_point(
      size = 2.5,
      position = ggplot2::position_dodge(width = 0.4)
    ) +
    ggplot2::coord_cartesian(xlim = x_lim) +
    ggplot2::labs(
      x      = x_label,
      y      = NULL,
      colour = "Threshold"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.margin     = ggplot2::margin(0, 5, 5, 5)
    )

  if (!is.null(palette)) {
    p_bot <- p_bot + ggplot2::scale_colour_manual(values = palette)
  } else {
    p_bot <- p_bot + ggplot2::scale_colour_viridis_d(end = 0.9)
  }

  # ================================================================
  # COMBINE
  # ================================================================
  combined <- p_top / p_mid / p_bot +
    patchwork::plot_layout(heights = height_ratios)

  combined
}


# ── Internal helper: extract threshold data ──────────────────────
#' @keywords internal
.extract_threshold_data <- function(draws, model, unique_items,
                                    item_name, person_name,
                                    is_ordinal, lower_prob, upper_prob) {

  result_list <- list()

  if (is_ordinal) {
    has_grouped <- any(grepl("^b_Intercept\\[.+,\\d+\\]$", names(draws)))

    for (item_label in unique_items) {
      if (has_grouped) {
        thresh_pattern <- paste0(
          "^b_Intercept\\[",
          gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
          ","
        )
        thresh_cols <- grep(thresh_pattern, names(draws), value = TRUE)
      } else {
        thresh_cols <- grep("^b_Intercept\\[\\d+\\]$", names(draws),
                            value = TRUE)
      }

      if (length(thresh_cols) == 0) next

      thresh_nums <- as.numeric(
        gsub(".*,(\\d+)\\]$|.*\\[(\\d+)\\]$", "\\1\\2", thresh_cols)
      )
      thresh_cols <- thresh_cols[order(thresh_nums)]

      for (idx in seq_along(thresh_cols)) {
        vals <- draws[[thresh_cols[idx]]]
        result_list[[length(result_list) + 1]] <- data.frame(
          item     = item_label,
          category = as.character(idx),
          estimate = mean(vals),
          lower    = stats::quantile(vals, probs = lower_prob),
          upper    = stats::quantile(vals, probs = upper_prob),
          stringsAsFactors = FALSE
        )
      }
    }
  } else {
    intercept_col <- grep("^b_Intercept$", names(draws), value = TRUE)
    if (length(intercept_col) == 0) {
      stop("Could not find intercept parameter.", call. = FALSE)
    }

    for (item_label in unique_items) {
      re_col <- paste0("r_", item_name, "[",  item_label, ",Intercept]")
      if (!re_col %in% names(draws)) {
        re_col_alt <- grep(
          paste0("^r_", item_name, "\\[",
                 gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
                 ",Intercept\\]$"),
          names(draws), value = TRUE
        )
        if (length(re_col_alt) == 1) {
          re_col <- re_col_alt
        } else {
          warning("Could not find random effect for item '",
                  item_label, "'. Skipping.", call. = FALSE)
          next
        }
      }

      item_location <- -(draws[[intercept_col]] + draws[[re_col]])

      result_list[[length(result_list) + 1]] <- data.frame(
        item     = item_label,
        category = "1",
        estimate = mean(item_location),
        lower    = stats::quantile(item_location, probs = lower_prob),
        upper    = stats::quantile(item_location, probs = upper_prob),
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, result_list)
  rownames(result) <- NULL
  result
}
