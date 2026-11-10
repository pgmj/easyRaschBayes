#' Item Category Probability Function Curves for Polytomous IRT Models
#'
#' Plots item category probability functions (ICPFs) for polytomous
#' Bayesian IRT models fitted with \pkg{brms}. For each item, the
#' probability of endorsing each response category is plotted as a
#' function of the latent variable (theta), with separate colored
#' curves per category. All items are displayed in a combined faceted
#' plot.
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
#' \code{cumulative}, \code{sratio}, and \code{cratio} families, with
#' the model's link function (for example logit, probit or cloglog).
#' The \code{acat} family is supported with the logit link only, and
#' models that estimate a discrimination parameter (\code{disc}) stop
#' with an error. For rating scale models (common thresholds with item
#' effects) each item's thresholds are shifted by its item effect.
#'
#' Posterior uncertainty in the thresholds propagates into credible
#' interval ribbons around the category probability curves.
#'
#' @references
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with
#' brms and Stan. \emph{Journal of Statistical Software}, \emph{100},
#' 1--54. \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link[brms]{posterior_epred}},
#' \code{\link[brms]{conditional_effects}},
#'
#' @examples
#' \donttest{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(ggplot2)
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
#' # Plot all items
#' plot_ipf(fit_pcm, item_var = item, person_var = id)
#'
#' # Plot a subset of items
#' plot_ipf(fit_pcm, item_var = item, person_var = id,
#'          items = c("q1", "q2", "q3"))
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

  linkinv <- .ordinal_linkinv(model, "plot_ipf")

  # Use a subsample of draws for speed
  max_draws <- min(n_draws_total, 500)
  draw_sample <- sample(seq_len(n_draws_total), max_draws)

  lower_prob <- (1 - prob) / 2
  upper_prob <- 1 - lower_prob

  # --- Compute category probabilities per item ---
  plot_data_list <- list()

  tau_all <- .item_tau_draws(draws, .item_levels(model, item_name),
                             item_name, ordinal = TRUE)

  for (item_label in unique_items) {
    thresh_sub <- tau_all[[as.character(item_label)]][draw_sample, ,
                                                      drop = FALSE]
    n_thresholds <- ncol(thresh_sub)
    n_cat <- n_thresholds + 1

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
          cum_probs <- linkinv(tau_s - theta)
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
            p_stop <- linkinv(tau_s[cc] - theta)
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
            p_cont <- linkinv(theta - tau_s[cc])
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
    ggplot2::geom_line(ggplot2::aes(linetype = .data$category),
                       linewidth = line_size) +
    ggplot2::facet_wrap(~ item, ncol = ncol) +
    ggplot2::labs(
      x        = expression(theta),
      y        = "Probability",
      colour   = "Category",
      fill     = "Category",
      linetype = "Category"
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
#' Produces a three-panel targeting plot with a shared logit x-axis:
#'
#' \enumerate{
#'   \item \strong{Top}: Histogram of person locations, with a reference
#'     line for the mean (or median) and shading for ±1 SD (or ±1 MAD).
#'   \item \strong{Middle}: Inverted histogram of item threshold
#'     locations, with the same summary annotations.
#'   \item \strong{Bottom}: One bar per item, either partitioned into
#'     response-category bands (\code{panel = "categories"}, the
#'     default) or drawn as a dot-and-whisker plot of the individual
#'     thresholds (\code{panel = "thresholds"}).
#' }
#'
#' Together, the top and middle panels form a back-to-back histogram
#' that shows whether the test is well targeted to the sample. The
#' bottom panel places the items on the same scale, so the category
#' bands show which response is the most likely one at the locations
#' where the persons are. The layout follows
#' \code{easyRasch2::RMtargeting()}.
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
#'   recentred so that the mean of all item threshold locations is
#'   zero. Person locations are shifted by the same constant. If
#'   \code{FALSE}, the raw brms parameterisation is used.
#' @param sort_items Character. How to order items in the bottom
#'   panel. \code{"data"} (the default) keeps the order in which items
#'   first appear in the model data, with the first item at the top.
#'   \code{"location"} sorts items by their mean threshold location,
#'   easiest at the top.
#' @param bins Integer. Number of bins for both histograms. Default is
#'   half the maximum possible sum score, rounded up, and at least 11.
#' @param prob Numeric in \eqn{(0, 1)}, or \code{NULL}. Width of the
#'   credible intervals for the item thresholds and of the intervals in
#'   the returned \code{categories} table. Default is 0.95. \code{NULL}
#'   hides the threshold intervals in the plot.
#' @param palette An optional character vector of colours, one per
#'   response category. Overrides the viridis colours in both panel
#'   types.
#' @param person_fill Fill colour for the person histogram. Default is
#'   \code{"#0072B2"} (blue).
#' @param threshold_fill Fill colour for the threshold histogram, and
#'   for the dichotomous dot-and-whisker panel. Default is
#'   \code{"#D55E00"} (vermillion).
#' @param height_ratios Numeric vector of length 3 giving the relative
#'   heights of the top, middle and bottom panels. Default is
#'   \code{c(3, 2, 5)}.
#' @param panel Character. \code{"categories"} (the default) draws each
#'   item as a bar partitioned into response-category bands, with the
#'   threshold estimates and their credible intervals below it.
#'   \code{"thresholds"} draws a dot-and-whisker plot of the thresholds.
#' @param person_estimate Character. \code{"eap"} (the default) uses
#'   the posterior means of the person random effects.
#'   \code{"wle"} uses Warm's weighted likelihood estimates from
#'   \code{\link{person_parameters}}. EAPs are shrunk toward the mean,
#'   so their histogram is narrower than that of the WLEs.
#' @param xlim Numeric vector of length 2. Initial limits of the shared
#'   x-axis, widened to whole numbers if any person location, threshold
#'   or interval falls outside. Default is \code{c(-4, 4)}.
#' @param category_labels Optional character vector of labels for the
#'   response categories, in ascending order and one per category. Used
#'   in the legend of the \code{"categories"} panel. Default \code{NULL}
#'   uses the category scores, starting at 0.
#' @param viridis_option Character. Viridis palette option for the
#'   category colours. Default is \code{"G"} (mako).
#' @param viridis_begin,viridis_end Numeric in \eqn{[0, 1]}. Start and end
#'   of the viridis palette. Defaults \code{0.9} and \code{0.2} run from
#'   light to dark, so higher categories are darker.
#' @param row_gap Numeric. Vertical spacing between item rows in the
#'   \code{"categories"} panel. Default \code{NULL} uses 1, widened when
#'   a category collapses or thresholds are disordered, to make room for
#'   their labels.
#' @param output Character. \code{"patchwork"} (the default) returns the
#'   combined plot. \code{"list"} returns the three panels and the
#'   posterior summaries of the category structure.
#'
#' @return With \code{output = "patchwork"}, a \code{patchwork} object.
#'   With \code{output = "list"}, a list with:
#' \describe{
#'   \item{p1, p2, p3}{The person histogram, threshold histogram and
#'     bottom panel as \code{ggplot} objects.}
#'   \item{categories}{A \code{\link[tibble]{tibble}} with one row per
#'     item and response category: \code{p_never_modal}, the posterior
#'     probability that the category is never the most likely response
#'     at any location, and \code{lower_bound}, \code{lower_bound_lo},
#'     \code{lower_bound_hi}, the posterior median and \code{prob}
#'     interval of the location where the category starts to be the most
#'     likely response, over the draws where it is (\code{NA} for
#'     category 0).}
#'   \item{disordered}{A \code{\link[tibble]{tibble}} with one row per
#'     item and adjacent threshold pair: \code{p_disordered}, the
#'     posterior probability that the higher threshold lies below the
#'     lower one.}
#' }
#'
#' @details
#' \strong{Item thresholds} come from the posterior draws (see
#' \code{\link{item_parameters}}). Points are posterior means and
#' intervals are equal-tailed credible intervals.
#'
#' \strong{Category bands.} With \code{panel = "categories"}, each band
#' spans the locations at which its response category is the most
#' likely response, computed from the posterior mean thresholds. When an
#' item's thresholds are ordered these boundaries are the thresholds
#' themselves. When they are not, the disordered run is pooled by
#' averaging, and the categories it skips over, which are never the
#' most likely response at any location, collapse to a red tick
#' labelled with the category number. Red arrows below the bar mark each
#' pair of disordered thresholds and give the distance between them in
#' logits. The two outer bands fade
#' towards the panel edge, since they are open-ended.
#'
#' \strong{Posterior probabilities.} The same category structure is
#' computed for every posterior draw. Each red tick is labelled with the
#' posterior probability that its category is never the most likely
#' response, and each pair of disordered thresholds with the posterior
#' probability that the pair is disordered. The full tables, including categories and pairs
#' without red marks, are returned with \code{output = "list"}.
#'
#' \strong{Person locations.} \code{person_estimate = "eap"} uses the
#' posterior means of the person random effects, shifted by the same
#' constant as the thresholds. \code{"wle"} uses the WLEs of
#' \code{\link{person_parameters}}, which are finite at extreme scores
#' and use the items each person answered.
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
#' \code{\link{item_parameters}}, \code{\link{person_parameters}}.
#'
#' @examples
#' \donttest{
#' library(brms)
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' library(ggplot2)
#' library(patchwork)
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
#' # Default: category bands, EAP person locations, items in data order
#' plot_targeting(fit_pcm)
#'
#' # WLE person locations, items sorted by location, labelled categories
#' plot_targeting(fit_pcm, person_estimate = "wle", sort_items = "location",
#'                category_labels = c("Not at all", "Several days",
#'                                   "More than half the days",
#'                                   "Nearly every day"))
#'
#' # Dot-and-whisker panel
#' plot_targeting(fit_pcm, panel = "thresholds")
#'
#' # Panels and posterior summaries of the category structure
#' res <- plot_targeting(fit_pcm, output = "list")
#' res$categories
#' res$disordered
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
    bins,
    prob = 0.95,
    palette = NULL,
    person_fill = "#0072B2",
    threshold_fill = "#D55E00",
    height_ratios = c(3, 2, 5),
    panel = c("categories", "thresholds"),
    person_estimate = c("eap", "wle"),
    xlim = c(-4, 4),
    category_labels = NULL,
    viridis_option = "G",
    viridis_begin = 0.9,
    viridis_end = 0.2,
    row_gap = NULL,
    output = c("patchwork", "list")
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

  sort_items      <- match.arg(sort_items)
  panel           <- match.arg(panel)
  person_estimate <- match.arg(person_estimate)
  output          <- match.arg(output)
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
  show_ci <- !is.null(prob)
  ci_prob <- if (show_ci) prob else 0.95

  # ================================================================
  # ITEM THRESHOLDS
  # ================================================================
  tau <- .model_tau_draws(model, item_name)
  thr_means <- lapply(tau, colMeans)
  shift <- if (center) mean(unlist(thr_means)) else 0
  tau <- lapply(tau, function(m) m - shift)

  # ================================================================
  # PERSON LOCATIONS
  # ================================================================
  if (person_estimate == "eap") {
    person_re <- brms::ranef(model)[[person_name]]
    if (is.null(person_re)) {
      stop("No random effects found for '", person_name,
           "'. Check your model formula.", call. = FALSE)
    }
    person_theta <- person_re[, "Estimate", "Intercept"] - shift
    person_label <- "EAP"
  } else {
    pp <- person_parameters(model, item_var = !!rlang::sym(item_name),
                            person_var = !!rlang::sym(person_name),
                            center = center)
    person_theta <- pp$person_estimates$wle
    person_label <- "WLE"
  }
  person_theta <- person_theta[is.finite(person_theta)]
  n_persons <- length(unique(model$data[[person_name]]))

  # ================================================================
  # PANELS
  # ================================================================
  built <- .targeting_build(
    tau = tau, person_theta = person_theta, person_label = person_label,
    n_persons = n_persons,
    data_order = as.character(unique(model$data[[item_name]])),
    center = center, robust = robust, sort_items = sort_items,
    bins = if (missing(bins)) NULL else bins,
    show_ci = show_ci, ci_prob = ci_prob, palette = palette,
    person_fill = person_fill, threshold_fill = threshold_fill,
    panel = panel, xlim = xlim, category_labels = category_labels,
    viridis_option = viridis_option, viridis_begin = viridis_begin,
    viridis_end = viridis_end, row_gap = row_gap
  )

  if (output == "list") return(built[c("p1", "p2", "p3", "categories",
                                       "disordered")])
  built$p1 / built$p2 / built$p3 +
    patchwork::plot_layout(heights = height_ratios)
}



# ── Internal: build the three targeting panels ───────────────────
#
# Shared by plot_targeting() and plot_targeting_hpcm(). tau is a named list
# of draws x thresholds matrices, already centred. person_theta is on the
# same scale. data_order gives the items in data order. cat_offset is the
# score of the lowest category (0, or 1 for the severity part of the hurdle
# PCM). extra_caption is appended to the caption. Returns list(p1, p2, p3,
# categories, disordered).
#' @keywords internal
#' @noRd
.targeting_build <- function(tau, person_theta, person_label, n_persons,
                             data_order, center, robust, sort_items, bins,
                             show_ci, ci_prob, palette, person_fill,
                             threshold_fill, panel, xlim, category_labels,
                             viridis_option, viridis_begin, viridis_end,
                             row_gap, cat_offset = 0L,
                             extra_caption = NULL) {
  lower_prob <- (1 - ci_prob) / 2
  item_levels <- names(tau)
  thr <- do.call(rbind, lapply(item_levels, function(it) {
    m <- tau[[it]]
    data.frame(
      Item      = it,
      k         = seq_len(ncol(m)),
      Threshold = paste0("T", seq_len(ncol(m))),
      Location  = colMeans(m),
      CI_low    = apply(m, 2, stats::quantile, probs = lower_prob),
      CI_high   = apply(m, 2, stats::quantile, probs = 1 - lower_prob),
      stringsAsFactors = FALSE
    )
  }))
  rownames(thr) <- NULL
  n_cat <- max(thr$k) + 1L
  max_score <- nrow(thr)

  # ================================================================
  # AXIS, BINS AND SUMMARIES
  # ================================================================
  if (is.null(bins)) bins <- max(11L, ceiling(max_score / 2))

  all_values <- c(person_theta, thr$Location)
  if (show_ci) all_values <- c(all_values, thr$CI_low, thr$CI_high)
  if (max(all_values, na.rm = TRUE) > xlim[2]) {
    xlim[2] <- ceiling(max(all_values, na.rm = TRUE))
  }
  if (min(all_values, na.rm = TRUE) < xlim[1]) {
    xlim[1] <- floor(min(all_values, na.rm = TRUE))
  }
  x_breaks <- seq(xlim[1], xlim[2], by = 1)
  x_label <- if (center) "Centred location (logit scale)" else
    "Location (logit scale)"

  if (robust) {
    p_center <- stats::median(person_theta)
    p_spread <- stats::mad(person_theta)
    t_center <- stats::median(thr$Location)
    t_spread <- stats::mad(thr$Location)
    center_label <- "Median"
    spread_label <- "MAD"
  } else {
    p_center <- mean(person_theta)
    p_spread <- stats::sd(person_theta)
    t_center <- mean(thr$Location)
    t_spread <- stats::sd(thr$Location)
    center_label <- "Mean"
    spread_label <- "SD"
  }
  annot <- function(cn, sp) {
    paste0(center_label, " = ", round(cn, 2), ", ", spread_label, " = ",
           round(sp, 2))
  }

  # ================================================================
  # TOP PANEL: person histogram
  # ================================================================
  p1 <- ggplot2::ggplot(data.frame(theta = person_theta),
                        ggplot2::aes(x = .data$theta)) +
    ggplot2::geom_histogram(bins = bins, fill = person_fill,
                            colour = "white", alpha = 0.85) +
    ggplot2::annotate("rect", xmin = p_center - p_spread,
                      xmax = p_center + p_spread, ymin = -Inf, ymax = Inf,
                      fill = person_fill, alpha = 0.12) +
    ggplot2::geom_vline(xintercept = p_center, linewidth = 0.8,
                        linetype = "dashed", colour = "grey20") +
    ggplot2::annotate("text", x = p_center, y = Inf, vjust = -0.5,
                      label = annot(p_center, p_spread), size = 3.2,
                      colour = "grey20") +
    ggplot2::scale_y_continuous(breaks = function(lim) {
      seq(0, floor(lim[2]), by = max(1, round(lim[2] / 6)))
    }) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::coord_cartesian(xlim = xlim, clip = "off") +
    ggplot2::labs(x = NULL, y = "Persons") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(5, 5, 0, 5)) +
    .erb_caption_theme()

  # ================================================================
  # MIDDLE PANEL: inverted threshold histogram
  # ================================================================
  p2 <- ggplot2::ggplot(data.frame(location = thr$Location),
                        ggplot2::aes(x = .data$location)) +
    ggplot2::geom_histogram(bins = bins, fill = threshold_fill,
                            colour = "white", alpha = 0.85) +
    ggplot2::annotate("rect", xmin = t_center - t_spread,
                      xmax = t_center + t_spread, ymin = -Inf, ymax = Inf,
                      fill = threshold_fill, alpha = 0.12) +
    ggplot2::geom_vline(xintercept = t_center, linewidth = 0.8,
                        linetype = "dashed", colour = "grey20") +
    ggplot2::annotate("text", x = t_center, y = -Inf, vjust = 1.5,
                      label = annot(t_center, t_spread), size = 3.2,
                      colour = "grey20") +
    ggplot2::scale_y_reverse(minor_breaks = NULL, breaks = function(lim) {
      max_val <- abs(floor(lim[1]))
      seq(0, max_val, by = max(1, round(max_val / 4)))
    }) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::coord_cartesian(xlim = xlim, clip = "off") +
    ggplot2::labs(x = NULL, y = "Thresholds") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0, 5, 0, 5)) +
    .erb_caption_theme()

  # ================================================================
  # BOTTOM PANEL
  # ================================================================
  # ggplot draws the first factor level at the bottom, so the order is
  # reversed to put the first (or easiest) item at the top
  if (sort_items == "location") {
    item_means <- tapply(thr$Location, thr$Item, mean)
    item_order <- rev(names(item_means)[order(item_means)])
  } else {
    item_order <- rev(data_order)
  }

  if (is.null(category_labels)) {
    category_labels <- as.character(seq_len(n_cat) - 1L + cat_offset)
  } else if (length(category_labels) != n_cat) {
    stop("`category_labels` must have one label per response category. ",
         "The model has ", n_cat, " categories and ",
         length(category_labels), " labels were supplied.", call. = FALSE)
  }
  fills <- if (!is.null(palette)) {
    if (length(palette) < n_cat) {
      stop("`palette` needs at least ", n_cat, " colours.", call. = FALSE)
    }
    palette[seq_len(n_cat)]
  } else {
    scales::viridis_pal(option = viridis_option, begin = viridis_begin,
                        end = viridis_end)(n_cat)
  }

  cat_post <- .category_posterior(tau, ci_prob)

  caption_body <- paste0(
    "Person location (", person_label, ") ", tolower(center_label), ": ",
    round(p_center, 2), " (", spread_label, " ", round(p_spread, 2),
    "). Item threshold location ", tolower(center_label), ": ",
    round(t_center, 2), " (", spread_label, " ", round(t_spread, 2),
    "). n = ", n_persons, " persons."
  )
  if (!is.null(extra_caption)) {
    caption_body <- paste(caption_body, extra_caption)
  }

  if (panel == "categories") {
    bands_preview <- .category_bands(thr, item_order, xlim)
    revs_preview <- .threshold_reversals(thr, item_order)
    if (is.null(row_gap)) {
      # Room for tick labels above and reversal labels below the bars
      n_rev_max <- if (is.null(revs_preview)) 0L else
        max(table(revs_preview$Item))
      # Tick labels above the bar are only needed for collapses that are
      # not explained by a reversal (those are labelled on the arrow)
      drops_above <- !is.null(bands_preview$drops) && (
        is.null(revs_preview) ||
          !all(paste(bands_preview$drops$Item, bands_preview$drops$cat) %in%
                 paste(revs_preview$Item, revs_preview$k)))
      row_gap <- 1 + 0.3 * drops_above + 0.36 * n_rev_max
    }
    band_caption <- paste(
      "Bands span the locations at which each response category is the",
      "most likely response, from the posterior mean thresholds. Outer",
      "bands fade because they are open-ended."
    )
    if (show_ci) {
      band_caption <- paste(band_caption, sprintf(paste(
        "Points and intervals below each bar are the thresholds with %d%%",
        "credible intervals, coloured by the category entered."
      ), round(ci_prob * 100)))
    }
    if (!is.null(bands_preview$drops)) {
      band_caption <- paste(
        band_caption,
        "A red tick marks a category that is never the most likely",
        "response, with its number and the posterior probability of that."
      )
    }
    if (!is.null(revs_preview)) {
      band_caption <- paste(
        band_caption,
        "Red arrows mark disordered thresholds, with the distance between",
        "them in logits and the posterior probability that the pair is",
        "disordered. When disordered thresholds make a category collapse,",
        "the arrow is labelled 'category: distance (P)' instead, with the",
        "probability that the category is never the most likely response."
      )
    }
    p3 <- .targeting_category_panel(
      thr = thr, item_order = item_order, xlim = xlim,
      category_labels = category_labels, fills = fills, show_ci = show_ci,
      row_gap = row_gap,
      caption_text = .erb_caption(paste(caption_body, band_caption)),
      x_label = x_label,
      drop_prob = cat_post$categories, rev_prob = cat_post$disordered,
      cat_offset = cat_offset
    )
  } else {
    thr$Item <- factor(thr$Item, levels = item_order)
    ci_caption <- if (show_ci) {
      sprintf("Error bars show %d%% credible intervals.",
              round(ci_prob * 100))
    } else ""
    caption_text <- .erb_caption(paste(caption_body, ci_caption))
    if (n_cat == 2L) {
      p3 <- ggplot2::ggplot(thr, ggplot2::aes(x = .data$Location,
                                              y = .data$Item)) +
        ggplot2::geom_point(size = 3, colour = threshold_fill)
      if (show_ci) {
        p3 <- p3 + ggplot2::geom_errorbar(
          ggplot2::aes(xmin = .data$CI_low, xmax = .data$CI_high),
          width = 0.25, linewidth = 0.5, colour = threshold_fill,
          orientation = "y"
        )
      }
    } else {
      thr$Threshold <- factor(thr$Threshold,
                              levels = paste0("T", seq_len(n_cat - 1L)))
      thr_cols <- .ci_contrast(fills[-1])
      shape_vals <- rep(c(16, 17, 15, 18, 8, 3, 4, 1, 2, 0),
                        length.out = n_cat - 1L)
      p3 <- ggplot2::ggplot(thr, ggplot2::aes(
        x = .data$Location, y = .data$Item,
        colour = .data$Threshold, shape = .data$Threshold
      ))
      if (show_ci) {
        p3 <- p3 + ggplot2::geom_errorbar(
          ggplot2::aes(xmin = .data$CI_low, xmax = .data$CI_high),
          width = 0.25, linewidth = 0.5,
          position = ggplot2::position_dodge(width = 0.4),
          orientation = "y"
        )
      }
      p3 <- p3 +
        ggplot2::geom_point(size = 2.5,
                            position = ggplot2::position_dodge(width = 0.4)) +
        ggplot2::scale_colour_manual(values = thr_cols, name = "Threshold") +
        ggplot2::scale_shape_manual(values = shape_vals, name = "Threshold")
    }
    p3 <- p3 +
      ggplot2::scale_x_continuous(breaks = x_breaks) +
      ggplot2::coord_cartesian(xlim = xlim) +
      ggplot2::labs(x = x_label, y = NULL, caption = caption_text) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom",
                     plot.margin = ggplot2::margin(0, 5, 5, 5)) +
      .erb_caption_theme()
  }

  cat_post$categories$category <- cat_post$categories$category + cat_offset
  list(p1 = p1, p2 = p2, p3 = p3,
       categories = cat_post$categories,
       disordered = cat_post$disordered)
}

# ── Internal helper: extract threshold data ──────────────────────
#' @keywords internal
.extract_threshold_data <- function(draws, model, unique_items,
                                    item_name, person_name,
                                    is_ordinal, lower_prob, upper_prob) {
  tau <- .item_tau_draws(draws, .item_levels(model, item_name), item_name,
                         is_ordinal)
  result_list <- list()
  for (item_label in as.character(unique_items)) {
    m <- tau[[item_label]]
    for (idx in seq_len(ncol(m))) {
      vals <- m[, idx]
      result_list[[length(result_list) + 1]] <- data.frame(
        item     = item_label,
        category = as.character(idx),
        estimate = mean(vals),
        lower    = unname(stats::quantile(vals, probs = lower_prob)),
        upper    = unname(stats::quantile(vals, probs = upper_prob)),
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, result_list)
  rownames(result) <- NULL
  result
}
