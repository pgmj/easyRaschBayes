#' Item Characteristic Curves with Class Intervals
#'
#' Plots Item Characteristic Curves (ICCs) for Bayesian Rasch-family
#' models fitted with \pkg{brms}. Each item panel shows the
#' model-expected item score curve (with a credible interval ribbon)
#' overlaid with observed average item scores computed within class
#' intervals. Optionally, observed scores can be split by a grouping
#' variable to visually assess differential item functioning (DIF).
#'
#' This is the Bayesian analogue of \code{ICCplot()} from the
#' \pkg{iarm} package, using the class interval method.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from
#'   an ordinal IRT model (e.g., \code{family = acat}) or a
#'   dichotomous model (\code{family = bernoulli()}).
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the
#'   person grouping variable in the model data. Default is
#'   \code{id}.
#' @param items An optional character vector of item names to plot.
#'   If \code{NULL} (the default), all items are plotted.
#' @param n_intervals Integer. The number of class intervals into
#'   which persons are binned along the sum score. Default is 5.
#' @param theta_range A numeric vector of length 2 specifying the
#'   range of theta for the expected curve. Default is
#'   \code{c(-4, 4)}.
#' @param n_points Integer. Number of evenly spaced theta values
#'   for computing the expected curve. Default is 200.
#' @param center Logical. If \code{TRUE} (the default), the scale
#'   is recentered so the grand mean of item thresholds = 0,
#'   consistent with \code{\link{plot_targeting}}.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the credible
#'   interval ribbon around the expected curve. Default is 0.95.
#' @param ncol Integer. Number of columns in the faceted layout.
#'   If \code{NULL}, chosen automatically.
#' @param line_size Numeric. Line width for the expected curve.
#'   Default is 0.8.
#' @param ribbon_alpha Numeric in \eqn{[0, 1]}. Transparency of
#'   the credible interval ribbon. Default is 0.3.
#' @param point_size Numeric. Size of observed score points.
#'   Default is 2.5.
#' @param dif_var An optional unquoted variable name for a grouping
#'   variable to assess DIF visually. If supplied, observed scores
#'   are computed separately per group and coloured accordingly.
#' @param dif_data An optional data frame containing the DIF
#'   variable. Required when \code{dif_var} is specified and the
#'   variable was not part of the model formula (since brms drops
#'   unused variables from \code{model$data}). Must have the same
#'   rows and row order as the original model data.
#' @param dif_labels An optional character vector of labels for
#'   the DIF groups. If \code{NULL}, factor levels are used.
#' @param dif_stats Logical. If \code{TRUE} (the default when
#'   \code{dif_var} is specified), the partial gamma coefficient
#'   for each item is annotated in the plot panel. The partial
#'   gamma measures the strength of DIF, stratified by total score.
#' @param min_n Integer. Minimum number of persons required in a
#'   class interval for the observed mean to be plotted. Intervals
#'   with fewer persons are dropped to avoid unstable estimates.
#'   Default is 5.
#' @param palette An optional character vector of colors. If
#'   \code{NULL}, viridis colors are used.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @details
#' \strong{Expected curve:} For each item, the expected item score
#' \eqn{E_i(\theta)} is computed for a dense grid of theta values
#' using the posterior draws of item thresholds. For the adjacent
#' category (PCM/acat) family:
#' \deqn{E_i(\theta) = \sum_{c=0}^{K} c \cdot P_{ic}(\theta)}
#' where \eqn{P_{ic}} are the category probabilities. For
#' dichotomous items, \eqn{E_i(\theta) = P_i(\theta)}.
#'
#' The posterior mean of \eqn{E_i(\theta)} is plotted as a line,
#' with a credible interval ribbon showing the \code{prob} interval
#' across posterior draws.
#'
#' \strong{Class intervals:} Following the approach in
#' \code{iarm::ICCplot()}, persons are binned into class intervals
#' by their ordinal \emph{sum score} (the sufficient statistic in
#' Rasch models), not by their EAP theta estimate. Within each
#' interval, the mean observed item response is computed and
#' plotted at the theta value corresponding to the mean sum score
#' in that interval, obtained by inverting the posterior-mean
#' expected total score function. Persons with extreme sum scores
#' (0 and maximum) are excluded from the class intervals, as their
#' theta positions cannot be estimated.
#'
#' \strong{Uncertainty around observed points:} For each class
#' interval, the standard error of the mean observed response is
#' computed and displayed as error bars showing \eqn{\pm 1.96} SE.
#' These reflect sampling variability of the observed mean within
#' each bin.
#'
#' \strong{DIF overlay:} When \code{dif_var} is specified, observed
#' scores are computed separately for each level of the DIF
#' variable, producing group-specific points connected by lines.
#' Deviations between groups that track alongside each other (but
#' offset from the expected curve) suggest uniform DIF. Crossing
#' group lines suggest non-uniform DIF.
#'
#' \strong{Partial gamma:} When \code{dif_stats = TRUE} (default
#' when DIF is requested), the Goodman-Kruskal partial gamma
#' coefficient is displayed in each item panel. This measures the
#' ordinal association between item response and DIF group,
#' stratified by total score. Values near 0 indicate no DIF;
#' values near \eqn{\pm 1} indicate strong DIF. The computation
#' reuses the concordant/discordant pair counting algorithm from
#' \code{\link{gk_gamma}}.
#'
#' @seealso
#' \code{\link{plot_ipf}} for category probability curves,
#' \code{\link{plot_targeting}} for person-item maps,
#' \code{\link{dif_statistic}} for formal Bayesian DIF testing,
#' \code{\link{gk_gamma}} for the underlying gamma algorithm.
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
#' df_pcm <- eRm::pcmdat2 %>%
#'   mutate(across(everything(), ~ .x + 1)) %>%
#'   rownames_to_column("id") %>%
#'   pivot_longer(!id, names_to = "item", values_to = "response")
#'
#' fit_pcm <- brm(
#'   response | thres(gr = item) ~ 1 + (1 | id),
#'   data   = df_pcm,
#'   family = acat,
#'   chains = 4, cores = 2, iter = 1000 # use more iter and cores
#' )
#'
#' # Basic ICC plot with 5 class intervals
#' plot_icc(fit_pcm)
#'
#' # Select specific items, more intervals
#' plot_icc(fit_pcm, items = c("I1", "I2"), n_intervals = 5)
#'
#' # With DIF overlay and partial gamma annotation
#' 
#' # same dataset, adding a `gender` DIF variable
#' df_pcm <- eRm::pcmdat2 %>%
#'   mutate(across(everything(), ~ .x + 1)) %>%
#'   rownames_to_column("id") %>%
#'   mutate(gender = sample(c("M", "F"), nrow(.), TRUE)) %>% 
#'   pivot_longer(!c(id,gender), names_to = "item", values_to = "response")
#'   
#' plot_icc(fit_pcm, dif_var = gender, dif_data = df_pcm,
#'          dif_labels = c("Female", "Male"), n_intervals = 5)
#' }
#'
#' @importFrom brms as_draws_df ranef
#' @importFrom rlang enquo as_name .data quo_is_null
#' @importFrom stats formula quantile plogis family sd aggregate
#'   as.formula complete.cases
#' @importFrom tibble as_tibble
#' @export
plot_icc <- function(
    model,
    item_var     = item,
    person_var   = id,
    items        = NULL,
    n_intervals  = 5,
    theta_range  = c(-4, 4),
    n_points     = 200,
    center       = TRUE,
    prob         = 0.95,
    ncol         = NULL,
    line_size    = 0.8,
    ribbon_alpha = 0.3,
    point_size   = 2.5,
    dif_var      = NULL,
    dif_data     = NULL,
    dif_labels   = NULL,
    dif_stats    = TRUE,
    min_n        = 5,
    palette      = NULL
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }
  
  item_name   <- rlang::as_name(rlang::enquo(item_var))
  person_name <- rlang::as_name(rlang::enquo(person_var))
  
  # --- DIF variable handling ---
  dif_quo  <- rlang::enquo(dif_var)
  has_dif  <- !rlang::quo_is_null(dif_quo)
  dif_name <- if (has_dif) rlang::as_name(dif_quo) else NULL
  
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
  
  model_data   <- model$data
  unique_items <- sort(unique(model_data[[item_name]]))
  if (!is.null(items)) {
    invalid <- setdiff(items, unique_items)
    if (length(invalid) > 0) {
      stop("Items not found in model data: ",
           paste(invalid, collapse = ", "), call. = FALSE)
    }
    unique_items <- unique_items[unique_items %in% items]
  }
  all_unique_items <- sort(unique(model_data[[item_name]]))
  
  # ================================================================
  # EXTRACT POSTERIOR DRAWS & DETECT MODEL TYPE
  # ================================================================
  draws       <- tibble::as_tibble(brms::as_draws_df(model))
  family_name <- stats::family(model)$family
  is_ordinal  <- grepl("acat|cumul|sratio|cratio",
                       family_name, ignore.case = TRUE)
  is_acat     <- grepl("acat", family_name, ignore.case = TRUE)
  is_cumul    <- grepl("cumul", family_name, ignore.case = TRUE)
  is_sratio   <- grepl("sratio", family_name, ignore.case = TRUE)
  is_cratio   <- grepl("cratio", family_name, ignore.case = TRUE)
  is_binary   <- grepl("bernoulli|binomial",
                       family_name, ignore.case = TRUE)
  
  if (!is_ordinal && !is_binary) {
    stop("Unsupported model family '", family_name,
         "'. Expected acat, cumulative, sratio, cratio, or ",
         "bernoulli.", call. = FALSE)
  }
  
  # ================================================================
  # CENTERING SHIFT
  # ================================================================
  lower_prob <- (1 - prob) / 2
  upper_prob <- 1 - lower_prob
  
  if (center) {
    if (is_ordinal) {
      all_thresh_means <- c()
      for (il in all_unique_items) {
        tc <- .find_thresh_cols(draws, il)
        if (length(tc) > 0) {
          all_thresh_means <- c(
            all_thresh_means,
            vapply(tc, function(col) mean(draws[[col]]), numeric(1))
          )
        }
      }
      shift <- mean(all_thresh_means)
    } else {
      item_locs <- vapply(all_unique_items, function(il) {
        loc_draws <- .get_item_location_draws_all(
          draws, il, item_name
        )
        mean(loc_draws)
      }, numeric(1))
      shift <- mean(item_locs)
    }
  } else {
    shift <- 0
  }
  
  # ================================================================
  # RESPONSE ADJUSTMENT FOR ORDINAL MODELS
  # ================================================================
  min_resp <- min(model_data[[resp_var]], na.rm = TRUE)
  
  # ================================================================
  # COMPUTE EXPECTED TOTAL SCORE CURVE E(R|theta) FOR CLASS INTERVALS
  # ================================================================
  theta_grid_fine <- seq(theta_range[1], theta_range[2],
                         length.out = 1000)
  
  e_total_score <- rep(0, length(theta_grid_fine))
  
  for (item_label in all_unique_items) {
    if (is_ordinal) {
      thresh_cols <- .find_thresh_cols(draws, item_label)
      if (length(thresh_cols) == 0) next
      tau_mean <- vapply(
        thresh_cols, function(col) mean(draws[[col]]), numeric(1)
      )
      tau_mean <- tau_mean - shift
      n_cat <- length(tau_mean) + 1
      cats  <- seq(0, n_cat - 1)
      
      for (t in seq_along(theta_grid_fine)) {
        p_cat <- .compute_cat_probs(
          theta_grid_fine[t], tau_mean,
          is_acat, is_cumul, is_sratio, is_cratio
        )
        e_total_score[t] <- e_total_score[t] + sum(cats * p_cat)
      }
    } else {
      loc_mean <- mean(.get_item_location_draws_all(
        draws, item_label, item_name
      ))
      loc_mean <- loc_mean - shift
      e_total_score <- e_total_score +
        stats::plogis(theta_grid_fine - loc_mean)
    }
  }
  
  # ================================================================
  # COMPUTE SUM SCORES AND BUILD PERSON-LEVEL DATA
  # ================================================================
  sum_scores_df <- stats::aggregate(
    stats::as.formula(paste(resp_var, "~", person_name)),
    data = model_data,
    FUN  = sum, na.rm = TRUE
  )
  colnames(sum_scores_df) <- c("person_id", "sum_score")
  
  if (min_resp > 0) {
    n_obs_per <- stats::aggregate(
      stats::as.formula(paste(resp_var, "~", person_name)),
      data = model_data, FUN = length
    )
    colnames(n_obs_per) <- c("person_id", "n_obs")
    sum_scores_df <- merge(sum_scores_df, n_obs_per, by = "person_id")
    sum_scores_df$sum_score <- sum_scores_df$sum_score -
      (min_resp * sum_scores_df$n_obs)
    sum_scores_df$n_obs <- NULL
  }
  
  if (is_ordinal) {
    n_thresholds <- vapply(all_unique_items, function(il) {
      length(.find_thresh_cols(draws, il))
    }, integer(1))
    max_score <- sum(n_thresholds)
  } else {
    max_score <- length(all_unique_items)
  }
  
  # ================================================================
  # CLASS INTERVALS (following iarm approach)
  # ================================================================
  person_ci <- sum_scores_df[
    sum_scores_df$sum_score > 0 & sum_scores_df$sum_score < max_score,
  ]
  
  if (nrow(person_ci) == 0) {
    stop("No persons with non-extreme sum scores. ",
         "Cannot create class intervals.", call. = FALSE)
  }
  
  breaks <- stats::quantile(
    person_ci$sum_score,
    probs = seq(0, 1, length.out = n_intervals + 1)
  )
  breaks <- unique(breaks)
  n_actual_intervals <- length(breaks) - 1
  
  if (n_actual_intervals < 2) {
    warning("Could not create enough distinct class intervals. ",
            "Try reducing n_intervals.", call. = FALSE)
  }
  
  person_ci$interval <- cut(
    person_ci$sum_score, breaks = breaks,
    include.lowest = TRUE, labels = FALSE
  )
  
  interval_mean_score <- stats::aggregate(
    sum_score ~ interval, data = person_ci, FUN = mean
  )
  
  interval_mean_score$theta <- vapply(
    interval_mean_score$sum_score,
    function(s) {
      idx <- which.min(abs(e_total_score - s))
      theta_grid_fine[idx]
    },
    numeric(1)
  )
  
  person_ci <- merge(
    person_ci,
    interval_mean_score[, c("interval", "theta")],
    by = "interval"
  )
  
  # ================================================================
  # DIF GROUPING VARIABLE
  # ================================================================
  if (has_dif) {
    if (!is.null(dif_data)) {
      dif_person <- data.frame(
        person_id  = as.character(dif_data[[person_name]]),
        .dif_group = as.character(dif_data[[dif_name]]),
        stringsAsFactors = FALSE
      )
      dif_person <- dif_person[!duplicated(dif_person$person_id), ]
    } else if (dif_name %in% names(model_data)) {
      dif_person <- data.frame(
        person_id  = as.character(model_data[[person_name]]),
        .dif_group = as.character(model_data[[dif_name]]),
        stringsAsFactors = FALSE
      )
      dif_person <- dif_person[!duplicated(dif_person$person_id), ]
    } else {
      stop("DIF variable '", dif_name, "' not found. ",
           "Supply 'dif_data' containing the variable.",
           call. = FALSE)
    }
    
    person_ci <- merge(person_ci, dif_person,
                       by = "person_id", all.x = TRUE)
    
    if (!is.null(dif_labels)) {
      old_levels <- sort(unique(person_ci$.dif_group))
      if (length(dif_labels) != length(old_levels)) {
        warning("Length of dif_labels (", length(dif_labels),
                ") doesn't match number of DIF groups (",
                length(old_levels), "). Using original labels.",
                call. = FALSE)
      } else {
        person_ci$.dif_group <- factor(
          person_ci$.dif_group,
          levels = old_levels, labels = dif_labels
        )
      }
    }
  }
  
  # ================================================================
  # COMPUTE PARTIAL GAMMA FOR DIF (if requested)
  # ================================================================
  gamma_annotations <- NULL
  
  if (has_dif && dif_stats) {
    wide_responses <- .build_wide_responses(
      model_data, resp_var, item_name, person_name,
      all_unique_items, min_resp, is_ordinal
    )
    
    wide_responses <- merge(
      wide_responses, dif_person, by = "person_id", all.x = TRUE
    )
    wide_responses <- wide_responses[
      stats::complete.cases(
        wide_responses[, c(".dif_group", "sum_score")]
      ),
    ]
    
    gamma_df <- data.frame(
      item  = character(),
      gamma = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (item_label in unique_items) {
      if (!item_label %in% names(wide_responses)) next
      
      item_resp  <- wide_responses[[item_label]]
      dif_group  <- as.numeric(as.factor(wide_responses$.dif_group))
      tot_score  <- wide_responses$sum_score
      
      ok <- stats::complete.cases(item_resp, dif_group, tot_score)
      if (sum(ok) < 10) {
        gamma_val <- NA_real_
      } else {
        gamma_val <- .partial_gk_gamma(
          item_resp[ok], dif_group[ok], tot_score[ok]
        )
      }
      
      gamma_df <- rbind(gamma_df, data.frame(
        item  = item_label,
        gamma = gamma_val,
        stringsAsFactors = FALSE
      ))
    }
    
    gamma_annotations <- gamma_df
  }
  
  # ================================================================
  # COMPUTE EXPECTED ITEM SCORE CURVES PER ITEM
  # ================================================================
  n_draws_total <- nrow(draws)
  max_draws     <- min(n_draws_total, 500)
  draw_sample   <- sample(seq_len(n_draws_total), max_draws)
  
  theta_grid <- seq(theta_range[1], theta_range[2],
                    length.out = n_points)
  
  expected_data_list <- list()
  
  for (item_label in unique_items) {
    if (is_ordinal) {
      thresh_cols <- .find_thresh_cols(draws, item_label)
      if (length(thresh_cols) == 0) {
        warning("Could not find threshold parameters for item '",
                item_label, "'. Skipping.", call. = FALSE)
        next
      }
      
      thresh_sub <- as.matrix(
        draws[draw_sample, thresh_cols, drop = FALSE]
      )
      thresh_sub <- thresh_sub - shift
      
      n_thresh <- ncol(thresh_sub)
      n_cat    <- n_thresh + 1
      cats     <- seq(0, n_cat - 1)
      
      escore_draws <- matrix(NA_real_, nrow = max_draws,
                             ncol = n_points)
      
      for (t in seq_along(theta_grid)) {
        theta <- theta_grid[t]
        for (s in seq_len(max_draws)) {
          tau_s <- thresh_sub[s, ]
          p_cat <- .compute_cat_probs(
            theta, tau_s, is_acat, is_cumul, is_sratio, is_cratio
          )
          escore_draws[s, t] <- sum(cats * p_cat)
        }
      }
      
    } else {
      item_loc_draws <- .get_item_location_draws(
        draws, item_label, item_name, draw_sample
      )
      item_loc_draws <- item_loc_draws - shift
      
      escore_draws <- matrix(NA_real_, nrow = max_draws,
                             ncol = n_points)
      for (t in seq_along(theta_grid)) {
        escore_draws[, t] <- stats::plogis(
          theta_grid[t] - item_loc_draws
        )
      }
    }
    
    e_mean  <- colMeans(escore_draws)
    e_lower <- apply(escore_draws, 2, stats::quantile,
                     probs = lower_prob)
    e_upper <- apply(escore_draws, 2, stats::quantile,
                     probs = upper_prob)
    
    expected_data_list[[length(expected_data_list) + 1]] <-
      data.frame(
        item    = item_label,
        theta   = theta_grid,
        e_score = e_mean,
        e_lower = e_lower,
        e_upper = e_upper,
        stringsAsFactors = FALSE
      )
  }
  
  expected_data <- do.call(rbind, expected_data_list)
  
  # ================================================================
  # COMPUTE OBSERVED SCORES BY CLASS INTERVAL PER ITEM
  # ================================================================
  observed_data_list <- list()
  
  for (item_label in unique_items) {
    item_rows      <- which(model_data[[item_name]] == item_label)
    item_responses <- as.numeric(model_data[[resp_var]][item_rows])
    item_persons   <- as.character(
      model_data[[person_name]][item_rows]
    )
    
    if (is_ordinal) {
      item_responses <- item_responses - min_resp
    }
    
    obs_df <- data.frame(
      person_id = item_persons,
      response  = item_responses,
      stringsAsFactors = FALSE
    )
    
    obs_df <- merge(obs_df, person_ci, by = "person_id",
                    all.x = FALSE)
    
    if (nrow(obs_df) == 0) next
    
    if (has_dif) {
      obs_agg <- stats::aggregate(
        response ~ interval + theta + .dif_group,
        data = obs_df, FUN = mean, na.rm = TRUE
      )
      obs_n <- stats::aggregate(
        response ~ interval + theta + .dif_group,
        data = obs_df, FUN = length
      )
      obs_se <- stats::aggregate(
        response ~ interval + theta + .dif_group,
        data = obs_df,
        FUN = function(x) {
          stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
        }
      )
      names(obs_n)[names(obs_n) == "response"]  <- "n"
      names(obs_se)[names(obs_se) == "response"] <- "se"
      obs_agg <- merge(obs_agg, obs_n,
                       by = c("interval", "theta", ".dif_group"))
      obs_agg <- merge(obs_agg, obs_se,
                       by = c("interval", "theta", ".dif_group"))
    } else {
      obs_agg <- stats::aggregate(
        response ~ interval + theta,
        data = obs_df, FUN = mean, na.rm = TRUE
      )
      obs_n <- stats::aggregate(
        response ~ interval + theta,
        data = obs_df, FUN = length
      )
      obs_se <- stats::aggregate(
        response ~ interval + theta,
        data = obs_df,
        FUN = function(x) {
          stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
        }
      )
      names(obs_n)[names(obs_n) == "response"]  <- "n"
      names(obs_se)[names(obs_se) == "response"] <- "se"
      obs_agg <- merge(obs_agg, obs_n, by = c("interval", "theta"))
      obs_agg <- merge(obs_agg, obs_se, by = c("interval", "theta"))
    }
    
    obs_agg <- obs_agg[obs_agg$n >= min_n, ]
    obs_agg$item <- item_label
    
    observed_data_list[[length(observed_data_list) + 1]] <- obs_agg
  }
  
  observed_data <- do.call(rbind, observed_data_list)
  
  # ================================================================
  # BUILD PLOT
  # ================================================================
  expected_data$item <- factor(expected_data$item,
                               levels = unique_items)
  observed_data$item <- factor(observed_data$item,
                               levels = unique_items)
  
  n_items_plot <- length(unique_items)
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_items_plot))
  }
  
  p <- ggplot2::ggplot()
  
  # Expected score curve with ribbon
  if (ribbon_alpha > 0) {
    p <- p + ggplot2::geom_ribbon(
      data    = expected_data,
      mapping = ggplot2::aes(
        x = .data$theta, ymin = .data$e_lower, ymax = .data$e_upper
      ),
      fill = "grey70", alpha = ribbon_alpha
    )
  }
  
  p <- p + ggplot2::geom_line(
    data    = expected_data,
    mapping = ggplot2::aes(x = .data$theta, y = .data$e_score),
    colour  = "grey40", linewidth = line_size
  )
  
  # Observed score points with error bars
  if (has_dif) {
    p <- p +
      ggplot2::geom_errorbar(
        data    = observed_data,
        mapping = ggplot2::aes(
          x      = .data$theta,
          ymin   = pmax(.data$response - (1.96 * .data$se), 0),
          ymax   = .data$response + (1.96 * .data$se),
          colour = .data$.dif_group
        ),
        width = 0.05, linewidth = 0.4
      ) +
      ggplot2::geom_line(
        data    = observed_data,
        mapping = ggplot2::aes(
          x      = .data$theta,
          y      = .data$response,
          colour = .data$.dif_group
        ),
        linewidth = 0.5
      ) +
      ggplot2::geom_point(
        data    = observed_data,
        mapping = ggplot2::aes(
          x      = .data$theta,
          y      = .data$response,
          colour = .data$.dif_group
        ),
        size = point_size
      )
    
    if (!is.null(palette)) {
      p <- p +
        ggplot2::scale_colour_manual(values = palette, name = "")
    } else {
      p <- p +
        ggplot2::scale_colour_viridis_d(end = 0.8, name = "")
    }
  } else {
    p <- p +
      ggplot2::geom_errorbar(
        data    = observed_data,
        mapping = ggplot2::aes(
          x    = .data$theta,
          ymin = pmax(.data$response - (1.96 * .data$se), 0),
          ymax = .data$response + (1.96 * .data$se)
        ),
        colour = "#1F78B4", width = 0.05, linewidth = 0.4
      ) +
      ggplot2::geom_point(
        data    = observed_data,
        mapping = ggplot2::aes(
          x = .data$theta, y = .data$response
        ),
        colour = "#1F78B4", size = point_size
      )
  }
  
  # Add partial gamma annotation if DIF is active
  if (has_dif && dif_stats && !is.null(gamma_annotations)) {
    gamma_annotations$item <- factor(gamma_annotations$item,
                                     levels = unique_items)
    gamma_annotations$label <- ifelse(
      is.na(gamma_annotations$gamma),
      "gamma == NA",
      paste0("gamma == ",
             format(round(gamma_annotations$gamma, 2), nsmall = 2))
    )
    
    p <- p +
      ggplot2::geom_text(
        data    = gamma_annotations,
        mapping = ggplot2::aes(
          x     = -Inf,
          y     = Inf,
          label = .data$label
        ),
        hjust   = -0.1,
        vjust   = 1.5,
        size    = 3.5,
        colour  = "black",
        parse   = TRUE,
        inherit.aes = FALSE
      )
  }
  
  p <- p +
    ggplot2::facet_wrap(~ item, ncol = ncol, scales = "free_y") +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::labs(
      x = expression(theta),
      y = "Item score"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold")
    )
  
  p
}


# ══════════════════════════════════════════════════════════════════
# INTERNAL HELPERS
# ══════════════════════════════════════════════════════════════════

# ── Internal: find threshold columns for an item ────────────────
#' @keywords internal
.find_thresh_cols <- function(draws, item_label) {
  thresh_pattern <- paste0(
    "^b_Intercept\\[",
    gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
    ","
  )
  thresh_cols <- grep(thresh_pattern, names(draws), value = TRUE)
  
  if (length(thresh_cols) == 0) {
    thresh_cols <- grep("^b_Intercept\\[\\d+\\]$", names(draws),
                        value = TRUE)
  }
  
  if (length(thresh_cols) == 0) return(character(0))
  
  thresh_nums <- as.numeric(
    gsub(".*,(\\d+)\\]$|.*\\[(\\d+)\\]$", "\\1\\2", thresh_cols)
  )
  thresh_cols[order(thresh_nums)]
}


# ── Internal: compute category probabilities ────────────────────
#' @keywords internal
.compute_cat_probs <- function(theta, tau,
                               is_acat, is_cumul,
                               is_sratio, is_cratio) {
  n_cat <- length(tau) + 1
  
  if (is_acat) {
    cum_eta   <- c(0, cumsum(theta - tau))
    max_eta   <- max(cum_eta)
    log_denom <- max_eta + log(sum(exp(cum_eta - max_eta)))
    probs     <- exp(cum_eta - log_denom)
    
  } else if (is_cumul) {
    cum_probs <- stats::plogis(tau - theta)
    probs     <- diff(c(0, cum_probs, 1))
    probs     <- pmax(probs, 0)
    probs     <- probs / sum(probs)
    
  } else if (is_sratio) {
    probs     <- numeric(n_cat)
    remaining <- 1
    for (cc in seq_along(tau)) {
      p_stop       <- stats::plogis(tau[cc] - theta)
      probs[cc]    <- remaining * p_stop
      remaining    <- remaining * (1 - p_stop)
    }
    probs[n_cat] <- remaining
    
  } else if (is_cratio) {
    probs     <- numeric(n_cat)
    remaining <- 1
    for (cc in seq_along(tau)) {
      p_cont       <- stats::plogis(theta - tau[cc])
      probs[cc]    <- remaining * (1 - p_cont)
      remaining    <- remaining * p_cont
    }
    probs[n_cat] <- remaining
  }
  
  probs
}


# ── Internal: get dichotomous item location draws (subsampled) ──
#' @keywords internal
.get_item_location_draws <- function(draws, item_label, item_name,
                                     draw_sample) {
  col_pattern <- paste0("^b_", item_name,
                        gsub("([.|()\\^{}+$*?])", "\\\\\\1",
                             item_label), "$")
  fe_col <- grep(col_pattern, names(draws), value = TRUE)
  
  if (length(fe_col) > 0) {
    return(-as.numeric(draws[[fe_col[1]]][draw_sample]))
  }
  
  intercept_col <- grep("^b_Intercept$", names(draws), value = TRUE)
  intercept_draws <- if (length(intercept_col) == 1) {
    as.numeric(draws[[intercept_col]][draw_sample])
  } else {
    rep(0, length(draw_sample))
  }
  
  re_col <- paste0("r_", item_name, "[", item_label, ",Intercept]")
  if (!re_col %in% names(draws)) {
    re_col_alt <- grep(
      paste0("^r_", item_name, "\\[",
             gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
             ",Intercept\\]$"),
      names(draws), value = TRUE
    )
    if (length(re_col_alt) == 1) re_col <- re_col_alt
    else stop("Could not find parameters for item '", item_label,
              "'.", call. = FALSE)
  }
  
  -(intercept_draws + as.numeric(draws[[re_col]][draw_sample]))
}


# ── Internal: get ALL dichotomous item location draws ────────────
#' @keywords internal
.get_item_location_draws_all <- function(draws, item_label,
                                         item_name) {
  col_pattern <- paste0("^b_", item_name,
                        gsub("([.|()\\^{}+$*?])", "\\\\\\1",
                             item_label), "$")
  fe_col <- grep(col_pattern, names(draws), value = TRUE)
  
  if (length(fe_col) > 0) {
    return(-as.numeric(draws[[fe_col[1]]]))
  }
  
  intercept_col <- grep("^b_Intercept$", names(draws), value = TRUE)
  intercept_draws <- if (length(intercept_col) == 1) {
    as.numeric(draws[[intercept_col]])
  } else {
    rep(0, nrow(draws))
  }
  
  re_col <- paste0("r_", item_name, "[", item_label, ",Intercept]")
  if (!re_col %in% names(draws)) {
    re_col_alt <- grep(
      paste0("^r_", item_name, "\\[",
             gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
             ",Intercept\\]$"),
      names(draws), value = TRUE
    )
    if (length(re_col_alt) == 1) re_col <- re_col_alt
    else stop("Could not find parameters for item '", item_label,
              "'.", call. = FALSE)
  }
  
  -(intercept_draws + as.numeric(draws[[re_col]]))
}


# ── Internal: build wide-format response matrix per person ──────
#' @keywords internal
.build_wide_responses <- function(model_data, resp_var, item_name,
                                  person_name, all_items,
                                  min_resp, is_ordinal) {
  wide <- stats::reshape(
    model_data[, c(person_name, item_name, resp_var)],
    idvar     = person_name,
    timevar   = item_name,
    v.names   = resp_var,
    direction = "wide",
    sep       = "___"
  )
  
  colnames(wide) <- gsub(
    paste0("^", resp_var, "___"), "", colnames(wide)
  )
  colnames(wide)[1] <- "person_id"
  
  item_cols <- intersect(all_items, colnames(wide))
  if (is_ordinal && min_resp > 0) {
    for (col in item_cols) {
      wide[[col]] <- wide[[col]] - min_resp
    }
  }
  
  wide$sum_score <- rowSums(
    wide[, item_cols, drop = FALSE], na.rm = TRUE
  )
  
  wide
}


# ── Internal: Partial Goodman-Kruskal gamma ─────────────────────
# Computes partial gamma between x and group, stratified by score.
# Reuses the concordant/discordant pair counting algorithm from
# gk_gamma_helpers.R (.gk_gamma_from_tab engine), accumulating
# C and D across strata before computing the ratio. This is the
# same statistic as iarm::partgam_DIF().
#' @keywords internal
.partial_gk_gamma <- function(x, group, score) {
  C_total <- 0
  D_total <- 0
  
  unique_scores <- sort(unique(score))
  
  for (s in unique_scores) {
    idx <- which(score == s)
    if (length(idx) < 2L) next
    
    x_s <- x[idx]
    g_s <- group[idx]
    
    # Build item-response x DIF-group contingency table
    ux <- sort.int(unique.default(x_s), method = "quick")
    ug <- sort.int(unique.default(g_s), method = "quick")
    nr <- length(ux)
    nc <- length(ug)
    if (nr < 2L || nc < 2L) next
    
    xi <- match(x_s, ux)
    gi <- match(g_s, ug)
    
    tab <- matrix(0L, nrow = nr, ncol = nc)
    for (k in seq_along(xi)) {
      tab[xi[k], gi[k]] <- tab[xi[k], gi[k]] + 1L
    }
    
    # Count concordant/discordant pairs for this stratum
    # (same suffix/prefix sum algorithm as .gk_gamma_from_tab)
    suffix <- matrix(0, nrow = nr + 1L, ncol = nc + 1L)
    for (i in seq.int(nr, 1L)) {
      for (j in seq.int(nc, 1L)) {
        suffix[i, j] <- tab[i, j] +
          suffix[i + 1L, j] +
          suffix[i, j + 1L] -
          suffix[i + 1L, j + 1L]
      }
    }
    
    pcol <- matrix(0, nrow = nr + 1L, ncol = nc + 1L)
    for (i in seq.int(nr, 1L)) {
      for (j in seq.int(2L, nc + 1L)) {
        pcol[i, j] <- tab[i, j - 1L] +
          pcol[i + 1L, j] +
          pcol[i, j - 1L] -
          pcol[i + 1L, j - 1L]
      }
    }
    
    for (i in seq_len(nr - 1L)) {
      for (j in seq_len(nc)) {
        n_ij <- tab[i, j]
        if (n_ij == 0L) next
        if (j < nc) C_total <- C_total + n_ij * suffix[i + 1L, j + 1L]
        if (j > 1L) D_total <- D_total + n_ij * pcol[i + 1L, j]
      }
    }
  }
  
  denom <- C_total + D_total
  if (denom == 0) return(NA_real_)
  (C_total - D_total) / denom
}