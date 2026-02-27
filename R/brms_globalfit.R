#' Residual PCA Contrast Plot for Bayesian IRT Models
#'
#' Assesses dimensionality of Bayesian IRT models by performing a
#' principal component analysis (PCA) of standardized residuals for
#' each posterior draw. The item loadings on the first residual
#' contrast are plotted against item locations, with posterior
#' uncertainty displayed as 2D density contours, crosshairs, or
#' both. A posterior predictive p-value for the first-contrast
#' eigenvalue tests whether the observed residual structure exceeds
#' what the model predicts under unidimensionality.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object from an
#'   ordinal IRT model (e.g., \code{family = acat}) or a dichotomous
#'   model (\code{family = bernoulli()}).
#' @param item_var An unquoted variable name identifying the item
#'   grouping variable in the model data. Default is \code{item}.
#' @param person_var An unquoted variable name identifying the person
#'   grouping variable in the model data. Default is \code{id}.
#' @param center Logical. If \code{TRUE} (the default), item
#'   locations are mean-centered to zero, matching the convention
#'   used in \code{\link{plot_pim}}.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the credible
#'   intervals for both loading and location whiskers. Default is
#'   0.95.
#' @param ndraws_use Optional positive integer. Number of posterior
#'   draws to use. If \code{NULL} (the default), up to 500 draws
#'   are used.
#' @param style Character. Visual style for displaying uncertainty.
#'   \code{"density"} (the default) overlays filled 2D density
#'   contours per item computed from the draw-level location and
#'   loading values, showing the full joint posterior uncertainty.
#'   \code{"crosshair"} shows point estimates with horizontal and
#'   vertical credible interval bars.
#'   \code{"both"} displays density contours with crosshairs on top.
#' @param density_alpha Numeric in \eqn{[0, 1]}. Maximum opacity
#'   of the density contours when \code{style} is \code{"density"}
#'   or \code{"both"}. Default is 0.3.
#' @param density_bins Integer. Number of contour bins for
#'   \code{\link[ggplot2]{geom_density_2d_filled}}. Default is 6.
#' @param density_palette An optional character vector of colors for
#'   the density contour fills (from low to high density). If
#'   \code{NULL} (the default), a blue sequential ramp is used. The
#'   length should match \code{density_bins}; the lowest
#'   (background) level is always transparent.
#' @param label_items Logical. If \code{TRUE} (the default), item
#'   names are displayed next to points.
#' @param point_size Numeric. Size of the item points. Default is
#'   2.5.
#' @param point_color Color for the item points and error bars.
#'   Default is \code{"#0072B2"}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{plot}{A \code{\link[ggplot2]{ggplot}} object showing item
#'     loadings on the first residual contrast (y-axis) versus item
#'     locations (x-axis).}
#'   \item{data}{A \code{\link[tibble]{tibble}} with columns:
#'     \code{item}, \code{location} (posterior mean item location),
#'     \code{location_lower}, \code{location_upper} (location CI),
#'     \code{loading} (posterior mean loading on first contrast),
#'     \code{loading_lower}, \code{loading_upper} (loading CI).}
#'   \item{eigenvalue}{A \code{\link[tibble]{tibble}} with columns:
#'     \code{eigenvalue_obs} (posterior mean observed eigenvalue),
#'     \code{eigenvalue_rep} (posterior mean replicated eigenvalue),
#'     \code{eigenvalue_diff} (posterior mean difference),
#'     \code{ppp} (posterior predictive p-value),
#'     \code{var_explained_obs}, \code{var_explained_rep} (posterior
#'     mean proportions of residual variance explained).}
#' }
#'
#' @details
#' The procedure for each posterior draw \eqn{s} is:
#'
#' \enumerate{
#'   \item Obtain category probabilities from
#'     \code{\link[brms]{posterior_epred}}. Compute expected values
#'     \eqn{E^{(s)}_{vi}} and variances \eqn{Var^{(s)}_{vi}}.
#'   \item Compute standardized residuals for observed data:
#'     \deqn{z^{obs(s)}_{vi} = \frac{X_{vi} - E^{(s)}_{vi}}
#'       {\sqrt{Var^{(s)}_{vi}}}}
#'   \item Generate replicated data \eqn{Y^{rep(s)}} from
#'     \code{\link[brms]{posterior_predict}} and compute standardized
#'     residuals:
#'     \deqn{z^{rep(s)}_{vi} = \frac{Y^{rep(s)}_{vi} -
#'       E^{(s)}_{vi}}{\sqrt{Var^{(s)}_{vi}}}}
#'   \item Reshape both sets of residuals into person \eqn{\times}
#'     item matrices and perform SVD on each.
#'   \item Extract the first-contrast eigenvalue and item loadings
#'     from both observed and replicated SVDs.
#'   \item Compare eigenvalues across draws: the posterior predictive
#'     p-value \code{ppp = mean(eigenvalue_obs > eigenvalue_rep)}
#'     tests whether the observed residual structure is stronger than
#'     what the model produces under its own assumptions.
#' }
#'
#' When \code{style = "density"} or \code{style = "both"}, the
#' draw-level (location, loading) pairs for each item are used to
#' construct filled 2D kernel density contours via
#' \code{\link[ggplot2]{geom_density_2d_filled}}. The lowest
#' contour level (outside all contours) is set to transparent so
#' the white panel background shows through. Higher density regions
#' use progressively darker fills.
#'
#' Item loadings are aligned across draws using majority-sign
#' alignment to resolve the sign indeterminacy of eigenvectors.
#'
#' @references
#' Smith, E. V. (2002). Detecting and evaluating the impact of
#' multidimensionality using item fit statistics and principal
#' component analysis of residuals. \emph{Journal of Applied
#' Measurement}, \emph{3}, 205--231.
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with
#' brms and Stan. \emph{Journal of Statistical Software}, \emph{100},
#' 1--54. \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link{plot_pim}} for person-item maps,
#' \code{\link{plot_ipf}} for item category probability curves,
#' \code{\link{q3_statistic}} for Q3 residual correlations (another
#' local dependence / dimensionality diagnostic).
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
#' # 2D density contours (default)
#' result <- plot_residual_pca(fit_pcm)
#' result$plot
#'
#' # Crosshair style
#' result_c <- plot_residual_pca(fit_pcm, style = "crosshair")
#' result_c$plot
#'
#' # Both combined
#' result_b <- plot_residual_pca(fit_pcm, style = "both")
#' result_b$plot
#'
#' # Custom warm palette
#' result_w <- plot_residual_pca(
#'   fit_pcm,
#'   density_palette = c("#FEE8C8", "#FDBB84", "#E34A33",
#'                       "#B30000", "#7F0000", "#4A0000")
#'   point_color = "#B30000"
#' )
#' result_w$plot
#' }
#'
#' @importFrom brms posterior_epred posterior_predict ndraws as_draws_df
#' @importFrom rlang enquo as_name .data
#' @importFrom stats formula quantile complete.cases aggregate
#' @importFrom tibble tibble
#' @importFrom grDevices colorRampPalette
#' @importFrom tidyr as_tibble
#' @export
plot_residual_pca <- function(
    model,
    item_var = item,
    person_var = id,
    center = TRUE,
    prob = 0.95,
    ndraws_use = NULL,
    style = c("both", "density", "crosshair"),
    density_alpha = 0.3,
    density_bins = 6,
    density_palette = NULL,
    label_items = TRUE,
    point_size = 2.5,
    point_color = "#0072B2"
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be a brmsfit object.", call. = FALSE)
  }

  style <- match.arg(style)
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

  resp_var <- as.character(formula(model)$formula[[2]])
  if (length(resp_var) > 1) resp_var <- resp_var[2]

  lower_prob <- (1 - prob) / 2
  upper_prob <- 1 - lower_prob

  # --- Draw subset ---
  if (is.null(ndraws_use)) {
    ndraws_use <- min(brms::ndraws(model), 500)
  } else {
    ndraws_use <- min(as.integer(ndraws_use), brms::ndraws(model))
  }
  draw_ids <- sample(seq_len(brms::ndraws(model)), ndraws_use)

  # --- Posterior predictions ---
  epred_array <- brms::posterior_epred(model, draw_ids = draw_ids)
  yrep_mat    <- brms::posterior_predict(model, draw_ids = draw_ids)

  n_draws <- dim(epred_array)[1]
  n_obs   <- dim(epred_array)[2]
  obs_response <- model$data[[resp_var]]

  # --- E and Var per draw x observation ---
  if (length(dim(epred_array)) == 3) {
    n_cat <- dim(epred_array)[3]
    cat_values <- seq_len(n_cat)
    E_mat <- apply(epred_array, c(1, 2), function(p) sum(cat_values * p))
    Var_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
    for (s in seq_len(n_draws)) {
      for (n in seq_len(n_obs)) {
        Var_mat[s, n] <- sum((cat_values - E_mat[s, n])^2 *
                               epred_array[s, n, ])
      }
    }
  } else {
    E_mat   <- epred_array
    Var_mat <- epred_array * (1 - epred_array)
  }

  Var_mat[Var_mat < 1e-12] <- 1e-12
  sd_mat <- sqrt(Var_mat)

  # --- Standardized residuals ---
  obs_mat    <- matrix(obs_response, nrow = n_draws, ncol = n_obs,
                       byrow = TRUE)
  zresid_obs <- (obs_mat - E_mat) / sd_mat
  zresid_rep <- (yrep_mat - E_mat) / sd_mat

  # --- Person x item mapping ---
  items   <- model$data[[item_name]]
  persons <- model$data[[person_name]]
  unique_items   <- unique(items)
  unique_persons <- unique(persons)
  k <- length(unique_items)
  n_persons <- length(unique_persons)

  person_idx <- match(persons, unique_persons)
  item_idx   <- match(items, unique_items)

  # --- PCA per draw: observed AND replicated ---
  loadings_obs_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  eigen_obs_draws    <- rep(NA_real_, n_draws)
  eigen_rep_draws    <- rep(NA_real_, n_draws)
  varexp_obs_draws   <- rep(NA_real_, n_draws)
  varexp_rep_draws   <- rep(NA_real_, n_draws)

  for (s in seq_len(n_draws)) {
    resid_obs_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)
    resid_rep_wide <- matrix(NA_real_, nrow = n_persons, ncol = k)

    for (obs in seq_len(n_obs)) {
      resid_obs_wide[person_idx[obs], item_idx[obs]] <- zresid_obs[s, obs]
      resid_rep_wide[person_idx[obs], item_idx[obs]] <- zresid_rep[s, obs]
    }

    complete_obs <- complete.cases(resid_obs_wide)
    complete_rep <- complete.cases(resid_rep_wide)
    resid_obs_c <- resid_obs_wide[complete_obs, , drop = FALSE]
    resid_rep_c <- resid_rep_wide[complete_rep, , drop = FALSE]

    if (nrow(resid_obs_c) < 3 || nrow(resid_rep_c) < 3) next

    sv_obs <- svd(resid_obs_c, nu = 0, nv = min(k, 2))
    loadings_obs_draws[s, ] <- sv_obs$v[, 1]
    eigen_obs_draws[s] <- sv_obs$d[1]^2 / (nrow(resid_obs_c) - 1)
    total_var_obs <- sum(sv_obs$d^2) / (nrow(resid_obs_c) - 1)
    varexp_obs_draws[s] <- eigen_obs_draws[s] / total_var_obs

    sv_rep <- svd(resid_rep_c, nu = 0, nv = min(k, 2))
    eigen_rep_draws[s] <- sv_rep$d[1]^2 / (nrow(resid_rep_c) - 1)
    total_var_rep <- sum(sv_rep$d^2) / (nrow(resid_rep_c) - 1)
    varexp_rep_draws[s] <- eigen_rep_draws[s] / total_var_rep
  }

  # --- Resolve sign indeterminacy ---
  ref_idx <- which(!is.na(loadings_obs_draws[, 1]))[1]
  if (is.na(ref_idx)) {
    stop("Could not compute residual PCA for any draw.", call. = FALSE)
  }
  ref_sign <- sign(loadings_obs_draws[ref_idx, ])

  for (s in seq_len(n_draws)) {
    if (is.na(loadings_obs_draws[s, 1])) next
    agreement <- sum(sign(loadings_obs_draws[s, ]) == ref_sign,
                     na.rm = TRUE)
    if (agreement < k / 2) {
      loadings_obs_draws[s, ] <- -loadings_obs_draws[s, ]
    }
  }

  # --- Item locations with full posterior (draw-level) ---
  draws_df <- tidyr::as_tibble(brms::as_draws_df(model))
  family_name <- family(model)$family
  is_ordinal <- grepl("acat|cumul|sratio|cratio",
                      family_name, ignore.case = TRUE)

  loc_draws_list <- .extract_item_location_draws(
    draws_df, model, unique_items, item_name, person_name, is_ordinal
  )

  loc_draws_mat <- matrix(NA_real_, nrow = n_draws, ncol = k)
  for (i in seq_len(k)) {
    loc_draws_mat[, i] <- loc_draws_list[[i]][draw_ids]
  }

  if (center) {
    shift_per_draw <- rowMeans(loc_draws_mat, na.rm = TRUE)
    loc_draws_mat <- loc_draws_mat - shift_per_draw
  }

  # Summarize locations
  loc_mean  <- colMeans(loc_draws_mat, na.rm = TRUE)
  loc_lower <- apply(loc_draws_mat, 2, stats::quantile,
                     probs = lower_prob, na.rm = TRUE)
  loc_upper <- apply(loc_draws_mat, 2, stats::quantile,
                     probs = upper_prob, na.rm = TRUE)

  # --- Summarize loadings ---
  loading_mean  <- colMeans(loadings_obs_draws, na.rm = TRUE)
  loading_lower <- apply(loadings_obs_draws, 2, stats::quantile,
                         probs = lower_prob, na.rm = TRUE)
  loading_upper <- apply(loadings_obs_draws, 2, stats::quantile,
                         probs = upper_prob, na.rm = TRUE)

  # --- Eigenvalue posterior predictive comparison ---
  valid <- !is.na(eigen_obs_draws) & !is.na(eigen_rep_draws) &
    eigen_obs_draws > 0 & eigen_rep_draws > 0
  ppp_eigen <- mean(eigen_obs_draws[valid] > eigen_rep_draws[valid])

  eigenvalue_summary <- tibble::tibble(
    eigenvalue_obs    = mean(eigen_obs_draws[valid]),
    eigenvalue_rep    = mean(eigen_rep_draws[valid]),
    eigenvalue_diff   = mean(eigen_obs_draws[valid] -
                               eigen_rep_draws[valid]),
    ppp               = ppp_eigen,
    var_explained_obs = mean(varexp_obs_draws[valid]),
    var_explained_rep = mean(varexp_rep_draws[valid])
  )

  # --- Build summary output data ---
  plot_df <- tibble::tibble(
    item           = unique_items,
    location       = loc_mean,
    location_lower = loc_lower,
    location_upper = loc_upper,
    loading        = loading_mean,
    loading_lower  = loading_lower,
    loading_upper  = loading_upper
  )

  # --- Build draw-level long data for density plots ---
  valid_draws <- which(!is.na(loadings_obs_draws[, 1]))
  draws_long_list <- vector("list", k)
  for (i in seq_len(k)) {
    draws_long_list[[i]] <- data.frame(
      item     = unique_items[i],
      location = loc_draws_mat[valid_draws, i],
      loading  = loadings_obs_draws[valid_draws, i],
      stringsAsFactors = FALSE
    )
  }
  draws_long <- do.call(rbind, draws_long_list)

  # --- Density fill palette ---
  # First value is always transparent (region outside all contours)
  # Remaining values are the actual contour fills (low → high density)
  if (is.null(density_palette)) {
    contour_colors <- grDevices::colorRampPalette(
      c("#DCEEF8", "#88BEDC", "#3A8EC2", "#0A5C96", "#003560")
    )(density_bins)
  } else {
    if (length(density_palette) < density_bins) {
      contour_colors <- grDevices::colorRampPalette(
        density_palette
      )(density_bins)
    } else {
      contour_colors <- density_palette[seq_len(density_bins)]
    }
  }
  fill_values <- c("transparent", contour_colors)



  # --- Build plot ---
  subtitle_text <- paste0(
    "First contrast eigenvalue: observed = ",
    round(eigenvalue_summary$eigenvalue_obs, 2),
    ", replicated = ",
    round(eigenvalue_summary$eigenvalue_rep, 2),
    ", ppp = ",
    round(eigenvalue_summary$ppp, 3)
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$location,
                                             y = .data$loading))

  # Density layer
  if (style %in% c("density", "both")) {
    p <- p +
      ggplot2::geom_density_2d_filled(
        data = draws_long,
        ggplot2::aes(x = .data$location, y = .data$loading,
                     group = .data$item),
        alpha = density_alpha,
        bins = density_bins,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = fill_values,
        guide  = "none"
      )
  }

  # Zero line
  p <- p +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", colour = "grey50"
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", colour = "grey50"
    )

  # Crosshair layer
  if (style %in% c("crosshair", "both")) {
    crosshair_alpha <- 1 #if (style == "both") 0.7 else 0.5
    p <- p +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data$location_lower,
                     xmax = .data$location_upper),
        width = 0, linewidth = 0.4, colour = point_color,
        alpha = crosshair_alpha
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$loading_lower,
                     ymax = .data$loading_upper),
        width = 0, linewidth = 0.4, colour = point_color,
        alpha = crosshair_alpha
      )
  }

  # Points always shown
  p <- p +
    ggplot2::geom_point(size = point_size, colour = point_color)

  # Labels
  if (label_items) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = .data$item),
        size = 4, colour = "grey30",
        max.overlaps = Inf,
        segment.colour = "grey70"
      )
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$item),
        size = 4, colour = "grey30",
        vjust = -1, hjust = 0.5
      )
    }
  }

  p <- p +
    ggplot2::labs(
      x = if (center) {
        expression("Centered item location" ~ (delta))
      } else {
        expression("Item location" ~ (delta))
      },
      y        = "Loading on first residual contrast",
      subtitle = subtitle_text
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background  = ggplot2::element_rect(fill = "white",
                                               colour = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      panel.grid.minor = ggplot2::element_blank()
    )

  list(plot = p, data = plot_df, eigenvalue = eigenvalue_summary)
}


# ── Internal: extract draw-level item locations ──────────────────
#' @keywords internal
.extract_item_location_draws <- function(draws, model, unique_items,
                                         item_name, person_name,
                                         is_ordinal) {

  result <- vector("list", length(unique_items))
  names(result) <- unique_items

  if (is_ordinal) {
    has_grouped <- any(grepl("^b_Intercept\\[.+,\\d+\\]$", names(draws)))

    for (i in seq_along(unique_items)) {
      item_label <- unique_items[i]

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

      if (length(thresh_cols) == 0) {
        result[[i]] <- rep(NA_real_, nrow(draws))
        next
      }

      thresh_mat <- as.matrix(draws[, thresh_cols, drop = FALSE])
      result[[i]] <- rowMeans(thresh_mat)
    }
  } else {
    intercept_col <- grep("^b_Intercept$", names(draws), value = TRUE)
    intercept_draws <- draws[[intercept_col]]

    for (i in seq_along(unique_items)) {
      item_label <- unique_items[i]

      re_col <- paste0("r_", item_name, "[", item_label, ",Intercept]")
      if (!re_col %in% names(draws)) {
        re_col_alt <- grep(
          paste0("^r_", item_name, "\\[",
                 gsub("([.|()\\^{}+$*?])", "\\\\\\1", item_label),
                 ",Intercept\\]$"),
          names(draws), value = TRUE
        )
        if (length(re_col_alt) == 1) re_col <- re_col_alt
        else {
          result[[i]] <- rep(NA_real_, nrow(draws))
          next
        }
      }

      result[[i]] <- -(intercept_draws + draws[[re_col]])
    }
  }

  result
}


# ── Internal: extract threshold data (used by plot_pim) ──────────
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
      re_col <- paste0("r_", item_name, "[", item_label, ",Intercept]")
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
