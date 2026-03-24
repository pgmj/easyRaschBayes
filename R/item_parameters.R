#' Extract Item Parameters from a Bayesian Rasch Model
#'
#' Extracts item difficulty (threshold) parameters from a fitted
#' Bayesian Rasch model. Returns a simple location table in both long
#' and wide formats, a full summary with posterior SEs and HDCIs,
#' item-level information, threshold ordering diagnostics, and
#' optionally the full posterior draws matrix.
#'
#' @param model A fitted \code{\link[brms]{brmsfit}} object. Supported
#'   parameterisations:
#'   \describe{
#'     \item{Polytomous ordinal (PCM)}{e.g., \code{family = acat} with
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
#' @param draws Logical. If \code{TRUE}, a draws matrix of full
#'   posterior draws is included in the output. Default is
#'   \code{FALSE}.
#' @param center Logical. If \code{TRUE} (the default), item
#'   parameters are shifted so that the grand mean of all threshold
#'   locations is zero, matching the frequentist CML convention and
#'   the centering used in \code{\link{plot_targeting}}.
#' @param prob Numeric in \eqn{(0, 1)}. Width of the highest density
#'   continuous interval (HDCI) reported in the summary. Default is
#'   0.95.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{locations}{A \code{\link[tibble]{tibble}} in long format
#'     with one row per item (dichotomous) or per item-threshold
#'     (polytomous), containing \code{item}, \code{threshold}
#'     (integer, always 1 for dichotomous), and \code{location}
#'     (posterior mean on the logit scale).}
#'   \item{locations_wide}{A \code{\link[tibble]{tibble}} in wide
#'     format with one row per item, sorted by mean location. For
#'     polytomous models, threshold columns are named \code{t1},
#'     \code{t2}, etc., and a \code{location} column gives the mean
#'     across thresholds. For dichotomous models, only \code{item}
#'     and \code{location} columns are present.}
#'   \item{summary}{A \code{\link[tibble]{tibble}} extending the long
#'     \code{locations} with: \code{se} (posterior SD),
#'     \code{hdci_lower} and \code{hdci_upper} (highest density
#'     continuous interval bounds at the level specified by
#'     \code{prob}), and \code{n_eff} (effective sample size for the
#'     parameter).}
#'   \item{item_information}{A \code{\link[tibble]{tibble}} with one
#'     row per item containing: \code{item}, \code{location} (mean
#'     item location, i.e. mean of thresholds for polytomous items),
#'     \code{info_at_location} (Fisher information at the item's own
#'     location), and \code{max_info} (maximum Fisher information
#'     across theta). For Rasch dichotomous items, both are 0.25.
#'     For polytomous items, information depends on the number and
#'     spacing of thresholds.}
#'   \item{threshold_order}{(Polytomous models only) A
#'     \code{\link[tibble]{tibble}} with one row per item containing:
#'     \code{item}, \code{n_thresholds}, \code{ordered} (logical:
#'     are all thresholds in ascending order?), and
#'     \code{prob_disordered} (posterior probability that at least
#'     one pair of adjacent thresholds is disordered, i.e.
#'     \eqn{\tau_{k+1} \le \tau_k}). \code{NULL} for dichotomous
#'     models.}
#'   \item{person_sd}{A \code{\link[tibble]{tibble}} with one row
#'     containing: \code{mean}, \code{sd}, \code{hdci_lower},
#'     \code{hdci_upper} — the posterior summary of the person-level
#'     standard deviation parameter \eqn{\sigma_\theta}.}
#'   \item{draws_matrix}{(Only if \code{draws = TRUE}) A numeric
#'     matrix with rows = thresholds (named
#'     \code{"item[threshold]"}) and columns = posterior draws. For
#'     dichotomous models, row names are item labels only.}
#' }
#'
#' @details
#' \strong{Dichotomous models} with random item effects
#' (\code{response ~ 1 + (1 | item) + (1 | id)}) parameterise item
#' difficulty as \eqn{\delta_i = -(b_0 + r_i)} where \eqn{b_0} is
#' the global intercept and \eqn{r_i} is the item random effect.
#' Models with fixed item effects (\code{response ~ 0 + item +
#' (1 | id)}) parameterise difficulty as \eqn{\delta_i = -b_i}.
#'
#' \strong{Polytomous (acat/PCM) models} with grouped thresholds
#' (\code{thres(gr = item)}) directly estimate item-specific
#' threshold parameters. Each row in the long output represents one
#' threshold within one item; each row in the wide output represents
#' one item.
#'
#' \strong{Item information} is computed from the posterior mean
#' item parameters using the standard Rasch/PCM information
#' formulae:
#' \describe{
#'   \item{Dichotomous}{
#'     \eqn{I_i(\theta) = P_i(\theta) Q_i(\theta)} where
#'     \eqn{P_i = \text{logistic}(\theta - \delta_i)}.}
#'   \item{Polytomous (PCM)}{
#'     \eqn{I_i(\theta) = \sum_c (c - E_i)^2 P_{ic}(\theta)} where
#'     \eqn{E_i = \sum_c c \cdot P_{ic}} is the expected score for
#'     item \eqn{i}.}
#' }
#'
#' \strong{Threshold ordering}: In the partial credit model,
#' disordered thresholds (\eqn{\tau_{k+1} \le \tau_k}) indicate that
#' the probability of responding in the intermediate category never
#' exceeds both adjacent categories — the category is empirically
#' "absorbed". This does not necessarily indicate misfit (see
#' Adams et al., 2012), but may suggest response categories should
#' be collapsed. The \code{prob_disordered} column reports the
#' posterior probability of at least one disordered pair per item,
#' providing a Bayesian alternative to post-hoc threshold checks.
#'
#' @references
#' Adams, R. J., Wu, M. L., & Wilson, M. (2012). The Rasch rating
#' model and the disordered threshold controversy. \emph{Educational
#' and Psychological Measurement}, \emph{72}(4), 547--573.
#' \doi{10.1177/0013164411432166}
#'
#' Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R
#' with brms and Stan. \emph{Journal of Statistical Software},
#' \emph{100}, 1--54. \doi{10.18637/jss.v100.i05}
#'
#' @seealso
#' \code{\link{person_parameters}}, \code{\link{plot_targeting}},
#' \code{\link{plot_ipf}}, \code{\link{posterior_to_prior}}.
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
#'   chains = 4, cores = 2, iter = 1000 # use more iter and cores
#' )
#'
#' ip <- item_parameters(fit_pcm)
#'
#' # Long format: one row per threshold
#' ip$locations
#'
#' # Wide format: one row per item, easy to scan
#' ip$locations_wide
#'
#' # Full summary with SE and HDCI
#' ip$summary
#'
#' # Item information
#' ip$item_information
#'
#' # Threshold ordering diagnostic
#' ip$threshold_order
#'
#' # Person SD
#' ip$person_sd
#'
#' # With full posterior draws
#' ip_draws <- item_parameters(fit_pcm, draws = TRUE)
#' dim(ip_draws$draws_matrix)
#'
#' # --- Dichotomous Rasch ---
#'
#' df_rm <- eRm::raschdat3 %>%
#'   as.data.frame() %>%
#'   rownames_to_column("id") %>%
#'   pivot_longer(!id, names_to = "item", values_to = "response")
#'
#' fit_rm <- brm(
#'   response ~ 1 + (1 | item) + (1 | id),
#'   data   = df_rm,
#'   family = bernoulli(),
#'   chains = 4, cores = 2, iter = 1000 # use more iter and cores
#' )
#'
#' ip_rm <- item_parameters(fit_rm)
#' # Wide and long are equivalent for dichotomous:
#' ip_rm$locations
#' ip_rm$locations_wide
#' ip_rm$summary
#' }
#'
#' @importFrom brms as_draws_df ndraws
#' @importFrom ggdist mean_hdci
#' @importFrom rlang enquo as_name .data
#' @importFrom stats family sd plogis optimize var acf
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr arrange
#' @export
item_parameters <- function(
    model,
    item_var   = item,
    person_var = id,
    draws      = FALSE,
    center     = TRUE,
    prob       = 0.95
) {
  
  # ================================================================
  # INPUT VALIDATION
  # ================================================================
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
  
  all_draws   <- tibble::as_tibble(brms::as_draws_df(model))
  family_name <- stats::family(model)$family
  is_ordinal  <- grepl("acat|cumul|sratio|cratio",
                       family_name, ignore.case = TRUE)
  
  unique_items <- sort(unique(model$data[[item_name]]))
  
  # ================================================================
  # DETECT MODEL PARAMETERISATION
  # ================================================================
  fe_pattern <- paste0("^b_", item_name)
  fe_cols    <- grep(fe_pattern, names(all_draws), value = TRUE)
  fe_cols    <- setdiff(fe_cols, "b_Intercept")
  
  has_fixed_items  <- !is_ordinal && length(fe_cols) > 0
  has_random_items <- !is_ordinal && any(grepl(
    paste0("^r_", item_name, "\\["), names(all_draws)
  ))
  
  if (!is_ordinal && !has_fixed_items && !has_random_items) {
    stop("Could not detect item parameters in the model. ",
         "Supported parameterisations:\n",
         "  - response ~ 0 + item + (1 | id)    [fixed items]\n",
         "  - response ~ 1 + (1 | item) + (1 | id)  [random items]\n",
         "  - response | thres(gr = item) ~ 1 + (1 | id) [PCM/acat]",
         call. = FALSE)
  }
  
  # ================================================================
  # EXTRACT DRAW-LEVEL ITEM PARAMETERS
  # ================================================================
  param_draws <- .extract_item_param_draws(
    all_draws, unique_items, item_name,
    is_ordinal, has_fixed_items
  )
  
  # ================================================================
  # CENTERING
  # ================================================================
  all_means <- vapply(
    unlist(param_draws, recursive = FALSE),
    mean, numeric(1)
  )
  shift <- if (center) mean(all_means) else 0
  
  param_draws <- lapply(param_draws, function(item_list) {
    lapply(item_list, function(draw_vec) draw_vec - shift)
  })
  
  # ================================================================
  # BUILD LOCATIONS TABLE — LONG FORMAT
  # ================================================================
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
  
  # ================================================================
  # BUILD LOCATIONS TABLE — WIDE FORMAT
  # ================================================================
  max_thresh <- max(vapply(param_draws, length, integer(1)))
  wide_rows  <- list()
  
  for (item_label in names(param_draws)) {
    thresh_list  <- param_draws[[item_label]]
    thresh_means <- vapply(thresh_list, mean, numeric(1))
    
    row <- list(item = item_label)
    
    if (max_thresh == 1) {
      # Dichotomous: single "location" column
      row[["location"]] <- round(thresh_means[1], 4)
    } else {
      # Polytomous: t1, t2, ... + mean location
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
  
  # ================================================================
  # BUILD SUMMARY TABLE (with SE, HDCI, n_eff)
  # ================================================================
  sum_rows  <- list()
  param_idx <- 0
  for (item_label in names(param_draws)) {
    thresh_list <- param_draws[[item_label]]
    for (k in seq_along(thresh_list)) {
      param_idx <- param_idx + 1
      vals <- thresh_list[[k]]
      
      hdci  <- ggdist::mean_hdci(vals, .width = prob)
      n_eff <- .compute_ess(vals)
      
      sum_rows[[param_idx]] <- data.frame(
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
  
  # ================================================================
  # ITEM INFORMATION
  # ================================================================
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
  
  # ================================================================
  # THRESHOLD ORDERING (polytomous only)
  # ================================================================
  threshold_order <- NULL
  if (is_ordinal) {
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
      
      n_draws_total    <- length(thresh_list[[1]])
      disordered_count <- 0
      for (s in seq_len(n_draws_total)) {
        draw_thresholds <- vapply(
          thresh_list, function(x) x[s], numeric(1)
        )
        if (any(diff(draw_thresholds) <= 0)) {
          disordered_count <- disordered_count + 1
        }
      }
      prob_disord <- disordered_count / n_draws_total
      
      ord_rows[[length(ord_rows) + 1]] <- data.frame(
        item            = item_label,
        n_thresholds    = n_thresh,
        ordered         = is_ordered,
        prob_disordered = round(prob_disord, 4),
        stringsAsFactors = FALSE
      )
    }
    threshold_order <- tibble::as_tibble(do.call(rbind, ord_rows))
  }
  
  # ================================================================
  # PERSON SD
  # ================================================================
  person_sd_col <- grep(
    paste0("^sd_", person_name, "__Intercept$"),
    names(all_draws), value = TRUE
  )
  if (length(person_sd_col) == 1) {
    sd_vals <- all_draws[[person_sd_col]]
    sd_hdci <- ggdist::mean_hdci(sd_vals, .width = prob)
    person_sd <- tibble::tibble(
      mean       = round(mean(sd_vals), 4),
      sd         = round(stats::sd(sd_vals), 4),
      hdci_lower = round(sd_hdci$ymin, 4),
      hdci_upper = round(sd_hdci$ymax, 4)
    )
  } else {
    warning("Could not find person-level SD parameter '",
            paste0("sd_", person_name, "__Intercept"),
            "' in the posterior draws.", call. = FALSE)
    person_sd <- NULL
  }
  
  # ================================================================
  # FULL POSTERIOR DRAWS (optional)
  # ================================================================
  out <- list(
    locations        = locations_long,
    locations_wide   = locations_wide,
    summary          = summary_tbl,
    item_information = item_info,
    threshold_order  = threshold_order,
    person_sd        = person_sd
  )
  
  if (draws) {
    draw_list <- list()
    row_names <- character(0)
    for (item_label in names(param_draws)) {
      thresh_list <- param_draws[[item_label]]
      for (k in seq_along(thresh_list)) {
        draw_list[[length(draw_list) + 1]] <- thresh_list[[k]]
        if (length(thresh_list) == 1) {
          row_names <- c(row_names, item_label)
        } else {
          row_names <- c(row_names,
                         paste0(item_label, "[", k, "]"))
        }
      }
    }
    draws_mat <- do.call(rbind, draw_list)
    rownames(draws_mat) <- row_names
    out$draws_matrix <- draws_mat
  }
  
  out
}


# ── Internal: extract draw-level item threshold parameters ───────
#' @keywords internal
.extract_item_param_draws <- function(draws, unique_items, item_name,
                                      is_ordinal, has_fixed_items) {
  
  result <- vector("list", length(unique_items))
  names(result) <- unique_items
  
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
      
      if (length(thresh_cols) == 0) {
        warning("No threshold parameters found for item '",
                item_label, "'.", call. = FALSE)
        result[[item_label]] <- list()
        next
      }
      
      thresh_nums <- as.numeric(
        gsub(".*,(\\d+)\\]$|.*\\[(\\d+)\\]$", "\\1\\2", thresh_cols)
      )
      thresh_cols <- thresh_cols[order(thresh_nums)]
      
      result[[item_label]] <- lapply(
        thresh_cols,
        function(col) as.numeric(draws[[col]])
      )
    }
    
  } else if (has_fixed_items) {
    for (item_label in unique_items) {
      col_pattern <- paste0("^b_", item_name,
                            gsub("([.|()\\^{}+$*?])", "\\\\\\1",
                                 item_label), "$")
      col_name <- grep(col_pattern, names(draws), value = TRUE)
      
      if (length(col_name) == 0) {
        warning("Could not find fixed effect for item '",
                item_label, "'.", call. = FALSE)
        result[[item_label]] <- list(rep(NA_real_, nrow(draws)))
        next
      }
      
      result[[item_label]] <- list(-as.numeric(draws[[col_name[1]]]))
    }
    
  } else {
    intercept_col <- grep("^b_Intercept$", names(draws), value = TRUE)
    intercept_draws <- if (length(intercept_col) == 1) {
      as.numeric(draws[[intercept_col]])
    } else {
      rep(0, nrow(draws))
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
                  item_label, "'.", call. = FALSE)
          result[[item_label]] <- list(rep(NA_real_, nrow(draws)))
          next
        }
      }
      
      result[[item_label]] <- list(
        -(intercept_draws + as.numeric(draws[[re_col]]))
      )
    }
  }
  
  result
}


# ── Internal: PCM item information at a given theta ──────────────
#' @keywords internal
.pcm_item_info <- function(theta, thresholds) {
  n_cat <- length(thresholds) + 1
  cats  <- seq(0, n_cat - 1)
  
  cum_eta   <- c(0, cumsum(theta - thresholds))
  max_eta   <- max(cum_eta)
  log_denom <- max_eta + log(sum(exp(cum_eta - max_eta)))
  p_cat     <- exp(cum_eta - log_denom)
  
  E_i <- sum(cats * p_cat)
  V_i <- sum((cats - E_i)^2 * p_cat)
  V_i
}


# ── Internal: compute effective sample size ──────────────────────
#' @keywords internal
.compute_ess <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) return(NA_real_)
  
  v <- stats::var(x)
  if (v < .Machine$double.eps) return(n)
  
  acf_vals <- stats::acf(x, lag.max = n - 1, plot = FALSE,
                         demean = TRUE)$acf[, 1, 1]
  
  sum_rho <- 0
  lag <- 2
  while (lag < length(acf_vals)) {
    pair_sum <- acf_vals[lag] + acf_vals[lag + 1]
    if (pair_sum < 0) break
    sum_rho <- sum_rho + pair_sum
    lag <- lag + 2
  }
  
  ess <- n / (1 + 2 * sum_rho)
  max(1, ess)
}