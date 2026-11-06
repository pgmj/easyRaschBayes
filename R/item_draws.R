# Shared extraction of item threshold draws from a brms Rasch-family model.
#
# Every function that needs item locations or thresholds goes through
# .item_tau_draws(). It returns, for each item, a draws x thresholds matrix
# on the scale where the person location theta enters as theta - tau
# (a dichotomous item has one "threshold", its difficulty).
#
# Supported parameterisations:
#   Ordinal, grouped thresholds:   y | thres(gr = item) ~ 1 + (1 | id)
#   Ordinal, common thresholds (rating scale model) with item effects:
#                                  y ~ 1 + (1 | item) + (1 | id), or + item
#   Bernoulli, random items:       y ~ 1 + (1 | item) + (1 | id)
#   Bernoulli, fixed items:        y ~ 0 + item + (1 | id), or 1 + item
#
# brms orders grouping levels as levels(factor(x)) and alters labels in
# parameter names (spaces, punctuation), so parameters are matched to item
# labels by position. Labels made only of letters and digits are unchanged
# by brms, and they are checked against the parameter names to catch any
# ordering mismatch.

# Ordered levels of the item variable, as brms uses them.
.item_levels <- function(model, item_name) {
  levels(factor(model$data[[item_name]]))
}

# Check a positional match of parameter labels to item levels.
.check_label_map <- function(param_labels, item_levels, what) {
  if (length(param_labels) != length(item_levels)) {
    stop("Found ", length(param_labels), " ", what, " but ",
         length(item_levels), " items.", call. = FALSE)
  }
  plain <- grepl("^[A-Za-z0-9]+$", item_levels)
  if (any(param_labels[plain] != item_levels[plain])) {
    stop("Could not match ", what, " to item labels by position. ",
         "Please report this with your model formula.", call. = FALSE)
  }
  invisible(TRUE)
}

# Item effect draws (draws x items) added to the linear predictor, from
# random item intercepts or fixed item effects. NULL if the model has none.
.item_effect_draws <- function(draws, item_levels, item_name) {
  nm <- names(draws)
  re_pat <- paste0("^r_", item_name, "\\[(.+),Intercept\\]$")
  re_cols <- grep(re_pat, nm, value = TRUE)
  if (length(re_cols)) {
    .check_label_map(sub(re_pat, "\\1", re_cols), item_levels,
                     "item random effects")
    return(as.matrix(draws[, re_cols, drop = FALSE]))
  }
  fe_pat <- paste0("^b_", item_name, "(.+)$")
  fe_cols <- grep(fe_pat, nm, value = TRUE)
  if (length(fe_cols)) {
    k <- length(item_levels)
    labs <- sub(fe_pat, "\\1", fe_cols)
    out <- matrix(0, nrow = nrow(draws), ncol = k)
    if (length(fe_cols) == k) {
      # No intercept: one coefficient per item
      .check_label_map(labs, item_levels, "item fixed effects")
      out[] <- as.matrix(draws[, fe_cols, drop = FALSE])
    } else if (length(fe_cols) == k - 1L) {
      # Treatment coding: the first level is the reference
      .check_label_map(labs, item_levels[-1], "item fixed effects")
      out[, -1] <- as.matrix(draws[, fe_cols, drop = FALSE])
    } else {
      stop("Found ", length(fe_cols), " item fixed effects for ", k,
           " items.", call. = FALSE)
    }
    return(out)
  }
  NULL
}

# Threshold draws per item.
#
# draws: data frame of posterior draws (as_draws_df).
# item_levels: item labels in brms order, see .item_levels().
# ordinal: TRUE for ordinal families, FALSE for bernoulli.
# Returns a named list (names = item_levels) of draws x thresholds matrices.
.item_tau_draws <- function(draws, item_levels, item_name, ordinal) {
  draws <- as.data.frame(draws)
  nm <- names(draws)
  k <- length(item_levels)
  n_draws <- nrow(draws)
  eff <- .item_effect_draws(draws, item_levels, item_name)

  if (ordinal) {
    grp_pat <- "^b_Intercept\\[(.+),(\\d+)\\]$"
    grp_cols <- grep(grp_pat, nm, value = TRUE)
    if (length(grp_cols)) {
      grp <- sub(grp_pat, "\\1", grp_cols)
      thr <- as.integer(sub(grp_pat, "\\2", grp_cols))
      groups <- unique(grp)
      .check_label_map(groups, item_levels, "threshold groups")
      out <- lapply(seq_len(k), function(i) {
        cols <- grp_cols[grp == groups[i]][order(thr[grp == groups[i]])]
        m <- as.matrix(draws[, cols, drop = FALSE])
        if (!is.null(eff)) m <- m - eff[, i]
        unname(m)
      })
    } else {
      com_cols <- grep("^b_Intercept\\[\\d+\\]$", nm, value = TRUE)
      if (!length(com_cols)) {
        stop("Could not find threshold parameters in the model.",
             call. = FALSE)
      }
      com_cols <- com_cols[order(as.integer(gsub("\\D", "", com_cols)))]
      tau <- as.matrix(draws[, com_cols, drop = FALSE])
      if (is.null(eff)) {
        stop("The model has common thresholds but no item effects, ",
             "so all items are identical. Use thres(gr = item) or add ",
             "(1 | item).", call. = FALSE)
      }
      out <- lapply(seq_len(k), function(i) unname(tau - eff[, i]))
    }
  } else {
    if (is.null(eff)) {
      stop("Could not find item parameters in the model. Supported: ",
           "response ~ 1 + (1 | item) + (1 | id) or ",
           "response ~ 0 + item + (1 | id).", call. = FALSE)
    }
    b0 <- if ("b_Intercept" %in% nm) draws[["b_Intercept"]] else 0
    out <- lapply(seq_len(k), function(i) matrix(-(b0 + eff[, i]), ncol = 1))
  }
  names(out) <- item_levels
  out
}

# Convenience wrapper for a brmsfit.
.model_tau_draws <- function(model, item_name,
                             draws = brms::as_draws_df(model)) {
  family_name <- stats::family(model)$family
  ordinal <- grepl("acat|cumul|sratio|cratio|hurdle_acat", family_name,
                   ignore.case = TRUE)
  .item_tau_draws(draws, .item_levels(model, item_name), item_name, ordinal)
}

# Inverse link for category probabilities, after checking that the model
# has unit discrimination and a link the formulas support. acat curves are
# written for the logit link. cumulative, sratio and cratio work with any
# link that stats::make.link() knows.
.ordinal_linkinv <- function(model, fun_name) {
  fam <- stats::family(model)
  vars <- brms::variables(model)
  if ("disc" %in% names(model$formula$pforms) ||
      ("disc" %in% vars &&
       any(abs(as.numeric(brms::as_draws_matrix(model, variable = "disc")) - 1) >
             1e-8))) {
    stop("`", fun_name, "()` assumes equal item discrimination (disc = 1), ",
         "but the model estimates `disc`.", call. = FALSE)
  }
  if (grepl("acat", fam$family) && !identical(fam$link, "logit")) {
    stop("`", fun_name, "()` supports the acat family with the logit link ",
         "only. The model uses '", fam$link, "'.", call. = FALSE)
  }
  if (identical(fam$link, "logit")) return(stats::plogis)
  linkinv <- tryCatch(stats::make.link(fam$link)$linkinv,
                      error = function(e) NULL)
  if (is.null(linkinv)) {
    stop("`", fun_name, "()` does not support the '", fam$link, "' link.",
         call. = FALSE)
  }
  linkinv
}
