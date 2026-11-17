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
# Fixed item effects with an intercept can use any contrasts (treatment,
# contr.sum, contr_sum_zero(), ...). .item_levels() attaches what is needed
# to rebuild the contrasts of the item factor as attribute "contrast_spec".
# Only .item_effect_draws() builds the matrix, and only for fixed item
# effects with an intercept, so models with random items or item
# thresholds never evaluate contrasts. Without the attribute (item levels
# built by hand) treatment coding is assumed. When the factor stores no
# contrasts, brms used the contrasts option of the fitting session, which
# is not saved, so R's default and the current option are both tried.
#
# brms orders grouping levels as levels(factor(x)) and alters labels in
# parameter names (spaces, punctuation), so parameters are matched to item
# labels by position. Labels made only of letters and digits are unchanged
# by brms, and they are checked against the parameter names to catch any
# ordering mismatch.

# Ordered levels of the item variable, as brms uses them, with the
# contrast specification of the item factor attached as attribute
# "contrast_spec" (see .contrast_spec()).
.item_levels <- function(model, item_name) {
  x <- model$data[[item_name]]
  lev <- levels(factor(x))
  attr(lev, "contrast_spec") <- .contrast_spec(x)
  lev
}

# Candidate contrast matrices of the item factor: the stored contrasts,
# or, when none are stored, R's default (treatment for unordered,
# polynomial for ordered factors) and the current contrasts option if it
# differs. Treatment coding when no specification is attached.
.item_contrast_candidates <- function(item_levels) {
  lev <- as.character(item_levels)
  spec <- attr(item_levels, "contrast_spec")
  if (is.null(spec)) {
    return(list(.contrast_from_spec(list(stored = "contr.treatment"), lev)))
  }
  if (!is.null(spec$stored)) return(list(.contrast_from_spec(spec, lev)))
  pick <- if (isTRUE(spec$ordered)) 2L else 1L
  default <- c("contr.treatment", "contr.poly")[[pick]]
  session <- getOption("contrasts", c("contr.treatment", "contr.poly"))[[pick]]
  cands <- list(default)
  if (!identical(session, default)) cands <- c(cands, list(session))
  lapply(cands, function(cm) .contrast_from_spec(list(stored = cm), lev))
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

# Names as brms writes them in parameter names: the replacements of the
# internal brms:::rename() (brms 2.21-2.23), applied to model matrix
# column names such as "itemcrit 2".
.brms_rename <- function(x) {
  pattern <- c(" ", "(", ")", "[", "]", ",", "\"", "'", "?", "+", "-",
               "*", "/", "^", "=", "$")
  replacement <- c(rep("", 9), "P", "M", "MU", "D", "E", "EQ", "USD")
  for (i in seq_along(pattern)) {
    x <- gsub(pattern[i], replacement[i], x, fixed = TRUE)
  }
  x
}

# Item effect draws (draws x items) added to the linear predictor, from
# random item intercepts or fixed item effects. NULL if the model has none.
# other_vars: names of the model's variables (names(model$data)), so that
# coefficients of a variable whose name starts with item_name (rater_exp
# for rater) are not taken for this factor's in the error check below.
.item_effect_draws <- function(draws, item_levels, item_name,
                               other_vars = NULL) {
  nm <- names(draws)
  re_pat <- paste0("^r_", item_name, "\\[(.+),Intercept\\]$")
  re_cols <- grep(re_pat, nm, value = TRUE)
  if (length(re_cols)) {
    .check_label_map(sub(re_pat, "\\1", re_cols), item_levels,
                     "item random effects")
    return(as.matrix(draws[, re_cols, drop = FALSE]))
  }
  # Fixed effects are matched by their exact names, which brms builds from
  # the model matrix column names (variable name + level or contrast
  # column name) with .brms_rename(). A prefix match would also pick up
  # other variables that start with the same name (rater, rater_exp).
  k <- length(item_levels)
  fe_k <- paste0("b_", .brms_rename(paste0(item_name, item_levels)))
  if (all(fe_k %in% nm)) {
    # No intercept: one coefficient per item
    return(as.matrix(draws[, fe_k, drop = FALSE]))
  }
  # With an intercept: item effects = coefficients x t(contrasts), for
  # each candidate contrast matrix whose column names match.
  matched <- list()
  for (cm in .item_contrast_candidates(item_levels)) {
    cn <- colnames(cm)
    if (is.null(cn)) cn <- as.character(seq_len(ncol(cm)))
    fe_k1 <- paste0("b_", .brms_rename(paste0(item_name, cn)))
    if (ncol(cm) == k - 1L && all(fe_k1 %in% nm)) {
      matched <- c(matched,
                   list(as.matrix(draws[, fe_k1, drop = FALSE]) %*% t(cm)))
    }
  }
  if (length(matched) == 1L) return(matched[[1L]])
  how <- paste0("Set options(contrasts = ...) to the contrasts used when ",
                "fitting, or set them on the factor with contrasts(data$",
                item_name, ") <- ... and refit.")
  if (length(matched) > 1L) {
    # Labels such as 1, 2, 3 give the same names under treatment and
    # numbered contrasts (contr.sum), so the names cannot tell them apart.
    stop("The fixed effects of '", item_name, "' match both treatment ",
         "coding and the current contrasts option. ", how, call. = FALSE)
  }
  # Coefficients that look like this factor's (k or k - 1 of them with its
  # name as prefix) but match no naming: most likely the contrasts differ
  # from those at fitting. Fewer or more are other variables.
  pref <- nm[startsWith(nm, paste0("b_", item_name))]
  other_vars <- as.character(other_vars)
  longer <- other_vars[startsWith(other_vars, item_name) &
                         other_vars != item_name]
  for (v in longer) {
    pref <- pref[!startsWith(pref, paste0("b_", .brms_rename(v)))]
  }
  n_prefix <- length(pref)
  if (n_prefix %in% c(k, k - 1L)) {
    stop("Could not match the fixed effects of '", item_name, "' to its ",
         "levels and contrasts. ", how, call. = FALSE)
  }
  NULL
}

# Threshold draws per item.
#
# draws: data frame of posterior draws (as_draws_df).
# item_levels: item labels in brms order, see .item_levels().
# ordinal: TRUE for ordinal families, FALSE for bernoulli.
# other_vars: passed to .item_effect_draws().
# Returns a named list (names = item_levels) of draws x thresholds matrices.
.item_tau_draws <- function(draws, item_levels, item_name, ordinal,
                            other_vars = NULL) {
  draws <- as.data.frame(draws)
  nm <- names(draws)
  k <- length(item_levels)
  n_draws <- nrow(draws)
  eff <- .item_effect_draws(draws, item_levels, item_name, other_vars)

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
  .item_tau_draws(draws, .item_levels(model, item_name), item_name, ordinal,
                  other_vars = names(model$data))
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
