# Many-facet Rasch models: contrasts and facet handling.
#
# A many-facet model adds a third variable, the facet (for example raters),
# to the person x item layout of the standard models. The fit functions
# treat each item x facet combination ("cell") as an item. .resolve_facet()
# collects what they need to know about the facet.

#' Orthonormal sum-to-zero contrasts
#'
#' Contrasts for a factor whose effects should sum to zero, with the same
#' prior variance for every level. Useful for small facets, such as a few
#' items or raters, modelled as fixed effects in a many-facet Rasch model.
#'
#' With sum-to-zero contrasts the model estimates the effect of each level
#' relative to the average level, instead of relative to a reference level.
#' R's [stats::contr.sum()] does this by setting the last level to minus the
#' sum of the others. Normal priors on the coefficients then give the last
#' level a larger prior variance than the others: \eqn{K - 1} times larger
#' for \eqn{K} levels, three times larger with four items. With little data
#' the prior is then not the same for all levels.
#'
#' `contr_sum_zero()` returns orthonormal contrasts (a normalised Helmert
#' basis). They describe the same sum-to-zero effects, but independent
#' normal priors on the coefficients give every level the same prior
#' variance. This is the coding used by Stan's `sum_to_zero_vector` type
#' (Morris, 2025), which brms does not offer directly.
#'
#' A normal(0, s) prior on the coefficients gives each level effect a prior
#' SD of \eqn{s \sqrt{(K - 1) / K}}. For a prior SD of \eqn{\sigma} per level,
#' use \eqn{s = \sigma \sqrt{K / (K - 1)}}. For example, with four levels and
#' \eqn{\sigma = 1}, use normal(0, 1.15).
#'
#' The item and facet extractors in this package read the contrasts stored
#' with the model data, so the level effects are reported on the original
#' scale whatever contrasts were used.
#'
#' @param n The number of levels, or a vector of level labels (for example
#'   `levels(x)`), as in [stats::contr.sum()].
#' @param contrasts If `FALSE`, one indicator column per level is returned
#'   instead, as for [stats::contr.sum()]. R uses this for factors without
#'   an intercept to pair with.
#' @param sparse Only `FALSE` is supported. The argument is there so that
#'   the function can be named in `options(contrasts = ...)` or set with
#'   `contrasts(x) <- "contr_sum_zero"`.
#'
#' @return A matrix with one row per level and \eqn{K - 1} columns. Columns
#'   are orthonormal and every column sums to zero. With
#'   `contrasts = FALSE`, the \eqn{K x K} identity matrix.
#'
#' @references
#' Morris, M. (2025). *The sum-to-zero constraint in Stan*. Stan case study.
#' <https://mc-stan.org/learn-stan/case-studies/sum_to_zero_vector.html>
#'
#' @seealso [stats::contr.sum()], [stats::contrasts()]
#'
#' @examples
#' contr_sum_zero(4)
#'
#' # The implied prior covariance of the level effects is proportional to
#' # diag(K) - 1/K: equal variances, equal negative correlations
#' H <- contr_sum_zero(c("A", "B", "C", "D"))
#' round(H %*% t(H), 3)
#'
#' \dontrun{
#' d <- data.frame(...)  # long data: id, item, rater, response
#' contrasts(d$item) <- contr_sum_zero(levels(d$item))
#' fit <- brms::brm(
#'   response ~ 1 + item + (1 | rater) + (1 | id),
#'   data = d, family = brms::acat(),
#'   prior = brms::prior(normal(0, 1.15), class = b)
#' )
#' }
#' @export
contr_sum_zero <- function(n, contrasts = TRUE, sparse = FALSE) {
  if (isTRUE(sparse)) {
    stop("`contr_sum_zero()` does not return sparse contrasts.",
         call. = FALSE)
  }
  if (is.numeric(n) && length(n) == 1L) {
    if (n != round(n)) {
      stop("'n' must be a whole number or a vector of level labels.",
           call. = FALSE)
    }
    lev <- as.character(seq_len(n))
  } else {
    lev <- as.character(n)
  }
  k <- length(lev)
  if (k < 2L) {
    stop("Contrasts need at least two levels.", call. = FALSE)
  }
  if (!contrasts) {
    m <- diag(k)
    dimnames(m) <- list(lev, lev)
    return(m)
  }
  h <- stats::contr.helmert(k)
  h <- sweep(h, 2L, sqrt(colSums(h^2)), "/")
  dimnames(h) <- list(lev, NULL)
  h
}

# What is needed to build the contrast matrix later: the stored contrasts
# attribute (NULL, a matrix, a function or its name) and whether the factor
# is ordered. Cheap to keep, since nothing is evaluated.
.contrast_spec <- function(x) {
  list(stored = if (is.factor(x)) attr(x, "contrasts") else NULL,
       ordered = is.ordered(x))
}

# Contrast matrix (levels x columns) from a specification: the stored
# contrasts, or without them the current contrasts option.
# .item_contrast_candidates() decides which specifications to try.
.contrast_from_spec <- function(spec, lev) {
  cm <- spec$stored
  if (is.null(cm)) {
    opt <- getOption("contrasts", c("contr.treatment", "contr.poly"))
    cm <- if (isTRUE(spec$ordered)) opt[[2L]] else opt[[1L]]
  }
  if (is.character(cm) || is.function(cm)) {
    cm <- match.fun(cm)(lev)
  }
  cm <- as.matrix(cm)
  if (nrow(cm) != length(lev)) {
    stop("The contrasts stored for the factor have ", nrow(cm),
         " rows, but it has ", length(lev), " levels.", call. = FALSE)
  }
  rownames(cm) <- lev
  cm
}

# Is `var` in the model formula as a fixed effect, a grouping factor, or a
# threshold group? Read from the formula, not the parameter names, because
# parameter name prefixes are ambiguous (b_r would match b_rater).
.formula_roles <- function(model, var) {
  f <- stats::formula(model)$formula
  rhs <- stats::delete.response(stats::terms(f))
  labs <- attr(rhs, "term.labels")
  bar <- grepl("|", labs, fixed = TRUE)
  groups <- .group_factors(trimws(sub(".*\\|+", "", labs[bar])))
  fixed_terms <- labs[!bar]
  fixed <- any(vapply(strsplit(fixed_terms, ":", fixed = TRUE),
                      function(p) length(p) == 1L && trimws(p) == var,
                      logical(1)))
  random <- var %in% groups
  lhs <- paste(deparse(f[[2L]]), collapse = "")
  # The name is escaped, and must end at `,` or `)`: `rater` must not
  # match `thres(gr = rater.grp)`
  var_re <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", var)
  thres <- grepl(paste0("thres\\(.*gr\\s*=\\s*", var_re, "\\s*[,)]"), lhs)
  list(fixed = fixed, random = random, thres = thres, groups = groups)
}

# Grouping factors of the random-effect terms, as brms expands them.
# `gr(rater, by = g)` gives `rater`, and nesting `id/rater` gives `id` and
# `id:rater`. Spaces are removed, so `id : rater` becomes `id:rater`.
.group_factors <- function(g) {
  g <- gsub("\\s+", "", g)
  g <- sub("^gr\\(([^,)]+).*$", "\\1", g)
  unlist(lapply(strsplit(g, "/", fixed = TRUE), function(p) {
    Reduce(function(a, b) paste0(a, ":", b), p, accumulate = TRUE)
  }), use.names = FALSE)
}

# Rank deficiency of the additive person + item + facet design.
#
# The item and facet effects are separated from the person locations only
# if the design connects them. The check eliminates the person effects from
# the normal equations of the additive design (a Schur complement built from
# count tables) and compares the rank of what remains, for the item and
# facet effects, with its full value (n_items + n_facet - 2). Returns the
# number of missing constraints, 0 for a connected design.
.design_rank_deficiency <- function(person, item, facet) {
  person <- factor(person)
  item <- factor(item)
  facet <- factor(facet)
  n_p <- as.numeric(table(person))
  n_ip <- unclass(table(item, person))
  n_fp <- unclass(table(facet, person))
  n_if <- unclass(table(item, facet))
  a <- rbind(cbind(diag(rowSums(n_if), nrow(n_if)), n_if),
             cbind(t(n_if), diag(colSums(n_if), ncol(n_if))))
  b <- rbind(n_ip, n_fp)
  m <- a - b %*% (t(b) / n_p)
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  tol <- max(dim(m)) * max(abs(ev)) * sqrt(.Machine$double.eps)
  full <- nlevels(item) + nlevels(facet) - 2L
  max(0L, full - sum(ev > tol))
}

# Connected groups of persons and facet elements (who was rated by whom),
# by label propagation on the person-facet graph. Returns a list of facet
# element labels per group.
.person_facet_groups <- function(person, facet) {
  pe <- as.integer(factor(person))
  fe <- as.integer(factor(facet))
  f_lab <- levels(factor(facet))
  p_comp <- seq_len(max(pe))
  repeat {
    f_comp <- as.numeric(tapply(p_comp[pe], fe, min))
    new_p <- pmin(p_comp, as.numeric(tapply(f_comp[fe], pe, min)))
    if (identical(new_p, p_comp)) break
    p_comp <- new_p
  }
  split(f_lab, factor(f_comp, levels = unique(f_comp)))
}

# Collect what the many-facet functions need to know about the facet.
#
# model: brmsfit. *_name: column names in model$data. Returns a list:
#   item_levels, facet_levels: item and facet labels in brms order, used to
#     order and type the item and facet columns of the output.
#   cell: factor of item x facet cells, one value per row of model$data,
#     labelled "item:facet", levels item-major in brms item order.
#   cell_item, cell_facet: item and facet label of each cell level.
#   roles: list(fixed, random, thres) for the facet in the formula.
#   halo: TRUE if the formula has a person x facet grouping term.
#   halo_term: that grouping term as brms expands it (e.g. "id:rater"),
#     NULL without one.
#   rank_deficiency: see .design_rank_deficiency(), NA if the facet is not
#     in the model (nothing to separate).
# Warns when the facet is in the model and the design does not separate
# the item and facet effects from the person locations.
.resolve_facet <- function(model, person_name, item_name, facet_name,
                           fun_name) {
  dat <- model$data
  for (v in c(person_name, item_name, facet_name)) {
    if (!v %in% names(dat)) {
      stop("`", fun_name, "()`: variable '", v, "' not found in the model ",
           "data. Only variables used in the model formula are kept there.",
           call. = FALSE)
    }
  }
  if (length(unique(c(person_name, item_name, facet_name))) < 3L) {
    stop("`", fun_name, "()`: `person_var`, `item_var` and `facet_var` ",
         "must be three different variables.", call. = FALSE)
  }
  # The output has a `cell` column next to the item and facet columns
  if ("cell" %in% c(item_name, facet_name)) {
    stop("`", fun_name, "()`: the output has a column `cell` for the ",
         "item x facet combinations, so `item_var` and `facet_var` cannot ",
         "be called 'cell'. Please rename the variable and refit.",
         call. = FALSE)
  }

  item_lev <- .item_levels(model, item_name)
  facet_lev <- levels(factor(dat[[facet_name]]))
  item_chr <- as.character(dat[[item_name]])
  facet_chr <- as.character(dat[[facet_name]])
  grid <- expand.grid(facet = facet_lev, item = as.character(item_lev),
                      stringsAsFactors = FALSE)
  grid_lab <- paste(grid$item, grid$facet, sep = ":")
  if (anyDuplicated(grid_lab)) {
    dup <- grid_lab[duplicated(grid_lab)][1L]
    stop("`", fun_name, "()`: ", item_name, " and ", facet_name, " labels ",
         "containing ':' give the same cell label '", dup, "' to two ",
         "different combinations. Please rename the levels.", call. = FALSE)
  }
  row_lab <- paste(item_chr, facet_chr, sep = ":")
  used <- grid_lab %in% row_lab
  cell <- factor(row_lab, levels = grid_lab[used])

  roles <- .formula_roles(model, facet_name)
  halo_terms <- intersect(roles$groups,
                          c(paste0(facet_name, ":", person_name),
                            paste0(person_name, ":", facet_name)))
  halo <- length(halo_terms) > 0L
  in_model <- roles$fixed || roles$random || roles$thres

  deficiency <- NA_integer_
  if (in_model) {
    deficiency <- .design_rank_deficiency(dat[[person_name]], item_chr,
                                          facet_chr)
    if (deficiency > 0L) {
      groups <- .person_facet_groups(dat[[person_name]], facet_chr)
      detail <- if (length(groups) > 1L) {
        paste0(" The ", facet_name, " elements fall into ", length(groups),
               " groups that never rate the same persons: ",
               paste0("{", vapply(groups, paste, character(1),
                                  collapse = ", "), "}", collapse = ", "),
               ".")
      } else {
        paste0(" Some ", item_name, " and ", facet_name, " elements are ",
               "only observed together, so their effects cannot be told ",
               "apart.")
      }
      warning("`", fun_name, "()`: the rating design is not connected (",
              deficiency, " missing link", if (deficiency > 1L) "s", ").",
              detail, " Measures are then compared across groups only ",
              "through the prior.", call. = FALSE)
    }
  }

  list(
    item_levels = as.character(item_lev),
    facet_levels = facet_lev,
    cell = cell,
    cell_item = grid$item[used],
    cell_facet = grid$facet[used],
    roles = roles[c("fixed", "random", "thres")],
    halo = halo,
    halo_term = if (halo) halo_terms[[1L]] else NULL,
    rank_deficiency = deficiency
  )
}

# Name of an optional tidy-eval argument, NULL when not given.
.optional_name <- function(quo) {
  if (rlang::quo_is_null(quo)) NULL else rlang::as_name(quo)
}

# Internal column names added to a copy of model$data. The leading dot
# keeps them apart from user variables. brms ignores columns it does not
# use, so predictions are unchanged.
.cell_col <- ".cell"
.unit_col <- ".unit"

# Prepare a many-facet call: resolve the facet, check one response per
# person x item x facet, and add the cell column to a copy of the model.
# Returns list(model, facet) with facet from .resolve_facet().
.facet_setup <- function(model, person_name, item_name, facet_name,
                         fun_name) {
  fr <- .resolve_facet(model, person_name, item_name, facet_name, fun_name)
  .check_unique_person_item(
    model$data[[person_name]], fr$cell, person_name,
    paste(item_name, "x", facet_name), fun_name
  )
  model$data[[.cell_col]] <- fr$cell
  list(model = model, facet = fr)
}

# Message when the model has a person x facet term (halo). The
# conditional method treats each draw's halo effects as known, and the
# theta method's expected values from posterior_epred() include them, so
# with either method fit is assessed after the halo is accounted for.
# `extra` is appended, for a caveat that applies to the calling statistic
# and method only (see .infit_halo_note()).
.facet_halo_message <- function(fr, fun_name, extra = "",
                                method = "conditional") {
  if (isTRUE(fr$halo)) {
    how <- if (method == "conditional") {
      "The conditional method treats its effects as known in each draw"
    } else {
      "The expected values from the posterior include its effects"
    }
    message(
      "`", fun_name, "()`: the model has a (1 | ", fr$halo_term, ") term. ",
      how, ", so the results describe fit after this dependence is ",
      "accounted for.", extra
    )
  }
  invisible(NULL)
}

# Caveat for conditional infit with the halo term modelled. Null run
# 2026-11-15 (dev/mfrm-design-plan.md, step 5), rejection rate at .05 for
# fitting units: about .02 per cell, .009-.019 pooled by facet, .004-.013
# pooled by item. Restscore and Q3 stayed close to .05 and get no caveat.
.infit_halo_note <- function(by, item_name, facet_name) {
  what <- switch(by,
                 cell = "about 2% of fitting cells",
                 facet = paste0("about 1-2% of fitting levels of ",
                                facet_name),
                 item = paste0("about 1% of fitting levels of ", item_name))
  paste0(" In simulations infit then flagged misfit less often than the ",
         "nominal level (", what, " at .05 instead of 5%), so some misfit ",
         "may be missed.")
}

# Replace the internal cell column of a statistic's output with columns
# `cell` (label "item:facet"), item_name and facet_name, in that order and
# in place of the cell column. All three are factors in brms order.
.label_cells <- function(df, fr, item_name, facet_name, col = .cell_col) {
  pos <- match(col, names(df))
  if (is.na(pos)) return(df)
  lab <- as.character(df[[col]])
  idx <- match(lab, levels(fr$cell))
  cell <- factor(lab, levels = levels(fr$cell))
  extra <- tibble::tibble(
    cell = cell,
    i = factor(fr$cell_item[idx], levels = fr$item_levels),
    f = factor(fr$cell_facet[idx], levels = fr$facet_levels)
  )
  names(extra)[2:3] <- c(item_name, facet_name)
  grouped <- dplyr::is_grouped_df(df)
  df <- dplyr::ungroup(df)
  out <- dplyr::bind_cols(df[seq_len(pos - 1L)], extra,
                          df[-seq_len(pos)])
  if (grouped) out <- dplyr::group_by(out, .data$cell)
  out
}

# Many-facet branch of infit_statistic().
.infit_facets <- function(model, resp_var, item_name, person_name,
                          facet_name, ndraws_use, outfit, method, by) {
  fs <- .facet_setup(model, person_name, item_name, facet_name,
                     "infit_statistic")
  fr <- fs$facet
  out_name <- switch(by, cell = .cell_col, item = item_name,
                     facet = facet_name)
  # The infit caveat comes from a null run with the conditional method
  .facet_halo_message(fr, "infit_statistic",
                      if (method == "conditional")
                        .infit_halo_note(by, item_name, facet_name) else "",
                      method = method)
  if (method == "conditional") {
    group_of <- switch(
      by,
      cell = NULL,
      item = function(l) fr$cell_item[match(as.character(l), levels(fr$cell))],
      facet = function(l) fr$cell_facet[match(as.character(l), levels(fr$cell))]
    )
    res <- .infit_conditional(fs$model, resp_var, .cell_col, person_name,
                              ndraws_use, outfit, group_of = group_of,
                              out_name = out_name)
  } else {
    # Theta infit is a sum over rows, so pooling over cells is the same as
    # grouping the rows by item or facet directly.
    res <- infit_statistic(fs$model, item_var = !!rlang::sym(out_name),
                           person_var = !!rlang::sym(person_name),
                           ndraws_use = ndraws_use, outfit = outfit,
                           method = "theta")
  }
  if (by == "cell") return(.label_cells(res, fr, item_name, facet_name))
  # Same type and order from both methods: a factor in brms order (the
  # conditional engine returns sorted character labels, the theta branch
  # the type of the data column)
  lev <- if (by == "item") fr$item_levels else fr$facet_levels
  res <- dplyr::ungroup(res)
  res[[out_name]] <- factor(as.character(res[[out_name]]), levels = lev)
  res <- dplyr::arrange(res, .data[[out_name]], .data$draw)
  dplyr::group_by(res, .data[[out_name]])
}

# Add a pair_type column to many-facet Q3 output, after item_2. Cell labels
# are mapped back to their item and facet element.
.q3_pair_types <- function(out, fr, item_name, facet_name) {
  lev <- levels(fr$cell)
  i1 <- match(as.character(out$item_1), lev)
  i2 <- match(as.character(out$item_2), lev)
  same_facet <- fr$cell_facet[i1] == fr$cell_facet[i2]
  same_item <- fr$cell_item[i1] == fr$cell_item[i2]
  pair_type <- ifelse(same_facet, paste("same", facet_name),
                      ifelse(same_item, paste("same", item_name), "neither"))
  pos <- match("item_2", names(out))
  dplyr::bind_cols(out[seq_len(pos)], tibble::tibble(pair_type = pair_type),
                   out[-seq_len(pos)])
}
