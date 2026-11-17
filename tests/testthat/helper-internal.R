# Internal helpers used only by the tests. They expose intermediate
# results of internal functions in a form that is easy to check.

# The same groups as .cond_weight_groups(), with W as a full U x k x C
# array (missing items get weights (1, 0, ..., 0), which force a score of
# 0). Used in tests, where k is small.
.cond_item_weights <- function(prob, person_idx, item_idx, n_persons, k) {
  cg <- .cond_weight_groups(prob, person_idx, item_idx, n_persons, k)
  n_cat <- ncol(prob)
  W <- array(0, dim = c(length(cg$W), k, n_cat))
  W[, , 1L] <- 1
  for (u in seq_along(cg$W)) {
    W[u, cg$items[[u]], ] <- cg$W[[u]]
  }
  list(key = cg$key, W = W)
}

# Contrast matrix (levels x columns) that model.matrix() used for a factor:
# the "contrasts" attribute if set, otherwise the contrasts option
# (treatment for unordered, polynomial for ordered factors). The option is
# read now, so a model fitted under a different option is only covered
# when the contrasts were stored on the factor.
.factor_contrast_matrix <- function(x, lev) {
  .contrast_from_spec(.contrast_spec(x), lev)
}
