# Internal input checks shared across functions.

# Stop if any person x item combination occurs more than once in the model
# data. Functions that place residuals in a person x item matrix (Q3,
# residual PCA) would otherwise silently keep only the last response of
# each duplicated cell. This happens, for example, with many-facet data
# where each rater x ratee combination has several criterion ratings.
.check_unique_person_item <- function(persons, items, person_name, item_name,
                                      fun_name) {
  key <- paste(as.character(persons), as.character(items), sep = "\r")
  dup <- duplicated(key)
  if (any(dup)) {
    n_cells <- length(unique(key[dup]))
    first <- which(key == key[which(dup)[1]])
    stop(
      "`", fun_name, "()` requires one response per ", person_name, " x ",
      item_name, " combination, but ", n_cells,
      " combination(s) occur more than once (for example ", person_name,
      " '", persons[first[1]], "' and ", item_name, " '", items[first[1]],
      "' appear ", length(first), " times). ",
      "Check `person_var` and `item_var`. With many-facet data, ",
      "`person_var` must identify a single response per item.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
