#' Tile Plot of Item Response Distributions
#'
#' Creates a tile (heat map) plot showing the distribution of
#' responses across all items and response categories. Each cell
#' displays the count (or percentage) of responses, with optional
#' conditional highlighting for cells with low counts. This is a
#' descriptive data visualization tool intended for use before
#' model fitting.
#'
#' @param data A data frame in wide format containing only the item
#'   response columns. Each column is one item, each row is one
#'   person. All columns must be numeric (integer-valued). Response
#'   categories may be coded starting from 0 or 1. Do not include
#'   person IDs, grouping variables, or other non-item columns.
#' @param cutoff Integer. Cells with counts below this value are
#'   highlighted (when \code{highlight = TRUE}). Default is 10.
#' @param highlight Logical. If \code{TRUE} (the default), cell
#'   labels with counts below \code{cutoff} are displayed in red.
#'   This includes cells with zero responses (empty categories),
#'   which is useful for identifying gaps in the response
#'   distribution.
#' @param percent Logical. If \code{TRUE}, cell labels show
#'   percentages instead of raw counts. Default is \code{FALSE}.
#' @param text_color Character. Color for cell label text (when
#'   not highlighted). Default is \code{"orange"}.
#' @param item_labels An optional character vector of descriptive
#'   labels for the items (y-axis). Must be the same length as
#'   \code{ncol(data)}. If \code{NULL} (the default), column names
#'   are used.
#' @param category_labels An optional character vector of labels
#'   for the response categories (x-axis). Must be the same length
#'   as the number of response categories spanning from the minimum
#'   to the maximum observed value. If \code{NULL} (the default),
#'   numeric category values are used.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @details
#' The plot displays items on the y-axis (in the same order as the
#' columns in \code{data}, from top to bottom) and response
#' categories on the x-axis. Cell shading represents the count
#' of responses (darker = more responses). Cell labels show either
#' raw counts or percentages.
#'
#' Categories with zero responses are explicitly shown (as cells
#' with \code{n = 0}), which helps identify gaps in the response
#' distribution — one of the primary purposes of this plot.
#'
#' \strong{Input requirements:}
#' \itemize{
#'   \item All columns must be numeric (integer-valued).
#'   \item The data frame must contain at least 2 columns (items)
#'     and at least 1 row (person).
#' }
#'
#' @examples
#' {
#' library(ggplot2)
#'
#' # Basic tile plot
#' plot_tile(eRm::pcmdat2)
#'
#' # With custom item labels
#' plot_tile(
#'   eRm::pcmdat2,
#'   item_labels = c("Mood", "Sleep", "Appetite", "Energy")
#' )
#'
#' # With custom category labels and percentages
#' plot_tile(
#'   eRm::pcmdat2,
#'   category_labels = c("Never", "Sometimes", "Often"),
#'   percent = TRUE
#' )
#'
#' # Adjust cutoff for highlighting
#' plot_tile(eRm::pcmdat2, cutoff = 20, highlight = TRUE)
#' }
#'
#' @importFrom rlang .data
#' @export
plot_tile <- function(
    data,
    cutoff          = 10,
    highlight       = TRUE,
    percent         = FALSE,
    text_color      = "orange",
    item_labels     = NULL,
    category_labels = NULL
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  
  # ================================================================
  # INPUT VALIDATION
  # ================================================================
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame in wide format ",
         "(one column per item, one row per person).",
         call. = FALSE)
  }
  
  if (ncol(data) < 2) {
    stop("'data' must contain at least 2 columns (items). ",
         "Found ", ncol(data), " column(s).", call. = FALSE)
  }
  
  if (nrow(data) < 1) {
    stop("'data' must contain at least 1 row (person). ",
         "Found 0 rows.", call. = FALSE)
  }
  
  # Check all columns are numeric
  non_numeric <- names(data)[!vapply(data, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("All columns must be numeric. Non-numeric column(s): ",
         paste(non_numeric, collapse = ", "), ". ",
         "Remove non-item columns (e.g., person IDs, grouping ",
         "variables) before passing to plot_tile().",
         call. = FALSE)
  }
  
  # Check values are integer-valued
  all_vals <- unlist(data, use.names = FALSE)
  all_vals <- all_vals[!is.na(all_vals)]
  if (length(all_vals) == 0) {
    stop("'data' contains only NA values.", call. = FALSE)
  }
  if (!all(all_vals == round(all_vals))) {
    stop("All response values must be integers. ",
         "Found non-integer values in the data.", call. = FALSE)
  }
  
  # Validate item_labels
  if (!is.null(item_labels)) {
    if (length(item_labels) != ncol(data)) {
      stop("'item_labels' must have the same length as the number ",
           "of columns in 'data'. Expected ", ncol(data),
           ", got ", length(item_labels), ".", call. = FALSE)
    }
  }
  
  # Determine category range across all items
  min_val <- min(all_vals)
  max_val <- max(all_vals)
  all_categories <- seq(min_val, max_val)
  n_categories <- length(all_categories)
  
  # Validate category_labels
  if (!is.null(category_labels)) {
    if (length(category_labels) != n_categories) {
      stop("'category_labels' must have the same length as the ",
           "number of response categories (", min_val, " to ",
           max_val, " = ", n_categories, " categories). Got ",
           length(category_labels), ".", call. = FALSE)
    }
  }
  
  # ================================================================
  # PREPARE DATA
  # ================================================================
  item_names <- names(data)
  if (!is.null(item_labels)) {
    label_map <- stats::setNames(item_labels, item_names)
  } else {
    label_map <- stats::setNames(item_names, item_names)
  }
  
  # Count responses per item x category, including zeros
  count_list <- list()
  for (item in item_names) {
    item_vals <- data[[item]]
    item_vals <- item_vals[!is.na(item_vals)]
    tab <- table(factor(item_vals, levels = all_categories))
    count_list[[item]] <- data.frame(
      item     = item,
      category = as.integer(names(tab)),
      n        = as.integer(tab),
      stringsAsFactors = FALSE
    )
  }
  count_df <- do.call(rbind, count_list)
  rownames(count_df) <- NULL
  
  # Compute percentage within each item
  item_totals <- stats::aggregate(n ~ item, data = count_df, FUN = sum)
  colnames(item_totals)[2] <- "total"
  count_df <- merge(count_df, item_totals, by = "item")
  count_df$percentage <- round(count_df$n / count_df$total * 100, 1)
  
  # Apply item labels
  count_df$item_label <- label_map[count_df$item]
  
  # Set factor levels: items in reverse column order (so first item
  # appears at the top of the y-axis)
  count_df$item_label <- factor(
    count_df$item_label,
    levels = rev(label_map[item_names])
  )
  
  # ================================================================
  # BUILD PLOT
  # ================================================================
  p <- ggplot2::ggplot(
    count_df,
    ggplot2::aes(
      x = .data$category,
      y = .data$item_label,
      fill = .data$n
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(
      expression(italic(n)),
      limits = c(0, NA)
    ) +
    ggplot2::labs(y = "Items") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8),
      panel.grid  = ggplot2::element_blank()
    )
  
  # X-axis: category labels or numeric breaks
  if (!is.null(category_labels)) {
    p <- p + ggplot2::scale_x_continuous(
      "Response category",
      expand = c(0, 0),
      breaks = all_categories,
      labels = category_labels
    )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      "Response category",
      expand = c(0, 0),
      breaks = all_categories
    )
  }
  
  # Cell labels: count or percentage, with optional highlighting
  if (percent) {
    if (highlight) {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(
            label = paste0(.data$percentage, "%"),
            color = ifelse(.data$n < cutoff, "red", text_color)
          )
        ) +
        ggplot2::guides(color = "none") +
        ggplot2::scale_color_identity()
    } else {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = paste0(.data$percentage, "%")),
          color = text_color
        )
    }
  } else {
    if (highlight) {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(
            label = .data$n,
            color = ifelse(.data$n < cutoff, "red", text_color)
          )
        ) +
        ggplot2::guides(color = "none") +
        ggplot2::scale_color_identity()
    } else {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = .data$n),
          color = text_color
        )
    }
  }
  
  p
}

#' Item Response Distribution Bar Chart
#'
#' Creates a faceted bar chart showing the response distribution for
#' each item, with counts and percentages displayed on each bar.
#' Each item gets its own panel, with response categories on the
#' x-axis and percentage of responses on the y-axis. This is a
#' descriptive data visualization tool intended for use before
#' model fitting.
#'
#' @param data A data frame in wide format containing only the item
#'   response columns. Each column is one item, each row is one
#'   person. All columns must be numeric (integer-valued). Response
#'   categories may be coded starting from 0 or 1. Do not include
#'   person IDs, grouping variables, or other non-item columns.
#' @param item_labels An optional character vector of descriptive
#'   labels for the items (facet strips). Must be the same length
#'   as \code{ncol(data)}. If \code{NULL} (the default), column
#'   names are used. Labels are displayed as
#'   \code{"column_name - label"}.
#' @param category_labels An optional character vector of labels
#'   for the response categories (x-axis). Must be the same length
#'   as the number of response categories spanning from the minimum
#'   to the maximum observed value. If \code{NULL} (the default),
#'   numeric category values are used.
#' @param ncol Integer. Number of columns in the faceted layout.
#'   Default is 1.
#' @param label_wrap Integer. Number of characters per line in
#'   facet strip labels before wrapping. Default is 25.
#' @param text_y Numeric. Vertical position (in percent units) for
#'   the count labels on each bar. Adjust upward if bars are tall.
#'   Default is 6.
#' @param viridis_option Character. Viridis palette option for the
#'   count text color. One of \code{"A"} through \code{"H"}.
#'   Default is \code{"A"}.
#' @param viridis_end Numeric in \eqn{[0, 1]}. End point of the
#'   viridis color scale for count text. Adjust if text is hard to
#'   read against the bar colors. Default is 0.9.
#' @param font Character. Font family for all text. Default is
#'   \code{"sans"}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @details
#' Each item is displayed as a separate facet panel with the item
#' label in the strip on the left side. Bars are colored by response
#' category using the viridis palette. Each bar shows the count
#' (\code{n = X}) as text.
#'
#' \strong{Input requirements:}
#' \itemize{
#'   \item All columns must be numeric (integer-valued).
#'   \item The data frame must contain at least 2 columns (items)
#'     and at least 1 row (person).
#' }
#'
#' @examples
#' {
#' library(ggplot2)
#'
#' # Basic response distribution plot
#' plot_bars(eRm::pcmdat2)
#'
#' # With custom item labels
#' plot_bars(
#'   eRm::pcmdat2,
#'   item_labels = c("Mood", "Sleep", "Appetite", "Energy")
#' )
#'
#' # Two-column layout with wrapped labels
#' plot_bars(
#'   eRm::pcmdat2,
#'   item_labels = c(
#'     "General mood and emotional wellbeing",
#'     "Quality of sleep at night",
#'     "Appetite and eating habits",
#'     "Overall energy level during the day"
#'   ),
#'   ncol = 2, label_wrap = 20
#' )
#'
#' # With custom category labels
#' plot_bars(
#'   eRm::pcmdat2,
#'   category_labels = c("Never", "Sometimes", "Often")
#' )
#' }
#'
#' @importFrom rlang .data
#' @export
plot_bars <- function(
    data,
    item_labels     = NULL,
    category_labels = NULL,
    ncol            = 1,
    label_wrap      = 25,
    text_y          = 6,
    viridis_option  = "A",
    viridis_end     = 0.9,
    font            = "sans"
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  
  # ================================================================
  # INPUT VALIDATION
  # ================================================================
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame in wide format ",
         "(one column per item, one row per person).",
         call. = FALSE)
  }
  
  if (ncol(data) < 2) {
    stop("'data' must contain at least 2 columns (items). ",
         "Found ", ncol(data), " column(s).", call. = FALSE)
  }
  
  if (nrow(data) < 1) {
    stop("'data' must contain at least 1 row (person). ",
         "Found 0 rows.", call. = FALSE)
  }
  
  non_numeric <- names(data)[!vapply(data, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("All columns must be numeric. Non-numeric column(s): ",
         paste(non_numeric, collapse = ", "), ". ",
         "Remove non-item columns (e.g., person IDs, grouping ",
         "variables) before passing to plot_item_responses().",
         call. = FALSE)
  }
  
  all_vals <- unlist(data, use.names = FALSE)
  all_vals <- all_vals[!is.na(all_vals)]
  if (length(all_vals) == 0) {
    stop("'data' contains only NA values.", call. = FALSE)
  }
  if (!all(all_vals == round(all_vals))) {
    stop("All response values must be integers. ",
         "Found non-integer values in the data.", call. = FALSE)
  }
  
  item_names <- names(data)
  
  if (!is.null(item_labels)) {
    if (length(item_labels) != length(item_names)) {
      stop("'item_labels' must have the same length as the number ",
           "of columns in 'data'. Expected ", length(item_names),
           ", got ", length(item_labels), ".", call. = FALSE)
    }
  }
  
  min_val <- min(all_vals)
  max_val <- max(all_vals)
  all_categories <- seq(min_val, max_val)
  n_categories <- length(all_categories)
  
  if (!is.null(category_labels)) {
    if (length(category_labels) != n_categories) {
      stop("'category_labels' must have the same length as the ",
           "number of response categories (", min_val, " to ",
           max_val, " = ", n_categories, " categories). Got ",
           length(category_labels), ".", call. = FALSE)
    }
  }
  
  # ================================================================
  # PREPARE DATA
  # ================================================================
  # Build facet labels: "itemnr - description" or just "itemnr"
  if (!is.null(item_labels)) {
    facet_labels <- paste0(item_names, " - ", item_labels)
  } else {
    facet_labels <- item_names
  }
  
  # Count responses per item x category, including zeros
  count_list <- list()
  for (i in seq_along(item_names)) {
    item <- item_names[i]
    item_vals <- data[[item]]
    item_vals <- item_vals[!is.na(item_vals)]
    tab <- table(factor(item_vals, levels = all_categories))
    total <- sum(tab)
    
    count_list[[item]] <- data.frame(
      item_name   = item,
      facet_label = facet_labels[i],
      category    = as.integer(names(tab)),
      n           = as.integer(tab),
      percent     = if (total > 0) {
        round(as.integer(tab) / total * 100, 1)
      } else {
        rep(0, length(tab))
      },
      stringsAsFactors = FALSE
    )
  }
  plot_df <- do.call(rbind, count_list)
  rownames(plot_df) <- NULL
  
  # Factor levels preserve column order
  plot_df$facet_label <- factor(
    plot_df$facet_label,
    levels = facet_labels
  )
  plot_df$category_fct <- factor(
    plot_df$category,
    levels = all_categories
  )
  
  # Apply category labels to the factor if provided
  if (!is.null(category_labels)) {
    levels(plot_df$category_fct) <- category_labels
  }
  
  # ================================================================
  # BUILD PLOT
  # ================================================================
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$category_fct, y = .data$percent)
  ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$category_fct)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0("n = ", .data$n),
        color = .data$category_fct
      ),
      y = text_y
    ) +
    ggplot2::facet_wrap(
      ~ facet_label,
      ncol           = ncol,
      strip.position = "left",
      labeller       = ggplot2::labeller(
        facet_label = ggplot2::label_wrap_gen(label_wrap)
      )
    ) +
    ggplot2::scale_fill_viridis_d(guide = "none") +
    ggplot2::scale_color_viridis_d(
      guide     = "none",
      direction = -1,
      option    = viridis_option,
      end       = viridis_end
    ) +
    ggplot2::scale_y_continuous(position = "right") +
    ggplot2::labs(
      x = "Response category",
      y = "% of responses"
    ) +
    ggplot2::theme_bw(base_family = font) +
    ggplot2::theme(
      legend.position    = "none",
      strip.text.y.left  = ggplot2::element_text(angle = 0)
    )
  
  p
}

#' Stacked Bar Chart of Item Response Distributions
#'
#' Creates a horizontal stacked bar chart showing the response
#' distribution for all items. Each bar represents one item, with
#' segments colored by response category. Counts are displayed as
#' text labels within each segment. This is a descriptive data
#' visualization tool intended for use before model fitting.
#'
#' @param data A data frame in wide format containing only the item
#'   response columns. Each column is one item, each row is one
#'   person. All columns must be numeric (integer-valued). Response
#'   categories may be coded starting from 0 or 1. Do not include
#'   person IDs, grouping variables, or other non-item columns.
#' @param item_labels An optional character vector of descriptive
#'   labels for the items (y-axis). Must be the same length as
#'   \code{ncol(data)}. If \code{NULL} (the default), column names
#'   are used.
#' @param category_labels An optional character vector of labels
#'   for the response categories (legend). Must be the same length
#'   as the number of response categories spanning from the minimum
#'   to the maximum observed value, ordered from lowest to highest
#'   category. If \code{NULL} (the default), numeric category
#'   values are used.
#' @param show_n Logical. If \code{TRUE} (the default), the count
#'   of responses is displayed as a text label inside each bar
#'   segment.
#' @param show_percent Logical. If \code{TRUE}, the percentage of
#'   responses is displayed instead of (or in addition to) counts.
#'   Default is \code{FALSE}.
#' @param text_color Character. Color for the count/percentage
#'   labels. Default is \code{"sienna1"}.
#' @param text_size Numeric. Size of the count/percentage labels.
#'   Default is 3.
#' @param min_label_n Integer. Minimum count required for a label
#'   to be displayed within a bar segment. Segments with fewer
#'   responses are left unlabelled to avoid clutter. Default is 0
#'   (all segments labelled).
#' @param viridis_option Character. Viridis palette option. One of
#'   \code{"A"} through \code{"H"}. Default is \code{"D"}.
#' @param viridis_end Numeric in \eqn{[0, 1]}. End point of the
#'   viridis color scale. Default is 0.99.
#' @param title Character. Plot title. Default is
#'   \code{"Item responses"}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @details
#' Items are displayed on the y-axis in the same order as the
#' columns in \code{data} (first column at the top). Each bar is
#' divided into segments representing response categories, with
#' the lowest category on the left and the highest on the right.
#' The total bar length equals the number of non-missing responses
#' for that item.
#'
#' Categories with zero responses still appear in the legend but
#' produce no visible bar segment, which helps identify gaps in
#' the response distribution.
#'
#' \strong{Input requirements:}
#' \itemize{
#'   \item All columns must be numeric (integer-valued).
#'   \item The data frame must contain at least 2 columns (items)
#'     and at least 1 row (person).
#' }
#'
#' @examples
#' {
#' library(ggplot2)
#'
#' # Basic stacked bar chart
#' plot_stackedbars(eRm::pcmdat2)
#'
#' # With custom item and category labels
#' plot_stackedbars(
#'   eRm::pcmdat2,
#'   item_labels = c("Mood", "Sleep", "Appetite", "Energy"),
#'   category_labels = c("Never", "Sometimes", "Often")
#' )
#'
#' # Show percentages, suppress small segments
#' plot_stackedbars(
#'   eRm::pcmdat2,
#'   show_percent = TRUE,
#'   show_n       = FALSE,
#'   min_label_n  = 5
#' )
#' }
#'
#' @importFrom rlang .data
#' @export
plot_stackedbars <- function(
    data,
    item_labels     = NULL,
    category_labels = NULL,
    show_n          = TRUE,
    show_percent    = FALSE,
    text_color      = "sienna1",
    text_size       = 3,
    min_label_n     = 0,
    viridis_option  = "D",
    viridis_end     = 0.99,
    title           = "Item responses"
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  
  # ================================================================
  # INPUT VALIDATION
  # ================================================================
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame in wide format ",
         "(one column per item, one row per person).",
         call. = FALSE)
  }
  
  if (ncol(data) < 2) {
    stop("'data' must contain at least 2 columns (items). ",
         "Found ", ncol(data), " column(s).", call. = FALSE)
  }
  
  if (nrow(data) < 1) {
    stop("'data' must contain at least 1 row (person). ",
         "Found 0 rows.", call. = FALSE)
  }
  
  non_numeric <- names(data)[!vapply(data, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("All columns must be numeric. Non-numeric column(s): ",
         paste(non_numeric, collapse = ", "), ". ",
         "Remove non-item columns (e.g., person IDs, grouping ",
         "variables) before passing to plot_item_stacked().",
         call. = FALSE)
  }
  
  all_vals <- unlist(data, use.names = FALSE)
  all_vals <- all_vals[!is.na(all_vals)]
  if (length(all_vals) == 0) {
    stop("'data' contains only NA values.", call. = FALSE)
  }
  if (!all(all_vals == round(all_vals))) {
    stop("All response values must be integers. ",
         "Found non-integer values in the data.", call. = FALSE)
  }
  
  item_names <- names(data)
  
  if (!is.null(item_labels)) {
    if (length(item_labels) != length(item_names)) {
      stop("'item_labels' must have the same length as the number ",
           "of columns in 'data'. Expected ", length(item_names),
           ", got ", length(item_labels), ".", call. = FALSE)
    }
  }
  
  min_val <- min(all_vals)
  max_val <- max(all_vals)
  all_categories <- seq(min_val, max_val)
  n_categories <- length(all_categories)
  
  if (!is.null(category_labels)) {
    if (length(category_labels) != n_categories) {
      stop("'category_labels' must have the same length as the ",
           "number of response categories (", min_val, " to ",
           max_val, " = ", n_categories, " categories). Got ",
           length(category_labels), ".", call. = FALSE)
    }
  }
  
  # ================================================================
  # PREPARE DATA
  # ================================================================
  # Build y-axis labels
  if (!is.null(item_labels)) {
    y_labels <- stats::setNames(item_labels, item_names)
  } else {
    y_labels <- stats::setNames(item_names, item_names)
  }
  
  # Count responses per item x category, including zeros
  count_list <- list()
  for (item in item_names) {
    item_vals <- data[[item]]
    item_vals <- item_vals[!is.na(item_vals)]
    tab <- table(factor(item_vals, levels = all_categories))
    total <- sum(tab)
    
    count_list[[item]] <- data.frame(
      item_name = item,
      category  = as.integer(names(tab)),
      n         = as.integer(tab),
      percent   = if (total > 0) {
        round(as.integer(tab) / total * 100, 1)
      } else {
        rep(0, length(tab))
      },
      stringsAsFactors = FALSE
    )
  }
  plot_df <- do.call(rbind, count_list)
  rownames(plot_df) <- NULL
  
  # Item factor: reversed so first column appears at top of y-axis
  plot_df$item_fct <- factor(
    plot_df$item_name,
    levels = rev(item_names),
    labels = rev(y_labels[item_names])
  )
  
  # Category factor: reversed for stacking order (lowest category
  # on the left = bottom of stack, which means it must be the
  # last factor level in geom_col's default stacking)
  plot_df$category_fct <- factor(
    plot_df$category,
    levels = rev(all_categories)
  )
  
  # Build text labels
  plot_df$label <- ""
  show_mask <- plot_df$n >= min_label_n
  
  if (show_n && show_percent) {
    plot_df$label[show_mask] <- paste0(
      plot_df$n[show_mask], " (", plot_df$percent[show_mask], "%)"
    )
  } else if (show_n) {
    plot_df$label[show_mask] <- as.character(plot_df$n[show_mask])
  } else if (show_percent) {
    plot_df$label[show_mask] <- paste0(
      plot_df$percent[show_mask], "%"
    )
  }
  
  # Suppress labels for zero-count segments (no visible bar)
  plot_df$label[plot_df$n == 0] <- ""
  
  # ================================================================
  # BUILD LEGEND LABELS
  # ================================================================
  # Legend should show categories from lowest to highest (top to
  # bottom), which is the reverse of the stacking factor
  if (!is.null(category_labels)) {
    legend_labels <- stats::setNames(
      rev(category_labels),
      rev(as.character(all_categories))
    )
  } else {
    legend_labels <- stats::setNames(
      rev(as.character(all_categories)),
      rev(as.character(all_categories))
    )
  }
  
  # ================================================================
  # BUILD PLOT
  # ================================================================
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x    = .data$n,
      y    = .data$item_fct,
      fill = .data$category_fct
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_viridis_d(
      "Response\ncategory",
      direction = -1,
      option    = viridis_option,
      end       = viridis_end,
      labels    = legend_labels
    ) +
    ggplot2::labs(
      title = title,
      x     = "Number of responses",
      y     = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right"
    )
  
  # Add text labels if any content
  if (show_n || show_percent) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        color    = text_color,
        size     = text_size,
        position = ggplot2::position_stack(vjust = 0.5)
      )
  }
  
  p
}