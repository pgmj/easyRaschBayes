# Internal helpers for plot_targeting().
#
# The category-band geometry is ported from easyRasch2::RMtargeting()
# (R/targeting_plot.R, easyRasch2 1.4.0.9000) so the two packages draw the
# same figure. Additions here: posterior probabilities for collapsed
# categories and disordered thresholds, computed draw by draw, and plain-text
# captions (no ggtext dependency).

# ── Category geometry ───────────────────────────────────────────

# Locations where each category is the most likely response, for one item.
# The boundaries are the slopes of the lower convex hull of
# {(k, S_k)}, S_k the cumulative thresholds. Ordered thresholds give the
# thresholds themselves; a disordered run is pooled by averaging, and the
# categories it skips are never the most likely response.
# Returns list(cats, dropped, bounds).
#' @keywords internal
#' @noRd
.modal_boundaries <- function(tau) {
  s_cum <- c(0, cumsum(tau))
  k_idx <- seq_along(s_cum) - 1L
  hull <- 1L
  for (i in 2:length(s_cum)) {
    while (length(hull) >= 2L) {
      a <- hull[length(hull) - 1L]
      b <- hull[length(hull)]
      above <- (s_cum[b] - s_cum[a]) * (k_idx[i] - k_idx[a]) >=
        (s_cum[i] - s_cum[a]) * (k_idx[b] - k_idx[a])
      if (above) hull <- hull[-length(hull)] else break
    }
    hull <- c(hull, i)
  }
  vert <- k_idx[hull]
  list(
    cats    = vert,
    dropped = setdiff(k_idx, vert),
    bounds  = diff(s_cum[hull]) / diff(k_idx[hull])
  )
}

# Category bands and collapsed categories for every item.
# thr: data.frame with Item, k (threshold index) and Location.
#' @keywords internal
#' @noRd
.category_bands <- function(thr, item_levels, xlim) {
  ord <- order(thr$Item, thr$k)
  taus <- split(thr$Location[ord],
                factor(thr$Item[ord], levels = item_levels))
  bands <- list()
  drops <- list()
  for (it in item_levels) {
    tau <- as.numeric(taus[[it]])
    mb <- .modal_boundaries(tau)
    bands[[it]] <- data.frame(
      Item = it,
      cat  = mb$cats,
      xmin = c(xlim[1], mb$bounds),
      xmax = c(mb$bounds, xlim[2]),
      stringsAsFactors = FALSE
    )
    if (length(mb$dropped)) {
      s_cum <- c(0, cumsum(tau))
      drops[[it]] <- data.frame(
        Item = it,
        cat  = mb$dropped,
        x    = vapply(mb$dropped, function(cc) {
          lo <- max(mb$cats[mb$cats < cc])
          hi <- min(mb$cats[mb$cats > cc])
          (s_cum[hi + 1L] - s_cum[lo + 1L]) / (hi - lo)
        }, numeric(1)),
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    bands = do.call(rbind, bands),
    drops = if (length(drops)) do.call(rbind, drops) else NULL
  )
}

# Disordered adjacent threshold pairs (point estimates).
#' @keywords internal
#' @noRd
.threshold_reversals <- function(thr, item_levels) {
  out <- lapply(item_levels, function(it) {
    tau <- thr$Location[thr$Item == it][order(thr$k[thr$Item == it])]
    d <- which(diff(tau) < 0)
    if (!length(d)) return(NULL)
    data.frame(
      Item = it, k = d, lo = tau[d + 1L], hi = tau[d],
      gap = tau[d] - tau[d + 1L],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out) || !nrow(out)) NULL else out
}

# ── Posterior summaries of category structure ───────────────────

# Per-draw category structure for every item.
# tau: named list of draws x thresholds matrices (already centred).
# Returns list(categories, disordered):
#   categories: item, category, p_never_modal, and the lower boundary of the
#     category's modal region (median and interval over the draws where the
#     category is modal; NA for category 0, which has no lower boundary).
#   disordered: item, thresholds ("k-(k+1)"), p_disordered.
#' @keywords internal
#' @noRd
.category_posterior <- function(tau, prob) {
  lo_q <- (1 - prob) / 2
  cats_out <- list()
  revs_out <- list()
  for (it in names(tau)) {
    m <- tau[[it]]
    n_cat <- ncol(m) + 1L
    dropped <- matrix(FALSE, nrow = nrow(m), ncol = n_cat)
    lower <- matrix(NA_real_, nrow = nrow(m), ncol = n_cat)
    for (s in seq_len(nrow(m))) {
      mb <- .modal_boundaries(m[s, ])
      dropped[s, mb$dropped + 1L] <- TRUE
      lower[s, mb$cats[-1] + 1L] <- mb$bounds
    }
    q <- function(x, p) {
      x <- x[!is.na(x)]
      if (length(x)) unname(stats::quantile(x, p)) else NA_real_
    }
    cats_out[[it]] <- tibble::tibble(
      item          = it,
      category      = seq_len(n_cat) - 1L,
      p_never_modal = colMeans(dropped),
      lower_bound   = apply(lower, 2, q, p = 0.5),
      lower_bound_lo = apply(lower, 2, q, p = lo_q),
      lower_bound_hi = apply(lower, 2, q, p = 1 - lo_q)
    )
    if (ncol(m) > 1L) {
      rev_draws <- m[, -1, drop = FALSE] < m[, -ncol(m), drop = FALSE]
      revs_out[[it]] <- tibble::tibble(
        item       = it,
        thresholds = paste0(seq_len(ncol(m) - 1L), "-", seq_len(ncol(m) - 1L) + 1L),
        p_disordered = colMeans(rev_draws)
      )
    }
  }
  list(
    categories = dplyr::bind_rows(cats_out),
    disordered = if (length(revs_out)) dplyr::bind_rows(revs_out) else
      tibble::tibble(item = character(), thresholds = character(),
                     p_disordered = numeric())
  )
}

# ".83", "1.00", "<.01" style probability labels.
#' @keywords internal
#' @noRd
.fmt_prob <- function(p) {
  out <- sub("^0\\.", ".", sprintf("%.2f", p))
  out[p > 0 & p < 0.005] <- "<.01"
  out
}

# ── Colour and caption helpers ──────────────────────────────────

# Darken pale colours until they read as a thin line on white.
#' @keywords internal
#' @noRd
.ci_contrast <- function(cols, floor = 130) {
  vapply(cols, function(cl) {
    v <- as.numeric(grDevices::col2rgb(cl)[, 1])
    while (sum(v * c(0.299, 0.587, 0.114)) > floor) v <- v * 0.78
    grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
  }, character(1), USE.NAMES = FALSE)
}

# Readable label colour for text on a filled band.
#' @keywords internal
#' @noRd
.band_label_colour <- function(fills) {
  lum <- apply(grDevices::col2rgb(fills), 2, function(v) {
    sum(v * c(0.299, 0.587, 0.114))
  })
  ifelse(lum > 145, "grey15", "white")
}

# Wrapped plain-text figure note.
#' @keywords internal
#' @noRd
.erb_caption <- function(text, width = 95L) {
  paste(strwrap(paste0("Note. ", text), width = width), collapse = "\n")
}

#' @keywords internal
#' @noRd
.erb_caption_theme <- function() {
  ggplot2::theme(
    plot.caption = ggplot2::element_text(hjust = 0, size = 9),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12))
  )
}

# ── Category band panel ─────────────────────────────────────────

# thr: Item, k, Location and, with show_ci, CI_low, CI_high.
# drop_prob: data.frame item, category, p_never_modal (or NULL).
# rev_prob: data.frame item, thresholds, p_disordered (or NULL).
# cat_offset: added to category numbers in labels (1 when the lowest
#   category shown is 1, as in the severity part of the hurdle PCM).
#' @keywords internal
#' @noRd
.targeting_category_panel <- function(thr, item_order, xlim, category_labels,
                                      fills, show_ci, row_gap, caption_text,
                                      x_label, drop_prob = NULL,
                                      rev_prob = NULL, cat_offset = 0L) {
  rev_col <- "#B2182B"
  bar_h <- 0.52
  hh <- bar_h / 2
  yof <- function(x) match(as.character(x), item_order) * row_gap

  bd <- .category_bands(thr, item_order, xlim)
  b <- bd$bands
  b$y <- yof(b$Item)
  b$category <- factor(category_labels[b$cat + 1L], levels = category_labels)
  b$fill <- fills[b$cat + 1L]

  # Open-ended outer bands fade towards the panel edge
  is_out <- b$cat == 0L | b$cat == max(b$cat)
  n_step <- 40L
  fade <- do.call(rbind, lapply(which(is_out), function(i) {
    r <- b[i, ]
    br <- seq(r$xmin, r$xmax, length.out = n_step + 1L)
    data.frame(
      y = r$y,
      category = r$category,
      alpha = if (r$cat == 0L) {
        seq(0.22, 1, length.out = n_step)
      } else {
        seq(1, 0.22, length.out = n_step)
      },
      xmin = br[-length(br)],
      xmax = br[-1],
      stringsAsFactors = FALSE
    )
  }))

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = fade,
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$y - hh, ymax = .data$y + hh,
        fill = .data$category, alpha = .data$alpha
      )
    ) +
    ggplot2::geom_rect(
      data = b[!is_out, ],
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$y - hh, ymax = .data$y + hh,
        fill = .data$category
      )
    )

  joins <- b[b$xmin > xlim[1], ]
  p <- p +
    ggplot2::geom_segment(
      data = joins,
      ggplot2::aes(
        x = .data$xmin, xend = .data$xmin,
        y = .data$y - hh, yend = .data$y + hh
      ),
      colour = "white", linewidth = 0.6, show.legend = FALSE
    )

  lab <- b
  lab$xc <- (pmax(lab$xmin, xlim[1]) + pmin(lab$xmax, xlim[2])) / 2
  keep <- (pmin(lab$xmax, xlim[2]) - pmax(lab$xmin, xlim[1])) >
    0.30 * diff(xlim) / 8
  lab <- lab[keep, ]
  if (nrow(lab)) {
    p <- p +
      ggplot2::geom_text(
        data = lab,
        ggplot2::aes(x = .data$xc, y = .data$y,
                     label = .data$cat + cat_offset),
        colour = .band_label_colour(lab$fill),
        size = 3.1, show.legend = FALSE
      )
  }

  # Threshold reversals. A reversal of thresholds k and k + 1 is what
  # makes category k collapse, and the two are then labelled together on
  # the reversal arrow below the bar.
  revs <- .threshold_reversals(thr, item_order)

  # Collapsed categories: red tick, labelled with the category number and
  # the posterior probability that it is never the most likely response.
  # Collapses explained by a reversal get their label on the arrow instead.
  drops <- bd$drops
  if (!is.null(drops)) {
    drops <- drops[order(drops$Item, drops$cat), ]
    drops$p_never <- NA_real_
    if (!is.null(drop_prob)) {
      pm <- match(paste(drops$Item, drops$cat),
                  paste(drop_prob$item, drop_prob$category))
      drops$p_never <- drop_prob$p_never_modal[pm]
      drops$tag <- paste0(drops$cat + cat_offset, " (P = ",
                          .fmt_prob(drops$p_never), ")")
    } else {
      drops$tag <- as.character(drops$cat + cat_offset)
    }
    drops$on_arrow <- if (is.null(revs)) {
      rep(FALSE, nrow(drops))
    } else {
      paste(drops$Item, drops$cat) %in% paste(revs$Item, revs$k)
    }
    drops$tag[drops$on_arrow] <- ""
    key <- paste(drops$Item, round(drops$x, 8))
    drops$label <- vapply(key, function(kk) {
      tg <- drops$tag[key == kk]
      paste(tg[nzchar(tg)], collapse = ", ")
    }, character(1), USE.NAMES = FALSE)
    drops <- drops[!duplicated(key), ]
    drops$y <- yof(drops$Item)
    p <- p +
      ggplot2::geom_segment(
        data = drops,
        ggplot2::aes(
          x = .data$x, xend = .data$x,
          y = .data$y - hh - 0.06, yend = .data$y + hh + 0.06
        ),
        colour = "white", linewidth = 2.8, lineend = "butt", show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        data = drops,
        ggplot2::aes(
          x = .data$x, xend = .data$x,
          y = .data$y - hh - 0.04, yend = .data$y + hh + 0.04
        ),
        colour = rev_col, linewidth = 1.7, lineend = "butt", show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = drops[nzchar(drops$label), , drop = FALSE],
        ggplot2::aes(x = .data$x, y = .data$y + hh + 0.14, label = .data$label),
        colour = rev_col, size = 2.7, fontface = "bold",
        hjust = 0.5, vjust = 0, show.legend = FALSE
      )
  }

  # Threshold estimates and intervals below the bar, coloured by the
  # category the threshold gives entry to
  if (show_ci) {
    ci <- thr
    ci$y <- yof(ci$Item) - hh - 0.17
    ci$col <- .ci_contrast(fills)[ci$k + 1L]
    p <- p +
      ggplot2::geom_errorbar(
        data = ci,
        ggplot2::aes(y = .data$y, xmin = .data$CI_low, xmax = .data$CI_high),
        width = 0.17, linewidth = 2.1, colour = "white",
        orientation = "y", show.legend = FALSE
      ) +
      ggplot2::geom_errorbar(
        data = ci,
        ggplot2::aes(y = .data$y, xmin = .data$CI_low, xmax = .data$CI_high),
        width = 0.17, linewidth = 0.8, colour = ci$col,
        orientation = "y", show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = ci,
        ggplot2::aes(x = .data$Location, y = .data$y),
        colour = ci$col, size = 1.15, show.legend = FALSE
      )
  }

  # Reversal spans with their size and posterior probability
  n_rev <- 0L
  if (!is.null(revs)) {
    revs <- revs[order(revs$Item, revs$k), ]
    revs$row <- stats::ave(revs$k, revs$Item, FUN = seq_along)
    n_rev <- max(revs$row)
    revs$y <- yof(revs$Item) - hh - 0.32 - (revs$row - 1L) * 0.36
    revs$label <- sprintf("%.2f", revs$gap)
    if (!is.null(rev_prob)) {
      pm <- match(paste(revs$Item, paste0(revs$k, "-", revs$k + 1L)),
                  paste(rev_prob$item, rev_prob$thresholds))
      revs$label <- paste0(revs$label, " (P = ",
                           .fmt_prob(rev_prob$p_disordered[pm]), ")")
    }
    # A reversal that makes category k collapse: name the category and
    # give the probability that it is never the most likely response
    if (!is.null(drops)) {
      dm <- match(paste(revs$Item, revs$k),
                  paste(drops$Item, drops$cat))
      hit <- !is.na(dm)
      revs$label[hit] <- paste0(
        drops$cat[dm[hit]] + cat_offset, ": ", sprintf("%.2f", revs$gap[hit]),
        ifelse(is.na(drops$p_never[dm[hit]]), "",
               paste0(" (P = ", .fmt_prob(drops$p_never[dm[hit]]), ")"))
      )
    }
    p <- p +
      ggplot2::geom_segment(
        data = revs,
        ggplot2::aes(x = .data$lo, xend = .data$hi, y = .data$y, yend = .data$y),
        colour = rev_col, linewidth = 0.8, show.legend = FALSE,
        arrow = ggplot2::arrow(
          ends = "both", length = ggplot2::unit(3, "pt"), type = "closed"
        )
      ) +
      # Label centred under the arrow, clear of the threshold intervals
      ggplot2::geom_text(
        data = revs,
        ggplot2::aes(x = (.data$lo + .data$hi) / 2, y = .data$y - 0.04,
                     label = .data$label),
        colour = rev_col, size = 2.5, hjust = 0.5, vjust = 1,
        show.legend = FALSE
      )
  }

  y_bottom <- max(0.55, 0.32 + 0.36 * n_rev + 0.15)

  p +
    ggplot2::scale_fill_manual(
      values = stats::setNames(fills, category_labels),
      name = NULL,
      drop = FALSE,
      guide = ggplot2::guide_legend(nrow = 1, override.aes = list(alpha = 1))
    ) +
    ggplot2::scale_alpha_identity(guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = seq_along(item_order) * row_gap,
      labels = item_order,
      expand = ggplot2::expansion(add = c(y_bottom, max(0.55, hh + 0.34)))
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(xlim[1], xlim[2], by = 1),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = x_label, y = NULL, caption = caption_text) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 5, 5, 5)
    ) +
    .erb_caption_theme()
}
