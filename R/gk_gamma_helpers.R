#' Goodman-Kruskal Gamma from Two Vectors (No table() Overhead)
#'
#' Computes the Goodman-Kruskal gamma coefficient directly from two
#' numeric vectors, bypassing \code{table()} for speed. Uses manual
#' tabulation with integer indexing and cumulative sums for efficient
#' concordant/discordant pair counting.
#'
#' @param x Integer or numeric vector (e.g., item scores).
#' @param y Integer or numeric vector of the same length (e.g., rest scores).
#'
#' @return A scalar: the gamma coefficient in \eqn{[-1, 1]}, or
#'   \code{NA_real_} if the table has fewer than 2 rows or 2 columns.
#'
#' @references
#' Goodman, L. A. & Kruskal, W. H. (1954). Measures of association for
#' cross classifications. \emph{Journal of the American Statistical
#' Association}, \emph{49}(268), 732--764.
#'
#' @keywords internal
gk_gamma_vec <- function(x, y) {
  # Fast manual tabulation
  ux <- sort.int(unique.default(x), method = "quick")
  uy <- sort.int(unique.default(y), method = "quick")
  nr <- length(ux)
  nc <- length(uy)
  if (nr < 2L || nc < 2L) return(NA_real_)
  
  xi <- match(x, ux)
  yi <- match(y, uy)
  
  # Build contingency table via integer indexing
  tab <- matrix(0L, nrow = nr, ncol = nc)
  for (idx in seq_along(xi)) {
    tab[xi[idx], yi[idx]] <- tab[xi[idx], yi[idx]] + 1L
  }
  
  .gk_gamma_from_tab(tab, nr, nc)
}


#' Goodman-Kruskal Gamma from a Contingency Table
#'
#' Computes the Goodman-Kruskal gamma coefficient for an ordinal
#' contingency table. This is an internal helper function.
#'
#' @param tab A matrix or table representing a contingency table with
#'   rows and columns in ordinal order.
#'
#' @return A scalar: the gamma coefficient in \eqn{[-1, 1]}.
#'
#' @details
#' Gamma is defined as \eqn{(C - D) / (C + D)}, where \eqn{C} is the
#' number of concordant pairs and \eqn{D} is the number of discordant
#' pairs. Uses cumulative suffix sums for \eqn{O(nr \times nc)}
#' complexity.
#'
#' @references
#' Goodman, L. A. & Kruskal, W. H. (1954). Measures of association for
#' cross classifications. \emph{Journal of the American Statistical
#' Association}, \emph{49}(268), 732--764.
#'
#' @keywords internal
gk_gamma <- function(tab) {
  tab <- as.matrix(tab)
  nr <- nrow(tab)
  nc <- ncol(tab)
  if (nr < 2L || nc < 2L) return(NA_real_)
  
  .gk_gamma_from_tab(tab, nr, nc)
}


#' Compute Gamma from a Pre-built Contingency Table
#'
#' Shared engine for \code{gk_gamma()} and \code{gk_gamma_vec()}.
#' Uses cumulative suffix sums and prefix column sums for
#' \eqn{O(nr \times nc)} concordant/discordant pair counting.
#'
#' @param tab Integer or numeric matrix (contingency table).
#' @param nr Number of rows.
#' @param nc Number of columns.
#'
#' @return A scalar: the gamma coefficient, or \code{NA_real_}.
#'
#' @keywords internal
.gk_gamma_from_tab <- function(tab, nr, nc) {
  
  
  # --- Suffix sum: suffix[i,j] = sum of tab[i:nr, j:nc] ---
  # Padded to (nr+1) x (nc+1) so out-of-bounds reads give 0
  suffix <- matrix(0, nrow = nr + 1L, ncol = nc + 1L)
  for (i in seq.int(nr, 1L)) {
    for (j in seq.int(nc, 1L)) {
      suffix[i, j] <- tab[i, j] +
        suffix[i + 1L, j] +
        suffix[i, j + 1L] -
        suffix[i + 1L, j + 1L]
    }
  }
  
  # --- Prefix-column sum: pcol[i,j] = sum of tab[i:nr, 1:j] ---
  # Padded to (nr+1) x (nc+1); column 0 stays 0 as the base case
  pcol <- matrix(0, nrow = nr + 1L, ncol = nc + 1L)
  for (i in seq.int(nr, 1L)) {
    # j runs over columns 1..nc, stored in columns 2..(nc+1)
    # so "column j-1" is always >= 1 and never hits index 0
    for (j in seq.int(2L, nc + 1L)) {
      pcol[i, j] <- tab[i, j - 1L] +
        pcol[i + 1L, j] +
        pcol[i, j - 1L] -
        pcol[i + 1L, j - 1L]
    }
  }
  
  C <- 0
  D <- 0
  for (i in seq_len(nr - 1L)) {
    for (j in seq_len(nc)) {
      n_ij <- tab[i, j]
      if (n_ij == 0L) next
      # Concordant: cells strictly below-right
      if (j < nc) {
        C <- C + n_ij * suffix[i + 1L, j + 1L]
      }
      # Discordant: cells strictly below-left
      # pcol column index is shifted by +1 (column j in pcol = index j+1)
      if (j > 1L) {
        D <- D + n_ij * pcol[i + 1L, j]
      }
    }
  }
  
  denom <- C + D
  if (denom == 0) return(NA_real_)
  (C - D) / denom
}