# Goodman-Kruskal Gamma from Two Vectors (No table() Overhead)

Computes the Goodman-Kruskal gamma coefficient directly from two numeric
vectors, bypassing [`table()`](https://rdrr.io/r/base/table.html) for
speed. Uses manual tabulation with integer indexing and cumulative sums
for efficient concordant/discordant pair counting.

## Usage

``` r
gk_gamma_vec(x, y)
```

## Arguments

- x:

  Integer or numeric vector (e.g., item scores).

- y:

  Integer or numeric vector of the same length (e.g., rest scores).

## Value

A scalar: the gamma coefficient in \\\[-1, 1\]\\, or `NA_real_` if the
table has fewer than 2 rows or 2 columns.

## References

Goodman, L. A. & Kruskal, W. H. (1954). Measures of association for
cross classifications. *Journal of the American Statistical
Association*, *49*(268), 732–764.
