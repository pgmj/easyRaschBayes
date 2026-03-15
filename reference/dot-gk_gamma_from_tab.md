# Compute Gamma from a Pre-built Contingency Table

Shared engine for
[`gk_gamma()`](https://pgmj.github.io/easyRaschBayes/reference/gk_gamma.md)
and
[`gk_gamma_vec()`](https://pgmj.github.io/easyRaschBayes/reference/gk_gamma_vec.md).
Uses cumulative suffix sums and prefix column sums for \\O(nr \times
nc)\\ concordant/discordant pair counting.

## Usage

``` r
.gk_gamma_from_tab(tab, nr, nc)
```

## Arguments

- tab:

  Integer or numeric matrix (contingency table).

- nr:

  Number of rows.

- nc:

  Number of columns.

## Value

A scalar: the gamma coefficient, or `NA_real_`.
