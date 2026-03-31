# Goodman-Kruskal Gamma from a Contingency Table

Computes the Goodman-Kruskal gamma coefficient for an ordinal
contingency table. This is an internal helper function.

## Usage

``` r
gk_gamma(tab)
```

## Arguments

- tab:

  A matrix or table representing a contingency table with rows and
  columns in ordinal order.

## Value

A scalar: the gamma coefficient in \\\[-1, 1\]\\.

## Details

Gamma is defined as \\(C - D) / (C + D)\\, where \\C\\ is the number of
concordant pairs and \\D\\ is the number of discordant pairs. Uses
cumulative suffix sums for \\O(nr \times nc)\\ complexity.

## References

Goodman, L. A. & Kruskal, W. H. (1954). Measures of association for
cross classifications. *Journal of the American Statistical
Association*, *49*(268), 732–764.
