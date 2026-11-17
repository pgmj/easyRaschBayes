# Orthonormal sum-to-zero contrasts

Contrasts for a factor whose effects should sum to zero, with the same
prior variance for every level. Useful for small facets, such as a few
items or raters, modelled as fixed effects in a many-facet Rasch model.

## Usage

``` r
contr_sum_zero(n, contrasts = TRUE, sparse = FALSE)
```

## Arguments

- n:

  The number of levels, or a vector of level labels (for example
  `levels(x)`), as in
  [`stats::contr.sum()`](https://rdrr.io/r/stats/contrast.html).

- contrasts:

  If `FALSE`, one indicator column per level is returned instead, as for
  [`stats::contr.sum()`](https://rdrr.io/r/stats/contrast.html). R uses
  this for factors without an intercept to pair with.

- sparse:

  Only `FALSE` is supported. The argument is there so that the function
  can be named in `options(contrasts = ...)` or set with
  `contrasts(x) <- "contr_sum_zero"`.

## Value

A matrix with one row per level and \\K - 1\\ columns. Columns are
orthonormal and every column sums to zero. With `contrasts = FALSE`, the
\\K x K\\ identity matrix.

## Details

With sum-to-zero contrasts the model estimates the effect of each level
relative to the average level, instead of relative to a reference level.
R's [`stats::contr.sum()`](https://rdrr.io/r/stats/contrast.html) does
this by setting the last level to minus the sum of the others. Normal
priors on the coefficients then give the last level a larger prior
variance than the others: \\K - 1\\ times larger for \\K\\ levels, three
times larger with four items. With little data the prior is then not the
same for all levels.

`contr_sum_zero()` returns orthonormal contrasts (a normalised Helmert
basis). They describe the same sum-to-zero effects, but independent
normal priors on the coefficients give every level the same prior
variance. This is the coding used by Stan's `sum_to_zero_vector` type
(Morris, 2025), which brms does not offer directly.

A normal(0, s) prior on the coefficients gives each level effect a prior
SD of \\s \sqrt{(K - 1) / K}\\. For a prior SD of \\\sigma\\ per level,
use \\s = \sigma \sqrt{K / (K - 1)}\\. For example, with four levels and
\\\sigma = 1\\, use normal(0, 1.15).

The item and facet extractors in this package read the contrasts stored
with the model data, so the level effects are reported on the original
scale whatever contrasts were used.

## References

Morris, M. (2025). *The sum-to-zero constraint in Stan*. Stan case
study.
<https://mc-stan.org/learn-stan/case-studies/sum_to_zero_vector.html>

## See also

[`stats::contr.sum()`](https://rdrr.io/r/stats/contrast.html),
[`stats::contrasts()`](https://rdrr.io/r/stats/contrasts.html)

## Examples

``` r
contr_sum_zero(4)
#>         [,1]       [,2]       [,3]
#> 1 -0.7071068 -0.4082483 -0.2886751
#> 2  0.7071068 -0.4082483 -0.2886751
#> 3  0.0000000  0.8164966 -0.2886751
#> 4  0.0000000  0.0000000  0.8660254

# The implied prior covariance of the level effects is proportional to
# diag(K) - 1/K: equal variances, equal negative correlations
H <- contr_sum_zero(c("A", "B", "C", "D"))
round(H %*% t(H), 3)
#>       A     B     C     D
#> A  0.75 -0.25 -0.25 -0.25
#> B -0.25  0.75 -0.25 -0.25
#> C -0.25 -0.25  0.75 -0.25
#> D -0.25 -0.25 -0.25  0.75

if (FALSE) { # \dontrun{
d <- data.frame(...)  # long data: id, item, rater, response
contrasts(d$item) <- contr_sum_zero(levels(d$item))
fit <- brms::brm(
  response ~ 1 + item + (1 | rater) + (1 | id),
  data = d, family = brms::acat(),
  prior = brms::prior(normal(0, 1.15), class = b)
)
} # }
```
