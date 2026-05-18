# Posterior Predictive Infit Statistic for the Hurdle Partial Credit Model

Computes conditional item infit statistics separately for the two
submodels of a hurdle partial credit model fitted with brms using the
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
custom family: (i) the *hurdle* submodel for \\P(Y \> 0)\\ (Bernoulli)
and (ii) the *partial credit* severity submodel for \\P(Y = k \mid Y \>
0)\\ on the positive categories. For each posterior draw, expected
values and variances are derived from the submodel-specific category
probabilities, and variance-weighted standardised residuals are computed
for both observed and replicated data.

## Usage

``` r
infit_statistic_hpcm(
  model,
  item_var = item,
  person_var = id,
  ndraws_use = NULL,
  outfit = FALSE
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object using the
  [`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
  custom family (i.e., `posterior_epred` returns an `S x N x K_total`
  array whose first category is the hurdle / zero probability).

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data (e.g., `item`).

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data (e.g., `id`).

- ndraws_use:

  Optional positive integer. If specified, a random subset of posterior
  draws of this size is used. If `NULL` (the default), all draws are
  used.

- outfit:

  Logical. If `TRUE`, outfit statistics are computed alongside infit.
  Default is `FALSE`.

## Value

A list with two elements, each a
[`tibble`](https://tibble.tidyverse.org/reference/tibble.html) in the
same format as the output of
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
(and directly compatible with
[`infit_post`](https://pgmj.github.io/easyRaschBayes/reference/infit_post.md)):

- `hurdle`:

  Item infit for the Bernoulli hurdle submodel on \\1\[Y \> 0\]\\,
  evaluated on *all* observations.

- `pcm`:

  Item infit for the partial credit severity submodel on \\P(Y = k \mid
  Y \> 0)\\, evaluated only on the observations with \\Y\_{obs} \> 0\\.

## Details

The hurdle PCM splits the generative process into:

1.  A Bernoulli hurdle with \\hu = P(Y = 0)\\.

2.  A partial credit / acat-logit severity process over the positive
    categories \\1, \ldots, K - 1\\, applied only when the hurdle is
    crossed.

`posterior_epred` for the `hurdle_acat` family returns an
`S x N x K_total` array whose first category is \\hu\\ and whose
remaining categories are \\(1 - hu) \cdot P\_{sev}(k)\\. The two
submodel infits are computed as follows:

**Hurdle submodel.** All observations contribute. The Bernoulli moments
are \\E\_{hurdle} = 1 - hu\\ and \\Var\_{hurdle} = hu \cdot (1 - hu)\\.
Observed and replicated scores are \\1\[Y\_{obs} \> 0\]\\ and
\\1\[Y\_{rep} \> 0\]\\ with \\Y\_{rep}\\ obtained from
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

**Partial credit submodel.** Only observations with \\Y\_{obs} \> 0\\
contribute. Conditional probabilities are \$\$P(Y = k \mid Y \> 0) =
epred\[, , k+1\] / (1 - hu), \quad k = 1, \ldots, K - 1.\$\$ The
conditional expected value and variance use category scores \\1, \ldots,
K - 1\\. Replicated severity values are drawn independently for each
(draw, observation) from these conditional probabilities via inverse-CDF
sampling, so the partial credit PPC is not contaminated by hurdle
misfit.

Within each submodel the infit per item is \$\$Infit_i^{(s)} = \sum_v
(X\_{vi} - E\_{vi}^{(s)})^2 / \sum_v Var\_{vi}^{(s)},\$\$ with the sum
restricted to the rows the submodel applies to (all rows for the hurdle;
rows with \\Y\_{obs} \> 0\\ for partial credit).

## References

Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013). *Rasch
Models in Health*. Iste and Wiley, pp. 86–90.

Kreiner, S. & Christensen, K. B. (2011). Exact evaluation of bias in
Rasch model residuals. *Advances in Mathematics Research*, *12*, 19–40.

Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, *27*(2), 261-279.
[doi:10.1037/met0000395](https://doi.org/10.1037/met0000395)

## See also

[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
for the single-submodel version,
[`infit_post`](https://pgmj.github.io/easyRaschBayes/reference/infit_post.md)
for summarising and plotting the draws,
[`q3_statistic_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic_hpcm.md)
for hPCM Q3 residual correlations.
