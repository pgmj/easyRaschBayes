# Posterior Predictive Q3 Residual Correlations for the Hurdle PCM

Computes Yen's Q3 residual correlations separately for the two submodels
of a hurdle partial credit model fitted with brms using the
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
custom family: (i) the *hurdle* submodel for \\P(Y \> 0)\\ (Bernoulli)
and (ii) the *partial credit* severity submodel for \\P(Y = k \mid Y \>
0)\\ on the positive categories. For each posterior draw, residuals are
computed from each submodel's expected values and Pearson correlations
between item-pair residuals are returned for both observed and
replicated data.

## Usage

``` r
q3_statistic_hpcm(model, item_var = item, person_var = id, ndraws_use = NULL)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object using the
  [`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
  custom family.

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data. Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data. Default is `id`.

- ndraws_use:

  Optional positive integer. If specified, a random subset of posterior
  draws of this size is used. If `NULL` (the default), all draws are
  used.

## Value

A list with two elements, each a
[`tibble`](https://tibble.tidyverse.org/reference/tibble.html) in the
same long format as
[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)
(and directly compatible with
[`q3_post`](https://pgmj.github.io/easyRaschBayes/reference/q3_post.md)):

- `hurdle`:

  Q3 correlations of hurdle residuals \\1\[Y \> 0\] - (1 - hu)\\, using
  *all* observations.

- `pcm`:

  Q3 correlations of partial credit residuals \\Y - E\_{pcm}\\, using
  only observations with \\Y\_{obs} \> 0\\. Rows with \\Y\_{obs} = 0\\
  are set to `NA` and excluded pairwise.

## Details

The hurdle PCM has a Bernoulli hurdle \\hu = P(Y = 0)\\ and a partial
credit severity process on the positive categories \\1, \ldots, K - 1\\.
`posterior_epred` for
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
returns an `S x N x K_total` array whose first category is \\hu\\ and
whose remaining categories are \\(1 - hu) \cdot P\_{sev}(k)\\.

**Hurdle residuals.** For each (draw, observation): \$\$d\_{hurdle} =
1\[Y \> 0\] - (1 - hu).\$\$ All observations contribute. Replicated
residuals use \\1\[Y\_{rep} \> 0\]\\ from
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

**Partial credit residuals.** For each (draw, observation) with
\\Y\_{obs} \> 0\\: \$\$d\_{pcm} = Y - \sum\_{k=1}^{K-1} k \cdot P(Y = k
\mid Y \> 0),\$\$ with conditional probabilities computed as \\epred\[,
, k+1\] / (1 - hu)\\. Replicated partial credit values are drawn
independently per (draw, observation) from these conditional
probabilities via inverse-CDF sampling, so the partial credit PPC is not
contaminated by hurdle misfit. Rows with \\Y\_{obs} = 0\\ are set to
`NA` in the wide residual matrix and excluded via
`use = "pairwise.complete.obs"` when correlations are computed.

## References

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement*, *8*(2), 125–145.

Christensen, K. B., Makransky, G. & Horton, M. (2017). Critical values
for Yen's Q3: Identification of local dependence in the Rasch model
using residual correlations. *Applied Psychological Measurement*,
*41*(3), 178–194.

Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, *27*(2), 261-279.
[doi:10.1037/met0000395](https://doi.org/10.1037/met0000395)

## See also

[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)
for the single-submodel version,
[`q3_post`](https://pgmj.github.io/easyRaschBayes/reference/q3_post.md)
for summarising and plotting the draws,
[`infit_statistic_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic_hpcm.md)
for hPCM infit.
