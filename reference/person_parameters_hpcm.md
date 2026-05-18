# Extract Person Parameters from a Hurdle Partial Credit Model

Extracts person trait estimates from a fitted hurdle Rasch partial
credit model
(`family = `[`hurdle_acat()`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)),
returning a per-submodel breakdown: a presence trait
\\\theta\_{hurdle}\\ (susceptibility) on the Bernoulli gate and a
severity trait \\\theta\_{pcm}\\ on the partial credit submodel. Each
submodel's output mirrors the shape of
[`person_parameters`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters.md)
so `res$hurdle` and `res$pcm` can be passed to existing downstream
functions
([`RMUreliability`](https://pgmj.github.io/easyRaschBayes/reference/RMUreliability.md),
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md),
etc.).

## Usage

``` r
person_parameters_hpcm(
  model,
  item_var = item,
  person_var = id,
  draws = FALSE,
  center = TRUE,
  theta_range = c(-7, 7)
)
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

- draws:

  Logical. If `TRUE`, a matrix of full posterior draws (persons x draws)
  is included in each submodel's output. Default is `FALSE`.

- center:

  Logical. If `TRUE` (the default), item parameters and person
  parameters are shifted within each submodel so that the mean item
  parameter is zero. Person traits are shifted by the same constant as
  the items in their submodel.

- theta_range:

  A numeric vector of length 2 giving the range for the Newton-Raphson
  WLE search on the hurdle submodel. Default is `c(-7, 7)`.

## Value

A list with three elements:

- `hurdle`:

  A list with the same structure as
  [`person_parameters`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters.md)
  for a dichotomous Rasch model: `person_estimates` (one row per person;
  `sum_score` is the number of items with \\Y \> 0\\; both EAP and
  Warm's WLE are provided), `score_table`, and optionally
  `draws_matrix`. The hurdle EAP is the sign-flipped brms random effect
  on `hu`, so higher values mean greater presence / susceptibility.

- `pcm`:

  A list with `person_estimates` (columns `id`, `sum_score`, `n_active`,
  `eap`, `eap_se`), `score_table`, and optionally `draws_matrix`.
  `sum_score` is the sum of \\Y\\ across items with \\Y\_{vi} \> 0\\ for
  person \\v\\; `n_active` is the count of such items. WLE columns are
  omitted intentionally (see Details).

- `correlation`:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  summarising the posterior of \\\rho(\theta\_{hurdle},
  \theta\_{pcm})\\, identical to the `correlation` element returned by
  [`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters_hpcm.md).

## Details

**Hurdle person trait.** Defined as \\\theta\_{hurdle, v} = -r\_{id}(v,
\texttt{hu\\Intercept})\\, i.e., the brms random effect on `hu` with its
sign flipped. Under this convention, the hurdle submodel reads \\P(Y \>
0) = \mathrm{plogis}(\theta\_{hurdle, v} - \delta\_{hurdle, i})\\, so
higher \\\theta\_{hurdle}\\ means greater presence / susceptibility. WLE
is computed by treating the gate as a dichotomous Rasch model on \\1\[Y
\> 0\]\\ with item difficulties \\\delta\_{hurdle, i}\\ from
[`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters_hpcm.md).

**Partial credit person trait.** Defined as \\\theta\_{pcm, v} =
r\_{id}(v, \texttt{Intercept})\\. Higher values mean greater severity
given presence.

**Why no WLE on the partial credit submodel.** Conditional on the
hurdle, the severity submodel is a PCM, but the set of items
contributing to person \\v\\'s severity score depends on \\v\\'s own
pattern of \\1\[Y\_{vi} \> 0\]\\: only items with \\Y\_{vi} \> 0\\
provide PCM information about \\\theta\_{pcm, v}\\. The sum of partial
credit scores is therefore not a sufficient statistic across the sample,
and a standard score-table WLE is not well defined. The Bayesian EAP,
which integrates over the (correlated) multivariate prior on
\\(\theta\_{hurdle}, \theta\_{pcm})\\, remains well defined for all
persons (including those with all-zero patterns) and is therefore the
recommended point estimate; see Magnus and Garnier-Villarreal (2022, p.
272ff) for discussion in the closely-related MZI-GRM.

**Centering.** When `center = TRUE`, the same shifts that
[`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters_hpcm.md)
applies to the item parameters are applied here to the corresponding
person traits. The likelihood is unchanged.

**Row ordering.** `person_estimates` and `draws_matrix` preserve the
order of first appearance of each person ID in the model data, matching
[`person_parameters`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters.md).

## References

Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, *27*(2), 261-279.
[doi:10.1037/met0000395](https://doi.org/10.1037/met0000395)

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika*, *54*(3), 427–450.
[doi:10.1007/BF02294627](https://doi.org/10.1007/BF02294627)

## See also

[`person_parameters`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters.md)
for the single-submodel version,
[`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters_hpcm.md)
for the item-side counterpart,
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
for the custom brms family.
