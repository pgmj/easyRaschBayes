# Extract Item Parameters from a Hurdle Partial Credit Model

Extracts item parameters from a fitted hurdle Rasch partial credit model
(`family = `[`hurdle_acat()`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)),
returning a per-submodel breakdown: hurdle item difficulties (Bernoulli
logit on \\P(Y \> 0)\\) and partial credit thresholds. Each submodel's
output mirrors the shape of
[`item_parameters`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters.md)
so existing plotting and downstream functions can be applied directly to
`res$hurdle` or `res$pcm`.

## Usage

``` r
item_parameters_hpcm(
  model,
  item_var = item,
  person_var = id,
  draws = FALSE,
  center = TRUE,
  prob = 0.95
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object using the
  [`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
  custom family. The recommended formula is


      bf(
        response | thres(gr = item) ~ 1 + (1 |g| id),
        hu ~ 0 + factor(item) + (1 |g| id)
      )

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data. Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data. Default is `id`.

- draws:

  Logical. If `TRUE`, a draws matrix of full posterior draws is included
  in each submodel's output. Default is `FALSE`.

- center:

  Logical. If `TRUE` (the default), item parameters are shifted within
  each submodel so their mean is zero, matching the convention in
  [`item_parameters`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters.md).
  Person parameters reported by
  [`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md)
  use the same shifts.

- prob:

  Numeric in \\(0, 1)\\. Width of the highest density continuous
  interval (HDCI) reported in the summary. Default is 0.95.

## Value

A list with three elements:

- `hurdle`:

  A list with the same structure as the output of
  [`item_parameters`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters.md)
  applied to a dichotomous Rasch model: `locations`, `locations_wide`,
  `summary`, `item_information`, `person_sd`, optionally `draws_matrix`.
  `location` is the brms posterior mean of `b_hu_factoritem<i>`; higher
  values mean more zeros (a harder hurdle to cross).

- `pcm`:

  A list with the same structure as
  [`item_parameters`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters.md)
  applied to a PCM: `locations` (long), `locations_wide` (t1, t2, ...,
  location), `summary`, `item_information`, `threshold_order`,
  `person_sd`, optionally `draws_matrix`. `location` columns are
  posterior means of `b_Intercept[<item>, <k>]`.

- `correlation`:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  `mean`, `sd`, `hdci_lower`, `hdci_upper` summarising the posterior of
  \\\rho(\theta\_{hurdle}, \theta\_{pcm})\\. This is the brms
  correlation `cor_id__Intercept__hu_Intercept` with its sign flipped to
  match the "higher = more presence" convention used for the hurdle
  person trait (see Details).

## Details

**Hurdle submodel.** The hurdle linear predictor is \\logit(hu) =
\delta\_{hurdle, i} + \tilde{r}\_v\\ with \\hu = P(Y = 0)\\; higher
\\\delta\_{hurdle, i}\\ means more zeros, so this is reported directly
as the item "location" (harder = higher value, standard Rasch convention
for the Bernoulli on \\P(Y \> 0)\\). Hurdle item information is the
Bernoulli variance \\p(1-p)\\ evaluated at \\\theta = \delta\_{hurdle,
i}\\, which is exactly \\1/4\\.

**Partial credit submodel.** Thresholds \\\tau\_{ik}\\ are the per-item
PCM threshold parameters. Item information uses the standard PCM
formula; threshold ordering diagnostics are the same as for
[`item_parameters`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters.md).

**Sign of the trait correlation.** The brms model reports
`cor_id__Intercept__hu_Intercept`, which is the correlation between the
brms random effects `r_id[, Intercept]` and `r_id[, hu_Intercept]`. The
latter has the opposite sign of the conventional "susceptibility" person
trait \\\theta\_{hurdle}\\ (because higher values of the brms random
effect mean more zeros, i.e., lower susceptibility). The correlation
reported here is therefore \\-\text{cor}\_{brms}\\, so that positive
values mean higher susceptibility goes with higher severity. The
marginal SDs are unchanged by the sign flip.

**Centering.** When `center = TRUE`, the hurdle item difficulties are
shifted by their mean and the PCM thresholds are shifted by the grand
mean of all PCM thresholds. The same shifts are applied to the
corresponding person traits in
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md),
preserving the underlying likelihood.

## References

Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, *27*(2), 261-279.
[doi:10.1037/met0000395](https://doi.org/10.1037/met0000395)

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`item_parameters`](https://pgmj.github.io/easyRaschBayes/reference/item_parameters.md)
for the single-submodel version,
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md)
for the person-side counterpart,
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
for the custom brms family.
