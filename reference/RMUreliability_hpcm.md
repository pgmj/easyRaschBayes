# Relative Measurement Uncertainty Reliability for the Hurdle PCM

Computes posterior reliability via Relative Measurement Uncertainty
(Bignardi, Kievit & Bürkner, 2025) separately for each of the two person
traits in a hurdle partial credit model fitted with the
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
custom family: the presence / susceptibility trait \\\theta\_{hurdle}\\
and the severity trait \\\theta\_{pcm}\\. Internally extracts
per-submodel person draws via
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md)
and applies
[`RMUreliability`](https://pgmj.github.io/easyRaschBayes/reference/RMUreliability.md)
to each.

## Usage

``` r
RMUreliability_hpcm(
  x,
  item_var = item,
  person_var = id,
  level = 0.95,
  verbose = FALSE,
  center = TRUE
)
```

## Arguments

- x:

  Either:

  - a fitted
    [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
    object using the
    [`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)
    family, or

  - the result of
    [`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md)`(model, draws = TRUE)`,
    i.e., a list with `$hurdle$draws_matrix` and `$pcm$draws_matrix`.

  Passing a pre-computed
  [`person_parameters_hpcm()`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md)
  result avoids re-extracting the posterior draws from `brms`.

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data. Used only when `x` is a `brmsfit`; ignored otherwise.
  Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data. Used only when `x` is a `brmsfit`; ignored otherwise.
  Default is `id`.

- level:

  Numeric in \\(0, 1)\\. Credibility level for the highest density
  continuous interval on each reliability. Default is `0.95`.

- verbose:

  Logical. If `TRUE`, print the number of subjects and posterior draws
  used for each submodel. Default is `FALSE`.

- center:

  Logical. If `TRUE` (the default), use centered person draws (mean item
  parameter = 0 in each submodel). The centering shift is a constant per
  draw and does not affect the correlation across draw halves, so
  reliability is invariant to this choice; the argument is exposed only
  for consistency with
  [`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md).

## Value

A list with two elements:

- `hurdle`:

  A data frame with one row containing `rmu_estimate`,
  `hdci_lowerbound`, and `hdci_upperbound` for reliability of
  \\\theta\_{hurdle}\\, as returned by
  [`RMUreliability`](https://pgmj.github.io/easyRaschBayes/reference/RMUreliability.md).

- `pcm`:

  Same format, for reliability of \\\theta\_{pcm}\\.

## Details

**Two reliabilities, two traits.** The hurdle PCM gives every person a
posterior on \\\theta\_{hurdle}\\ (presence / susceptibility) and on
\\\theta\_{pcm}\\ (severity given presence). Each posterior has its own
measurement uncertainty, and the RMU framework therefore yields one
reliability per trait. The two values address different questions:
\\\theta\_{hurdle}\\ reliability is about how well the items rank
persons by overall endorsement (a function of how many zeros appear and
how discriminating the items are at the gate), while \\\theta\_{pcm}\\
reliability is about how well the conditional severity is recovered (a
function of how many positive responses each person provides and how
spread out the PCM thresholds are). Reporting both is necessary, as
emphasised in Magnus and Garnier-Villarreal (2022, p. 277).

**What about all-zero responders?** For persons with \\Y = 0\\ on every
item, no items provide information about \\\theta\_{pcm}\\ directly.
Their posterior on \\\theta\_{pcm}\\ is shaped by the
multivariate-normal prior linking it to \\\theta\_{hurdle}\\ (and is
therefore wide but well-defined). These persons still contribute to the
RMU calculation; their contribution is appropriately downweighted
through their large posterior variance.

**Computational note.** When `x` is a `brmsfit`, this function calls
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md)`(..., draws = TRUE)`
internally, which extracts the full posterior draws of the person random
effects. For large models you may prefer to call
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md)
once and pass its result to `RMUreliability_hpcm`,
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md),
etc.

## References

Bignardi, G., Kievit, R., & Bürkner, P.-C. (2025). A general method for
estimating reliability using Bayesian Measurement Uncertainty.
*PsyArXiv*.
[doi:10.31234/osf.io/h54k8_v1](https://doi.org/10.31234/osf.io/h54k8_v1)

Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, *27*(2), 261-279.
[doi:10.1037/met0000395](https://doi.org/10.1037/met0000395)

## See also

[`RMUreliability`](https://pgmj.github.io/easyRaschBayes/reference/RMUreliability.md)
for the single-trait version,
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/reference/person_parameters_hpcm.md)
for the underlying person draws extraction.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a hurdle PCM with the hurdle_acat() family
fit <- brms::brm(
  brms::bf(
    response | brms::thres(gr = item) ~ 1 + (1 |g| id),
    hu ~ 0 + factor(item) + (1 |g| id)
  ),
  data     = dat,
  family   = hurdle_acat(),
  stanvars = hurdle_acat_stanvars()
)

# Option 1: pass the brmsfit directly
rel <- RMUreliability_hpcm(fit)
rel$hurdle
rel$pcm

# Option 2: reuse a cached person_parameters_hpcm() result
pp  <- person_parameters_hpcm(fit, draws = TRUE)
rel <- RMUreliability_hpcm(pp)
} # }
```
