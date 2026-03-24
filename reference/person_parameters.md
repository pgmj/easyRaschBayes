# Extract Person Parameters from a Bayesian Rasch Model

Extracts person ability estimates from a fitted Bayesian Rasch model.
Returns both Bayesian EAP (expected a posteriori) estimates with
posterior SDs and frequentist WLE (Warm's weighted likelihood) estimates
with standard errors, plus a lookup table mapping ordinal sum scores to
both scales.

## Usage

``` r
person_parameters(
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
  object. Supported parameterisations:

  Polytomous ordinal (PCM)

  :   e.g., `family = acat` with `thres(gr = item)`, producing
      item-specific thresholds.

  Dichotomous Rasch (random items)

  :   e.g., `response ~ 1 + (1 | item) + (1 | id)` with
      `family = bernoulli()`.

  Dichotomous 1PL (fixed items)

  :   e.g., `response ~ 0 + item + (1 | id)` with
      `family = bernoulli()`.

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data. Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data. Default is `id`.

- draws:

  Logical. If `TRUE`, a matrix of full posterior draws (persons x draws)
  is included in the output. Default is `FALSE`.

- center:

  Logical. If `TRUE` (the default), person parameters and item
  difficulties are recentered so that the mean item difficulty is zero,
  matching the convention in frequentist CML Rasch estimation. If
  `FALSE`, raw brms parameterisation is used.

- theta_range:

  A numeric vector of length 2 specifying the range for the
  Newton-Raphson WLE search. Default is `c(-7, 7)`.

## Value

A list with the following elements:

- person_estimates:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per person containing: the person ID, `sum_score`, `eap`
  (posterior mean of person random effect), `eap_se` (posterior SD),
  `wle` (Warm's weighted likelihood estimate), and `wle_se` (asymptotic
  SE of the WLE). Rows are ordered to match the original person order in
  the model data (i.e., the order of first appearance of each person
  ID).

- score_table:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  mapping each observed ordinal sum score to its mean EAP, mean EAP SE,
  WLE, and WLE SE. Extreme scores (0 and maximum) receive WLE estimates
  at the boundary of `theta_range` with `NA` standard errors.

- draws_matrix:

  (Only if `draws = TRUE`) A numeric matrix with rows = persons and
  columns = posterior draws. Row names are person IDs. Rows are ordered
  to match the original person order in the model data. Can be passed
  directly to
  [`RMUreliability`](https://pgmj.github.io/easyRaschBayes/reference/RMUreliability.md).

## Details

**EAP estimates** are extracted as the posterior means of the person
random effects (`r_id[j, Intercept]`) from
[`as_draws_df`](https://mc-stan.org/posterior/reference/draws_df.html),
with the posterior SD serving as the standard error. These are the
standard Bayesian point estimates and reflect shrinkage toward the
population mean.

**WLE estimates** (Warm, 1989) are computed from the posterior mean item
parameters using Newton-Raphson iteration with adaptive step damping on
the Warm-corrected likelihood. WLE adds a bias correction term
\\J(\theta) / (2 I(\theta))\\ to the score equations, where
\\I(\theta)\\ is the test information and \\J(\theta) = \sum_i \sum_c
(c - E_i)^3 P\_{ic}\\ is the sum of third central moments. This produces
estimates with reduced finite-sample bias compared to MLE, especially at
extreme scores (Warm, 1989).

The Newton-Raphson algorithm uses adaptive step damping following the
approach in iarm (Mueller): the maximum allowed step size shrinks by a
factor of 1.05 each iteration, preventing overshoot and ensuring
convergence for near-extreme scores.

For extreme scores (sum score = 0 or maximum possible), the WLE is not
well-defined. These cases are assigned the boundary values of
`theta_range` with `NA` standard errors.

When `center = TRUE` (the default), item difficulty parameters are
shifted so their mean is zero, and EAP person parameters are shifted by
the same constant. WLE is computed from the centered item parameters.
This matches the convention in frequentist CML Rasch estimation.

**Row ordering:** Both `person_estimates` and `draws_matrix` preserve
the original person order from the model data (order of first appearance
of each person ID). This allows direct row-binding with the source data
without re-matching.

## References

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika*, *54*(3), 427–450.
[doi:10.1007/BF02294627](https://doi.org/10.1007/BF02294627)

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013). *Rasch
Models in Health*. Iste and Wiley, pp. 63–70.

## See also

[`RMUreliability`](https://pgmj.github.io/easyRaschBayes/reference/RMUreliability.md),
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md),
[`ranef`](https://rdrr.io/pkg/nlme/man/random.effects.html),
[`as_draws_df`](https://mc-stan.org/posterior/reference/draws_df.html).

## Examples

``` r
# \donttest{
library(brms)
library(dplyr)
library(tidyr)
library(tibble)

# --- Dichotomous Rasch Model (random items) ---

df_rm <- eRm::raschdat3 %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

fit_rm <- brm(
  response ~ 1 + (1 | item) + (1 | id),
  data   = df_rm,
  family = bernoulli(),
  chains = 4, cores = 2, iter = 1000 # use more iter and cores
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

# Basic usage — rows match original person order
pp <- person_parameters(fit_rm)
#> Error: object 'fit_rm' not found
pp$person_estimates
#> Error: object 'pp' not found
pp$score_table
#> Error: object 'pp' not found

# With full posterior draws (e.g., for RMUreliability)
pp_draws <- person_parameters(fit_rm, draws = TRUE)
#> Error: object 'fit_rm' not found
RMUreliability(pp_draws$draws_matrix)
#> Error: object 'pp_draws' not found

# --- Polytomous PCM (acat) ---

df_pcm <- eRm::pcmdat2 %>%
  mutate(across(everything(), ~ .x + 1)) %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

fit_pcm <- brm(
  response | thres(gr = item) ~ 1 + (1 | id),
  data   = df_pcm,
  family = acat,
  chains = 4, cores = 2, iter = 1000 # use more iter and cores
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

pp_pcm <- person_parameters(fit_pcm)
#> Error: object 'fit_pcm' not found
pp_pcm$score_table
#> Error: object 'pp_pcm' not found
# }
```
