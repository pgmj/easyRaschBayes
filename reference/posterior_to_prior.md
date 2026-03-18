# Extract Informative Priors from a Fitted Bayesian IRT Model

Takes a fitted
[`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
object and constructs a
[`brmsprior`](https://paulbuerkner.com/brms/reference/set_prior.html)
object in which each item parameter receives a `normal(mean, sd)` prior
derived from its posterior distribution. The person-level random effect
SD prior is also updated. The returned prior can be passed directly to
`update` (or [`brm`](https://paulbuerkner.com/brms/reference/brm.html))
to refit the model with empirical Bayes / informative priors — useful
for anchoring scales, warm-starting a model on new data, or regularising
estimation with small samples.

## Usage

``` r
posterior_to_prior(model, item_var = item, person_var = id, mult = 1)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object. Supported parameterisations:

  Polytomous ordinal

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

- mult:

  Numeric multiplier applied to each posterior SD before it is used as
  the prior SD. Values \> 1 widen the priors (less informative); values
  \< 1 tighten them. Default is 1 (use posterior SD directly).

## Value

A [`brmsprior`](https://paulbuerkner.com/brms/reference/set_prior.html)
object that can be supplied to the `prior` argument of
[`brm`](https://paulbuerkner.com/brms/reference/brm.html) or `update`.

## Details

The function extracts all posterior draws via
[`as_draws_df`](https://mc-stan.org/posterior/reference/draws_df.html),
computes the mean and SD of each parameter's marginal posterior, and
constructs `normal(mean, sd * mult)` priors.

**Polytomous ordinal models** with grouped thresholds
(`thres(gr = item)`): each threshold receives its own prior via
`brms::set_prior("normal(...)", class = "Intercept", group = item, coef = threshold_index)`.

**Dichotomous Rasch models** parameterised as
`response ~ 1 + (1 | item) + (1 | id)`: priors are set on the global
intercept (`class = "Intercept"`), the item-level SD
(`class = "sd", group = item_var`), and the person-level SD.

**Dichotomous 1PL models** parameterised as
`response ~ 0 + item + (1 | id)`: each item-specific fixed effect (e.g.,
`b_itemI1`) receives its own `normal(mean, sd)` prior via
`brms::set_prior(..., class = "b", coef = "itemI1")`.

In all cases the person-level SD receives a `normal(mean, sd * mult)`
prior (brms applies the lower bound of zero automatically for SD
parameters).

## Examples

``` r
# \donttest{
library(brms)
library(dplyr)
library(tidyr)
library(tibble)

# --- Partial Credit Model ---

df_pcm <- eRm::pcmdat2 %>%
  mutate(across(everything(), ~ .x + 1)) %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

fit_pcm <- brm(
  response | thres(gr = item) ~ 1 + (1 | id),
  data   = df_pcm,
  family = acat,
  chains = 4,
  cores  = 4,
  iter   = 2000
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

# Extract posterior-informed priors
new_priors <- posterior_to_prior(fit_pcm)
#> Error: object 'fit_pcm' not found
new_priors
#> Error: object 'new_priors' not found

# Widen the priors by a factor of 2
wide_priors <- posterior_to_prior(fit_pcm, mult = 2)
#> Error: object 'fit_pcm' not found

# --- Dichotomous 1PL (fixed item effects) ---

df_rm <- eRm::raschdat3 %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

fit_1pl <- brm(
  response ~ 0 + item + (1 | id),
  data   = df_rm,
  family = bernoulli(),
  chains = 4,
  cores  = 4,
  iter   = 2000
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

priors_1pl <- posterior_to_prior(fit_1pl)
#> Error: object 'fit_1pl' not found
priors_1pl
#> Error: object 'priors_1pl' not found
# }
```
