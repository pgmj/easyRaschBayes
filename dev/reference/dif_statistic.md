# Differential Item Functioning (DIF) Analysis for Bayesian IRT Models

Tests for differential item functioning (DIF) in Bayesian Rasch-family
models fitted with brms by comparing item parameters across subgroups
defined by an exogenous variable. The function fits a DIF model that
includes group-by-item interactions and summarizes the posterior
distribution of the DIF effects.

## Usage

``` r
dif_statistic(
  model,
  group_var,
  item_var = item,
  person_var = id,
  data = NULL,
  dif_type = c("uniform", "non-uniform"),
  prob = 0.95,
  rope = 0.5,
  refit = TRUE,
  ...
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object from the baseline (no-DIF) model.

- group_var:

  An unquoted variable name identifying the grouping variable for DIF
  testing (e.g., `gender`). Must be a factor or character variable with
  exactly 2 levels in the current implementation.

- item_var:

  An unquoted variable name identifying the item grouping variable.
  Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable.
  Default is `id`.

- data:

  An optional data frame containing all variables needed for the DIF
  model, including the group variable. If `NULL` (the default), the
  function attempts to use `model$data`. Since the baseline model
  formula typically does not include the group variable, brms will have
  dropped it from the stored model data. In that case, you must supply
  the original data frame here.

- dif_type:

  Character. For polytomous ordinal models only. `"uniform"` (the
  default) tests for a uniform location shift per item via a
  `group:item` fixed-effect interaction. `"non-uniform"` fits
  group-specific thresholds per item and computes per-threshold DIF
  effects as the difference between groups. Ignored for dichotomous
  models.

- prob:

  Numeric in \\(0, 1)\\. Width of the credible intervals. Default is
  0.95.

- rope:

  Numeric. Half-width of the Region of Practical Equivalence (ROPE)
  around zero for DIF effects, on the logit scale. Default is 0.5,
  corresponding to a practically negligible DIF effect. Set to 0 to skip
  ROPE analysis.

- refit:

  Logical. If `TRUE` (the default), the DIF model is fitted
  automatically by updating the baseline model via
  [`update`](https://rdrr.io/r/stats/update.html), which reuses the
  compiled Stan code for faster sampling. If `FALSE`, only the DIF model
  formula is returned (useful for manual fitting with custom settings).

- ...:

  Additional arguments passed to
  [`update.brmsfit`](https://paulbuerkner.com/brms/reference/update.brmsfit.html)
  when refitting the DIF model (e.g., `cores`, `control`).

## Value

A list with the following elements:

- summary:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item (for uniform DIF) or per item × threshold (for
  non-uniform DIF) containing: `item`, optionally `threshold`,
  `dif_estimate` (posterior mean), `dif_lower`, `dif_upper` (credible
  interval), `dif_sd` (posterior SD), `pd` (probability of direction),
  `rope_percentage` (proportion inside ROPE), and `flag`
  (classification).

- dif_draws:

  A matrix of posterior draws for the DIF effects (draws × effects), for
  further analysis.

- dif_model:

  The fitted DIF `brmsfit` object (if `refit = TRUE`), or `NULL`.

- dif_formula:

  The `brmsformula` used for the DIF model.

- baseline_model:

  The original baseline model.

- plot:

  A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  forest plot of DIF effects with credible intervals and ROPE.

## Details

For polytomous models, two types of DIF can be tested:

- Uniform DIF (`dif_type = "uniform"`, default):

  A single location shift per item across groups, modelled as a
  `group:item` fixed-effect interaction. This tests whether the average
  item difficulty differs between groups.

- Non-uniform / threshold-level DIF (`dif_type = "non-uniform"`):

  Each item receives group-specific thresholds via
  `thres(gr = interaction(item, group))`. DIF effects are computed as
  the difference in each threshold between groups, revealing whether DIF
  affects specific response categories.

The function constructs a DIF model by adding a group-by-item
interaction to the baseline model:

- **Dichotomous models** (`family = bernoulli()`): The baseline
  `response ~ 1 + (1 | item) + (1 | id)` becomes
  `response ~ 1 + group + (1 + group | item) + (1 | id)`, where the
  group slope varying by item captures item-specific DIF.

- **Polytomous uniform DIF** (`dif_type = "uniform"`): The baseline
  `response | thres(gr = item) ~ 1 + (1 | id)` becomes
  `response | thres(gr = item) ~ 1 + group:item + (1 | id)`.

- **Polytomous non-uniform DIF** (`dif_type = "non-uniform"`): The
  baseline becomes `response | thres(gr = item_group) ~ 1 + (1 | id)`,
  where `item_group = interaction(item, group)`. Each item × group
  combination gets its own thresholds. DIF effects are the differences
  between group-specific thresholds for each item, computed draw-by-draw
  from the posterior.

DIF effects are summarized using:

- Probability of Direction (pd):

  The proportion of the posterior on the dominant side of zero. Values
  \> 0.975 indicate strong directional evidence.

- ROPE:

  The Region of Practical Equivalence (Kruschke, 2018). If \> 95\\ the
  DIF effect is practically negligible. If \> 95\\ outside, the effect
  is practically significant.

- Credible Interval:

  If the CI excludes zero, there is evidence of DIF at the specified
  credibility level.

## References

Kruschke, J. K. (2018). Rejecting or accepting parameter values in
Bayesian estimation. *Advances in Methods and Practices in Psychological
Science*, *1*(2), 270–280.

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
for item fit,
[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
for local dependence,
[`brm`](https://paulbuerkner.com/brms/reference/brm.html),
[`hypothesis`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html).

## Examples

``` r
# \donttest{
library(brms)
library(dplyr)
library(tidyr)
library(tibble)

# --- Dichotomous Rasch with DIF testing ---

set.seed(123)
df <- expand.grid(id = 1:200, item = paste0("I", 1:10)) %>%
  mutate(
    gender = rep(sample(c("M", "F"), 200, TRUE), each = 10),
    theta  = rep(rnorm(200), each = 10),
    delta  = rep(seq(-2, 2, length.out = 10), 200),
    dif    = ifelse(item == "I3" & gender == "F", 1.0,
             ifelse(item == "I7" & gender == "F", -0.8, 0)),
    p      = plogis(theta - delta - dif),
    response = rbinom(n(), 1, p)
  )

fit_base <- brm(
  response ~ 1 + (1 | item) + (1 | id),
  data   = df,
  family = bernoulli(),
  chains = 4, 
  cores  = 2, # use more cores if you have
  iter   = 1000 # use at least 2000 
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

dif_result <- dif_statistic(
  model     = fit_base,
  group_var = gender,
  data      = df
)
#> Error: object 'fit_base' not found

dif_result$summary
#> Error: object 'dif_result' not found
dif_result$plot
#> Error: object 'dif_result' not found

# --- Partial Credit Model: uniform DIF ---

df_pcm <- eRm::pcmdat2 %>%
  mutate(across(everything(), ~ .x + 1)) %>%
  rownames_to_column("id") %>%
  mutate(gender = sample(c("M", "F"), n(), TRUE)) %>%
  pivot_longer(!c(id, gender),
               names_to = "item", values_to = "response")

fit_pcm <- brm(
  response | thres(gr = item) ~ 1 + (1 | id),
  data   = df_pcm,
  family = acat,
  chains = 4, 
  cores  = 2, # use more cores if you have
  iter   = 1000 # use at least 2000 
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

# Uniform DIF (default): one shift per item
dif_uni <- dif_statistic(fit_pcm, group_var = gender, data = df_pcm)
#> Error: object 'fit_pcm' not found
dif_uni$plot
#> Error: object 'dif_uni' not found

# Non-uniform DIF: threshold-level effects
dif_nu <- dif_statistic(fit_pcm, group_var = gender, data = df_pcm,
                         dif_type = "non-uniform")
#> Error: object 'fit_pcm' not found
dif_nu$summary
#> Error: object 'dif_nu' not found
dif_nu$plot
#> Error: object 'dif_nu' not found
# }
```
