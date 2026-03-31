# Posterior Predictive Infit Statistic for Bayesian IRT Models

Computes a Bayesian analogue of the conditional item infit statistic (as
described in Christensen, Kreiner & Mesbah, 2013) for Rasch-family
models fitted with brms. For each posterior draw, expected values and
variances are derived from the category probabilities returned by
[`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
and variance-weighted standardised residuals are computed for both
observed and replicated data. The result can be summarised into
posterior predictive p-values to assess item fit.

## Usage

``` r
infit_statistic(
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
  object from an ordinal IRT model (e.g., `family = acat` for a partial
  credit model or `family = bernoulli()` for a dichotomous Rasch model).

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
  Default is `FALSE` (infit only), since outfit is highly sensitive to
  outliers and rarely recommended for Rasch diagnostics.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
the following columns:

- item:

  The item identifier.

- draw:

  Integer index of the posterior draw.

- infit:

  The observed infit statistic for that item and draw.

- infit_rep:

  The replicated infit statistic (based on posterior predicted data) for
  that item and draw.

- outfit:

  (Only if `outfit = TRUE`) The observed outfit statistic for that item
  and draw.

- outfit_rep:

  (Only if `outfit = TRUE`) The replicated outfit statistic for that
  item and draw.

The output is grouped by the item variable. Posterior predictive
p-values can be obtained by computing, e.g., `mean(infit_rep > infit)`
within each item.

## Details

The procedure adapts the conditional infit/outfit statistics
(Christensen et al., 2013; Kreiner & Christensen, 2011; Müller, 2020) to
the Bayesian framework:

1.  For each posterior draw \\s\\, category probabilities
    \\P^{(s)}(X\_{vi} = c)\\ are obtained from
    [`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html).

2.  The conditional expected value and variance for each observation are
    computed as: \$\$E^{(s)}\_{vi} = \sum_c c \cdot P^{(s)}(X\_{vi} =
    c)\$\$ \$\$Var^{(s)}\_{vi} = \sum_c (c - E^{(s)}\_{vi})^2 \cdot
    P^{(s)}(X\_{vi} = c)\$\$

3.  Standardised squared residuals are: \$\$Z^{2(s)}\_{vi} = (X\_{vi} -
    E^{(s)}\_{vi})^2 / Var^{(s)}\_{vi}\$\$

4.  The infit statistic for item \\i\\ is the variance-weighted mean of
    \\Z^2\\ across persons: \$\$Infit_i^{(s)} = \frac{\sum_v
    Var\_{vi}^{(s)} Z^{2(s)}\_{vi}} {\sum_v Var\_{vi}^{(s)}}\$\$

5.  If requested, the outfit is the unweighted mean of \\Z^2\\.

## References

Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013). *Rasch
Models in Health*. Iste and Wiley, pp. 86–90.

Kreiner, S. & Christensen, K. B. (2011). Exact evaluation of Bias in
Rasch model residuals. *Advances in Mathematics Research*, 12, 19–40.

Müller, M. (2020). Item fit statistics for Rasch analysis: can we trust
them? *Journal of Statistical Distributions and Applications*, *7*(1).
[doi:10.1186/s40488-020-00108-7](https://doi.org/10.1186/s40488-020-00108-7)

## See also

[`fit_statistic_pcm`](https://pgmj.github.io/easyRaschBayes/dev/reference/fit_statistic_pcm.md)
for a general-purpose posterior predictive fit statistic with
user-supplied criterion functions,
[`fit_statistic_rm`](https://pgmj.github.io/easyRaschBayes/dev/reference/fit_statistic_rm.md)
for a general-purpose posterior predictive fit statistic with
user-supplied criterion functions,
[`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

## Examples

``` r
# \donttest{
library(brms)
library(dplyr)
library(tidyr)
library(tibble)

# --- Partial Credit Model (polytomous) ---

df_pcm <- eRm::pcmdat2 %>%
  mutate(across(everything(), ~ .x + 1)) %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

fit_pcm <- brm(
  response | thres(gr = item) ~ 1 + (1 | id),
  data   = df_pcm,
  family = acat,
  chains = 4,
  cores  = 1, # use more cores if you have
  iter   = 500 # use at least 2000 
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

# Compute infit per item
item_infit <- infit_statistic(
  model      = fit_pcm,
  ndraws_use = 100 # use at least 500
)
#> Error: object 'fit_pcm' not found

# Post-process draws
infit_results <- infit_post(item_infit)
#> Error: object 'item_infit' not found
infit_results$summary
#> Error: object 'infit_results' not found
infit_results$hdi
#> Error: object 'infit_results' not found
infit_results$plot
#> Error: object 'infit_results' not found

# --- Dichotomous Rasch Model ---

df_rm <- eRm::raschdat3 %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

fit_rm <- brm(
  response ~ 1 + (1 | item) + (1 | id),
  data   = df_rm,
  family = bernoulli(),
  chains = 4,
  cores  = 1, # use more cores if you have
  iter   = 500 # use at least 2000 
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

item_infit_rm <- infit_statistic(
  model      = fit_rm,
  ndraws_use = 100 # use at least 500
)
#> Error: object 'fit_rm' not found

# Post-process draws
infit_results <- infit_post(item_infit_rm)
#> Error: object 'item_infit_rm' not found
infit_results$summary
#> Error: object 'infit_results' not found
infit_results$hdi
#> Error: object 'infit_results' not found
infit_results$plot
#> Error: object 'infit_results' not found

# }
```
