# Posterior Predictive Infit Statistic for Bayesian IRT Models

Computes item infit (and optionally outfit) for each posterior draw of a
Rasch-family model fitted with brms, for the observed data and for
replicated data. By default (`method = "conditional"`) the expected
values and variances are conditional on each person's total score, as in
the conditional infit of Christensen et al. (2013) and Kreiner &
Christensen (2011), and the replicated data keep each person's observed
total score. With `method = "theta"` they are computed given each draw's
person locations instead.

## Usage

``` r
infit_statistic(
  model,
  item_var = item,
  person_var = id,
  ndraws_use = NULL,
  outfit = FALSE,
  method = c("conditional", "theta")
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object. For `method = "conditional"` it must be a Rasch-family model
  in which the total score is sufficient for the person location:
  `family = acat` (partial credit or rating scale model) or
  `family = bernoulli()` (dichotomous Rasch model), both with the logit
  link, a person random intercept `(1 | id)`, and no `disc` term.
  `method = "theta"` accepts any ordinal or Bernoulli model.

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

- method:

  Either `"conditional"` (the default), with expectations conditional on
  the total score, or `"theta"`, with expectations given each draw's
  person locations. `"theta"` gives the same computation as
  [`infit_statistic_hpcm`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic_hpcm.md)
  and as `infit_statistic()` in easyRaschBayes 0.3.0.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
the following columns:

- item:

  The item identifier, named after `item_var`.

- draw:

  Integer index of the posterior draw.

- infit:

  The observed infit statistic for that item and draw.

- infit_rep:

  The replicated infit statistic for that item and draw.

- outfit:

  (Only if `outfit = TRUE`) The observed outfit statistic for that item
  and draw.

- outfit_rep:

  (Only if `outfit = TRUE`) The replicated outfit statistic for that
  item and draw.

The output is grouped by the item variable. Use
[`infit_post`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
for posterior predictive p-values and a plot.

## Details

Both methods use the infit and outfit mean squares (Christensen et al.,
2013; Müller, 2020) as discrepancy measures in a posterior predictive
check. For each posterior draw \\s\\, each observation gets an expected
value \\E^{(s)}\_{vi}\\ and a variance \\Var^{(s)}\_{vi}\\, and
\$\$Infit_i^{(s)} = \frac{\sum_v (X\_{vi} - E^{(s)}\_{vi})^2} {\sum_v
Var^{(s)}\_{vi}}, \qquad Outfit_i^{(s)} = \frac{1}{N_i} \sum_v
\frac{(X\_{vi} - E^{(s)}\_{vi})^2}{Var^{(s)}\_{vi}}.\$\$ The replicated
statistics use the same \\E\\ and \\Var\\ with replicated responses in
place of \\X\\.

**Conditional method.** Under a Rasch-family model the total score
\\R_v\\ is sufficient for the person location. The function takes the
category probabilities from
[`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
removes the person location from them, and computes \\E^{(s)}\_{vi} =
E(X\_{vi} \mid R_v)\\ and \\Var^{(s)}\_{vi} = Var(X\_{vi} \mid R_v)\\
with elementary symmetric functions of that draw's item parameters.
Persons with the minimum or maximum total score have zero variance on
every item and are left out, as in
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html).
Replicated responses are drawn conditional on each person's observed
total score, so the reference distribution does not depend on the person
locations. Given the same item parameters, the observed values equal the
conditional infit and outfit of
[`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html). The
observed infit varies across draws only through the uncertainty in the
item parameters. With missing responses, the total score is the sum over
the items a person answered, scored from 0, and the conditional
distribution is computed for each person's set of answered items.

**Theta method.** \\E\\ and \\Var\\ are computed from the category
probabilities given the person location \\\theta_v^{(s)}\\ of each draw,
and replicated responses come from
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html).
The person locations are estimated from the same data, which makes
posterior predictive p-values of this kind conservative. They tend to
cluster around 0.5, so values near 0 or 1 are rarer than p-values of the
same size would be (Robins et al., 2000).

The conditional method is better calibrated than the theta method, but
its p-values are also conservative, because the item parameters are
estimated from the same data. In simulations under a fitting Rasch
model, two-sided p-values below .05 occurred for 0.5% to 2% of items,
with the lowest rates for polytomous items. Misfit can therefore go
undetected. It is suggested to also evaluate item fit with
[`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md),
whose p-values were close to uniform in the same simulations.

With either method, the posterior predictive p-value is not a measure of
the size of misfit. For a fixed amount of misfit it moves toward 0 or 1
as the sample size grows.

Replication uses random numbers, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible results.

## References

Christensen, K. B., Kreiner, S. & Mesbah, M. (Eds.) (2013). *Rasch
Models in Health*. Iste and Wiley, pp. 86–90.

Kreiner, S. & Christensen, K. B. (2011). Exact evaluation of Bias in
Rasch model residuals. *Advances in Mathematics Research*, 12, 19–40.

Müller, M. (2020). Item fit statistics for Rasch analysis: can we trust
them? *Journal of Statistical Distributions and Applications*, *7*(1).
[doi:10.1186/s40488-020-00108-7](https://doi.org/10.1186/s40488-020-00108-7)

Robins, J. M., van der Vaart, A. & Ventura, V. (2000). Asymptotic
distribution of P values in composite null models. *Journal of the
American Statistical Association*, *95*(452), 1143–1156.
[doi:10.1080/01621459.2000.10474310](https://doi.org/10.1080/01621459.2000.10474310)

## See also

[`infit_post`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
for posterior predictive p-values and a plot,
[`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
for the item-restscore association,
[`fit_statistic_pcm`](https://pgmj.github.io/easyRaschBayes/dev/reference/fit_statistic_pcm.md)
and
[`fit_statistic_rm`](https://pgmj.github.io/easyRaschBayes/dev/reference/fit_statistic_rm.md)
for general-purpose posterior predictive fit statistics with
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
