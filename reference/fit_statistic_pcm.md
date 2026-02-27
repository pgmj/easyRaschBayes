# Posterior Predictive Item Fit Statistic for Bayesian IRT Models

Computes posterior predictive item (or person) fit statistics for
Bayesian IRT models fitted with brms. For each posterior draw, observed
and replicated data are compared via a user-supplied criterion function,
grouped by item, person, or any other variable. Posterior predictive
p-values can then be derived from the output to assess fit.

## Usage

``` r
fit_statistic_pcm(model, criterion, group, ndraws_use = NULL)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object.

- criterion:

  A function with signature `function(y, p)` that computes a pointwise
  fit criterion. For ordinal and categorical models, `y` is the observed
  (or replicated) response category and `p` is the model-predicted
  probability of that category. For binary models, `y` is the binary
  response and `p` is the predicted probability of success.

- group:

  An unquoted variable name (e.g., `item` or `id`) indicating the
  grouping variable over which the fit statistic is aggregated.
  Typically `item` for item fit or `id` for person fit.

- ndraws_use:

  Optional positive integer. If specified, a random subset of posterior
  draws of this size is used, which can speed up computation for large
  models. If `NULL` (the default), all draws are used.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
the following columns:

- `group`:

  The grouping variable (e.g., item name or person id).

- draw:

  Integer index of the posterior draw.

- crit:

  The observed fit statistic (criterion applied to observed data) summed
  within each group and draw.

- crit_rep:

  The replicated fit statistic (criterion applied to posterior predicted
  data) summed within each group and draw.

- crit_diff:

  The difference `crit_rep - crit`.

The output is grouped by the grouping variable. Posterior predictive
p-values can be obtained by computing `mean(crit_rep > crit)` within
each group.

## Details

The function implements the posterior predictive checking approach for
item fit described in Bürkner (2020). The procedure works as follows:

1.  Draw posterior expected category probabilities via
    [`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
    and posterior predicted responses via
    [`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

2.  For ordinal or categorical models (3D array output from
    `posterior_epred`), extract the probability assigned to the observed
    response category and to the replicated response category for each
    draw and observation.

3.  Apply the user-supplied `criterion` function to compute pointwise
    fit values for both observed and replicated data.

4.  Aggregate (sum) the criterion values within each level of `group`
    and each posterior draw.

A common choice for ordinal IRT models is the categorical log-likelihood
criterion `function(y, p) log(p)`. For binary (e.g., dichotomous Rasch)
models, the Bernoulli log-likelihood
`function(y, p) y * log(p) + (1 - y) * log(1 - p)` may be used instead.

## References

Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS)
with Bayesian Item Response Models. *Journal of Intelligence*, *8*(1).
[doi:10.3390/jintelligence8010005](https://doi.org/10.3390/jintelligence8010005)

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`fit_statistic_rm`](https://pgmj.github.io/easyRaschBayes/reference/fit_statistic_rm.md)
for dichotomous Rasch models,
[`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
for expected predictions,
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
for posterior predictive samples,
[`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html) for
graphical posterior predictive checks.

## Examples

``` r
if (FALSE) { # \dontrun{
library(brms)
library(dplyr)
library(tidyr)
library(tibble)

# --- Polytomous Rasch (Partial Credit Model) ---

# Prepare data in long format
df_pcm <- eRm::pcmdat2 %>%
  mutate(across(everything(), ~ .x + 1)) %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

# Fit a Partial Credit Model using the adjacent category family
fit_pcm <- brm(
  response | thres(gr = item) ~ 1 + (1 | id),
  data    = df_pcm,
  family  = acat,
  chains  = 4,
  cores   = 4,
  iter    = 2000
)

# Categorical log-likelihood criterion (for polytomous models)
ll_categorical <- function(y, p) log(p)

# Compute item fit statistics
item_fit <- fit_statistic_pcm(
  model      = fit_pcm,
  criterion  = ll_categorical,
  group      = item,
  ndraws_use = 500
)

# Summarise: posterior predictive p-values per item
item_fit %>%
  group_by(item) %>%
  summarise(
    observed   = mean(crit),
    replicated = mean(crit_rep),
    ppp        = mean(crit_rep > crit)
  )

# Use ggplot2 to make a histogram
library(ggplot2)
item_fit %>%
  ggplot(aes(crit_diff)) +
  geom_histogram(aes(fill = ifelse(crit_diff > 0, "above","below"))) +
  facet_wrap("item") +
  theme_bw() +
  theme(legend.position = "none")

# Compute person fit statistics
person_fit <- fit_statistic_pcm(
  model      = fit_pcm,
  criterion  = ll_categorical,
  group      = id,
  ndraws_use = 500
)
} # }
```
