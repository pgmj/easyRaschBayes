# Posterior Predictive Item Fit Statistic for Binary Bayesian IRT Models

Computes posterior predictive item (or person) fit statistics for
dichotomous Bayesian IRT models fitted with brms. For each posterior
draw, observed and replicated data are compared via a user-supplied
criterion function, grouped by item, person, or any other variable.
Posterior predictive p-values can then be derived from the output to
assess fit.

## Usage

``` r
fit_statistic_rm(model, criterion, group, ndraws_use = NULL)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object with a binary response (e.g., `family = bernoulli()`).

- criterion:

  A function with signature `function(y, p)` that computes a pointwise
  fit criterion, where `y` is the binary response (0 or 1) and `p` is
  the predicted probability of success. A common choice is the Bernoulli
  log-likelihood: `function(y, p) y * log(p) + (1 - y) * log(1 - p)`.

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

This function is the binary-response counterpart of
[`fit_statistic_pcm`](https://pgmj.github.io/easyRaschBayes/reference/fit_statistic_pcm.md),
which handles polytomous (ordinal / categorical) models. For dichotomous
models,
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
returns a 2D matrix (S x N) of success probabilities, so the criterion
function receives the observed binary response and the corresponding
probability directly.

The procedure follows the posterior predictive checking approach
described in Bürkner (2020):

1.  Draw posterior expected success probabilities via
    [`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
    and posterior predicted binary responses via
    [`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

2.  Apply the user-supplied `criterion` function pointwise to both
    observed and replicated data paired with the predicted
    probabilities.

3.  Aggregate (sum) the criterion values within each level of `group`
    and each posterior draw.

The standard criterion for binary models is the Bernoulli
log-likelihood: \$\$\ell(y, p) = y \log(p) + (1 - y) \log(1 - p).\$\$

## References

Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS)
with Bayesian Item Response Models. *Journal of Intelligence*, *8*(1).
[doi:10.3390/jintelligence8010005](https://doi.org/10.3390/jintelligence8010005)

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`fit_statistic_pcm`](https://pgmj.github.io/easyRaschBayes/reference/fit_statistic_pcm.md)
for polytomous (ordinal/categorical) models,
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

# --- Dichotomous Rasch Model ---

# Prepare binary response data in long format
df_rm <- eRm::raschdat3 %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

# Fit a dichotomous Rasch model
fit_rm <- brm(
  response ~ 1 + (1 | item) + (1 | id),
  data   = df_rm,
  family = bernoulli(),
  chains = 4,
  cores  = 4,
  iter   = 2000
)

# Bernoulli log-likelihood criterion
ll_bernoulli <- function(y, p) y * log(p) + (1 - y) * log(1 - p)

# Compute item fit statistics
item_fit <- fit_statistic_rm(
  model      = fit_rm,
  criterion  = ll_bernoulli,
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
person_fit <- fit_statistic_rm(
  model      = fit_rm,
  criterion  = ll_bernoulli,
  group      = id,
  ndraws_use = 500
)

person_fit %>%
  group_by(id) %>%
  summarise(
    observed   = mean(crit),
    replicated = mean(crit_rep),
    ppp        = mean(crit_rep > crit)
  )

# --- 1PL model with item-specific intercepts ---

# Alternative parameterisation with fixed item effects
fit_1pl <- brm(
  response ~ 0 + item + (1 | id),
  data   = df_rm,
  family = bernoulli(),
  chains = 4,
  cores  = 4,
  iter   = 2000
)

item_fit_1pl <- fit_statistic_rm(
  model      = fit_1pl,
  criterion  = ll_bernoulli,
  group      = item,
  ndraws_use = 500
)

item_fit_1pl %>%
  group_by(item) %>%
  summarise(
    observed   = mean(crit),
    replicated = mean(crit_rep),
    ppp        = mean(crit_rep > crit)
  )
} # }
```
