# Item-Restscore Association Conditional on the Total Score

Computes the item-restscore association (Kreiner, 2011) for Rasch-family
models fitted with brms, with a posterior predictive check that
conditions on each person's total score. The observed Goodman-Kruskal
gamma between each item and its rest-score (total score minus that item)
is compared with (i) its expected value under the model, conditional on
the total scores, and (ii) the gammas of datasets replicated conditional
on the total scores. Both are computed for each posterior draw of the
item parameters.

## Usage

``` r
item_restscore_statistic(
  model,
  item_var = item,
  person_var = id,
  ndraws_use = NULL
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object from a Rasch-family model in which the total score is
  sufficient for the person location: `family = acat` (partial credit or
  rating scale model) or `family = bernoulli()` (dichotomous Rasch
  model), both with the logit link, a person random intercept
  `(1 | id)`, and no `disc` term. Other models stop with an error.

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

## Value

A list with two elements:

- result:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item, sorted by `gamma_diff` in decreasing order, and the
  columns:

  item

  :   The item identifier, named after `item_var`.

  gamma_obs

  :   The observed gamma between the item and the rest-score.

  gamma_exp

  :   The expected gamma conditional on the total scores, as in Kreiner
      (2011), averaged over posterior draws.

  gamma_diff

  :   `gamma_obs - gamma_exp`.

  ppp

  :   Mid-p posterior predictive p-value, the proportion of draws where
      `gamma_obs` exceeds the replicated gamma, counting ties as one
      half. Values near 1 indicate stronger association than the model
      implies, values near 0 weaker.

  ppp_adj

  :   Two-sided posterior predictive p-value adjusted for the number of
      items by a single-step max-T procedure. See Details.

  gamma_rep

  :   Mean of the replicated gamma across draws.

  gamma_diff_q025, gamma_diff_q975

  :   95% credible interval for `gamma_obs` minus the expected gamma. It
      reflects uncertainty in the item parameters only. See Details.

  gamma_diff_q005, gamma_diff_q995

  :   99% credible interval for the same difference.

  gamma_rep_q025, gamma_rep_q975

  :   95% posterior predictive interval for the gamma of a replicated
      dataset.

  gamma_rep_q005, gamma_rep_q995

  :   99% posterior predictive interval for the same.

- draws:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) in
  long format with one row per item and draw, and the columns `item`
  (named after `item_var`), `draw`, `gamma_rep` (replicated gamma) and
  `gamma_exp` (expected gamma).

## Details

Observed gamma above the expected value means the item discriminates
more than the equal discrimination of the Rasch model implies, and
observed gamma below it means the item discriminates less.

Under a Rasch-family model the total score is sufficient for the person
location. Given the total score, the distribution of a person's
responses depends on the item parameters only. For each posterior draw
\\s\\ the function:

1.  Takes the category probabilities from
    [`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
    and removes the person location from them, which leaves item weights
    that depend on that draw's item parameters only.

2.  Computes, with elementary symmetric functions, the probability of
    each item score given each person's total score. Summed over
    persons, these give the expected item \\\times\\ rest-score table
    and its gamma, \\\gamma^{exp(s)}\_i\\.

3.  Draws a replicated dataset in which every person keeps their
    observed total score, and computes the replicated gamma
    \\\gamma^{rep(s)}\_i\\.

The observed gamma \\\gamma^{obs}\_i\\ is computed once from the data.
Persons with the minimum or maximum total score carry no information
about fit and their responses are reproduced as observed in every
replicated dataset.

Because replication conditions on the total scores, the reference
distribution does not depend on the person locations. The data are still
used to estimate the item parameters, but these are few relative to the
number of persons. This removes the main reason why posterior predictive
p-values based on person parameters drawn from the same data are
conservative (Robins et al., 2000).

The `gamma_diff` interval describes how far the observed sample departs
from the model expectation, given the uncertainty in the item
parameters. It contains no sampling variation, so it is not a test, and
it is not a credible interval for misfit in the population. Use `ppp` or
`ppp_adj` to judge whether a departure is larger than chance would
produce. As with any test, they will flag smaller departures as the
sample size grows.

`ppp_adj` standardises each item's observed gamma by the mean and
standard deviation of its replicated gammas. Its reference is the
largest absolute standardised deviation across items within each
replicated dataset, so it accounts for testing all items at once.

With missing responses, the total and rest-scores are sums over the
items a person answered, scored from 0. The conditional distribution is
computed for each person's set of answered items.

Replication uses random numbers, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible results.

## References

Kreiner, S. (2011). A note on item-restscore association in Rasch
models. *Applied Psychological Measurement*, *35*(7), 557–561.

Goodman, L. A. & Kruskal, W. H. (1954). Measures of association for
cross classifications. *Journal of the American Statistical
Association*, *49*(268), 732–764.

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

Robins, J. M., van der Vaart, A. & Ventura, V. (2000). Asymptotic
distribution of P values in composite null models. *Journal of the
American Statistical Association*, *95*(452), 1143–1156.
[doi:10.1080/01621459.2000.10474310](https://doi.org/10.1080/01621459.2000.10474310)

## See also

[`item_restscore_post`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_post.md)
for a summary table and plot,
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
for posterior predictive infit/outfit,
[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
for Bayesian Q3 residual correlations,
[`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html).

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
  cores  = 1, # use more cores if you have
  iter   = 500 # use at least 2000 
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

# Item-restscore association
irs <- item_restscore_statistic(
  model      = fit_pcm,
  ndraws_use = 100 # use at least 500
)
#> Error: object 'fit_pcm' not found

# Post-process draws
irs_results <- item_restscore_post(irs)
#> Error: object 'irs' not found
irs_results$summary
#> Error: object 'irs_results' not found
irs_results$plot
#> Error: object 'irs_results' not found

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

irs_rm <- item_restscore_statistic(
  model      = fit_rm,
  ndraws_use = 100 # use at least 500
)
#> Error: object 'fit_rm' not found

# Post-process draws
irs_results <- item_restscore_post(irs_rm)
#> Error: object 'irs_rm' not found
irs_results$summary
#> Error: object 'irs_results' not found
irs_results$plot
#> Error: object 'irs_results' not found
# }
```
