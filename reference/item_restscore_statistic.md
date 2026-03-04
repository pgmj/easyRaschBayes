# Posterior Predictive Item-Restscore Association for Bayesian IRT Models

Computes a Bayesian analogue of the item-restscore association test
(Kreiner, 2011) for Rasch-family models fitted with brms. For each
posterior draw, the Goodman-Kruskal gamma coefficient between each
item's score and the rest-score (total score minus that item) is
computed for both observed and replicated data. Posterior predictive
p-values indicate whether the observed association is stronger than the
model predicts, which signals violations of local independence or
unidimensionality.

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
  object from an ordinal IRT model (e.g., `family = acat` for a partial
  credit model) or a dichotomous model (`family = bernoulli()`).

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

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
the following columns:

- item:

  The item identifier.

- gamma_obs:

  Posterior mean of the observed Goodman-Kruskal gamma between this item
  and the rest-score.

- gamma_rep:

  Posterior mean of the replicated gamma.

- gamma_diff:

  Posterior mean of `gamma_obs - gamma_rep`. Positive values indicate
  the observed item-restscore association is stronger than the model
  expects.

- ppp:

  Posterior predictive p-value: `mean(gamma_obs > gamma_rep)` across
  draws. Values close to 1 indicate the item discriminates more than the
  model predicts (too high discrimination). Values close to 0 indicate
  the item discriminates less than expected (too low discrimination,
  e.g., noise or miskeyed item).

- gamma_obs_q025, gamma_obs_q975:

  95\\ the observed gamma.

- gamma_obs_q005, gamma_obs_q995:

  99\\ the observed gamma.

- gamma_diff_q025, gamma_diff_q975:

  95\\ the gamma difference.

- gamma_diff_q005, gamma_diff_q995:

  99\\ the gamma difference.

## Details

The item-restscore association is a key diagnostic in Rasch measurement.
Under the Rasch model, each item should relate to the latent trait (and
hence the rest-score) only through the modelled relationship.
Goodman-Kruskal's gamma is a rank-based measure of association for
ordinal cross-tabulations that is well-suited for this purpose (Kreiner,
2011).

The procedure for each posterior draw \\s\\ is:

1.  Obtain replicated responses \\Y^{rep(s)}\\ from
    [`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

2.  For each item \\i\\ and each person \\v\\, compute the rest-score:
    \\R^{obs}\_{vi} = \sum\_{j \neq i} X\_{vj}\\ for observed data and
    \\R^{rep(s)}\_{vi} = \sum\_{j \neq i} Y^{rep(s)}\_{vj}\\ for
    replicated data.

3.  Cross-tabulate item score \\\times\\ rest-score and compute the
    Goodman-Kruskal gamma for both observed and replicated data.

4.  Compare the two gammas across draws.

Items with `ppp` close to 1 have observed item-restscore association
that is consistently stronger than the model predicts. This typically
indicates that the item discriminates more than assumed under the
equal-discrimination Rasch model (i.e., a violation of the Rasch
assumption). Items with `ppp` close to 0 discriminate less than
expected.

## References

Kreiner, S. (2011). A note on item-restscore association in Rasch
models. *Applied Psychological Measurement*, *35*(7), 557–561.

Goodman, L. A. & Kruskal, W. H. (1954). Measures of association for
cross classifications. *Journal of the American Statistical
Association*, *49*(268), 732–764.

Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS)
with Bayesian Item Response Models. *Journal of Intelligence*, *8*(1).
[doi:10.3390/jintelligence8010005](https://doi.org/10.3390/jintelligence8010005)

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`fit_statistic_pcm`](https://pgmj.github.io/easyRaschBayes/reference/fit_statistic_pcm.md)
for posterior predictive fit statistics,
[`fit_statistic_rm`](https://pgmj.github.io/easyRaschBayes/reference/fit_statistic_rm.md)
for posterior predictive fit statistics,
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
for Bayesian infit/outfit,
[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)
for Bayesian Q3 residual correlations,
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

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

# Item-restscore association
irs <- item_restscore_statistic(
  model      = fit_pcm,
  item_var   = item,
  person_var = id,
  ndraws_use = 500
)
#> Error: object 'fit_pcm' not found

# Flag items with too-strong discrimination (ppp > 0.95)
irs %>% filter(ppp > 0.95)
#> Error: object 'irs' not found

# Flag items with too-weak discrimination (ppp < 0.05)
irs %>% filter(ppp < 0.05)
#> Error: object 'irs' not found

# --- Dichotomous Rasch Model ---

df_rm <- eRm::rainger %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")
#> Error: 'rainger' is not an exported object from 'namespace:eRm'

fit_rm <- brm(
  response ~ 1 + (1 | item) + (1 | id),
  data   = df_rm,
  family = bernoulli(),
  chains = 4,
  cores  = 4,
  iter   = 2000
)
#> Error: object 'df_rm' not found

irs_rm <- item_restscore_statistic(
  model      = fit_rm,
  item_var   = item,
  person_var = id,
  ndraws_use = 500
)
#> Error: object 'fit_rm' not found

irs_rm %>%
  arrange(ppp)
#> Error: object 'irs_rm' not found
# }
```
