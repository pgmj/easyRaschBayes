# Posterior Predictive Q3 Residual Correlations for Bayesian IRT Models

Computes a Bayesian analogue of Yen's Q3 statistic (Yen, 1984) for
detecting local dependence between item pairs in Rasch-family models
fitted with brms. For each posterior draw, residual correlations are
computed for both observed and replicated data, yielding a posterior
predictive p-value for each item pair that is automatically calibrated
without requiring knowledge of the sampling distribution.

## Usage

``` r
q3_statistic(model, item_var = item, person_var = id, ndraws_use = NULL)
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

- item_1:

  First item in the pair.

- item_2:

  Second item in the pair.

- q3_obs:

  Posterior mean of the observed Q3 residual correlation.

- q3_rep:

  Posterior mean of the replicated Q3 residual correlation.

- q3_diff:

  Posterior mean of `q3_obs - q3_rep`. Large positive values indicate
  that the observed residual correlation exceeds what the model expects.

- ppp:

  Posterior predictive p-value: `mean(q3_obs > q3_rep)` across draws.
  Values close to 1 indicate local dependence (observed correlation
  systematically higher than replicated).

- q3_obs_q025, q3_obs_q975:

  2.5\\ (95\\ observed Q3.

- q3_obs_q005, q3_obs_q995:

  0.5\\ (99\\ observed Q3.

- q3_diff_q025, q3_diff_q975:

  2.5\\ (95\\ Q3 differences.

- q3_diff_q005, q3_diff_q995:

  0.5\\ (99\\ Q3 differences.

## Details

The procedure works as follows for each posterior draw \\s\\:

1.  Compute expected values \\E^{(s)}\_{vi}\\ from the category
    probabilities returned by
    [`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html).
    For ordinal models: \\E^{(s)}\_{vi} = \sum_c c \cdot P^{(s)}(X\_{vi}
    = c)\\. For binary models: \\E^{(s)}\_{vi} = P^{(s)}(X\_{vi} = 1)\\.

2.  Compute observed residuals: \\d^{(s)}\_{vi} = X\_{vi} -
    E^{(s)}\_{vi}\\.

3.  Compute replicated residuals: \\d^{rep(s)}\_{vi} =
    Y^{rep(s)}\_{vi} - E^{(s)}\_{vi}\\, where \\Y^{rep}\\ is drawn via
    [`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

4.  For each item pair \\(i, j)\\, compute Q3 as the Pearson correlation
    of residuals across all persons who responded to both items.

5.  Aggregate across draws to obtain posterior means, quantiles, and
    posterior predictive p-values.

The key advantage over parametric bootstrapping is that the reference
distribution is obtained directly from the posterior, automatically
accounting for parameter uncertainty and the negative bias inherent in
Q3 (which depends on test length and person ability distribution).

## References

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement*, *8*(2), 125–145.

Christensen, K. B., Makransky, G. & Horton, M. (2017). Critical values
for Yen's Q3: Identification of local dependence in the Rasch model
using residual correlations. *Applied Psychological Measurement*,
*41*(3), 178–194.
[doi:10.1177/0146621616677520](https://doi.org/10.1177/0146621616677520)

Bürkner, P.-C. (2020). Analysing Standard Progressive Matrices (SPM-LS)
with Bayesian Item Response Models. *Journal of Intelligence*, *8*(1).
[doi:10.3390/jintelligence8010005](https://doi.org/10.3390/jintelligence8010005)

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`fit_statistic_pcm`](https://pgmj.github.io/easyRaschBayes/dev/reference/fit_statistic_pcm.md)
for posterior predictive fit statistics with user-supplied criterion
functions,
[`fit_statistic_rm`](https://pgmj.github.io/easyRaschBayes/dev/reference/fit_statistic_rm.md)
for posterior predictive fit statistics with user-supplied criterion
functions,
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
for Bayesian infit/outfit,
[`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
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
  cores  = 1, # use more cores if you have
  iter   = 500 # use at least 2000 
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

# Q3 residual correlations
q3_results <- q3_statistic(
  model      = fit_pcm,
  ndraws_use = 100 # use at least 500
)
#> Error: object 'fit_pcm' not found

# Flag item pairs with ppp > 0.95 as locally dependent
q3_results %>%
  filter(ppp > 0.95) %>%
  arrange(desc(q3_diff))
#> Error: object 'q3_results' not found

# Inspect 99% credible intervals for Q3 differences
q3_results %>%
  filter(q3_diff_q005 > 0)
#> Error: object 'q3_results' not found

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

q3_rm <- q3_statistic(
  model      = fit_rm,
  ndraws_use = 100 # use at least 500
)
#> Error: object 'fit_rm' not found

q3_rm %>%
  filter(ppp > 0.95)
#> Error: object 'q3_rm' not found
# }
```
