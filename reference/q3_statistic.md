# Posterior Predictive Q3 Residual Correlations for Bayesian IRT Models

Computes a Bayesian analogue of Yen's Q3 statistic (Yen, 1984) for
detecting local dependence between item pairs in Rasch-family models
fitted with brms. For each posterior draw, residual correlations are
computed for both observed and replicated data, yielding draw-level Q3
values that can be summarized and visualized via
[`q3_post`](https://pgmj.github.io/easyRaschBayes/reference/q3_post.md).
By default (`method = "conditional"`) the residuals are taken from the
expected responses given each person's total score, and the replicated
data keep each person's observed total score. With `method = "theta"`
they are computed given each draw's person locations instead.

## Usage

``` r
q3_statistic(
  model,
  item_var = item,
  person_var = id,
  ndraws_use = NULL,
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
  the model data. Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data. Default is `id`.

- ndraws_use:

  Optional positive integer. If specified, a random subset of posterior
  draws of this size is used. If `NULL` (the default), all draws are
  used.

- method:

  Either `"conditional"` (the default), with expected responses
  conditional on the total score, or `"theta"`, with expected responses
  given each draw's person locations, as in easyRaschBayes 0.3.0.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) in long
format with one row per draw per item pair, containing:

- draw:

  Integer draw index.

- item_pair:

  Character label of the item pair (`"item1 : item2"`).

- item_1:

  First item in the pair.

- item_2:

  Second item in the pair.

- q3:

  Observed Q3 residual correlation for this draw.

- q3_rep:

  Replicated Q3 residual correlation for this draw.

This long-format output parallels the structure of
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
and can be passed directly to
[`q3_post`](https://pgmj.github.io/easyRaschBayes/reference/q3_post.md)
for summary tables and plots.

## Details

**Conditional method.** Under a Rasch-family model the total score is
sufficient for the person location, so the distribution of a person's
responses given the total does not depend on the location. For each
posterior draw \\s\\, the expected response \\E^{(s)}(X\_{vi} \mid
R_v)\\ given the total score \\R_v\\ is computed from the item
parameters of that draw with elementary symmetric functions. The
observed residuals are \\X\_{vi} - E^{(s)}(X\_{vi} \mid R_v)\\.
Replicated responses are drawn given each person's total score, and
their residuals use the same expected values. Persons with the minimum
or maximum total score have no residual variation and do not contribute.

Given the total score, a person's residuals sum to 0, so the observed Q3
values centre on about \\-1/(k-1)\\ with \\k\\ items, not on 0
(Christensen et al., 2017). The replicated Q3 values have the same
centre, so the comparison of observed and replicated Q3 is not affected,
and `q3_diff` in
[`q3_post`](https://pgmj.github.io/easyRaschBayes/reference/q3_post.md)
centres on 0.

**Theta method.** For each draw, the expected responses
\\E^{(s)}\_{vi}\\ are computed from the category probabilities of
[`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
given the draw's person locations. The observed residuals are
\\X\_{vi} - E^{(s)}\_{vi}\\, and the replicated residuals use responses
from
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
for the same draw. Because the person locations are posterior draws
rather than point estimates, the observed Q3 values centre on 0 under
the model, as the replicated values do.

With both methods, Q3 for an item pair is the Pearson correlation of the
residuals across all persons who responded to both items.

**Calibration.** In simulations under a fitting Rasch model (dichotomous
and partial credit models, 300 and 1000 persons, 8 or 10 items), the
conditional ppp from
[`q3_post`](https://pgmj.github.io/easyRaschBayes/reference/q3_post.md)
was close to uniform: 4.2\\ theta method was conservative, with 1.9\\
locations are drawn given the observed responses. With many item pairs,
some pair will usually exceed .95 by chance. The adjusted `q3_p_adj` in
[`q3_post`](https://pgmj.github.io/easyRaschBayes/reference/q3_post.md)
accounts for this. Under the conditional method, the proportion of
simulated datasets with any pair flagged by it at .05 was 3.5\\
condition, pooled).

Replication uses random numbers, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible results.

## References

Yen, W. M. (1984). Effects of local item dependence on the fit and
equating performance of the three-parameter logistic model. *Applied
Psychological Measurement*, *8*(2), 125–145.
[doi:10.1177/014662168400800201](https://doi.org/10.1177/014662168400800201)

Christensen, K. B., Makransky, G. & Horton, M. (2017). Critical values
for Yen's Q3: Identification of local dependence in the Rasch model
using residual correlations. *Applied Psychological Measurement*,
*41*(3), 178–194.
[doi:10.1177/0146621616677520](https://doi.org/10.1177/0146621616677520)

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`q3_post`](https://pgmj.github.io/easyRaschBayes/reference/q3_post.md)
for postprocessing summaries and plots,
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
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

df_pcm <- phq9[1:200, 1:9] %>%
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
q3_draws <- q3_statistic(fit_pcm, ndraws_use = 500)
#> Error: object 'fit_pcm' not found

# Postprocess
result <- q3_post(q3_draws)
#> Error: object 'q3_draws' not found
result$summary
#> Error: object 'result' not found
result$hdi
#> Error: object 'result' not found
result$plot
#> Error: object 'result' not found

# --- Dichotomous Rasch Model ---

# PHQ-9 items dichotomised (any symptom vs none), for illustration only
df_rm <- phq9[1:200, 1:9] %>%
  mutate(across(everything(), ~ as.integer(.x > 0))) %>%
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

q3_draws <- q3_statistic(fit_rm, ndraws_use = 500)
#> Error: object 'fit_rm' not found

# Postprocess
result <- q3_post(q3_draws)
#> Error: object 'q3_draws' not found
result$summary
#> Error: object 'result' not found
result$hdi
#> Error: object 'result' not found
result$plot
#> Error: object 'result' not found
# }
```
