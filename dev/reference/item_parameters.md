# Extract Item Parameters from a Bayesian Rasch Model

Extracts item difficulty (threshold) parameters from a fitted Bayesian
Rasch model. Returns a simple location table in both long and wide
formats, a full summary with posterior SEs and HDCIs, item-level
information, threshold ordering diagnostics, and optionally the full
posterior draws matrix.

## Usage

``` r
item_parameters(
  model,
  item_var = item,
  person_var = id,
  draws = FALSE,
  center = TRUE,
  prob = 0.95
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object. Supported parameterisations:

  Polytomous ordinal (PCM)

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

- draws:

  Logical. If `TRUE`, a draws matrix of full posterior draws is included
  in the output. Default is `FALSE`.

- center:

  Logical. If `TRUE` (the default), item parameters are shifted so that
  the grand mean of all threshold locations is zero, matching the
  frequentist CML convention and the centering used in
  [`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).

- prob:

  Numeric in \\(0, 1)\\. Width of the highest density continuous
  interval (HDCI) reported in the summary. Default is 0.95.

## Value

A list with the following elements:

- locations:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) in
  long format with one row per item (dichotomous) or per item-threshold
  (polytomous), containing `item`, `threshold` (integer, always 1 for
  dichotomous), and `location` (posterior mean on the logit scale).

- locations_wide:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) in
  wide format with one row per item, sorted by mean location. For
  polytomous models, threshold columns are named `t1`, `t2`, etc., and a
  `location` column gives the mean across thresholds. For dichotomous
  models, only `item` and `location` columns are present.

- summary:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  extending the long `locations` with: `se` (posterior SD), `hdci_lower`
  and `hdci_upper` (highest density continuous interval bounds at the
  level specified by `prob`), and `n_eff` (effective sample size for the
  parameter).

- item_information:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item containing: `item`, `location` (mean item location,
  i.e. mean of thresholds for polytomous items), `info_at_location`
  (Fisher information at the item's own location), and `max_info`
  (maximum Fisher information across theta). For Rasch dichotomous
  items, both are 0.25. For polytomous items, information depends on the
  number and spacing of thresholds.

- threshold_order:

  (Polytomous models only) A
  [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item containing: `item`, `n_thresholds`, `ordered`
  (logical: are all thresholds in ascending order?), and
  `prob_disordered` (posterior probability that at least one pair of
  adjacent thresholds is disordered, i.e. \\\tau\_{k+1} \le \tau_k\\).
  `NULL` for dichotomous models.

- person_sd:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row containing: `mean`, `sd`, `hdci_lower`, `hdci_upper` — the
  posterior summary of the person-level standard deviation parameter
  \\\sigma\_\theta\\.

- draws_matrix:

  (Only if `draws = TRUE`) A numeric matrix with rows = thresholds
  (named `"item[threshold]"`) and columns = posterior draws. For
  dichotomous models, row names are item labels only.

## Details

**Dichotomous models** with random item effects
(`response ~ 1 + (1 | item) + (1 | id)`) parameterise item difficulty as
\\\delta_i = -(b_0 + r_i)\\ where \\b_0\\ is the global intercept and
\\r_i\\ is the item random effect. Models with fixed item effects
(`response ~ 0 + item + (1 | id)`) parameterise difficulty as \\\delta_i
= -b_i\\.

**Polytomous (acat/PCM) models** with grouped thresholds
(`thres(gr = item)`) directly estimate item-specific threshold
parameters. Each row in the long output represents one threshold within
one item; each row in the wide output represents one item.

**Item information** is computed from the posterior mean item parameters
using the standard Rasch/PCM information formulae:

- Dichotomous:

  \\I_i(\theta) = P_i(\theta) Q_i(\theta)\\ where \\P_i =
  \text{logistic}(\theta - \delta_i)\\.

- Polytomous (PCM):

  \\I_i(\theta) = \sum_c (c - E_i)^2 P\_{ic}(\theta)\\ where \\E_i =
  \sum_c c \cdot P\_{ic}\\ is the expected score for item \\i\\.

**Threshold ordering**: In the partial credit model, disordered
thresholds (\\\tau\_{k+1} \le \tau_k\\) indicate that the probability of
responding in the intermediate category never exceeds both adjacent
categories — the category is empirically "absorbed". This does not
necessarily indicate misfit (see Adams et al., 2012), but may suggest
response categories should be collapsed. The `prob_disordered` column
reports the posterior probability of at least one disordered pair per
item, providing a Bayesian alternative to post-hoc threshold checks.

## References

Adams, R. J., Wu, M. L., & Wilson, M. (2012). The Rasch rating model and
the disordered threshold controversy. *Educational and Psychological
Measurement*, *72*(4), 547–573.
[doi:10.1177/0013164411432166](https://doi.org/10.1177/0013164411432166)

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`person_parameters`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md),
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md),
[`plot_ipf`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_ipf.md),
[`posterior_to_prior`](https://pgmj.github.io/easyRaschBayes/dev/reference/posterior_to_prior.md).

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
  chains = 4, cores = 2, iter = 1000 # use more iter and cores
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

ip <- item_parameters(fit_pcm)
#> Error: object 'fit_pcm' not found

# Long format: one row per threshold
ip$locations
#> Error: object 'ip' not found

# Wide format: one row per item, easy to scan
ip$locations_wide
#> Error: object 'ip' not found

# Full summary with SE and HDCI
ip$summary
#> Error: object 'ip' not found

# Item information
ip$item_information
#> Error: object 'ip' not found

# Threshold ordering diagnostic
ip$threshold_order
#> Error: object 'ip' not found

# Person SD
ip$person_sd
#> Error: object 'ip' not found

# With full posterior draws
ip_draws <- item_parameters(fit_pcm, draws = TRUE)
#> Error: object 'fit_pcm' not found
dim(ip_draws$draws_matrix)
#> Error: object 'ip_draws' not found

# --- Dichotomous Rasch ---

df_rm <- eRm::raschdat3 %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

fit_rm <- brm(
  response ~ 1 + (1 | item) + (1 | id),
  data   = df_rm,
  family = bernoulli(),
  chains = 4, cores = 2, iter = 1000 # use more iter and cores
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

ip_rm <- item_parameters(fit_rm)
#> Error: object 'fit_rm' not found
# Wide and long are equivalent for dichotomous:
ip_rm$locations
#> Error: object 'ip_rm' not found
ip_rm$locations_wide
#> Error: object 'ip_rm' not found
ip_rm$summary
#> Error: object 'ip_rm' not found
# }
```
