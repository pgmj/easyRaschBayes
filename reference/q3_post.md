# Summarize and Plot Posterior Predictive Q3 Residual Correlations

Postprocesses the output of
[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)
to produce summary tables of posterior predictive Q3 statistics and a
combined slab + interval plot comparing observed Q3 residual
correlations to the posterior predictive distribution for each item
pair.

## Usage

``` r
q3_post(
  q3_draws,
  ci = 0.84,
  n_pairs = NULL,
  sort_by = c("q3_diff", "q3_obs", "ppp")
)
```

## Arguments

- q3_draws:

  A data frame (or tibble) as returned by
  [`q3_statistic`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)
  containing at minimum the columns `item_pair`, `q3` (observed Q3 per
  draw), and `q3_rep` (replicated Q3 per draw).

- ci:

  Numeric in \\(0, 1)\\. Width of the highest density continuous
  interval of the replicated (posterior predictive) Q3 distribution,
  used in the `hdi` table and the slab display. Default is 0.84.

- n_pairs:

  Integer. Maximum number of item pairs to display in the plot, selected
  by largest absolute `q3_diff`. If `NULL` (the default), all pairs are
  shown. Useful when the number of item pairs is large.

- sort_by:

  Character. How to order item pairs on the y-axis. `"q3_diff"` (the
  default) sorts by the posterior mean difference between observed and
  replicated Q3 (largest at top). `"q3_obs"` sorts by the posterior mean
  observed Q3. `"ppp"` sorts by the posterior predictive p-value.

## Value

A list with three elements:

- summary:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item pair containing: `item_pair`, `item_1`, `item_2`,
  `q3_obs` (posterior mean observed Q3), `q3_rep` (posterior mean
  replicated Q3), `q3_diff` (posterior mean difference), `q3_ppp`
  (posterior predictive p-value: proportion of draws where the observed
  Q3 exceeds the replicated Q3), and `q3_p_adj` (one-sided p-value for
  local dependence, adjusted for the number of item pairs; small values
  flag dependence).

- hdi:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item pair containing: `item_pair`, `item_1`, `item_2`,
  `ld` (proportion of draws where observed Q3 exceeds the upper HDI
  bound of the replicated distribution), and `lr` (proportion of draws
  where observed Q3 falls below the lower HDI bound). See Details for
  how to read these.

- plot:

  A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object showing the posterior predictive distribution of replicated Q3
  (grey filled slab) overlaid with the observed Q3 distribution
  (coloured slab + interval), with a dashed reference line at the mean
  replicated Q3 over all pairs (about \\-1/(k-1)\\ for the conditional
  method, 0 for the theta method).

## Details

Two complementary summary tables are provided:

- `summary`:

  Reports posterior mean observed and replicated Q3 values alongside the
  posterior predictive p-value (ppp). The ppp is the proportion of draws
  where the observed Q3 exceeds the replicated Q3. Under good fit (no
  local dependence), the ppp should be near 0.5. A ppp close to 1
  indicates the observed correlation is systematically higher than
  expected (local dependence); a ppp close to 0 indicates it is
  systematically lower (local repulsion, e.g., speed-accuracy
  tradeoffs). `q3_p_adj` is a single-step max-T adjustment over item
  pairs, as `p_adj` in
  [`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_statistic.md),
  but one-sided: each pair's mean observed Q3 is standardised by the
  mean and SD of its replicated Q3, and compared with the largest
  standardised replicated value over all pairs in each draw. Unlike
  `q3_ppp`, small values flag local dependence.

- `hdi`:

  Reports the proportion of draws where the observed Q3 falls above
  (local dependence) or below (local repulsion) the HDI of the
  replicated distribution. The HDI is computed over all draws, so each
  draw's observed value is compared to the same fixed interval. The
  observed Q3 varies across draws only through parameter uncertainty,
  while the replicated Q3 also varies through the sampling of new
  responses. The observed distribution is therefore narrower, and these
  proportions move quickly from near 0 to near 1 as the observed values
  pass an interval bound. They are exceedance proportions, not posterior
  probabilities of local dependence.

With `q3_statistic(method = "conditional")`, the default, `q3_ppp` was
close to uniform in simulations under a fitting Rasch model, and
`q3_p_adj` kept the proportion of datasets with any pair flagged at .05
near or below .05. With `method = "theta"`, `q3_ppp` is conservative,
and so is `q3_p_adj` (2\\ datasets with any pair flagged at .05). See
[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md).

**Reading `q3_p_adj`.** For the pair with the smallest value, `q3_p_adj`
is a test of local independence over all pairs, and it keeps the
familywise error rate at .05 when no pair is dependent. The values for
the other pairs do not have this guarantee when some pair is dependent,
since all pairs are computed from the same residuals (Ferreira, 2024,
Example 2). The suggested use is to act on the most deviant pair, for
example by removing an item, combining the two items into a superitem or
modelling the dependence, then refit the model and recompute Q3
(Christensen et al., 2017). The adjustment covers the Q3 pairs of one
call, not other fit statistics computed on the same data.

Each pair is compared with its own mean replicated Q3 rather than with
the mean over all pairs. Given the total score, expected Q3 differs
between pairs, and is highest for pairs of the easiest and hardest
items. In simulations, a reference at the mean over all pairs made such
a pair the most deviant one in 11\\ dependence, against about 2\\

## References

Christensen, K. B., Makransky, G. & Horton, M. (2017). Critical values
for Yen's Q3: Identification of local dependence in the Rasch model
using residual correlations. *Applied Psychological Measurement*,
*41*(3), 178–194.
[doi:10.1177/0146621616677520](https://doi.org/10.1177/0146621616677520)

Ferreira, J. A. (2024). Methods of testing a 'small' or 'moderate'
number of hypotheses simultaneously: An account focusing on the control
of the probability of at least one incorrect rejection and of the false
discovery rate. *Journal of Statistical Theory and Practice*, *19*, 6.
[doi:10.1007/s42519-024-00412-4](https://doi.org/10.1007/s42519-024-00412-4)

The plot uses two layers from the ggdist package:

- `stat_slab`:

  Displays the posterior predictive (replicated) Q3 distribution as a
  filled density slab per item pair, shaded by credible interval level.

- `stat_slabinterval`:

  Displays the observed Q3 distribution per item pair as a
  semi-transparent slab with point and interval summaries.

## See also

[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md),
[`infit_post`](https://pgmj.github.io/easyRaschBayes/reference/infit_post.md),
[`item_restscore_post`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_post.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(brms)
library(ggplot2)

# Assuming fit_pcm is a fitted brmsfit object
q3_draws <- q3_statistic(fit_pcm, ndraws_use = 500)

result <- q3_post(q3_draws)
result$summary
result$hdi
result$plot

# Show only top 10 pairs by Q3 difference
result_top <- q3_post(q3_draws, n_pairs = 10)
result_top$plot
} # }
```
