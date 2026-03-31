# Summarize and Plot Posterior Predictive Q3 Residual Correlations

Postprocesses the output of
[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
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
  [`q3_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
  containing at minimum the columns `item_pair`, `q3` (observed Q3 per
  draw), and `q3_rep` (replicated Q3 per draw).

- ci:

  Numeric in \\(0, 1)\\. Width of the credible interval used for the
  posterior predictive HDI and the slab display. Default is 0.84.

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
  replicated Q3), `q3_diff` (posterior mean difference), and `q3_ppp`
  (posterior predictive p-value: proportion of draws where the observed
  Q3 exceeds the replicated Q3).

- hdi:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item pair containing: `item_pair`, `item_1`, `item_2`,
  `ld` (local dependence probability: proportion of draws where observed
  Q3 exceeds upper HDI bound of replicated distribution), and `lr`
  (local repulsion probability: proportion of draws where observed Q3
  falls below lower HDI bound).

- plot:

  A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object showing the posterior predictive distribution of replicated Q3
  (grey filled slab) overlaid with the observed Q3 distribution
  (coloured slab + interval), with a dashed reference line at 0.

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
  tradeoffs).

- `hdi`:

  Reports the probability that the observed Q3 falls above (local
  dependence) or below (local repulsion) the HDI of the replicated
  distribution. This provides a more distributional assessment than the
  ppp alone.

The plot uses two layers from the ggdist package:

- `stat_slab`:

  Displays the posterior predictive (replicated) Q3 distribution as a
  filled density slab per item pair, shaded by credible interval level.

- `stat_slabinterval`:

  Displays the observed Q3 distribution per item pair as a
  semi-transparent slab with point and interval summaries.

## See also

[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md),
[`infit_post`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md),
[`item_restscore_post`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_post.md).

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
