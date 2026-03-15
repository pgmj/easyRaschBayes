# Summarize and Plot Posterior Predictive Infit Statistics

Postprocesses the output of
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
to produce summary tables of posterior predictive infit statistics and a
combined slab + interval plot comparing observed infit values to the
posterior predictive distribution.

## Usage

``` r
infit_post(infit_draws, ci = 0.84)
```

## Arguments

- infit_draws:

  A data frame (or tibble) as returned by
  [`infit_statistic`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
  containing at minimum the columns `item`, `infit` (observed infit per
  draw), and `infit_rep` (replicated infit per draw).

- ci:

  Numeric in \\(0, 1)\\. Width of the credible interval used for the
  posterior predictive HDI and the slab display. Default is 0.84.

## Value

A list with three elements:

- summary:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item containing: `item`, `infit_obs` (posterior mean of
  observed infit), `infit_rep` (posterior mean of replicated infit), and
  `infit_ppp` (posterior predictive p-value: proportion of draws where
  the replicated infit exceeds the observed infit). Values near 0.5
  indicate good fit; values near 0 suggest underfit; values near 1
  suggest overfit.

- hdi:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item containing: `item`, `underfit` (posterior probability
  that the observed infit exceeds the upper HDI bound of the replicated
  distribution), and `overfit` (posterior probability that the observed
  infit falls below the lower HDI bound).

- plot:

  A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object showing the posterior predictive distribution of replicated
  infit (grey filled slab) overlaid with the observed infit distribution
  (coloured slab + interval), with a dashed reference line at 1 (perfect
  fit).

## Details

Two complementary summary tables are provided:

- `summary`:

  Reports the posterior mean observed and replicated infit values
  alongside the posterior predictive p-value (ppp). The ppp is the
  proportion of draws where the replicated infit exceeds the observed
  infit. Under good fit, the ppp should be near 0.5. A ppp near 0
  indicates the observed infit is consistently larger than expected
  (underfit); a ppp near 1 indicates it is consistently smaller
  (overfit).

- `hdi`:

  Reports the probability that the observed infit falls above (underfit)
  or below (overfit) the HDI of the replicated distribution. This
  provides a more distributional assessment than the ppp alone.

The plot uses two layers from the ggdist package:

- `stat_slab`:

  Displays the posterior predictive (replicated) infit distribution as a
  filled density slab per item, shaded by credible interval level.

- `stat_slabinterval`:

  Displays the observed infit distribution per item as a
  semi-transparent slab with point and interval summaries.

Under good model fit, the observed infit distribution should overlap
substantially with the replicated distribution. Items where the observed
distribution sits systematically above the replicated HDI indicate
underfit (more variation than expected); items below indicate overfit
(less variation than expected).

## See also

[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md),
[`item_restscore_post`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_post.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(brms)
library(ggplot2)

# Assuming fit_pcm is a fitted brmsfit object
infit_draws <- infit_statistic(fit_pcm)

result <- infit_post(infit_draws)
result$summary
result$hdi
result$plot
} # }
```
