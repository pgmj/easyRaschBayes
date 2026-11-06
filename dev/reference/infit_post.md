# Summarize and Plot Posterior Predictive Infit Statistics

Postprocesses the output of
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
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
  [`infit_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
  or one element of the output of
  [`infit_statistic_hpcm`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic_hpcm.md),
  containing the columns `infit` (observed infit per draw) and
  `infit_rep` (replicated infit per draw). The item identifier is taken
  from a column named `item` if there is one, and otherwise from the
  first column.

- ci:

  Numeric in \\(0, 1)\\. Width of the highest density continuous
  interval of the replicated (posterior predictive) infit distribution,
  used in the `hdi` table and the slab display. Default is 0.84.

## Value

A list with three elements:

- summary:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item containing: `item`, `infit_obs` (posterior mean of
  observed infit), `infit_rep` (mean of replicated infit across draws),
  and `infit_ppp` (posterior predictive p-value: proportion of draws
  where the observed infit exceeds the replicated infit, counting ties
  as one half). Values near 0.5 indicate good fit, values near 1 suggest
  underfit and values near 0 suggest overfit.

- hdi:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item containing: `item`, `underfit` (proportion of draws
  where the observed infit exceeds the upper HDI bound of the replicated
  distribution), and `overfit` (proportion of draws where the observed
  infit falls below the lower HDI bound). See Details for how to read
  these.

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
  proportion of draws where the observed infit exceeds the replicated
  infit, with ties counted as one half. Under good fit, the ppp should
  be near 0.5. A ppp near 1 indicates the observed infit is consistently
  larger than expected (underfit), and a ppp near 0 that it is
  consistently smaller (overfit). The direction is the same as for the
  ppp of
  [`q3_post`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
  and
  [`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md).

- `hdi`:

  Reports the proportion of draws where the observed infit falls above
  (underfit) or below (overfit) the HDI of the replicated distribution.
  The HDI is computed over all draws, so each draw's observed value is
  compared to the same fixed interval. The observed infit varies across
  draws only through parameter uncertainty, while the replicated infit
  also varies through the sampling of new responses. The observed
  distribution is therefore narrower, and these proportions move quickly
  from near 0 to near 1 as the observed values pass an interval bound.
  They are exceedance proportions, not posterior probabilities of
  misfit.

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

The infit ppp is conservative. Values near 0 or 1 are rarer than for a
uniform p-value, and misfit can go undetected (see
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)).
It is suggested to also evaluate item fit with
[`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
and
[`item_restscore_post`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_post.md).

## See also

[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md),
[`item_restscore_post`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_post.md).

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
