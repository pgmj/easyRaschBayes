# Summarize and Plot Posterior Predictive Item-Restscore Associations

Postprocesses the output of
[`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_statistic.md)
to produce a summary table and a slab plot comparing observed
item-restscore gamma associations to the posterior predictive
distribution.

## Usage

``` r
item_restscore_post(item_restscore)
```

## Arguments

- item_restscore:

  A list as returned by
  [`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_statistic.md),
  containing at minimum:

  result

  :   A data frame with columns including `item` and summary statistics
      (first 5 columns are used).

  draws

  :   A data frame with columns `item`, `gamma` (observed gamma per
      draw), and `gamma_rep` (replicated gamma per draw).

## Value

A list with two elements:

- summary:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  the first 5 columns of the result table, rounded to 3 decimal places
  and sorted by item.

- plot:

  A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object showing the posterior predictive distribution of replicated
  gamma (grey filled slab) with the observed gamma values overlaid as
  orange diamond points.

## Details

The item-restscore gamma association measures the strength of the
relationship between each item's responses and the rest score (total
score excluding that item). Under good fit, the observed gamma should
fall within the posterior predictive distribution.

The plot displays:

- Grey slab:

  The posterior predictive distribution of replicated gamma values,
  shaded by 84\\ interval levels.

- Orange diamonds:

  The observed gamma values per draw, plotted as points on top of the
  slab.

Items where the observed gamma (orange) falls consistently outside the
replicated distribution (grey) indicate poor fit in terms of item
discrimination.

## See also

[`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_statistic.md),
[`infit_post`](https://pgmj.github.io/easyRaschBayes/reference/infit_post.md).

## Examples

``` r
if (FALSE) { # \dontrun{
library(brms)
library(ggplot2)

# Assuming fit_pcm is a fitted brmsfit object
irs <- item_restscore_statistic(fit_pcm)

result <- item_restscore_post(irs)
result$summary
result$plot
} # }
```
