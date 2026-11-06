# Summarize and Plot Item-Restscore Associations

Postprocesses the output of
[`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
to produce a summary table and a slab plot comparing observed
item-restscore gamma associations to the distribution of gamma in
datasets replicated conditional on the total scores.

## Usage

``` r
item_restscore_post(item_restscore)
```

## Arguments

- item_restscore:

  A list as returned by
  [`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md),
  with elements `result` and `draws`. The first column of both is the
  item identifier.

## Value

A list with two elements:

- summary:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item and the columns item, `gamma_obs`, `gamma_exp`,
  `gamma_diff`, `ppp` and `ppp_adj`, rounded to 3 decimals and sorted by
  item.

- plot:

  A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object showing the distribution of replicated gamma (grey slab) with
  the observed gamma overlaid as an orange diamond.

## Details

The item-restscore gamma association measures the strength of the
relationship between each item's responses and the rest score (total
score excluding that item). Under good fit, the observed gamma should
fall within the distribution of replicated gamma.

The plot displays:

- Grey slab:

  The distribution of replicated gamma values, shaded at the 84% and 95%
  interval levels.

- Orange diamond:

  The observed gamma.

Items where the observed gamma (orange) falls outside the replicated
distribution (grey) indicate poor fit in terms of item discrimination.
See
[`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
for the definitions of the summary columns.

## See also

[`item_restscore_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md),
[`infit_post`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md).

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
