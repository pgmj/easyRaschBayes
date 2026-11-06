# Person-Item Map (Targeting Plot) for Bayesian IRT Models

Produces a three-panel targeting plot with a shared logit x-axis:

## Usage

``` r
plot_targeting(
  model,
  item_var = item,
  person_var = id,
  robust = FALSE,
  center = TRUE,
  sort_items = c("data", "location"),
  bins,
  prob = 0.95,
  palette = NULL,
  person_fill = "#0072B2",
  threshold_fill = "#D55E00",
  height_ratios = c(3, 2, 5),
  panel = c("categories", "thresholds"),
  person_estimate = c("eap", "wle"),
  xlim = c(-4, 4),
  category_labels = NULL,
  viridis_option = "G",
  viridis_begin = 0.9,
  viridis_end = 0.2,
  row_gap = NULL,
  output = c("patchwork", "list")
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object from an ordinal IRT model (e.g., `family = acat`) or a
  dichotomous model (`family = bernoulli()`).

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data. Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data. Default is `id`.

- robust:

  Logical. If `FALSE` (the default), the histogram annotations use mean
  ± SD. If `TRUE`, median ± MAD is used instead.

- center:

  Logical. If `TRUE` (the default), the scale is recentred so that the
  mean of all item threshold locations is zero. Person locations are
  shifted by the same constant. If `FALSE`, the raw brms
  parameterisation is used.

- sort_items:

  Character. How to order items in the bottom panel. `"data"` (the
  default) keeps the order in which items first appear in the model
  data, with the first item at the top. `"location"` sorts items by
  their mean threshold location, easiest at the top.

- bins:

  Integer. Number of bins for both histograms. Default is half the
  maximum possible sum score, rounded up, and at least 11.

- prob:

  Numeric in \\(0, 1)\\, or `NULL`. Width of the credible intervals for
  the item thresholds and of the intervals in the returned `categories`
  table. Default is 0.95. `NULL` hides the threshold intervals in the
  plot.

- palette:

  An optional character vector of colours, one per response category.
  Overrides the viridis colours in both panel types.

- person_fill:

  Fill colour for the person histogram. Default is `"#0072B2"` (blue).

- threshold_fill:

  Fill colour for the threshold histogram, and for the dichotomous
  dot-and-whisker panel. Default is `"#D55E00"` (vermillion).

- height_ratios:

  Numeric vector of length 3 giving the relative heights of the top,
  middle and bottom panels. Default is `c(3, 2, 5)`.

- panel:

  Character. `"categories"` (the default) draws each item as a bar
  partitioned into response-category bands, with the threshold estimates
  and their credible intervals below it. `"thresholds"` draws a
  dot-and-whisker plot of the thresholds.

- person_estimate:

  Character. `"eap"` (the default) uses the posterior means of the
  person random effects. `"wle"` uses Warm's weighted likelihood
  estimates from
  [`person_parameters`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md).
  EAPs are shrunk toward the mean, so their histogram is narrower than
  that of the WLEs.

- xlim:

  Numeric vector of length 2. Initial limits of the shared x-axis,
  widened to whole numbers if any person location, threshold or interval
  falls outside. Default is `c(-4, 4)`.

- category_labels:

  Optional character vector of labels for the response categories, in
  ascending order and one per category. Used in the legend of the
  `"categories"` panel. Default `NULL` uses the category scores,
  starting at 0.

- viridis_option:

  Character. Viridis palette option for the category colours. Default is
  `"G"` (mako).

- viridis_begin, viridis_end:

  Numeric in \\\[0, 1\]\\. Start and end of the viridis palette.
  Defaults `0.9` and `0.2` run from light to dark, so higher categories
  are darker.

- row_gap:

  Numeric. Vertical spacing between item rows in the `"categories"`
  panel. Default `NULL` uses 1, widened when a category collapses or
  thresholds are reversed, to make room for their labels.

- output:

  Character. `"patchwork"` (the default) returns the combined plot.
  `"list"` returns the three panels and the posterior summaries of the
  category structure.

## Value

With `output = "patchwork"`, a `patchwork` object. With
`output = "list"`, a list with:

- p1, p2, p3:

  The person histogram, threshold histogram and bottom panel as `ggplot`
  objects.

- categories:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item and response category: `p_never_modal`, the posterior
  probability that the category is never the most likely response at any
  location, and `lower_bound`, `lower_bound_lo`, `lower_bound_hi`, the
  posterior median and `prob` interval of the location where the
  category starts to be the most likely response, over the draws where
  it is (`NA` for category 0).

- reversals:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per item and adjacent threshold pair: `p_reversed`, the
  posterior probability that the higher threshold lies below the lower
  one.

## Details

1.  **Top**: Histogram of person locations, with a reference line for
    the mean (or median) and shading for ±1 SD (or ±1 MAD).

2.  **Middle**: Inverted histogram of item threshold locations, with the
    same summary annotations.

3.  **Bottom**: One bar per item, either partitioned into
    response-category bands (`panel = "categories"`, the default) or
    drawn as a dot-and-whisker plot of the individual thresholds
    (`panel = "thresholds"`).

Together, the top and middle panels form a back-to-back histogram that
shows whether the test is well targeted to the sample. The bottom panel
places the items on the same scale, so the category bands show which
response is the most likely one at the locations where the persons are.
The layout follows `easyRasch2::RMtargeting()`.

**Item thresholds** come from the posterior draws (see
[`item_parameters`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters.md)).
Points are posterior means and intervals are equal-tailed credible
intervals.

**Category bands.** With `panel = "categories"`, each band spans the
locations at which its response category is the most likely response,
computed from the posterior mean thresholds. When an item's thresholds
are ordered these boundaries are the thresholds themselves. When they
are not, the disordered run is pooled by averaging, and the categories
it skips over, which are never the most likely response at any location,
collapse to a red tick labelled with the category number. Red arrows
below the bar give the size of each threshold reversal in logits. The
two outer bands fade towards the panel edge, since they are open-ended.

**Posterior probabilities.** The same category structure is computed for
every posterior draw. Each red tick is labelled with the posterior
probability that its category is never the most likely response, and
each reversal with the posterior probability that the pair is reversed.
The full tables, including categories and pairs without red marks, are
returned with `output = "list"`.

**Person locations.** `person_estimate = "eap"` uses the posterior means
of the person random effects, shifted by the same constant as the
thresholds. `"wle"` uses the WLEs of
[`person_parameters`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md),
which are finite at extreme scores and use the items each person
answered.

## References

Wright, B. D. & Stone, M. H. (1979). *Best Test Design*. MESA Press.

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`plot_ipf`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_ipf.md)
for item category probability curves,
[`item_parameters`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters.md),
[`person_parameters`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md).

## Examples

``` r
# \donttest{
library(brms)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

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

# Default: category bands, EAP person locations, items in data order
plot_targeting(fit_pcm)
#> Error: object 'fit_pcm' not found

# WLE person locations, items sorted by location, labelled categories
plot_targeting(fit_pcm, person_estimate = "wle", sort_items = "location",
               category_labels = c("Never", "Sometimes", "Often"))
#> Error: object 'fit_pcm' not found

# Dot-and-whisker panel
plot_targeting(fit_pcm, panel = "thresholds")
#> Error: object 'fit_pcm' not found

# Panels and posterior summaries of the category structure
res <- plot_targeting(fit_pcm, output = "list")
#> Error: object 'fit_pcm' not found
res$categories
#> Error: object 'res' not found
res$reversals
#> Error: object 'res' not found

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

plot_targeting(fit_rm, sort_items = "location")
#> Error: object 'fit_rm' not found
# }
```
