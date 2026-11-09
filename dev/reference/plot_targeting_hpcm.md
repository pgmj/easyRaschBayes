# Person-Item Targeting Plot for the Hurdle Partial Credit Model

Builds two person-item targeting plots, one for each part of a hurdle
partial credit model fitted with the
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/dev/reference/hurdle_acat.md)
custom family, with the same design as
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md):
a person histogram, an inverted histogram of the item locations, and a
bottom panel with response-category bands or a dot-and-whisker plot of
the thresholds. The presence part \\(\theta\_{hurdle},
\delta\_{hurdle})\\ and the severity part \\(\theta\_{pcm},
\tau\_{ik})\\ are different latent variables, so each gets its own plot
and axis.

## Usage

``` r
plot_targeting_hpcm(
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
  object using the
  [`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/dev/reference/hurdle_acat.md)
  custom family.

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data. Default is `item`.

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data. Default is `id`.

- robust:

  Logical. If `TRUE`, the histogram panels mark the median and MAD
  instead of the mean and SD. Default is `FALSE`.

- center:

  Logical. If `TRUE` (the default), each part is shifted so that its
  mean item location is zero, with the person locations shifted by the
  same constant, as in
  [`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters_hpcm.md)
  and
  [`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters_hpcm.md).

- sort_items:

  Character. `"data"` (the default) keeps the data order, first item at
  the top. `"location"` sorts by mean item location, easiest item at the
  top.

- bins:

  Integer. Number of histogram bins. Default as in
  [`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).

- prob:

  Numeric in \\(0, 1)\\. Width of the credible intervals of the item
  locations and of the intervals in the returned tables. Default is
  0.95. `NULL` hides the intervals in the plots.

- palette:

  An optional character vector of colours, one per response category.
  Overrides the viridis colours.

- person_fill, threshold_fill:

  Fill colours for the person and item location histograms.

- height_ratios:

  Numeric vector of length 3 giving the relative heights of the three
  panels. Default is `c(3, 2, 5)`.

- panel:

  Character. `"categories"` (the default) draws each item as a bar
  partitioned into response-category bands. `"thresholds"` draws a
  dot-and-whisker plot.

- xlim:

  Numeric vector of length 2. Initial limits of the x-axis, widened to
  whole numbers if needed. Default is `c(-4, 4)`.

- category_labels:

  Optional list with elements `hurdle` (two labels) and/or `pcm` (one
  label per severity category, starting at category 1). Defaults are
  `c("0", "1+")` for the hurdle part and the category scores for the
  severity part.

- viridis_option, viridis_begin, viridis_end:

  Viridis palette for the category colours, as in
  [`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).

- row_gap:

  Numeric. Vertical spacing between item rows in the `"categories"`
  panel. Default `NULL` chooses it from the number of collapsed
  categories and disordered thresholds.

- output:

  Character. `"patchwork"` (the default) returns the two combined plots.
  `"list"` returns, for each part, the three panels and the posterior
  summaries of the category structure.

## Value

With `output = "patchwork"`, a list with elements `hurdle` and `pcm`,
each a patchwork object. With `output = "list"`, the same two elements,
each a list with `p1`, `p2`, `p3`, `categories` and `disordered` as
described in
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).
In the severity part, categories are numbered from 1.

## Details

**Hurdle part.** The person location is the posterior mean of the brms
random effect on `hu` with its sign flipped, so that higher values mean
a symptom is more likely to be present. The item location is the hurdle
item effect, where higher values mean more zeros. Then \\P(Y\_{vi} \> 0)
= \mathrm{plogis}(\theta\_{hurdle, v} - \delta\_{hurdle, i})\\, a
dichotomous Rasch model, and each item has two bands: 0 and any response
above 0.

**Severity part.** The person location is the posterior mean of the
random intercept, and the thresholds are those of the partial credit
model for the responses above 0. The bands show the severity categories
1 and up. Persons who answered 0 on every item gave no severity
responses, so their severity location comes from the prior and the
correlation with the hurdle location only. They are left out of the
severity person histogram, and the caption gives their number.

Category bands, red marks for categories that are never the most likely
response and for disordered thresholds, and their posterior
probabilities are computed as in
[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).

**Combining the two plots.** Each element is a patchwork object:


    plots <- plot_targeting_hpcm(fit)
    patchwork::wrap_plots(plots$hurdle, plots$pcm, ncol = 2)

## References

Wright, B. D. & Stone, M. H. (1979). *Best Test Design*. MESA Press.

Magnus, B. E. & Garnier-Villarreal, M. (2022). A multidimensional
zero-inflated graded response model for ordinal symptom data.
*Psychological Methods*, *27*(2), 261-279.
[doi:10.1037/met0000395](https://doi.org/10.1037/met0000395)

## See also

[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md)
for the standard Rasch models,
[`item_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters_hpcm.md),
[`person_parameters_hpcm`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters_hpcm.md).
