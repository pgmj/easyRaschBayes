# Person-Item Map (Targeting Plot) for Bayesian IRT Models

Plots a person-item map (also known as a Wright map or targeting plot)
for Bayesian IRT models fitted with brms. The plot consists of three
vertically stacked panels sharing the same latent variable (theta /
logit) x-axis:

## Usage

``` r
plot_targeting(
  model,
  item_var = item,
  person_var = id,
  robust = FALSE,
  center = TRUE,
  sort_items = c("data", "location"),
  bins = 30,
  prob = 0.95,
  palette = NULL,
  person_fill = "#0072B2",
  threshold_fill = "#D55E00",
  height_ratios = c(3, 2, 5)
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

  Logical. If `TRUE` (the default), the scale is recentered so that the
  grand mean of all item threshold locations is zero, following the
  convention in frequentist Rasch analysis. Person estimates are shifted
  by the same constant. If `FALSE`, the raw brms parameterisation is
  used.

- sort_items:

  Character. How to order items on the y-axis of the bottom panel.
  `"data"` (the default) preserves the order in which items first appear
  in the model data, with the first item at the top. `"location"` sorts
  items by their mean threshold location (easiest at top, hardest at
  bottom).

- bins:

  Integer. Number of bins for both histograms. Default is 30.

- prob:

  Numeric in \\(0, 1)\\. Width of the credible intervals for the item
  threshold whiskers. Default is 0.95.

- palette:

  An optional character vector of colors for the response categories. If
  `NULL` (the default), the `viridis` discrete scale is used.

- person_fill:

  Fill color for the person histogram. Default is `"#0072B2"` (blue).

- threshold_fill:

  Fill color for the threshold histogram. Default is `"#D55E00"`
  (vermillion).

- height_ratios:

  Numeric vector of length 3 specifying the relative heights of the top
  (person), middle (threshold), and bottom (dot-whisker) panels. Default
  is `c(3, 2, 5)`.

## Value

A `patchwork` object (combined `ggplot`).

## Details

1.  **Top**: A histogram of person ability estimates, with a reference
    line for the mean (or median) and shading for ±1 SD (or ±1 MAD).

2.  **Middle**: An inverted histogram of item threshold locations, with
    a reference line for the mean (or median) and shading for ±1 SD (or
    ±1 MAD), mirroring the top panel to visualise the overlap between
    person abilities and item difficulties.

3.  **Bottom**: A dot-and-whisker plot of item thresholds by item, with
    credible intervals and color-coded response categories.

Together, the top and middle panels form a half-moon (or back-to-back
histogram) display that makes it easy to assess whether the test is
well-targeted to the sample.

**Person estimates** are obtained as the posterior means of the person
random effects from the fitted model via
[`ranef`](https://rdrr.io/pkg/nlme/man/random.effects.html).

**Item thresholds** are extracted from the posterior draws. For models
with grouped thresholds (`thres(gr = item)`), each item has its own set
of threshold parameters. For models with a single set of thresholds
(e.g., dichotomous Rasch with `(1 | item)`), the item random effects are
subtracted from the global thresholds to obtain item-specific locations.

When `center = TRUE` (the default), the grand mean of all item threshold
posterior means is computed and subtracted from every threshold
estimate, its credible interval bounds, and every person estimate. This
is a uniform translation of the entire scale that preserves all relative
distances and matches the zero-centered item difficulty convention used
in frequentist CML estimation.

## References

Wright, B. D. & Stone, M. H. (1979). *Best Test Design*. MESA Press.

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`plot_ipf`](https://pgmj.github.io/easyRaschBayes/reference/plot_ipf.md)
for item category probability curves,
[`ranef`](https://rdrr.io/pkg/nlme/man/random.effects.html),
[`as_draws_df`](https://mc-stan.org/posterior/reference/draws_df.html).

## Examples

``` r
if (FALSE) { # \dontrun{
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
  cores  = 4,
  iter   = 2000
)

# Default: centered, mean ± SD, items in data order
plot_targeting(fit_pcm)

# Uncentered (raw brms parameterisation)
plot_targeting(fit_pcm, center = FALSE)

# Robust: median ± MAD, items sorted by location
plot_targeting(fit_pcm, robust = TRUE, sort_items = "location")

# --- Dichotomous Rasch Model ---

df_rm <- eRm::rainger %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

fit_rm <- brm(
  response ~ 1 + (1 | item) + (1 | id),
  data   = df_rm,
  family = bernoulli(),
  chains = 4,
  cores  = 4,
  iter   = 2000
)

plot_targeting(fit_rm, sort_items = "location")
} # }
```
