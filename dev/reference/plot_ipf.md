# Item Category Probability Function Curves for Polytomous IRT Models

Plots item category probability functions (ICPFs) for polytomous
Bayesian IRT models fitted with brms. For each item, the probability of
endorsing each response category is plotted as a function of the latent
variable (theta), with separate colored curves per category. All items
are displayed in a combined faceted plot, similar to the trace plots
produced by `itemplot` in the mirt package.

## Usage

``` r
plot_ipf(
  model,
  item_var = item,
  person_var = id,
  items = NULL,
  theta_range = c(-4, 4),
  n_points = 100,
  ncol = NULL,
  line_size = 0.8,
  ribbon_alpha = 0.15,
  prob = 0.95,
  category_labels = NULL,
  palette = NULL
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object from a polytomous IRT model (e.g., `family = acat` for a
  partial credit model or `family = cumulative` for a graded response
  model).

- item_var:

  An unquoted variable name identifying the item grouping variable in
  the model data (e.g., `item`).

- person_var:

  An unquoted variable name identifying the person grouping variable in
  the model data (e.g., `id`).

- items:

  An optional character vector of item names to plot. If `NULL` (the
  default), all items in the model are plotted.

- theta_range:

  A numeric vector of length 2 specifying the range of the latent
  variable (theta) for the x-axis. Default is `c(-4, 4)`.

- n_points:

  Integer. Number of evenly spaced theta values at which to evaluate the
  category probabilities. Default is 100.

- ncol:

  Integer. Number of columns in the faceted plot layout. If `NULL` (the
  default), an appropriate number is chosen automatically.

- line_size:

  Numeric. Line width for the probability curves. Default is 0.8.

- ribbon_alpha:

  Numeric in \\\[0, 1\]\\. Transparency of the credible interval
  ribbons. Default is 0.15. Set to 0 to hide ribbons.

- prob:

  Numeric in \\(0, 1)\\. Width of the credible interval for the ribbons.
  Default is 0.95.

- category_labels:

  An optional character vector of labels for the response categories. If
  `NULL` (the default), categories are labelled as integers starting
  from 1.

- palette:

  An optional character vector of colors, one per response category. If
  `NULL` (the default), the `viridis` discrete scale from ggplot2 is
  used.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. The plot can be further customised using standard ggplot2
functions.

## Details

The function computes category probabilities directly from the posterior
draws of the item threshold parameters. For the brms `acat` (adjacent
category / partial credit) family with logit link, the density is:
\$\$P(Y = y \| \eta) = \frac{\exp\bigl(\sum\_{k=1}^{y}(\eta -
\tau_k)\bigr)}{\sum\_{k=0}^{K} \exp\bigl(\sum\_{j=1}^{k}(\eta -
\tau_j)\bigr)}\$\$ where \\\eta\\ is the linear predictor (i.e., theta
for a Rasch model with no additional fixed effects) and \\\tau_k\\ are
the item thresholds. Analogous formulas are used for the `cumulative`,
`sratio`, and `cratio` families.

Posterior uncertainty in the thresholds propagates into credible
interval ribbons around the category probability curves — a Bayesian
advantage over point-estimate-based plots from packages like mirt or
eRm.

## References

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
[`conditional_effects`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html),
`itemplot`.

## Examples

``` r
# \donttest{
library(brms)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

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

# Plot all items
plot_ipf(fit_pcm, item_var = item, person_var = id)
#> Error: object 'fit_pcm' not found

# Plot a subset of items
plot_ipf(fit_pcm, item_var = item, person_var = id,
         items = c("I1", "I2", "I3"))
#> Error: object 'fit_pcm' not found

# Customise appearance
plot_ipf(fit_pcm, item_var = item, person_var = id,
         theta_range = c(-6, 6), ncol = 3, prob = 0.90) +
  theme_minimal() +
  labs(title = "Item Category Probability Functions")
#> Error: object 'fit_pcm' not found
# }
```
