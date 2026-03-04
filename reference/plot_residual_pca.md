# Residual PCA Contrast Plot for Bayesian IRT Models

Assesses dimensionality of Bayesian IRT models by performing a principal
component analysis (PCA) of standardized residuals for each posterior
draw. The item loadings on the first residual contrast are plotted
against item locations, with posterior uncertainty displayed as 2D
density contours, crosshairs, or both. A posterior predictive p-value
for the first-contrast eigenvalue tests whether the observed residual
structure exceeds what the model predicts under unidimensionality.

## Usage

``` r
plot_residual_pca(
  model,
  item_var = item,
  person_var = id,
  center = TRUE,
  prob = 0.95,
  ndraws_use = NULL,
  style = c("both", "density", "crosshair"),
  density_alpha = 0.3,
  density_bins = 6,
  density_palette = NULL,
  label_items = TRUE,
  point_size = 2.5,
  point_color = "#0072B2"
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

- center:

  Logical. If `TRUE` (the default), item locations are mean-centered to
  zero, matching the convention used in
  [`plot_targeting`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md).

- prob:

  Numeric in \\(0, 1)\\. Width of the credible intervals for both
  loading and location whiskers. Default is 0.95.

- ndraws_use:

  Optional positive integer. Number of posterior draws to use. If `NULL`
  (the default), up to 500 draws are used.

- style:

  Character. Visual style for displaying uncertainty. `"density"` (the
  default) overlays filled 2D density contours per item computed from
  the draw-level location and loading values, showing the full joint
  posterior uncertainty. `"crosshair"` shows point estimates with
  horizontal and vertical credible interval bars. `"both"` displays
  density contours with crosshairs on top.

- density_alpha:

  Numeric in \\\[0, 1\]\\. Maximum opacity of the density contours when
  `style` is `"density"` or `"both"`. Default is 0.3.

- density_bins:

  Integer. Number of contour bins for
  [`geom_density_2d_filled`](https://ggplot2.tidyverse.org/reference/geom_density_2d.html).
  Default is 6.

- density_palette:

  An optional character vector of colors for the density contour fills
  (from low to high density). If `NULL` (the default), a blue sequential
  ramp is used. The length should match `density_bins`; the lowest
  (background) level is always transparent.

- label_items:

  Logical. If `TRUE` (the default), item names are displayed next to
  points.

- point_size:

  Numeric. Size of the item points. Default is 2.5.

- point_color:

  Color for the item points and error bars. Default is `"#0072B2"`.

## Value

A list with three elements:

- plot:

  A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object showing item loadings on the first residual contrast (y-axis)
  versus item locations (x-axis).

- data:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  columns: `item`, `location` (posterior mean item location),
  `location_lower`, `location_upper` (location CI), `loading` (posterior
  mean loading on first contrast), `loading_lower`, `loading_upper`
  (loading CI).

- eigenvalue:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  columns: `eigenvalue_obs` (posterior mean observed eigenvalue),
  `eigenvalue_rep` (posterior mean replicated eigenvalue),
  `eigenvalue_diff` (posterior mean difference), `ppp` (posterior
  predictive p-value), `var_explained_obs`, `var_explained_rep`
  (posterior mean proportions of residual variance explained).

## Details

The procedure for each posterior draw \\s\\ is:

1.  Obtain category probabilities from
    [`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html).
    Compute expected values \\E^{(s)}\_{vi}\\ and variances
    \\Var^{(s)}\_{vi}\\.

2.  Compute standardized residuals for observed data:
    \$\$z^{obs(s)}\_{vi} = \frac{X\_{vi} - E^{(s)}\_{vi}}
    {\sqrt{Var^{(s)}\_{vi}}}\$\$

3.  Generate replicated data \\Y^{rep(s)}\\ from
    [`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
    and compute standardized residuals: \$\$z^{rep(s)}\_{vi} =
    \frac{Y^{rep(s)}\_{vi} - E^{(s)}\_{vi}}{\sqrt{Var^{(s)}\_{vi}}}\$\$

4.  Reshape both sets of residuals into person \\\times\\ item matrices
    and perform SVD on each.

5.  Extract the first-contrast eigenvalue and item loadings from both
    observed and replicated SVDs.

6.  Compare eigenvalues across draws: the posterior predictive p-value
    `ppp = mean(eigenvalue_obs > eigenvalue_rep)` tests whether the
    observed residual structure is stronger than what the model produces
    under its own assumptions.

When `style = "density"` or `style = "both"`, the draw-level (location,
loading) pairs for each item are used to construct filled 2D kernel
density contours via
[`geom_density_2d_filled`](https://ggplot2.tidyverse.org/reference/geom_density_2d.html).
The lowest contour level (outside all contours) is set to transparent so
the white panel background shows through. Higher density regions use
progressively darker fills.

Item loadings are aligned across draws using majority-sign alignment to
resolve the sign indeterminacy of eigenvectors.

## References

Smith, E. V. (2002). Detecting and evaluating the impact of
multidimensionality using item fit statistics and principal component
analysis of residuals. *Journal of Applied Measurement*, *3*, 205–231.

Bürkner, P.-C. (2021). Bayesian Item Response Modeling in R with brms
and Stan. *Journal of Statistical Software*, *100*, 1–54.
[doi:10.18637/jss.v100.i05](https://doi.org/10.18637/jss.v100.i05)

## See also

[`plot_targeting`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md)
for person-item maps,
[`plot_ipf`](https://pgmj.github.io/easyRaschBayes/reference/plot_ipf.md)
for item category probability curves,
[`q3_statistic`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)
for Q3 residual correlations (another local dependence / dimensionality
diagnostic).

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

# 2D density contours (default)
result <- plot_residual_pca(fit_pcm)
#> Error: object 'fit_pcm' not found
result$plot
#> Error: object 'result' not found

# Crosshair style
result_c <- plot_residual_pca(fit_pcm, style = "crosshair")
#> Error: object 'fit_pcm' not found
result_c$plot
#> Error: object 'result_c' not found

# Both combined
result_b <- plot_residual_pca(fit_pcm, style = "both")
#> Error: object 'fit_pcm' not found
result_b$plot
#> Error: object 'result_b' not found

# Custom warm palette
result_w <- plot_residual_pca(
  fit_pcm,
  density_palette = c("#FEE8C8", "#FDBB84", "#E34A33",
                      "#B30000", "#7F0000", "#4A0000"),
  point_color = "#B30000"
)
#> Error: object 'fit_pcm' not found
result_w$plot
#> Error: object 'result_w' not found
# }
```
