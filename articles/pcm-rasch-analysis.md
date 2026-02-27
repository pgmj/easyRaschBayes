# Partial Credit Model Analysis with easyRaschBayes

## Overview

This vignette demonstrates a Partial Credit Model (PCM) Rasch analysis
workflow using `easyRaschBayes`. The analysis uses the
[`eRm::pcmdat2`](https://rdrr.io/pkg/eRm/man/eRm.data.html) dataset, a
small polytomous item response data set included in the `eRm` package.

All functions in this package work with a Bayesian `brms` model object
fitted with the `acat` (adjacent categories) family, which parameterises
the PCM. Dichotomous Rasch models can also be fit using `brms` and
analyzed with the functions in this package. A code example is available
[here](https://pgmj.github.io/reliability.html#rasch-dichotomous-model).

There is brief text explaining the output and interpretation in this
vignette. For a more extensive treatment of the Rasch analysis steps,
please see the [`easyRasch`
vignette](https://pgmj.github.io/raschrvignette/RaschRvign.html).

## Data Preparation

``` r
library(easyRaschBayes)
library(brms)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
```

[`eRm::pcmdat2`](https://rdrr.io/pkg/eRm/man/eRm.data.html) is in wide
format with item responses coded 0, 1, 2, …. The `brms` `acat` family
expects response categories starting at **1**, so we add 1 to every
response before reshaping to long format.

``` r
df_pcm <- eRm::pcmdat2 %>%
  mutate(across(everything(), ~ .x + 1)) %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response")

head(df_pcm)
#> # A tibble: 6 × 3
#>   id    item  response
#>   <chr> <chr>    <dbl>
#> 1 1     I1           2
#> 2 1     I2           2
#> 3 1     I3           2
#> 4 1     I4           2
#> 5 2     I1           1
#> 6 2     I2           1
```

## Fitting the PCM

The model is fitted once and saved to disk. The code chunk below shows
the fitting call (not evaluated during `R CMD check`). A pre-fitted
model stored at `fits/fit_pcm.rds` is loaded instead.

``` r
prior_pcm <- prior("normal(0, 5)", class = "Intercept") +
  prior("normal(0, 3)", class = "sd", group = "id")
```

``` r
fit_pcm <- brm(
  response | thres(gr = item) ~ 1 + (1 | id),
  data    = df_pcm,
  family  = acat,
  prior   = prior_pcm,
  chains  = 4,
  cores   = 4,
  iter    = 2000
)
saveRDS(fit_pcm, "fits/fit_pcm.rds")
```

``` r
fit_pcm <- readRDS("fits/fit_pcm.rds")
```

## Item Fit: Infit and Outfit Statistics

[`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
computes posterior predictive infit and outfit statistics for each item.
Values near 1.0 indicate good fit; values substantially above 1 suggest
underfit (unexpected responses), values below 1 suggest overfit (too
predictable).

The `ndraws_use` argument limits the number of posterior draws used,
which speeds up computation during exploration. For final reporting, use
all draws (set `ndraws_use = NULL` or omit it).

``` r
fit_stats <- infit_statistic(fit_pcm, ndraws_use = 500)

# Posterior predictive p-values summarised per item
fit_stats %>%
  group_by(item) %>%
  summarise(
    infit_obs = round(mean(infit),3),
    infit_rep = round(mean(infit_rep),3),
    infit_ppp = round(mean(infit_rep > infit),3)
  )
```

`infit_obs` indicates the observed conditional infit, which can be
compared to `infit_rep`, which is akin to the model expected value.
Posterior predictive p-values (`*_ppp`) close to 0.5 indicate that the
observed statistic falls near the centre of the posterior predictive
distribution, implying good fit. Values near 0 or 1 warrant further
investigation.

## Item–Rest Score Association

[`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_statistic.md)
computes Goodman-Kruskal’s gamma between each item’s observed responses
and the rest score (total score minus the focal item). In a well-fitting
Rasch model, gamma should be positive and moderate; very high values may
indicate redundancy, very low or negative values suggest the item does
not relate well to the latent trait.

``` r
rest_stats <- item_restscore_statistic(fit_pcm, ndraws_use = 500)

rest_stats[,1:5] %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))
```

The output is limited to columns 1 through 5 in the output above.

## Dimensionality: Residual PCA

[`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/reference/plot_residual_pca.md)
performs a principal components analysis on the person-item residuals
and plots the loadings on the first contrast together with the item
locations and the uncertainty of both. Substantial loadings on the first
contrast suggest multidimensionality.

``` r
pca <- plot_residual_pca(fit_pcm, ndraws_use = 500)
pca$plot
```

Items with positive loadings cluster on one end, negative loadings on
the other. If the observed largest eigenvalue is smaller than the
replicated, unidimensionality is supported. The ppp should not be close
to 1.

## Local Dependence: Q3 Residual Correlations

[`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)
computes Yen’s Q3 statistic — the correlation between person-item
residuals for every item pair. After conditioning on the latent trait,
residuals should be uncorrelated; elevated Q3 values indicate local
dependence (LD). Our primary metric here is the ppp, that should not be
close to 1. The output is filtered on ppp values above 0.95.

``` r
q3_stats <- q3_statistic(fit_pcm, ndraws_use = 500)

q3_stats %>% 
  filter(ppp > 0.95)
```

Pairs flagged as locally dependent should be examined for substantive
overlap in item content.

## Item Category Probabilities

This plot shows the probability of using a response category on the y
axis and the latent score on the x axis. The crossover points, where
lines meet, are the item category threshold locations. Uncertainty is
shown with the shaded area around each line.

``` r
plot_ipf(fit_pcm, theta_range = c(-6,5))
```

## Person–Item Targeting

[`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/reference/plot_targeting.md)
produces a Wright map (person–item targeting plot) showing the
distribution of person locations alongside the item threshold locations
on the same logit scale. Good targeting occurs when person and item
distributions overlap substantially.

``` r
plot_targeting(fit_pcm)
```

If the person distribution is systematically above or below the item
thresholds, the test may be too easy or too hard for the sample.

## Reliability: Relative Measurement Uncertainty

[`RMUreliability()`](https://pgmj.github.io/easyRaschBayes/reference/RMUreliability.md)
provides a Bayesian reliability estimate via Relative Measurement
Uncertainty (RMU, see Bignardi et al., 2025). It requires a matrix of
person location draws with dimensions $$persons \times draws$$. The
output is a point estimate and lower/upper 95% highest density
continuous intervals (HDCI).

``` r
person_draws <- fit_pcm %>%
  as_draws_df() %>%
  as_tibble() %>% 
  select(starts_with("r_id")) %>%
  t()
rmu <- RMUreliability(person_draws)
rmu
```

RMU values range from 0 to 1, with higher values indicating higher
reliability, similarly to traditional reliability metrics.
