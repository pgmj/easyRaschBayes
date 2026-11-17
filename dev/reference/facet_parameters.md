# Facet Parameters from a Many-Facet Rasch Model

Extracts the measures of a third facet, for example raters, from a
many-facet Rasch model fitted with brms, with posterior summaries, the
number of ratings, the observed average rating and a fair average rating
for each element.

## Usage

``` r
facet_parameters(
  model,
  facet_var,
  item_var = item,
  person_var = id,
  draws = FALSE,
  prob = 0.95
)
```

## Arguments

- model:

  A fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  object from a many-facet model with the facet as a random intercept, a
  fixed effect, or a threshold group. Examples:

  Rating scale, random raters

  :   `rating ~ 1 + item + (1 | rater) + (1 | id)`, `family = acat`.

  Partial credit by item

  :   `rating | thres(gr = item) ~ 1 + (1 | rater) + (1 | id)`.

  Rater-specific rating scales

  :   `rating | thres(gr = rater) ~ 1 + item + (1 | id)`.

  Dichotomous

  :   `y ~ 1 + item + (1 | rater) + (1 | id)`, `family = bernoulli()`.

- facet_var:

  An unquoted variable name identifying the facet, for example `rater`.

- item_var:

  An unquoted variable name identifying the items. Default is `item`.

- person_var:

  An unquoted variable name identifying the persons. Default is `id`.

- draws:

  Logical. If `TRUE`, the centred posterior draws of the facet locations
  are included as a matrix (elements x draws). Default is `FALSE`.

- prob:

  Numeric in \\(0, 1)\\. Width of the highest density continuous
  interval (HDCI). Default is 0.95.

## Value

A list with the following elements:

- summary:

  A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per facet element and the columns: the facet variable, `n`
  (number of ratings), `observed_avg` (mean observed rating),
  `fair_avg`, `fair_avg_lower` and `fair_avg_upper` (posterior mean and
  HDCI of the fair average, see Details), `location` (posterior mean
  location, positive = lower ratings), `se` (posterior SD),
  `hdci_lower`, `hdci_upper` and `n_eff` (bulk effective sample size).

- thresholds:

  (Only with `thres(gr = facet)`) A
  [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  one row per element and threshold: the facet variable, `threshold`,
  `location`, `se`, `hdci_lower` and `hdci_upper`, on the same scale as
  `summary`. `NULL` otherwise.

- facet_sd:

  (Only for a random facet) A
  [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
  the posterior mean, SD and HDCI of the facet standard deviation.
  `NULL` otherwise.

- draws_matrix:

  (Only if `draws = TRUE`) A numeric matrix of centred location draws
  with one row per element and one column per posterior draw.

## Details

**Locations.** The model is read as \\\theta_v - \delta_i - \lambda_j\\
on the logit scale, where \\\delta_i\\ is the item location and
\\\lambda_j\\ the location of facet element \\j\\. A positive location
means lower ratings, so for raters it is severity, on the same side as
item difficulty in
[`item_parameters`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters.md).
The locations are centred within each posterior draw so that they
average 0. The result is the same whether the facet is a random effect
or a fixed effect with any contrasts, including
[`contr_sum_zero`](https://pgmj.github.io/easyRaschBayes/dev/reference/contr_sum_zero.md).
With facet-specific thresholds (`thres(gr = rater)`) the location is the
mean of the element's thresholds, and the `thresholds` table shows each
threshold. brms sets the number of thresholds per element from the
categories that element used, so set it in the formula, as in
`thres(3, gr = rater)` for four categories. The function stops if the
elements have different numbers of thresholds. Widely spaced thresholds
for a rater mean the middle categories cover a wide range of locations
and are used often (central tendency). Closely spaced thresholds mean
the extreme categories are used often.

**Fair average.** The observed average rating of a rater depends on
which persons and items the rater happened to rate. The fair average is
the expected rating the element would give a person at the mean location
of the persons in the sample, averaged over the items, computed in each
posterior draw. It corresponds to the "Fair-M average" of the Facets
program. Person x facet (halo) effects are set to 0. Ratings are on the
scale brms uses: 0 and 1 for
[`bernoulli()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
categories 1 to \\C\\ for ordinal families.

**Reliability.** With `draws = TRUE`, the draws matrix can be passed to
[`RMUreliability`](https://pgmj.github.io/easyRaschBayes/dev/reference/RMUreliability.md).
For raters, a high value means that the raters differ reliably in
severity, which is a problem for the ratings, not a sign of their
quality. With few elements, for example 3 to 5 raters, the correlation
across elements that RMU is based on is very uncertain.

## References

Linacre, J. M. (1989). *Many-facet Rasch measurement*. MESA Press.

Eckes, T. (2015). *Introduction to many-facet Rasch measurement:
Analyzing and evaluating rater-mediated assessments* (2nd ed.). Peter
Lang.

## See also

[`item_parameters`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters.md),
[`contr_sum_zero`](https://pgmj.github.io/easyRaschBayes/dev/reference/contr_sum_zero.md),
[`infit_statistic`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
(with `facet_var`),
[`RMUreliability`](https://pgmj.github.io/easyRaschBayes/dev/reference/RMUreliability.md).

## Examples

``` r
# \donttest{
library(brms)

# Simulated ratings: 60 persons, 3 items, 4 raters, two raters per person
set.seed(1)
d <- do.call(rbind, lapply(1:60, function(p) {
  expand.grid(id = p, item = c("i1", "i2", "i3"),
              rater = paste0("R", c((p - 1) %% 4 + 1, p %% 4 + 1)),
              stringsAsFactors = FALSE)
}))
severity <- c(R1 = -0.5, R2 = 0, R3 = 0.2, R4 = 0.6)
theta <- rnorm(60)
d$y <- rbinom(nrow(d), 1, plogis(theta[d$id] - severity[d$rater]))

fit <- brm(y ~ 1 + item + (1 | rater) + (1 | id), data = d,
           family = bernoulli(), chains = 2, iter = 1000)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

fp <- facet_parameters(fit, facet_var = rater, draws = TRUE)
#> Error: object 'fit' not found
fp$summary
#> Error: object 'fp' not found
RMUreliability(fp$draws_matrix)
#> Error: object 'fp' not found
# }
```
