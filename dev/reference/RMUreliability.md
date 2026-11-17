# Estimate reliability (Relative Measurement Uncertainty) from Bayesian measurement models

This function measures reliability using posterior draws from a fitted
Bayesian model.

## Usage

``` r
RMUreliability(input_draws, verbose = FALSE, level = 0.95)
```

## Arguments

- input_draws:

  A matrix or data frame of posterior draws. Rows represent subjects and
  columns represent draws.

- verbose:

  Logical. Print the number of subjects and draws. Default is FALSE.

- level:

  Numeric. Credibility level for the highest density continuous
  interval. Default is 0.95.

## Value

A data frame from
[`ggdist::mean_hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html)
with one row and the columns `rmu_estimate` (mean of the reliability
draws), `hdci_lowerbound`, `hdci_upperbound`, `.width`, `.point` and
`.interval`.

## Details

To use this function, you will need to provide a matrix (input_draws)
that contains the posterior draws for the parameter you wish to
calculate reliability. The function assumes that rows of input_draws
represent subjects and columns represent posterior draws.

For an example of how to apply this function to calculate mean score
reliability using brms, see [this
tutorial](https://www.bignardi.co.uk/8_bayes_reliability/tutorial_rmu_sum_score_reliability.html).

For an example of how to apply this function to go/go-no task data using
brms, see [this
tutorial](https://www.bignardi.co.uk/8_bayes_reliability/tutorial_calculating_rmu_gonogo.html).

The draws (columns) are shuffled and split into two halves. Each column
of the first half is correlated with the matching column of the second
half across subjects, which gives `floor(K / 2)` reliability draws for
`K` posterior draws. Their mean and highest density continuous interval
are returned. The split is random, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible results. Missing values in `input_draws` give a warning and
an `NA` result.

## References

Bignardi, G., Kievit, R., & Bürkner, P. C. (2026). A general method for
estimating reliability using Bayesian Measurement Uncertainty. PsyArXiv.
<https://osf.io/preprints/psyarxiv/h54k8>

## Examples

``` r
# \donttest{
library(brms)
#> Loading required package: Rcpp
#> Loading 'brms' package (version 2.23.0). Useful instructions
#> can be found by typing help('brms'). A more detailed introduction
#> to the package is available through vignette('brms_overview').
#> 
#> Attaching package: ‘brms’
#> The following object is masked from ‘package:stats’:
#> 
#>     ar
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
library(tidyr)
library(tibble)

# --- Dichotomous Rasch Model ---

# PHQ-9 items dichotomised (any symptom vs none), for illustration only
df_rm <- phq9[1:200, 1:9] %>%
  mutate(across(everything(), ~ as.integer(.x > 0))) %>%
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

# Extract posterior draws from brms model

posterior_draws = fit_rm %>%
  as_draws_df() %>%
  as_tibble() %>% 
  select(starts_with("r_id")) %>%
  t()
#> Error: object 'fit_rm' not found

# Calculate RMU

RMUreliability(posterior_draws)
#> Error: object 'posterior_draws' not found
# }
```
