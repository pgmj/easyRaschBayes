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

  Logical. Print detailed information about the input data. Default is
  TRUE.

- level:

  Numeric. Credibility level for the highest density continuous
  interval. Default is 0.95.

## Value

A list containing:

- hdci: A data frame with a point-estimate (posterior mean) and highest
  density continuous interval for reliability, calculated using the
  ggdist::mean_hdci function

- reliability_posterior_draws: A numeric vector of posterior draws for
  reliability, of length K/2 (K = number of columns/draws in your
  input_draws matrix)

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

## References

Bignardi, G., Kievit, R., & Bürkner, P. C. (2025). A general method for
estimating reliability using Bayesian Measurement Uncertainty. PsyArXiv.
[doi:10.31234/osf.io/h54k8](https://osf.io/preprints/psyarxiv/h54k8)

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
