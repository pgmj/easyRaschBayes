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
# See https://www.bignardi.co.uk/8_bayes_reliability/tutorial_rmu_sum_score_reliability.html
# for more details on this example

# Simulate data

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
set.seed(1)
N                   = 5000 # number of subjects (mice)
J                   = 3    # number of measurements per subject
true_score_variance = 1
error_variance      = 10

df = expand.grid(j = 1:J, mouse = 1:N)

true_scores       = rnorm(N, mean = 10, sd = sqrt(true_score_variance))
measurement_error = rnorm(N*J, mean = 0, sd = sqrt(error_variance))

df$measurement = true_scores[df$mouse] + measurement_error

df_average_lengths = df %>%
  group_by(mouse) %>%
  summarise(average_measurement = mean(measurement))

# Reliability should equal this:

true_score_variance/(true_score_variance+error_variance/J)
#> [1] 0.2307692

# Approximately the same as:

cor(df_average_lengths$average_measurement, true_scores)^2
#> [1] 0.2414736

# Fit model and calculate RMU

brms_model = brm(
  measurement ~ 1 + (1 | mouse),
  data    = df, 
  iter = 500,
  warmup = 150
)
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')

# Extract posterior draws from brms model

posterior_draws = brms_model %>%
  as_draws_df() %>%
  as_tibble %>% 
  select(starts_with("r_mouse")) %>%
  t()
#> Error: object 'brms_model' not found

# Calculate RMU

RMUreliability(posterior_draws)$hdci
#> Error: object 'posterior_draws' not found
# }
```
