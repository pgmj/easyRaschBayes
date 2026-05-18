# Log-Likelihood Method for `hurdle_acat`

Not called directly. brms invokes this when
[`log_lik`](https://mc-stan.org/rstantools/reference/log_lik.html) is
called on a model fitted with
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md).

## Usage

``` r
log_lik_hurdle_acat(i, prep)
```

## Arguments

- i:

  Observation index (supplied by brms).

- prep:

  A brms prep object (supplied by brms).

## Value

A numeric vector of length `prep$ndraws`.

## Author

Kristoffer Magnusson
