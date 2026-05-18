# Posterior-Predict Method for `hurdle_acat`

Not called directly. brms invokes this when
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
is called on a model fitted with
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md).

## Usage

``` r
posterior_predict_hurdle_acat(i, prep, ...)
```

## Arguments

- i:

  Observation index (supplied by brms).

- prep:

  A brms prep object (supplied by brms).

- ...:

  Currently unused.

## Value

An integer vector of length `prep$ndraws` containing replicated category
labels in \\0, 1, \ldots, K - 1\\.

## Author

Kristoffer Magnusson
