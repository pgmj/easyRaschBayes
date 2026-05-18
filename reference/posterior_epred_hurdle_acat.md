# Posterior-Epred Method for `hurdle_acat`

Not called directly. brms invokes this when
[`posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
is called on a model fitted with
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md).
Returns an \\S \times N \times K\\ array of category probabilities,
where slot 1 is \\hu = P(Y = 0)\\ and slots \\2..K\\ are \\(1 - hu)
\cdot P\_{sev}(k)\\ for \\k = 1, \ldots, K - 1\\.

## Usage

``` r
posterior_epred_hurdle_acat(prep)
```

## Arguments

- prep:

  A brms prep object (supplied by brms).

## Value

A numeric \\S \times N \times K\_{total}\\ array of category
probabilities.

## Author

Kristoffer Magnusson
