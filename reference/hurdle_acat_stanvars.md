# Stan Function Block for `hurdle_acat`

Returns a
[`stanvar`](https://paulbuerkner.com/brms/reference/stanvar.html)
containing the Stan implementation of the hurdle-acat log-PMF used by
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md).
Pass it as the `stanvars` argument to
[`brm`](https://paulbuerkner.com/brms/reference/brm.html).

## Usage

``` r
hurdle_acat_stanvars()
```

## Value

A `stanvars` object adding the `hurdle_acat_merged_lpmf` (and helper
`acat_logit_lpmf`) to the `functions` block of the generated Stan model.

## See also

[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md)

## Author

Kristoffer Magnusson
