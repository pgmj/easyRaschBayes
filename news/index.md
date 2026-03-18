# Changelog

## easyRaschBayes 0.1.1

- [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
  is now faster and defaults to not output outfit statistics.
- [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_statistic.md)
  is now faster and has output similar to
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/reference/infit_post.md).
- [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)
  is now faster and has output similar to
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/reference/infit_post.md).
- Three post-processing helper functions added that output a list with
  result table(s) and a plot:
  - [`infit_post()`](https://pgmj.github.io/easyRaschBayes/reference/infit_post.md)
    for the output of
    [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/infit_statistic.md)
  - [`item_restscore_post()`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_post.md)
    for the output of
    [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/item_restscore_statistic.md)
  - [`q3_post()`](https://pgmj.github.io/easyRaschBayes/reference/q3_post.md)
    for the output of
    [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/reference/q3_statistic.md)
- [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/reference/plot_residual_pca.md)
  is now faster.
- New function
  [`posterior_to_prior()`](https://pgmj.github.io/easyRaschBayes/reference/posterior_to_prior.md)
  that extracts priors from a `brmsfit` model that may be useful in
  model fit assessment. An article will be added, exploring this
  further.

## easyRaschBayes 0.1.0

CRAN release: 2026-03-11

- Initial release.
