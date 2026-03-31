# Changelog

## easyRaschBayes 0.2.0.1 (dev)

- Fix for
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  theta values.

## easyRaschBayes 0.2.0

CRAN release: 2026-03-28

- Updates:
  - [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
    is now faster and defaults to not output outfit statistics, since
    these are of dubious value (see [Müller
    2020](https://doi.org/10.1186/s40488-020-00108-7) and [Johansson,
    2025](https://pgmj.github.io/rasch_itemfit/)).
  - [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
    is now faster and has output similar to
    [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md).
  - [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
    is now faster and has output similar to
    [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md).
  - [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
    is now faster.
- Three post-processing helper functions that output a list with result
  table(s) and a plot:
  - [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
    for the output of
    [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
  - [`item_restscore_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_post.md)
    for the output of
    [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
  - [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
    for the output of
    [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
- Other new functions:
  - [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
    can both produce a basic conditional Item Characteristic Curves plot
    and optionally a DIF ICC plot, the latter also reporting a partial
    gamma DIF magnitude coefficient. This is inspired by
    `iarm::ICCplot()`
  - [`item_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters.md)
    to retrieve item threshold locations
  - [`person_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md)
    to retrieve person locations (latent scores) using EAP and WLE, and
    also a transformation table from ordinal sum score to EAP/WLE.
  - [`posterior_to_prior()`](https://pgmj.github.io/easyRaschBayes/dev/reference/posterior_to_prior.md)
    that extracts priors from a `brmsfit` model that may be useful in
    model fit assessment. An article will be added later, exploring this
    further.
  - [`plot_tile()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_tile.md),
    [`plot_bars()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_bars.md),
    and
    [`plot_stackedbars()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_stackedbars.md)
    uses wide format item data to create plots to review response
    patterns.

## easyRaschBayes 0.1.0

CRAN release: 2026-03-11

- Initial release.
