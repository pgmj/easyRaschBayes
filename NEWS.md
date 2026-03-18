# easyRaschBayes 0.1.1

- `infit_statistic()` is now faster and defaults to not output outfit statistics.
- `item_restscore_statistic()` is now faster and has output similar to `infit_post()`.
- `q3_statistic()` is now faster and has output similar to `infit_post()`.
- Three post-processing helper functions added that output a list with result table(s) and a plot:
  - `infit_post()` for the output of `infit_statistic()`
  - `item_restscore_post()` for the output of `item_restscore_statistic()`
  - `q3_post()` for the output of `q3_statistic()`
- `plot_residual_pca()` is now faster.
- New function `posterior_to_prior()` that extracts priors from a `brmsfit` model 
that may be useful in model fit assessment. An article will be added, exploring this further.


# easyRaschBayes 0.1.0

- Initial release.
