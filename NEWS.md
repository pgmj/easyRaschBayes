# easyRaschBayes 0.2.0.1 (dev)

- Fix for `plot_icc()` theta values.

# easyRaschBayes 0.2.0

- Updates:
  - `infit_statistic()` is now faster and defaults to not output outfit statistics, 
since these are of dubious value (see [Müller 2020](https://doi.org/10.1186/s40488-020-00108-7) 
and [Johansson, 2025](https://pgmj.github.io/rasch_itemfit/)).
  - `item_restscore_statistic()` is now faster and has output similar to `infit_post()`.
  - `q3_statistic()` is now faster and has output similar to `infit_post()`.
  - `plot_residual_pca()` is now faster.

- Three post-processing helper functions that output a list with result table(s) and a plot:
  - `infit_post()` for the output of `infit_statistic()`
  - `item_restscore_post()` for the output of `item_restscore_statistic()`
  - `q3_post()` for the output of `q3_statistic()`

- Other new functions:
  - `plot_icc()` can both produce a basic conditional Item Characteristic Curves plot 
  and optionally a DIF ICC plot, the latter also reporting a partial gamma DIF 
  magnitude coefficient. This is inspired by `iarm::ICCplot()`
  - `item_parameters()` to retrieve item threshold locations
  - `person_parameters()` to retrieve person locations (latent scores) using EAP and WLE, 
  and also a transformation table from ordinal sum score to EAP/WLE.
  - `posterior_to_prior()` that extracts priors from a `brmsfit` model 
  that may be useful in model fit assessment. An article will be added later, 
  exploring this further.
  - `plot_tile()`, `plot_bars()`, and `plot_stackedbars()` uses wide format item data to create plots to review
  response patterns.

# easyRaschBayes 0.1.0

- Initial release.
