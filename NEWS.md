# easyRaschBayes 0.4.0

## Breaking changes

- `item_restscore_statistic()`, `infit_statistic()` and `q3_statistic()` now
  condition on each person's total score, both for the expected values and
  for the replicated data. They accept only Rasch-family models (`bernoulli`
  or `acat`, logit link, no `disc`). `method = "theta"` in
  `infit_statistic()` and `q3_statistic()` gives the 0.3.0 results.
  **Results change.**
  - `item_restscore_statistic()` adds Kreiner's expected gamma
    (`gamma_exp`), a mid-p `ppp` and a max-T adjusted `p_adj`. `gamma_diff`
    is now observed minus expected, and `gamma_obs_q*` is replaced by
    `gamma_rep_q*`.
  - Conditional infit matches `iarm::out_infit()` for the same item
    parameters.
  - `q3_post()` adds `q3_p_adj`, adjusted over item pairs, and its reference
    line is the mean replicated Q3.
- `infit_post()` `infit_ppp` and `q3_post()` `q3_ppp` compared each draw's
  observed value with the mean of the replicated values, not with the same
  draw's replicated value. `infit_ppp` is now also P(observed > replicated)
  with mid-p ties, so values near 1 indicate underfit, as for `q3_ppp`.
  **Results change.**
- `plot_targeting()` and `plot_targeting_hpcm()` follow
  `easyRasch2::RMtargeting()`. The bottom panel shows response-category
  bands, with categories that are never most likely and disordered
  thresholds marked in red and labelled with their posterior probabilities.
  `panel = "thresholds"` gives a dot-and-whisker plot. New arguments
  `panel`, `xlim`, `category_labels`, `viridis_*`, `row_gap` and
  `output = "list"`, plus `person_estimate` in `plot_targeting()`. The
  defaults for `bins` and the x-axis change, and `sort_items = "location"`
  puts the easiest item at the top. 'scales' is added to Imports.

## New features

- New datasets `phq9` (600 NHANES respondents, nine items with four
  categories, as in easyRasch2) and `phq9_nhanes` (1000 NHANES respondents,
  including zero responses, for the hurdle partial credit model).
- `plot_icc()` gains `method` and `score_breaks` for forming class
  intervals, as in `easyRasch2::RMitemICCPlot()`, and a caption stating the
  grouping. `n_intervals` is renamed `class_intervals` in `plot_icc()` and
  `plot_icc_hpcm()`, and the old name is deprecated.
- `plot_icc()` with `dif_var` separates groups by shape and line type as
  well as colour. `plot_ipf()` adds line types, and `plot_tile()` and
  `plot_stackedbars()` use `theme_bw()`.
- `plot_icc()`, `plot_icc_hpcm()` and `plot_residual_pca()` report how many
  persons with missing responses they leave out.
- `item_restscore_post()` and `infit_post()` accept any item column name.
- `item_parameters()` `n_eff` is now `posterior::ess_bulk()`. 'posterior' is
  added to Imports.

## Bug fixes

- `dif_statistic()` for polytomous models loaded differences in group means
  onto the items as DIF. Impact is now estimated separately (new `impact`
  element), and all DIF types share one sign: positive means easier for the
  focal group. **Results change for polytomous models.**
- Rating scale models gave every item the same thresholds.
  **Results change for rating scale models.** `0 + item` models no longer
  fail in `plot_residual_pca()` and `plot_targeting()`.
- `person_parameters()` and `person_parameters_hpcm()` now compute the WLE
  of persons with missing responses from the items they answered, with a
  new `n_items` column.
- `person_parameters()` and `plot_icc()` counted sum scores from the lowest
  observed category instead of category 1.
- `q3_statistic()`, `q3_statistic_hpcm()`, `plot_residual_pca()` and
  `item_restscore_statistic()` stop if a person and item combination occurs
  more than once. Previously only the last response was kept.
- `plot_icc_hpcm()` placed class intervals in the severity panel at their
  mean EAP, which made observed points look misfitting, and shifted
  responses down when there were no zeros. **Plots change.**
- `plot_ipf()` and `plot_icc()` use the model's link for cumulative, sratio
  and cratio models, and stop for non-logit `acat` models and models that
  estimate `disc`.
- `plot_residual_pca()` aligns eigenvector signs across draws by the dot
  product with the mean direction.
- `plot_icc()` partial gamma used the wrong column for numeric item labels.
- `posterior_epred_hurdle_acat()` failed when items had different numbers
  of categories.
- `infit_statistic_hpcm()` and `q3_statistic_hpcm()` failed with
  ordered-factor responses or a single draw, and with item labels that brms
  renames.

## Documentation

- `infit_statistic()` and `infit_post()` state that the infit ppp is
  conservative, also with the conditional method, and suggest
  `item_restscore_statistic()` as well.
- Corrected `RMUreliability()` return value, `plot_residual_pca()` `style`
  default and centring, and `plot_icc()` partial gamma sign. Added caveats
  to `posterior_to_prior()` and `plot_icc_hpcm()`. `infit_statistic_hpcm()`
  is no longer described as conditional infit.

# easyRaschBayes 0.3.0

- New custom brms family `hurdle_acat` for hurdle partial credit models (H-PCM).
  This is written by [Kristoffer Magnusson](https://rpsychologist.com/)
- Several adapted (experimental) functions for H-PCM, all using the `*_hpcm` suffix.
  - `infit_statistic_hpcm()`
  - `q3_statistic_hpcm()`
  - `item_parameters_hpcm()`
  - `person_parameters_hpcm()`
  - `RMUreliability_hpcm()`
  - `plot_icc_hpcm()`
  - `plot_targeting_hpcm()`
- Bug fix - changed post processing functions to use `median_hdci()` 
  instead of `median_hdi()` to avoid issues with multimodal posteriors.
- Bug fix - corrected the WLE person estimates in `person_parameters()` and 
  `person_parameters_hpcm()`. The polytomous Warm bias correction had the wrong 
  sign, biasing estimates outward; WLE now matches `catR`/`TAM`. Extreme (minimum 
  and maximum) scores now receive proper finite WLE estimates instead of being 
  clamped to the `theta_range` bounds.

# easyRaschBayes 0.2.0.1

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
