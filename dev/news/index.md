# Changelog

## easyRaschBayes (development version)

- Bug fix -
  [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md),
  [`q3_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic_hpcm.md),
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  and
  [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
  now stop with an error if a person x item combination occurs more than
  once in the model data. Previously only the last response of each
  duplicated combination was kept, without warning. Results for data
  with one response per person and item are unchanged.
- Breaking -
  [`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md)
  follows `easyRasch2::RMtargeting()`. The bottom panel now shows
  response-category bands (`panel = "categories"`), with red marks for
  categories that are never most likely and for disordered thresholds,
  labelled with their posterior probabilities. New arguments `panel`,
  `person_estimate` (EAP default, or WLE), `xlim`, `category_labels`,
  `viridis_*`, `row_gap` and `output = "list"`, which also returns the
  posterior category tables. Defaults for `bins` and the x-axis follow
  RMtargeting, `prob = NULL` hides intervals, and
  `sort_items = "location"` now puts the easiest item at the top.
  ‘scales’ is added to Imports.
- Breaking -
  [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
  now conditions on the total score. It adds Kreiner’s expected gamma
  (`gamma_exp`), replicates data given each person’s total score, uses a
  mid-p `ppp` and adds a max-T `p_adj`. `gamma_diff` is now observed
  minus expected. The zero-width `gamma_obs_q*` columns are replaced by
  `gamma_rep_q*`. Only Rasch-family models (`bernoulli`, `acat`, logit,
  no `disc`) are accepted. Results change.
  [`item_restscore_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_post.md)
  follows the new output and accepts any item column name.
- Breaking -
  [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
  gains `method = c("conditional", "theta")`. The default now computes
  conditional infit given the total score, matching
  [`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html) for
  the same item parameters, and replicates data given each person’s
  total score. Only Rasch-family models are accepted. `method = "theta"`
  gives the 0.3.0 results.
- Breaking -
  [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
  gains `method = c("conditional", "theta")`. The default now takes
  residuals from the expected responses given each person’s total score
  and replicates data given the total score, so observed and replicated
  Q3 centre on about -1/(k-1). Only Rasch-family models are accepted.
  `method = "theta"` gives the 0.3.0 results.
  [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
  adds `q3_p_adj`, a one-sided max-T adjustment over item pairs, and its
  reference line is the mean replicated Q3. Results change.
- Breaking -
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
  `infit_ppp` is now P(observed \> replicated) with mid-p ties, so
  values near 1 indicate underfit, as in
  [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md).
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
  accepts any item column name.
- New datasets `phq9`, the same PHQ-9 data as in easyRasch2 (600 NHANES
  respondents, nine items with four categories), and `phq9_nhanes`, a
  random sample of 1000 NHANES PHQ-9 respondents with the zero responses
  of a general population sample, for the hurdle partial credit model.
- [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  gains `method` (`"quantile"`, `"width"`, `"score"`, `"manual"`) and
  `score_breaks`, as in easyRasch2’s `RMitemICCPlot()`. `n_intervals` is
  renamed `class_intervals` in
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  and
  [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md)
  (the old name still works, with a deprecation warning). A caption
  states the grouping used, also when tied sum scores give fewer
  quantile groups than requested.
- [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md)
  severity panel places each class interval where the item curve equals
  the model-expected mean for the same responses, as
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  does with the total score. Placing intervals at their mean EAP had
  made the observed points look misfitting. Plots change.
- [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  with `dif_var` distinguishes groups by shape and line type as well as
  colour.
- Breaking -
  [`plot_targeting_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting_hpcm.md)
  uses the same design as
  [`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md):
  response-category bands with red marks and posterior probabilities for
  collapsed categories and disordered thresholds, or
  `panel = "thresholds"`. New arguments `panel`, `xlim`,
  `category_labels`, `viridis_*`, `row_gap` and `output = "list"`; the
  default `bins` follows
  [`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).
  Persons who answered 0 on every item are left out of the severity
  histogram, with a caption.
- [`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md)
  labels a category collapse caused by disordered thresholds once, on
  the arrow, and spaces rows more widely, so labels no longer overlap
  when many items have disordered thresholds.
- Docs -
  [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
  and
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
  state that the infit ppp is conservative, also with the conditional
  method, and suggest also using
  [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md).
- Bug fix -
  [`dif_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/dif_statistic.md)
  for polytomous models loaded any difference in group latent means onto
  every item as DIF, and the uniform DIF parameter was not identified.
  Both types now estimate impact separately (new `impact` element) with
  DIF relative to the average item. Uniform DIF uses item random slopes
  of a group indicator. All types share one sign: positive = easier for
  the focal (second) group. Item labels that brms renames no longer give
  silent zeros. Results change for polytomous models.
- Bug fix - rating scale models (`acat` with `(1 | item)` and common
  thresholds) gave every item the same thresholds, and `0 + item` models
  failed in
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  and
  [`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).
  All functions now share one item parameter extractor. Results change
  for rating scale models.
- Bug fix -
  [`person_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md)
  and
  [`person_parameters_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters_hpcm.md)
  gave persons with missing responses the full-test WLE for their sum
  score. The WLE now uses the items each person answered, and a new
  `n_items` column counts them. `score_table` covers persons who
  answered every item.
- Bug fix -
  [`person_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md)
  and
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  counted ordinal sum scores from the lowest category observed in the
  data, not from category 1, which shifted all scores when nobody used
  category 1.
- Bug fix -
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  partial gamma used the wrong column for numeric item labels.
- Bug fix -
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  aligned eigenvector signs across draws by counting signs, which could
  pull a small group of strongly loading items toward zero. Signs now
  follow the dot product with the mean direction.
- Bug fix -
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  for `hurdle_acat` failed when items had different numbers of
  categories.
- Bug fix -
  [`plot_ipf()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_ipf.md)
  and
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  always used the logistic curve. They now use the model’s link for
  cumulative, sratio and cratio models, and stop for non-logit acat
  models and models that estimate `disc`.
- Bug fix -
  [`infit_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic_hpcm.md)
  and
  [`q3_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic_hpcm.md)
  failed with ordered-factor responses or a single draw.
  [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md)
  no longer shifts responses down when there are no zeros. Hurdle item
  parameters are found for labels that brms renames.
- [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md),
  [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md)
  (hurdle panel) and
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  now report how many persons with missing responses they leave out.
- [`item_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters.md)
  `n_eff` is now
  [`posterior::ess_bulk()`](https://mc-stan.org/posterior/reference/ess_bulk.html),
  which accounts for disagreement between chains. ‘posterior’ is added
  to Imports.
- [`plot_tile()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_tile.md)
  and
  [`plot_stackedbars()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_stackedbars.md)
  use
  [`theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).
  The
  [`plot_stackedbars()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_stackedbars.md)
  legend now runs lowest to highest from the top.
  [`plot_ipf()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_ipf.md)
  adds line types to the category colours.
- Docs - corrected
  [`RMUreliability()`](https://pgmj.github.io/easyRaschBayes/dev/reference/RMUreliability.md)
  return value, the
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  `style` default and centring, the
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  partial gamma sign, and added caveats to
  [`posterior_to_prior()`](https://pgmj.github.io/easyRaschBayes/dev/reference/posterior_to_prior.md)
  (independent priors) and
  [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md)
  (EAP binning).
- Bug fix -
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
  `infit_ppp` and
  [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
  `q3_ppp` compared each draw’s observed value with the mean of the
  replicated values, not with the replicated value of the same draw.
  Results change.
- Docs -
  [`infit_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic_hpcm.md)
  is no longer described as conditional infit. The `hdi` tables of
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
  and
  [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
  are described as exceedance proportions. Results are unchanged.

## easyRaschBayes 0.3.0

CRAN release: 2026-06-18

- New custom brms family `hurdle_acat` for hurdle partial credit models
  (H-PCM). This is written by [Kristoffer
  Magnusson](https://rpsychologist.com/)
- Several adapted (experimental) functions for H-PCM, all using the
  `*_hpcm` suffix.
  - [`infit_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic_hpcm.md)
  - [`q3_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic_hpcm.md)
  - [`item_parameters_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters_hpcm.md)
  - [`person_parameters_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters_hpcm.md)
  - [`RMUreliability_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/RMUreliability_hpcm.md)
  - [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md)
  - [`plot_targeting_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting_hpcm.md)
- Bug fix - changed post processing functions to use `median_hdci()`
  instead of `median_hdi()` to avoid issues with multimodal posteriors.
- Bug fix - corrected the WLE person estimates in
  [`person_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md)
  and
  [`person_parameters_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters_hpcm.md).
  The polytomous Warm bias correction had the wrong sign, biasing
  estimates outward; WLE now matches `catR`/`TAM`. Extreme (minimum and
  maximum) scores now receive proper finite WLE estimates instead of
  being clamped to the `theta_range` bounds.

## easyRaschBayes 0.2.0.1

CRAN release: 2026-04-23

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
    [`iarm::ICCplot()`](https://rdrr.io/pkg/iarm/man/ICCplot.html)
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
