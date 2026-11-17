# Changelog

## easyRaschBayes 0.4.0.9001

### New features

- Many-facet Rasch models:
  [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md),
  [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md),
  [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
  and
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  get a `facet_var` argument (for example `facet_var = rater`). Each
  item x rater combination is treated as an item.
  `infit_statistic(by = "item")` or `by = "facet"` pools over raters or
  items. Q3 output gets a `pair_type` column and leaves out cell pairs
  that no person shares.
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
  and
  [`item_restscore_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_post.md)
  keep the item and facet columns, and
  [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
  keeps `pair_type`. A warning flags rating designs that do not connect
  all raters and items. A message notes when the model has a person x
  rater (halo) term, since the conditional methods then assess fit after
  the halo is accounted for.
- [`facet_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/facet_parameters.md):
  locations (severity) of the facet elements, for example raters, with
  the number of ratings, observed average and fair average rating. With
  `draws = TRUE` the draws can go to
  [`RMUreliability()`](https://pgmj.github.io/easyRaschBayes/dev/reference/RMUreliability.md).
- [`contr_sum_zero()`](https://pgmj.github.io/easyRaschBayes/dev/reference/contr_sum_zero.md):
  orthonormal sum-to-zero contrasts, which give every level of a fixed
  factor the same prior variance (unlike
  [`contr.sum()`](https://rdrr.io/r/stats/contrast.html)).
- New vignette “Many-Facet Rasch Models with easyRaschBayes”.
- The conditional methods of
  [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md),
  [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
  and
  [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
  are faster and use far less memory when persons answer only some items
  (planned missing data, many-facet designs). Results for complete data
  are unchanged. With missing responses the replicated data differ only
  by simulation noise.
- [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
  and the theta method of
  [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
  are faster with many item pairs.

### Bug fixes

- Item parameters from models with fixed item effects and an intercept
  now respect the contrasts of the item factor. Before, treatment coding
  was assumed. With other contrasts (for example
  [`contr.sum()`](https://rdrr.io/r/stats/contrast.html)) the functions
  stopped with an error, or, for item labels with spaces or punctuation,
  returned wrong item locations. Treatment-coded models are unchanged.
- Fixed item effects are matched by their exact parameter names. Before,
  a covariate whose name starts with the item variable name (for example
  `item_order` next to `item`) was read as an item effect: the functions
  stopped with an error, or, for item labels with spaces or punctuation,
  returned wrong item locations. For contrasts set by
  `options(contrasts = ...)` instead of on the factor, the same option
  must be set when the functions are run. If the names fit both that
  option and treatment coding, the functions stop with an error. Other
  results are unchanged.
- [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
  failed when a person with a total score of 0 had a pattern of missing
  responses that no other person had. Results for other data are
  unchanged.

### Documentation

- Help pages and vignettes now say that a p-value “flags misfit less
  often than its nominal level” instead of calling it “conservative”.
- Percent signs in the
  [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
  and
  [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
  help pages were shown with a backslash.

## easyRaschBayes 0.4.0

CRAN release: 2026-11-10

### Breaking changes

- [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md),
  [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
  and
  [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
  now condition on each person’s total score, both for the expected
  values and for the replicated data. They accept only Rasch-family
  models (`bernoulli` or `acat`, logit link, no `disc`).
  `method = "theta"` in
  [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
  and
  [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md)
  gives the 0.3.0 results. **Results change.**
  - [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
    adds Kreiner’s expected gamma (`gamma_exp`), a mid-p `ppp` and a
    max-T adjusted `p_adj`. `gamma_diff` is now observed minus expected,
    and `gamma_obs_q*` is replaced by `gamma_rep_q*`.
  - Conditional infit matches
    [`iarm::out_infit()`](https://rdrr.io/pkg/iarm/man/out_infit.html)
    for the same item parameters.
  - [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
    adds `q3_p_adj`, adjusted over item pairs, and its reference line is
    the mean replicated Q3.
- [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
  `infit_ppp` and
  [`q3_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_post.md)
  `q3_ppp` compared each draw’s observed value with the mean of the
  replicated values, not with the same draw’s replicated value.
  `infit_ppp` is now also P(observed \> replicated) with mid-p ties, so
  values near 1 indicate underfit, as for `q3_ppp`. **Results change.**
- [`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md)
  and
  [`plot_targeting_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting_hpcm.md)
  follow `easyRasch2::RMtargeting()`. The bottom panel shows
  response-category bands, with categories that are never most likely
  and disordered thresholds marked in red and labelled with their
  posterior probabilities. `panel = "thresholds"` gives a
  dot-and-whisker plot. New arguments `panel`, `xlim`,
  `category_labels`, `viridis_*`, `row_gap` and `output = "list"`, plus
  `person_estimate` in
  [`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).
  The defaults for `bins` and the x-axis change, and
  `sort_items = "location"` puts the easiest item at the top. ‘scales’
  is added to Imports.

### New features

- New datasets `phq9` (600 NHANES respondents, nine items with four
  categories, as in easyRasch2) and `phq9_nhanes` (1000 NHANES
  respondents, including zero responses, for the hurdle partial credit
  model).
- [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  gains `method` and `score_breaks` for forming class intervals, as in
  `easyRasch2::RMitemICCPlot()`, and a caption stating the grouping.
  `n_intervals` is renamed `class_intervals` in
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  and
  [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md),
  and the old name is deprecated.
- [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  with `dif_var` separates groups by shape and line type as well as
  colour.
  [`plot_ipf()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_ipf.md)
  adds line types, and
  [`plot_tile()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_tile.md)
  and
  [`plot_stackedbars()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_stackedbars.md)
  use
  [`theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).
- [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md),
  [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md)
  and
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  report how many persons with missing responses they leave out.
- [`item_restscore_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_post.md)
  and
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
  accept any item column name.
- [`item_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_parameters.md)
  `n_eff` is now
  [`posterior::ess_bulk()`](https://mc-stan.org/posterior/reference/ess_bulk.html).
  ‘posterior’ is added to Imports.

### Bug fixes

- [`dif_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/dif_statistic.md)
  for polytomous models loaded differences in group means onto the items
  as DIF. Impact is now estimated separately (new `impact` element), and
  all DIF types share one sign: positive means easier for the focal
  group. **Results change for polytomous models.**
- Rating scale models gave every item the same thresholds. **Results
  change for rating scale models.** `0 + item` models no longer fail in
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  and
  [`plot_targeting()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_targeting.md).
- [`person_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md)
  and
  [`person_parameters_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters_hpcm.md)
  now compute the WLE of persons with missing responses from the items
  they answered, with a new `n_items` column.
- [`person_parameters()`](https://pgmj.github.io/easyRaschBayes/dev/reference/person_parameters.md)
  and
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  counted sum scores from the lowest observed category instead of
  category 1.
- [`q3_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic.md),
  [`q3_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic_hpcm.md),
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  and
  [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
  stop if a person and item combination occurs more than once.
  Previously only the last response was kept.
- [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md)
  placed class intervals in the severity panel at their mean EAP, which
  made observed points look misfitting, and shifted responses down when
  there were no zeros. **Plots change.**
- [`plot_ipf()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_ipf.md)
  and
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  use the model’s link for cumulative, sratio and cratio models, and
  stop for non-logit `acat` models and models that estimate `disc`.
- [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  aligns eigenvector signs across draws by the dot product with the mean
  direction.
- [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  partial gamma used the wrong column for numeric item labels.
- [`posterior_epred_hurdle_acat()`](https://pgmj.github.io/easyRaschBayes/dev/reference/posterior_epred_hurdle_acat.md)
  failed when items had different numbers of categories.
- [`infit_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic_hpcm.md)
  and
  [`q3_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/q3_statistic_hpcm.md)
  failed with ordered-factor responses or a single draw, and with item
  labels that brms renames.

### Documentation

- [`infit_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic.md)
  and
  [`infit_post()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_post.md)
  state that the infit ppp is conservative, also with the conditional
  method, and suggest
  [`item_restscore_statistic()`](https://pgmj.github.io/easyRaschBayes/dev/reference/item_restscore_statistic.md)
  as well.
- Corrected
  [`RMUreliability()`](https://pgmj.github.io/easyRaschBayes/dev/reference/RMUreliability.md)
  return value,
  [`plot_residual_pca()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_residual_pca.md)
  `style` default and centring, and
  [`plot_icc()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc.md)
  partial gamma sign. Added caveats to
  [`posterior_to_prior()`](https://pgmj.github.io/easyRaschBayes/dev/reference/posterior_to_prior.md)
  and
  [`plot_icc_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/plot_icc_hpcm.md).
  [`infit_statistic_hpcm()`](https://pgmj.github.io/easyRaschBayes/dev/reference/infit_statistic_hpcm.md)
  is no longer described as conditional infit.

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
