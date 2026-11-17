# easyRaschBayes code review, 2026-11-05

Review of all code not already rewritten in 0.3.0.9000 (the conditional
restscore and infit work). Four parallel reviews, findings checked against the
code. "Verified" means reproduced by running code (synthetic data, small
fits or mocked inputs). Scratch scripts are in the session scratchpad and are
not kept. Nothing below is fixed yet.

## Critical

1. **FIXED 2026-11-05.** **Polytomous DIF confounds impact with DIF.** `R/brms_dif.R:286-309`.
   Neither the uniform (`... + group:item`) nor the non-uniform
   (`thres(gr = item___group)`) formula has a group mean term, and `(1 | id)`
   has mean 0 in both groups. A latent mean difference is loaded onto every
   item. Verified: simulated impact 0.8 with DIF on I3 only. Uniform contrast:
   the five clean items show about +0.8, I3 shows 0. Non-uniform: five clean
   items flagged, I3 not. The dichotomous formula
   (`1 + group + (1 + group | item)`) is fine and recovered I3 = -0.92
   (truth -1). Fix needs a group mean in the person distribution plus an
   identification constraint (anchor items or sum-to-zero DIF), or the
   dichotomous random-slope approach carried over to thresholds.
2. **FIXED 2026-11-05.** **Uniform polytomous DIF parameter is not identified.**
   `R/brms_dif.R:290-295, 606-614`. Without a group main effect brms codes
   `group:item` with full dummies, each pair collinear with the item's
   thresholds. `b_groupB:item<i>` is held only by the flat prior. Verified:
   posterior SD about 2, every item "Inconclusive".

## High

3. **FIXED 2026-11-05.** **Persons with missing responses get the full-test WLE for their partial
   sum score.** `R/person_parameters.R:319-352`, same in
   `person_parameters_hpcm()` (`R/brms_hpcm.R:1569-1590`) and the hurdle
   class intervals in `plot_icc_hpcm()`. Verified: 6 of 12 items answered,
   score 0 gets -3.52 (12-item value) instead of -2.97.
4. **FIXED 2026-11-05.** **Ordinal sum-score offset uses the lowest observed category, not
   category 1.** `R/person_parameters.R:302-314`. If no one uses category 1,
   every sum score drops by 2 x n_items. Verified with mocked data. Ordered
   factor responses error.
5. **FIXED 2026-11-05.** **DIF sign convention differs between modes and is undocumented.**
   Dichotomous and uniform: positive = easier for `levels[2]`. Non-uniform:
   positive = harder. Forest plot label says only "B vs. A". Verified.
6. **FIXED 2026-11-05.** **`plot_icc()` partial gamma shifted one column with numeric item labels.**
   `R/cicc_plot.R:467-470`, `wide_responses[[item_label]]` indexes by
   position. Verified on raschdat1. Fix: `as.character(item_label)`.
7. **FIXED 2026-11-05.** **Residual PCA sign alignment pulls real loadings toward zero.**
   `R/brms_globalfit.R:327-338`. Majority vote of signs against draw 1 lets
   many near-zero items outvote a small strongly loading cluster. Verified on
   synthetic eigenvectors: 86 of 500 draws mis-flipped, CI for the strong
   items crossed zero. Little effect on the two real fits. Fix: flip when the
   dot product with a running reference is negative.
8. **FIXED 2026-11-05.** **Rating scale models (`acat`, `(1 | item)`, no `thres(gr = item)`):
   every item gets identical locations.** `.extract_item_location_draws`
   (`R/brms_globalfit.R:563-566`), `.extract_threshold_data`
   (`R/brms_plots.R:837-840`), `plot_ipf()` (`R/brms_plots.R:203-207`).
   `r_item` is never added. Verified with an RSM fit.
9. **FIXED 2026-11-05.** **`0 + item` models crash** `plot_residual_pca()` and `plot_targeting()`
   (`b_item*` columns never looked for). Verified with mock draws.
10. **FIXED 2026-11-05.** **`posterior_epred_hurdle_acat()` errors when items have different
    numbers of categories** (`R/brms_hpcm.R:422-432`, no zero padding unlike
    brms). Triggered whenever an item's top category is unobserved, since
    brms sets thresholds per item from the observed maximum. Breaks
    `infit_statistic_hpcm()` and `q3_statistic_hpcm()`. Verified.

## Medium

11. FIXED 2026-11-05 (ranef labels, safe threshold labels, errors instead of zeros). DIF columns stay all zeros, flagged "Negligible", when a coefficient name
    is not found (e.g. group levels with spaces). `R/brms_dif.R:599, 626-628,
    726, 739-749`. Suspected.
12. PARTLY FIXED 2026-11-05 (now warns). Non-uniform DIF pairs threshold t with t across groups even when a
    category is missing in one group. `R/brms_dif.R:687-693`. Suspected.
13. FIXED (complete responders, message). `plot_icc()` with missing responses places persons too low on theta.
    `R/cicc_plot.R:322-355, 383-390`. Suspected, iarm shares it.
14. FIXED (link via make.link; acat logit only; disc stops). `plot_ipf()` ignores the link function and `disc`, always logistic.
    `R/brms_plots.R:224-284`. Suspected.
15. DOCUMENTED (design limit, not fixable in brms priors). `posterior_to_prior()` builds independent normals, dropping threshold
    correlations, which overstates precision of the overall level (SD of mean
    threshold 0.102 posterior vs 0.074 prior on fit_pcm). Undocumented.
    Verified. The 1.7 logit/probit rescaling is correct.
16. FIXED in the plot_targeting port. `plot_targeting()` `sort_items = "location"` puts the easiest item at the
    bottom, docs say top. `R/brms_plots.R:742-746`. Verified.
17. FIXED (message, docs). `plot_residual_pca()` drops every person with any missing response,
    silently. `R/brms_globalfit.R:300`. Suspected.
18. FIXED. person_parameters roxygen describes the old Newton-Raphson WLE and NA at
    extremes, code uses uniroot and finite extremes. Verified.
19. FIXED (docs match code). `RMUreliability()` docs promise a list with `hdci` and posterior draws,
    code returns only the `mean_hdci` data.frame; `verbose` default is FALSE
    in code, TRUE in docs. The method itself matches Bignardi. Verified.
20. FIXED (brms coding, single draw; checked with the mocked fit). hPCM: ordered-factor responses crash infit and Q3, and would shift scores
    in `person_parameters_hpcm()`. `ndraws_use = 1` crashes infit and Q3
    (`drop = TRUE`). Verified with a mocked fit.

## Low

21. FIXED (posterior::ess_bulk with chains). `.compute_ess()` runs Geyer on concatenated chains, ignores between-chain
    disagreement. `posterior::ess_bulk` would be safer.
22. FIXED. `.extract_threshold_data` is defined twice (`R/brms_globalfit.R:614`,
    `R/brms_plots.R:817`), identical apart from whitespace; the brms_plots.R
    copy wins. Delete one.
23. FIXED (docs). `plot_residual_pca()` `style` default is "both", docs say "density".
24. DOCUMENTED. PCA and targeting centre differently (item means vs all thresholds),
    docs say they match.
25. FIXED (positional matching everywhere, incl. hurdle and posterior_to_prior). Regex escaping in threshold and ranef patterns omits `[`, `]`, `\`.
    brms-renamed labels (spaces, punctuation) break extraction in several
    places, including `.extract_hurdle_item_draws`.
26. DOCUMENTED. `plot_icc_hpcm()` bins persons by their own EAP and plots E[f(theta)] at
    E[theta], flattening the observed points. Interpretation caveat.
27. FIXED. `plot_stackedbars()` legend runs highest to lowest despite the comment.
28. FIXED (plot_targeting in the port). House style: `plot_tile()` and `plot_stackedbars()` use `theme_minimal()`;
    `plot_ipf()` and `plot_targeting()` distinguish categories by colour
    only; `plot_targeting()` uses deprecated `geom_errorbarh()`.
29. FIXED. `gk_gamma_vec()` errors on NA input (callers filter first, safe now).

## Checked and fine

- Item locations (RM and PCM) match psychotools CML up to centring and prior
  shrinkage. WLE values match an independent maximiser at every score.
  Person rows align with persons despite ranef order. Q3 matches a
  hand computation. RMU follows Bignardi et al. Gamma pair counting correct.
- Dichotomous DIF recovers simulated DIF. ROPE, pd, CI as documented.
- `plot_icc()` partial gamma matches `iarm::partgam_DIF` with character item
  labels; model curves and class intervals match the docs.
- `fit_statistic_pcm/rm()` follow Bürkner (2020), per-draw pairing correct.
- Residual PCA residuals, SVD eigenvalue and ppp match a hand computation.
- hurdle_acat: Stan and R likelihoods agree, posterior_predict frequencies
  match exact probabilities, sign conventions and centring correct, hPCM
  infit and Q3 formulas correct.
- No other summarise-overwrite cases and no other `\%` in roxygen.
- Note: fitting `hurdle_acat` fails with CRAN brms 2.23
  (`stan_hurdle_ordinal_lpmf` assertion); the package documents a GitHub
  brms requirement.

## plot_targeting() port to the RMtargeting() style: decisions (2026-11-05)

1. Person estimates: new argument, EAP (posterior means of r_id) as the
   default, WLE as the alternative.
2. Adopt RMtargeting() defaults: category-band bottom panel
   (`panel = c("categories", "thresholds")`), bins = half the maximum score
   with a floor of 11, `xlim = c(-4, 4)` expanded, integer breaks,
   `ci_level`-style NULL to hide intervals, captions, `output = "list"`,
   `category_labels`, viridis arguments, `row_gap`.
3. Include the Bayesian additions: per-draw posterior probability that each
   category is never the most likely response and that each adjacent
   threshold pair is reversed, and intervals for band boundaries.
4. plot_targeting_hpcm() is not part of this round.
Also fix findings 8, 9, 16 and the colour-only categories while porting.

## DIF fix as built (2026-11-05)

- Uniform polytomous: `rhs + difgroup + (0 + difgroup | item)`, numeric 0/1
  focal indicator. Impact = `b_difgroup`, DIF = item slopes from
  `ranef(summary = FALSE)` (original labels).
- Non-uniform: thresholds per item x group with index labels `i<k>g<1|2>`.
  Per draw, focal minus reference threshold differences are centred on their
  mean over items (items weighted equally); impact = minus that mean, DIF =
  minus the centred differences.
- Dichotomous model unchanged, extraction moved to ranef/fixef.
- All types: positive DIF = easier for the focal group, positive impact =
  focal group higher. New `impact` element in the output.
- Verified on the review's simulated data (impact 0.8, DIF on I3): uniform
  I3 -0.61, impact 0.69; non-uniform I3 -0.89 / -0.51, impact 0.64, one of
  10 null thresholds flagged (I6 t1); dichotomous I3 -0.92, impact 0.74.
  Labels with spaces give identical results.

## High fixes as built (2026-11-05)

- New `R/item_draws.R`: `.item_tau_draws()` is the one extractor for item
  thresholds (grouped thresholds, rating scale with random or fixed item
  effects, bernoulli random or fixed items). Matches parameters to item
  labels by brms level order, checked against plain labels. Used by
  person_parameters, item_parameters (and the hpcm callers), the PCA, the
  targeting threshold table, plot_ipf and plot_icc. The duplicate
  `.extract_threshold_data` and the old per-function extractors are gone.
  `.find_thresh_cols` remains for plot_icc_hpcm only.
- WLE per set of answered items via `.wle_by_pattern()`, in
  person_parameters and person_parameters_hpcm (the hpcm change is not run
  against a fit, since CRAN brms cannot fit hurdle_acat).
  plot_icc_hpcm class intervals still assume complete responses.
- Offsets: ordinal base category 1, ordered factors converted, in
  person_parameters and plot_icc. plot_icc_hpcm keeps its own shift
  (finding 26/hpcm 6).
- Smoke run: item_parameters, person_parameters, plot_targeting,
  plot_residual_pca, plot_icc, plot_ipf all run on PCM, RM random items,
  RSM and RM `0 + item`. Complete-data person_parameters identical to the
  committed version. Partial gamma with numeric labels equals iarm.

## plot_targeting port as built (2026-11-05)

- Helpers in `R/targeting_helpers.R`, geometry ported from easyRasch2
  RMtargeting; captions plain text (no ggtext). New `.category_posterior()`
  runs the modal-category hull per draw: P(category never modal), lower
  boundary median and interval, P(adjacent pair reversed).
- Red tick labels "c (P = .xx)", reversal labels "gap (P = .xx)" centred
  under the arrow; row gap widens with ticks and reversals.
- Existing argument positions kept, new arguments appended. `bins` default
  now half the max score (min 11). `palette` overrides viridis in both
  panels. Thresholds panel uses colour + shape.
- Smoke: all four model types x five option sets. Tests: grid check of the
  hull boundaries against argmax of category probabilities.
