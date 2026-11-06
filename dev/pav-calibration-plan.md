# PAV-adjusted calibration checks for ordinal Rasch models (plan)

Status: plan only, 2026-11-05. Nothing implemented or run.

## Why

Säilynoja, Johnson, Martin & Vehtari (2025, arXiv 2503.01509, Section 6 and
Discussion) recommend PAV-adjusted calibration plots of the cumulative event
probabilities for ordinal models, one plot per event (M - 1 plots). Bar plots
of category frequencies are uninformative for ordinal models, since nearly any
model reproduces the marginal frequencies. For flexible models and small data
they recommend leave-one-out predictive distributions.

## What bayesplot has (checked in github-repos-cloned/bayesplot, commit 3b03f62a, 2026-10-12, 1.16.0.9000)

- `R/ppc-calibration.R`: `ppc_calibration()`, `_grouped()`, `_overlay()`,
  `ppc_loo_calibration()`, `ppc_loo_calibration_grouped()`,
  `ppc_calibration_data()`. All abort unless y is 0/1. NEWS: "for models with
  binary outcomes". No branch or commit on ordinal/categorical/cumulative
  calibration.
- Ordinal use is possible by hand only: binarise each cumulative event and loop.
- The `yrep` route (`.calibration_interval_data_from_yrep`) is the paper
  version: p = colMeans(yrep), PAV once on y ranked by p, band from PAV of each
  yrep draw in the same ranking. The `prep` route averages per-draw PAV curves
  rank by rank, which is the artefact found in the ordbeta null simulation
  (analyses/ordbeta_ppc). LOO resamples each observation's draws with PSIS
  weights and then runs the same steps.

## Starting point

`pav_calibration_core()` in `analyses/ordbeta_ppc/ppc_ordbeta.R`: takes an
S x N matrix of per-draw event probabilities, optional LOO weights, and
`band = "joint"` (one posterior draw per replicate, Säilynoja) or
`"independent"` (events drawn from p_bar, Dimitriadis et al. 2021). Only event
construction needs to change.

## Design

1. **Models.** `bernoulli` (1 event); `acat`, `cumulative`, `sratio`,
   `cratio` (K - 1 events from the S x N x K `posterior_epred` array);
   `hurdle_acat` (`posterior_epred_hurdle_acat` returns S x N x K). Unordered
   `categorical` (one-versus-others) out of scope.
2. **Events.** P(y >= k), k = 2..K ("at least category k"), the complement of
   the paper's P(y <= k). Same information, reads in the direction of
   thresholds and ICCs. Per-draw event probabilities by reverse cumulative sum
   over the category axis.
3. **Grouping.** Default: one panel per event pooled over items (Säilynoja
   Fig. 29). `by_item = TRUE`: item x event panels, a model-based item
   calibration check. Needs a cap or paging for long scales. Items with fewer
   categories get fewer events. Skip empty events.
4. **In-sample vs LOO.** In-sample default `band = "joint"`. LOO: PSIS weights
   from log_lik, force `band = "independent"` (each observation has its own LOO
   posterior, no joint draw to simulate from), as in the ordbeta prototype.
   Accept a user `psis_object` (moment matching, reloo). Return Pareto k with
   the data and message on high values. High k is likely with few items per
   person. `log_lik_hurdle_acat` is exported, so LOO works for hurdle models.
5. **Bands.** Ordbeta null results: bands are conservative in-sample and LOO.
   Only fresh y gave near-nominal exceedance. In-sample curves bend outward at
   the extremes under a correct model, LOO removes most of this. Document: read
   the curve shape, the band is a guide, the share outside is descriptive and
   not a test.
6. **API.** `plot_calibration(model, item_var = item, person_var = id,
   loo = FALSE, psis_object = NULL, band = c("joint", "independent"),
   by_item = FALSE, prob = 0.95, B = 1000, ndraws_use = NULL)`, argument style
   as in `infit_statistic()`/`plot_targeting()`. Data access: companion
   `calibration_data()` or a `return_data` argument (check what the other
   `plot_*` functions do first). theme_bw(), readable facet strips, colour +
   line type for curve/band/diagonal. Optional quantile dots for the
   distribution of predictions.
7. **Memory.** The epred array is S x N x K (N = 10,000, S = 4,000, K = 5 is
   about 1.6 GB). Use `ndraws_use`, build events one at a time without
   copying the array. PAV via `stats::isoreg` (or `monotone::monotone` if
   installed, as bayesplot does).
8. **Dependencies.** `loo` (psis, weights) in Imports, or Suggests with
   `requireNamespace()`. Already a brms dependency. bayesplot not needed.
9. **Tests.** Small simulated PCM (fit once, cache): output shape, events per
   item, unequal category counts, `bernoulli` one event, `hurdle_acat` event
   construction, LOO path with supplied `psis_object`, alignment when brms
   dropped NA responses.
10. **Validation (ask before running).** Reuse the ordbeta null design for a
    correctly specified PCM: band exceedance in-sample, LOO, fresh y. One
    misfit condition (an item with lower discrimination). Check whether the
    in-sample tail bend and the LOO correction replicate.
11. **Docs.** NEWS, `_pkgdown.yml` reference section, vignette model-checking
    section.

## Upstream option

A bayesplot issue proposing an ordinal wrapper over the binary core. The
`yrep` route has no independent band for LOO, so an easyRaschBayes version
may still be useful.

## Related

- dev/mfrm-design-plan.md (MFRM vignette uses this check).
- Manual binarised `ppc_loo_calibration()` loop for an MFRM, written
  earlier, the starting point for the ordinal version.

## Check of a conditional variant (2026-11-06, container, exploratory)

Idea: for Rasch models, use P(X_i >= k | total score) from the conditional
engine instead of theta-based predictions, with the band from replicates
drawn conditional on each person's total score. Script
`dev/conditional_calibration_check.R`, results in
`dev/results/conditional_calibration/`. Frequentist analogue only: CML item
parameters (eRm), comparator is ML theta plug-in with independent
replicates, not the brms posterior or LOO. PCM, 9 items, 4 categories,
200 reps, n = 500 and 200, item 1 with discrimination 1 (null), 0.5 or 2.

- Conditional curves have mean deviation exactly 0 per item and event (CML
  fits the category counts), so they show shape (slope) only, never
  over- or underprediction of an event overall.
- Pointwise band exceedance under the null: 2.1-2.7% (conditional),
  1.9-3.6% (theta). Both conservative. "Any point outside" happens in
  73-89% of null items, useless as a flag.
- Mean PAV slope under the null 1.01-1.03 for both, so no extra in-sample
  bend from theta here (ML theta depends only on the raw score).
- Item 1 with the highest exceedance, n = 500 / 200 (chance 0.11):
  low discrimination cond .995 / .865, theta .895 / .685;
  high discrimination cond .805 / .340, theta .985 / .650.
  Misfit spills less onto other items with the conditional version.
- Conclusion: works, but for Rasch models it repeats what item-restscore
  and plot_icc() already show (discrimination). It cannot be used for the
  case that motivated calibration (PCM vs hurdle PCM), since the hurdle
  model has no sufficient total score.
