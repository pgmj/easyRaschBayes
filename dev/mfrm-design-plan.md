# Many-facet Rasch models in easyRaschBayes: vignette and function changes (plan)

Status: plan only. Written 2026-11-05, revised 2026-11-06 after the
conditional item fit rewrite (0.3.0.9000). Nothing MFRM-specific is
implemented or run. Prompted by an MFRM application to an expert panel,
fitted with the brms model below.

## The model

A many-facet Rasch model (MFRM) is a rating scale or partial credit model
with extra additive facets. Data are long, one row per rating, with one
column per facet (for example ratee, rater, criterion). In brms:

```r
rating ~ criterion + (1 | rater) + (1 | ratee) + (1 | rater:ratee)
family = acat(link = "logit")
```

- `acat` with common thresholds = rating scale MFRM. `thres(gr = criterion)`
  gives criterion-specific thresholds (partial credit by criterion).
  `thres(gr = rater)` gives rater-specific scales, the Bayesian version of
  Facets' rater-specific rating scale, useful for central tendency and
  extremity effects.
- `(1 | rater:ratee)` is a testlet term for halo (one rater's view of one
  ratee carried across criteria). Its SD is the size of the dependence.
- Facets with few elements (4 criteria) as fixed effects with sum-to-zero
  contrasts, since an SD from four values is mostly prior.
- Rater discrimination (generalized rating scale MFRM, Guo & Wind 2021) via
  `bf(rating ~ ..., disc ~ (1 | rater))`. Needs identification (fix the disc
  intercept). The total score is then not sufficient, so the conditional
  fit functions refuse such models (by design, see below).
- Weakly informative priors on the logit scale: thresholds normal(0, 3),
  fixed facet effects normal(0, 1), SDs exponential(1). Prior sensitivity as
  standard.

The "person" in the existing functions is not always the ratee. Which column
plays which role depends on the question. That is the main design problem.

## What changed on 2026-11-05/06 that matters here

- **Duplicate-cell guard done.** `.check_unique_person_item()` in
  `R/input_checks.R` stops Q3, residual PCA and the conditional functions
  when a person x item combination occurs more than once (the old silent
  overwrite in `resid_obs_wide[lin_idx] <- ...`). Step 1 of the old order of
  work is therefore done.
- **Conditional fit engine** (`R/conditional_score.R`). Works row by row: it
  takes category probabilities from `posterior_epred`, removes the person
  location using the log-odds of the person's first answered item, groups
  persons by their item weights (hash), and conditions on the total score
  with elementary symmetric functions. Consequences for MFRM:
  - With `person_var = ratee` and the "item" a rater x criterion cell, each
    cell has its own location including rater severity. The ratee's total
    over all ratings is sufficient as long as the model is Rasch-family
    (acat or bernoulli, logit, no `disc`).
  - Incomplete designs (each ratee rated by a few raters) are missing
    responses, which the engine already handles (per answered-cell pattern).
  - `item_restscore_statistic()` therefore has an MFRM analogue after all:
    conditional gamma of each rater x criterion cell against the ratee's rest
    score. The old plan's "no natural analogue" verdict is wrong.
  - Rater-level (or criterion-level) conditional infit is an aggregate over
    rows: sum of squared conditional residuals over a rater's rows divided by
    the summed conditional variances. Cheap once E and Var exist per row.
- **Null calibration** (dev/results/null_restscore*, 2026-11-06, standard
  person x item designs): conditional restscore ppp close to uniform,
  conditional infit ppp conservative (two-sided rej05 0.5-2%). The infit
  docs now suggest also using item-restscore. Expect the same pattern for
  facets, to be checked with sparse rating designs.
- **Shared item extractor** `.item_tau_draws()` (`R/item_draws.R`). Handles
  one item variable: grouped thresholds, common thresholds with item random
  effects, or item fixed effects with treatment coding. Does not handle
  several random facets or sum-to-zero contrasts.
- **DIF** now separates impact from DIF and defines DIF relative to the
  average item. The same identification logic applies to facet interactions.
- **plot_targeting** ported from RMtargeting, with helpers in
  `R/targeting_helpers.R` (modal category bands, collapse and reversal
  probabilities per draw) that a facet map can reuse.

## How the existing functions map (revised 2026-11-06, from reading the code, not tested)

| Function | Works for MFRM now? | Change needed | Size |
|---|---|---|---|
| `infit_statistic()` | Only if the data has a single cell column (rater x criterion) passed as `item_var`. `model$data` holds only formula variables, so the user cannot pass `interaction(rater, criterion)` today | Accept `item_var = c(rater, criterion)` and build the cell internally. Add `by = rater` (or criterion) that aggregates conditional residuals and variances to facet-level infit/outfit. `method = "theta"` aggregates the same way | Small to medium |
| `item_restscore_statistic()` | Same as infit | Same cell `item_var`. No `by` (gamma per cell does not aggregate sensibly). Check the max-T adjustment with many cells | Small |
| `infit_post()`, `item_restscore_post()` | Yes | Labels and sorting for facet output, many cells per plot | Small |
| `q3_statistic()`, `plot_residual_pca()` | Yes if a cell id column exists | Accept `person_var = c(rater, ratee)` (cell), criterion as item: Q3 across criteria within rater x ratee = halo check | Small |
| `RMUreliability()` | Yes (user supplies draws) | `RMUreliability_facet(model, facet)` extracting r_facet draws. Document per facet: for raters, high RMU means raters differ, not good rating quality | Small |
| `plot_icc()` | No | Conditional ICC per rater (or cell), class intervals on the ratee total score | Medium |
| `.item_tau_draws()` | Partly | Broader facet extractor: several random facets, sum-to-zero fixed facets, `thres(gr = rater)` | Medium |
| `item_parameters()`, `person_parameters()` | Partly | Replaced for facets by `facet_parameters()` (below). Keep as is for standard models | n/a |
| `dif_statistic()` | No | Replaced for facets by `interaction_statistic()` (below) | n/a |
| `plot_targeting()` | No | Replaced for facets by `plot_facet_map()` (below), reusing the targeting helpers | n/a |
| `plot_ipf()` | Criterion or rater curves if thresholds are grouped by that facet | Check, doc only | Small |

General: a shared helper `.resolve_facet()` (columns exist in `model$data`,
cell construction, number of elements, random or fixed in the formula) to
avoid repeating validation in each function.

## New functions

1. **`facet_parameters(model, facet)`**. Element measures with credible
   intervals for random or sum-coded fixed facets, one documented sign
   convention per facet (positive = higher ratings). Includes the fair
   average: expected rating per element with other facets at their mean,
   via `posterior_epred(newdata)`, comparable to Facets' Fair Avg. Medium.
2. **`plot_facet_map(model, facets)`**. Facets-style variable map: one
   column per facet (posterior means with intervals) plus the category bands
   from `R/targeting_helpers.R`. theme_bw, colour + shape on every series.
   Medium.
3. **`interaction_statistic(model, facet_a, facet_b)`**. Bias analysis
   (rater x criterion, rater x ratee). Refits with `(1 | a:b)` or a fixed
   `a:b`, summarises element-pair effects relative to the average pair (as in
   the DIF fix). No default ROPE (no rule-of-thumb 0.5 logit). Require a
   user ROPE or report posteriors only. Medium, slow because of the refit.
4. **`rater_agreement_statistic(model, rater, ratee, criterion)`** (new
   2026-11-06). Exact agreement between rater pairs on the same ratee and
   criterion, observed vs posterior predictive replicates. Facets reports
   observed vs expected exact agreement and a Rasch-kappa. Agreement above
   the model's expectation suggests raters rating to a local norm rather than
   independently. Small to medium.
5. **Rater-rest association**: dropped. The conditional restscore with
   rater x criterion cells covers most of it.
6. **Calibration plot** with `by = facet` (see dev/pav-calibration-plan.md).
   Use a `by` argument from the start.

## Vignette: "Many-facet Rasch models with brms" (precompiled like pcm-rasch-analysis)

1. What an MFRM is, long data layout, choosing the object of measurement.
2. Data. Prefer a simulated dataset with known rater leniency, one
   inconsistent rater (low discrimination), a halo SD, and a ceiling, so each
   diagnostic can be checked against the truth. Optionally one public dataset.
   Candidates to verify: datasets in `immer`, `sirt` (data.ratings*), `TAM`.
   Not the reviewed paper's data (not public).
3. Model: rating scale MFRM, then the halo term, then rater-specific
   thresholds. Priors and prior predictive check.
4. Element measures, fair averages, facet map.
5. Fit per facet: conditional infit aggregated by rater, and conditional
   restscore per rater x criterion cell. Infit is conservative, so lead with
   restscore. The masking point from Guo & Wind (2021): one misfitting rater
   shifts the others.
6. Local dependence: Q3 between criteria within rater x ratee cells, and the
   halo SD. Compare SEs and RMU with and without the halo term.
7. Interactions (bias): rater x criterion. Rater agreement.
8. Reliability: RMU per facet with interpretation per facet.
9. Calibration: PAV cumulative-event plots, in-sample and LOO.
10. Sample size: small panels give wide intervals. Say what the model can
    and cannot show with 8 raters.

## Open questions

- **Halo term and conditional fit.** `(1 | rater:ratee)` passes the model
  check (only person-level terms other than the intercept are refused).
  Conditioning still works with the halo draws in the weights, but every
  ratee then becomes its own weight group (slower), and fit means fit after
  the halo is absorbed. Decide whether to allow, warn or refuse. Check on
  simulated data.
- Rater-specific thresholds and the conditional engine: it uses
  `posterior_epred`, so per-row thresholds should be read correctly. Verify
  against iarm with a cell-as-item layout and CML parameters.
- Calibration of conditional restscore and infit in sparse rating designs
  (few raters per ratee, many cells, many missing). Null run on the host.
- Max-T adjustment with many cells (rater x criterion can be 40+).
- Identification and priors for rater discrimination (only for the theta
  method and descriptive use, since conditional fit refuses `disc`).
- Fixed vs random for small facets: offer both, recommend fixed for < 5-6
  elements?

## Order of work (revised 2026-11-06)

1. ~~Q3 duplicate-cell check~~ (done, `R/input_checks.R`).
2. `.resolve_facet()`, cell `item_var`/`person_var` in the conditional
   functions, Q3 and residual PCA, and `by` aggregation in
   `infit_statistic()`. Test against iarm with cells as items.
3. Simulated MFRM data + model fits for the vignette, and a small null run
   for sparse designs (host scripts, ask before running). Settles the halo
   question.
4. `facet_parameters()` + fair average, broader facet extractor,
   `RMUreliability_facet()`.
5. `plot_facet_map()`.
6. `interaction_statistic()`, `rater_agreement_statistic()`.
7. Calibration function (separate plan).
8. Vignette.

## Effort (estimate 2026-11-06, in sessions like 2026-11-05)

- Minimal (steps 2 and 4 plus a vignette using existing functions for the
  rest): about 2 sessions.
- Full plan (steps 2-8): about 5 to 6 sessions, plus two small host runs
  (vignette fits, sparse-design null run).

Each new function: tests, NEWS (brief), `_pkgdown.yml`, README check
(easyrasch2-doc-review-checklist applies to this package too).
