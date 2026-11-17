# Many-facet Rasch models in easyRaschBayes: design plan

Status: plan only, nothing MFRM-specific implemented. First written
2026-11-05, revised 2026-11-06 after the conditional fit rewrite, revised
again 2026-11-13 (0.4.0.9001) with these decisions from the user:

- Keep `person_var` and `item_var` as in the standard models. The new facet
  is `facet_var` (neutral name). Examples and the vignette use raters, the
  usual reason for an MFRM, so `facet_var = rater` in the code shown to users.
- Dichotomous MFRMs are supported throughout, not only rating scales.
- Small fixed facets use orthonormal sum-to-zero contrasts (see below).
- Minimal scope first (steps 1-4 below plus a vignette), the rest later.
- The halo term is allowed in the conditional functions, with a message.

## The models

Data are long, one row per rating, columns `id` (person, the object of
measurement), `item` (criterion), `rater` (the extra facet), response.

```r
# Polytomous: rating scale MFRM
rating ~ 1 + item + (1 | rater) + (1 | id)                  # family = acat("logit")
rating | thres(gr = item)  ~ 1 + (1 | rater) + (1 | id)     # partial credit by item
rating | thres(gr = rater) ~ 1 + item + (1 | id)            # rater-specific scales

# Dichotomous
y ~ 1 + item + (1 | rater) + (1 | id)                       # family = bernoulli("logit")

# Halo (testlet) term, either model
... + (1 | rater:id)
```

- Rater discrimination via `disc ~ (1 | rater)` (Guo & Wind 2021) is
  descriptive only. The total score is then not sufficient and the
  conditional functions refuse the model, as they do now.
- Priors on the logit scale: thresholds normal(0, 3), SDs exponential(1),
  fixed facet coefficients see the next section.

## Coding of small fixed facets (sum-to-zero)

Source: Morris, "The Sum-to-Zero Constraint in Stan", Stan case study
(mc-stan.org/learn-stan/case-studies/sum_to_zero_vector.html). Read
2026-11-13. Two slips in it: the first prior snippet has the sqrt in the
wrong place (`sqrt((N / (N - 1.0)) * sigma)`, the later model code has
`s * sigma` with `s = sqrt(N / (N - 1))`), and the "sum_to_zero_vector, tiny
dataset" block refits the hard-constraint model.

- `contr.sum` is the case study's "hard" constraint: the last level is minus
  the sum of the others. With iid normal(0, s) on the K - 1 coefficients the
  last level has prior variance (K - 1) s^2, the others s^2. With 4 items the
  last one gets three times the prior variance. This matters where the prior
  carries weight, which is the small-panel case where fixed facets are
  recommended.
- Orthonormal contrasts (normalised Helmert, H with t(H) H = I and
  H t(H) = I - J/K) span the same sum-to-zero space as Stan's
  `sum_to_zero_vector`. iid normal(0, s) on the coefficients gives every
  level variance s^2 (K - 1) / K and correlation -1 / (K - 1). Prior
  normal(0, sigma * sqrt(K / (K - 1))) gives marginal SD sigma per level.
- Checked 2026-11-13 with brms 2.23.0: custom contrasts set with
  `contrasts(d$item) <- H` are used in the design matrix and kept in
  `attr(model$data$item, "contrasts")`. brms has no sum_to_zero_vector type
  (not in the 2.23 namespace), so contrasts are the route.
- New exported helper `contr_sum_zero(K)` (or taking a factor), returning
  H with level names as row names. Docs give the prior scale.
- Random facets (persons, raters with many elements) stay as brms random
  intercepts, identified softly by the mean-zero prior. Reported element
  effects are centred within each draw, as in the DIF fix. Summaries only,
  no change to the model.
- Interactions (facet x item bias): centre the interaction draws over both
  margins within each draw, so each cell is relative to its facet element's
  and its item's average. This is the definition of "relative to the average
  pair" for `interaction_statistic()`.

**Latent issue in the existing extractor (standard models too).**
`.item_effect_draws()` in `R/item_draws.R` assumes treatment coding whenever
there are K - 1 item coefficients. With any other contrasts (contr.sum,
Helmert, the new helper) the label check stops the function if the item
labels are plain letters and digits, but with other labels (spaces,
punctuation) the item locations come out silently wrong. Fix: read the
contrast matrix from `attr(model$data[[item_name]], "contrasts")` (or the
default `contr.treatment` when absent) and map coefficients back with
`b %*% t(C)`. This also removes the treatment special case. Do this in
step 2, with a test for each coding, and note it in NEWS if any released
behaviour changes.

## Disconnected designs

The case study needs a separate constraint per connected component of a
spatial graph. The rating-design analogue: groups of raters who never rate
the same persons. Measures are then not comparable across groups except
through the prior (Facets calls these subsets). `.resolve_facet()` builds the
person-facet bipartite graph and warns when it has more than one component,
naming the components' facet elements.

## Interface

| Argument | Role | Default |
|---|---|---|
| `person_var` | Object of measurement | `id` |
| `item_var` | Item or criterion | `item` |
| `facet_var` | The extra facet (rater, task, occasion) | `NULL` = standard model, unchanged behaviour |

- Output gets a `facet` column next to `item` (fixed names, as `item` is
  fixed now regardless of `item_var`'s name). Results per item x facet cell
  have one row per cell.
- `infit_statistic(..., by = c("cell", "item", "facet"))`: `"cell"` is the
  default when `facet_var` is set. Aggregation sums conditional squared
  residuals and conditional variances over the rows of each element. Same for
  `method = "theta"`.
- `item_restscore_statistic()`: per cell only (gamma per cell does not
  aggregate sensibly). Check the max-T adjustment with 40+ cells.
- `q3_statistic()`, `plot_residual_pca()`: with `facet_var`, residuals are
  taken within person x facet units and correlated between items. This is
  the halo check (Q3 between criteria within one rater's ratings of one
  person).
- Internal `.resolve_facet()`: columns exist in `model$data`, cell factor,
  number of elements, random or fixed in the formula, contrasts, connectivity.

## How the existing functions map

| Function | Change | Size |
|---|---|---|
| `infit_statistic()` | `facet_var`, cells as items for the conditional engine, `by` aggregation | Small to medium |
| `item_restscore_statistic()` | `facet_var`, cells as items | Small |
| `infit_post()`, `item_restscore_post()` | Labels and sorting for cells, many cells per plot | Small |
| `q3_statistic()`, `plot_residual_pca()` | `facet_var`, person x facet units | Small |
| `RMUreliability()` | `RMUreliability_facet(model, facet)` from centred element draws. Docs: for raters, high RMU means raters differ, not good rating | Small |
| `.item_tau_draws()` / `.item_effect_draws()` | Contrast-aware mapping (fixes the latent issue above), several random facets, `thres(gr = facet)` | Medium |
| `plot_icc()` | Conditional ICC per cell or facet element (later) | Medium |
| `plot_ipf()` | Doc only, works when thresholds are grouped by the facet | Small |
| `item_parameters()`, `person_parameters()` | Unchanged for standard models, facets via `facet_parameters()` | n/a |
| `dif_statistic()`, `plot_targeting()` | Facet versions are new functions (below) | n/a |

The conditional engine (`R/conditional_score.R`) needs no change in
principle: cells become items, the facet effect is part of each cell's
location, and incomplete designs are missing responses. The sufficiency
check passes `(1 | rater:id)` today by accident (the `sd_id__` prefix does
not match `sd_rater:id__`). Make that explicit and emit the halo message.

## New functions

1. `contr_sum_zero()`: orthonormal sum-to-zero contrasts. Small.
2. `facet_parameters(model, facet)`: element measures with credible
   intervals, centred within draws, positive = higher ratings. Fair average:
   expected rating (dichotomous: expected probability) per element with the
   other facets at their mean, via `posterior_epred(newdata)`. Medium.
3. `plot_facet_map()`: one column per facet plus the category bands from
   `R/targeting_helpers.R` (item locations only for dichotomous models).
   theme_bw, colour and shape on every series. Medium. (Later scope.)
4. `interaction_statistic(model, a, b)`: bias analysis by refit with an
   `a:b` term, double-centred. No default ROPE. Medium, slow. (Later.)
5. `rater_agreement_statistic()`: observed vs posterior predictive exact
   agreement between pairs of facet elements on the same person and item.
   Small to medium. (Later.)
6. Calibration plot with `by = facet`, see dev/pav-calibration-plan.md.
   (Later.)

## Simulated data and null run (host scripts, ask before running)

- Polytomous and dichotomous versions, same design: known rater leniency, one
  inconsistent rater, a halo SD, a ceiling, an incomplete (sparse but
  connected) rating design, and one disconnected variant for the warning.
- Null run: conditional restscore and infit ppp per cell and aggregated, max-T
  with many cells, with and without the halo term. Settles the halo message
  wording. Prior check: contr.sum vs contr_sum_zero for one small fixed
  facet (prior draws only, cheap).

## Vignette: "Many-facet Rasch models with brms" (precompiled)

1. What an MFRM is, long data, choosing the object of measurement.
2. Simulated data (polytomous main thread, dichotomous section).
3. Models: rating scale MFRM, halo term, facet-specific thresholds. Fixed vs
   random facets and `contr_sum_zero()` for small fixed facets. Priors and
   prior predictive check.
4. Element measures and fair averages.
5. Fit per cell and per facet element. Lead with restscore, since infit is
   conservative. One misfitting rater shifts the others (Guo & Wind 2021).
6. Local dependence: Q3 within person x rater, the halo SD, SEs and RMU with
   and without the halo term.
7. Reliability per facet with the interpretation per facet.
8. What a small panel can and cannot show.
(Later scope adds the facet map, interactions, agreement, calibration.)

## Order of work

1. ~~Q3 duplicate-cell check~~ (done, `R/input_checks.R`).
2. ~~`.resolve_facet()` (incl. connectivity), contrast-aware item extractor,
   `contr_sum_zero()`~~ (done 2026-11-13, `R/facets.R`, `R/item_draws.R`,
   tests in `tests/testthat/test-facets.R`). As built, two changes from the
   plan above: (a) connectivity is the rank deficiency of the additive
   person + item + facet design (persons eliminated by a Schur complement of
   count tables), which also catches items nested in raters, not only
   person-facet groups; the person-facet groups are named in the warning
   when they are the cause. (b) `.resolve_facet()` does not require the
   facet in the formula. It reports the roles (fixed, random, thres) and
   skips the connectivity check when the facet has no role in the model.
   Correction (same day): brms keeps only formula variables in
   `model$data`, so a facet absent from the formula cannot be found at
   all. "Fit without rater effects, then test rater misfit" would need a
   `data` argument. Not planned.
3. ~~`facet_var` in `infit_statistic()` (with `by`), `item_restscore_statistic()`,
   `q3_statistic()`, `plot_residual_pca()`, the post functions~~ (done
   2026-11-13). As built:
   - Wrappers add an internal `.cell` column to a copy of `model$data` and
     call the standard function with it as the item (brms ignores extra
     columns). Output: `cell`, item column, facet column, named after the
     user's variables (the item column is named after `item_var`, not fixed
     as I first said).
   - Conditional infit pools by `group` inside `.cond_infit_draw()` (sums of
     squared residuals, variances, z^2 and counts per group). Theta infit by
     item/facet is just theta infit grouped by that column.
   - **Q3 design changed.** Not "residuals within person x facet units":
     conditioning on the unit total removes a halo shift, so that version
     cannot see halo. Instead cells are items, Q3 for every cell pair, with
     `pair_type` = "same rater" / "same item" / "neither". Halo shows as
     positive Q3 for same-rater pairs. Pairs with no common person give NA
     (changed in 5b: they are now left out).
   - Residual PCA does use person x facet units (theta residuals), where a
     halo shows as a first component with same-sign loadings.
   - Validation: iarm only on a complete cell design (iarm refits complete
     cases with missing data, then errors), brute force for sparse cells and
     pooling, mocked-brms end-to-end tests. Smoke-tested on three small fits
     (RM, RSM, RM + halo) in the container.
   - Found and fixed a released bug: `cbind(integer(0), x)` gives a 1 x 1
     index, so restscore failed for a weight group with no admissible rest
     score (total 0 with a unique missing pattern; every person with halo).
     Same guard added in `.cond_sample()`.
4. `facet_parameters()` + fair average, `RMUreliability_facet()`.
   Design decisions 2026-11-13 (before coding):
   - **Sign: severity, like item difficulty.** The measure is the facet
     element's location on the same side as item difficulty: the model is
     theta - item - facet, so positive = lower ratings (a severe rater).
     This is the Facets convention for raters and matches
     `item_parameters()`. It reverses the earlier plan line "positive =
     higher ratings"; consistency with item locations decided it.
   - **Centred within each draw** (mean of the elements 0), so the result
     does not depend on the coding (random effects, treatment, contr.sum,
     contr_sum_zero) or on where brms puts the intercept.
   - **Supported layouts:** facet as random intercept or fixed effect (any
     contrasts, via the contrast-aware `.item_effect_draws()` applied to the
     facet), or facet-specific thresholds `thres(gr = facet)`, where the
     measure is the mean of the element's thresholds and a `thresholds`
     table is added (rater-specific scales: central tendency shows as
     narrow threshold spacing).
   - **Fair average** (Facets "Fair-M"): expected rating for the element
     with the person location at the sample mean of the person locations
     in that draw, averaged over items, on the response scale of the data
     (1..C for acat, 0/1 for bernoulli). Halo terms are left out (mean 0).
     Computed from the extracted parameters, not `posterior_epred(newdata)`,
     because brms cannot place a person at the sample mean of the person
     effects. Averaging expected ratings over items, rather than plugging in
     an "average item", keeps it defined for item-specific thresholds.
   - **No `RMUreliability_facet()`.** `facet_parameters(draws = TRUE)`
     returns the centred draws matrix (elements x draws), and
     `RMUreliability(fp$draws_matrix)` gives the RMU. One function fewer,
     and the docs can carry the caution: with 3-5 raters a correlation
     across elements is very noisy, and for raters high RMU means raters
     differ reliably, which is not good rating quality.
   As built 2026-11-13 (`R/facet_parameters.R`, tests in test-facets.R):
   - `facet_parameters(model, facet_var, item_var, person_var, draws,
     prob)`. Output list: `summary` (facet column named after the
     variable, n, observed_avg, fair_avg + HDCI, location, se, hdci_lower,
     hdci_upper, n_eff, same column names and rounding as
     `item_parameters()`), `thresholds` (thres(gr = facet) only),
     `facet_sd` (random facet only), `draws_matrix` (draws = TRUE).
   - Column is `location`, not "measure" or "severity", to match
     `item_parameters()`. The docs say positive = lower ratings.
   - Facet effects reuse `.item_effect_draws()` with the facet name, so
     random and fixed facets with any contrasts go through the same
     contrast-aware code as items. Tested: random coding and a
     contr_sum_zero fixed coding of the same effects give identical
     locations and fair averages.
   - Mean person location = row mean of the r_<person> draws. With no
     person effects found it is 0.
   - Fair average: acat and bernoulli vectorised, other ordinal families
     via `.compute_cat_probs()` per draw (slow but correct). Checked on the
     small RM fit against `posterior_epred(newdata, re_formula = ~ (1 |
     rater))` with persons at 0: equal to 3 decimals (the sample person
     mean is near 0 there; they differ in general by design).
   - Ratings scale: bernoulli 0/1, ordinal 1..C (brms category index),
     same for observed_avg, so the two columns are comparable.
   - `.facet_parameters_core()` takes the draws as an argument so tests can
     use fake draws on an empty brms shell.
5. Simulated data + fits (host), null run (host).
   Scripts written 2026-11-13, NOT RUN (user: write, do not run yet).
   Parse-checked only. Run the smoke tests first.
   - `vignettes/fit-mfrm-vignette-models.R` (added to .Rbuildignore).
     150 persons, 4 criteria, 6 raters, each person rated by 3 raters with
     cyclic offsets {0, 1, 3} mod 6 (connected, 12 of 24 cells per person,
     1,800 ratings). Truth saved in `vignettes/fits/mfrm_truth.rds`:
     severity (-0.6 ... 0.4), R6 inconsistent (discrimination 0.4), R3
     central tendency (thresholds -2.6, 0, 2.6 vs -1.5, 0, 1.5), halo SD
     0.6, 4 categories. Dichotomous `pass` from the same persons, design,
     halo and inconsistency. Fits (4 chains x 2000): rsm, rsm_halo,
     rsm_scales (thres(gr = rater) + halo), rm, rm_halo, plus two
     prior-only fits (contr.sum vs contr_sum_zero, same coefficient prior).
     Items use contr_sum_zero with prior normal(0, sqrt(4/3)) so the
     vignette demonstrates the helper. FIT=data|poly|dich|prior.
     Changes from the plan list: the "ceiling" feature is dropped (no
     diagnostic in scope targets it). The prior-only contrast check moved
     here from the null run, since it needs no replications and the
     vignette can show it.
   - `dev/null_mfrm_calibration.R` + `_summary.R`. 2 models (rm, rsm with
     4 categories) x 3 halo settings x n = 300, 600, 200 reps (2,400 fits,
     10 workers, resumable per rep). 5 items x 8 raters = 40 cells, each
     person rated by 2 adjacent raters (cyclic), 10 of 40 cells per person.
     Halo settings: none; modelled (SD 0.6 in data and model, a null
     condition for the halo question); unmodelled (in data only, a power
     check reported separately). Per rep: restscore ppp/p_adj per cell,
     infit ppp per cell, pooled by rater and by item with a max-T p_adj
     computed in the script, Q3 ppp/p_adj per cell pair by pair type,
     facet_parameters() coverage/bias/RMSE for severity against the centred
     truth. Summary: rej05/fwer05 tables, ECDF-minus-uniform figure (theme_bw,
     Okabe-Ito + linetype + shape; run the dataviz validator on the output).
     No run-time measurement yet: extrapolate from the smoke test.
   - Doc fix found while writing the vignette script: `facet_parameters()`
     help said closely spaced thresholds mean central tendency. Reversed:
     wide spacing = middle categories used more = central tendency.
   - 2026-11-13, first smoke test on the host printed rstan's "the number
     of chains is less than 1; sampling not done". Not an error: it comes
     from the `chains = 0` compile step. The compile step now captures that
     output and prints "compiled <model> in N s", and each job prints a
     progress line with its run time.
   - **Null run results (run on host, finished 2026-11-15, 2,400 of 2,400
     reps OK, NDRAWS = 1000, results in dev/results/null_mfrm.rds, figure
     null_mfrm_ecdf.png).** Rates are two-sided at .05 per unit (8,000
     cells, 1,600 rater units, 1,000 item units per condition), FWER over
     200 reps (MC SE about .015).
     - Convergence: median max Rhat 1.01 (max over all parameters, 2
       chains), above 1.05 in at most 1% of reps. Divergences in 5-21% of
       reps, median 1 when present. Rejection rates and coverage are the
       same in flagged and unflagged reps (restscore .048 vs .047, infit
       .031 vs .031, severity coverage within .01).
     - Restscore per cell: calibrated. rej05 .042-.054, ppp SD .28-.29.
       Max-T p_adj FWER .015-.055.
     - Infit per cell, no halo: rej05 .041-.047 (the known slight
       conservatism of conditional infit). Pooled by rater .029-.034, by
       item .010-.028.
     - **Infit with the halo term modelled: conservative.** Per cell
       rej05 .018-.021, by rater .009-.019, by item .004-.013. The ECDF
       shows both tails thin (p-values pulled to the middle). Likely
       mechanism: the halo effects are fitted to the same few ratings per
       person x rater and treated as known, which absorbs residual
       variation in the observed data. Errs toward missing misfit, not
       toward false alarms. Restscore is barely affected (.042).
     - Max-T p_adj for pooled infit: FWER .000-.045 in every null
       condition.
     - Q3 ppp: no halo .043-.053, halo modelled .026-.042. Max-T p_adj FWER
       .000-.055.
     - Halo in data, not in model (power check): mean same-rater Q3 -0.07
       (RM) and 0.00 (RSM) against -0.11 replicated, other pair types
       pushed down to -0.14/-0.20. So the halo is clearly visible in the Q3
       means. But per pair it is weak, and max-T over 780 pairs detects it
       (FWER) in 4-10% of reps for RM and 30-64% for RSM. Infit and
       restscore do not react to halo (rates as under the null), as
       expected for a person x rater shift. Candidate follow-up: one
       pooled halo statistic (mean same-rater Q3 against its replicated
       distribution), which would need its own null run.
     - Severity (facet_parameters, 95% HDCI): coverage .91-.95 with no
       halo or halo modelled (lowest RM n = 300 with halo, .912), RMSE
       .06-.16. **Halo not modelled: coverage .81-.90**, since the
       dependence makes the intervals too narrow. Point estimates stay
       unbiased. Supports recommending the halo term when raters score
       several items of the same person.
     - gamma_obs - gamma_exp: +.013 at n = 300, +.0065 at n = 600 (halves
       with n, the small-sample bias of observed gamma). The ppp compares
       against replicated gamma, so calibration is not affected.
     - Fair minus observed average per rater: at most .01 (RM) and .036
       (RSM) on the rating scale, as expected in a balanced rotation.
     - Cost: statistics 570-1,576 s per halo-modelled rep, 30-124 s
       otherwise. Fits 11-285 s.
5b. Code review fixes (2026-11-13, `dev/code-review/code-review-2026-11-13.md`, done
   before the smoke test was rerun). Status per finding is in that file.
   Design choices:
   - **ESF over answered items only** (finding 1). `.cond_draw_setup()`
     builds each group's ESF tables over the items its members answered.
     The compact arrays are padded to the largest answered set, with weights
     (1, 0, ..., 0) at padding positions, so `.cond_sample()` stays
     vectorised over persons unchanged. `pos_item` maps positions back to
     item columns, and `.cond_replicate(st)` returns the n x k matrix. For
     complete data the positions are the items and the random stream is the
     same as before, so results are unchanged. With missing responses the
     replicated data follow the same distribution but a different random
     stream. On the reviewer's example (10 items x 20 raters, 5
     categories, 2 raters per person, n = 1000, every person a group) the
     suffix array is 13 MB instead of about 1.3 GB per draw.
   - **Rest-score distributions only where needed.** After that fix,
     profiling showed most time in the k full convolutions of
     `.cond_expected_probs()`. It now fills the rest-score distribution
     only at t = r - x for the group's totals, one inner product per t for
     all items at once. Exact, tested against full convolution. Per draw:
     reviewer example 2.0 s to 1.0 s, null-run design (n = 600, 10 of 40
     cells, halo groups) 0.3 s to 0.2 s. A vectorised `.cond_convolve()`
     (rowsum, then Toeplitz product) gave nothing, since the cost was the
     number of calls, and was reverted.
   - **Q3** (findings 2, 7): `q3_post()` splits rows by pair once instead of
     scanning the table per pair, and the theta branch of `q3_statistic()`
     builds its output in one tibble. With `facet_var`, cell pairs that no
     person shares are dropped (they had NA in every draw). `q3_post()$hdi`
     now has `pair_type` too.
   - **Hurdle families** (finding 3): `facet_parameters()` matches the
     ordinal family names exactly, so hurdle_acat stops with an error.
   - **Lazy contrasts** (findings 4, 8): `.item_levels()` attaches only a
     contrast specification (stored attribute, ordered flag), and
     `.item_contrasts()` builds the matrix in the k - 1 fixed-effect branch
     of `.item_effect_draws()`. Random-item and threshold models never
     evaluate contrasts. `.resolve_facet()` no longer computes the facet
     contrasts, and drops `n_facet`.
   - **Output type and order** (finding 5): the cell, item and facet columns
     of all facet output are factors in brms order (`levels(factor(x))`),
     from both methods. brms order is alphabetical for character data
     (R1, R10, R2), as for items elsewhere in the package.
   - **Label collisions** (finding 6): `.resolve_facet()` stops when two
     item x facet combinations get the same "item:facet" label. Residual
     PCA units are built from integer codes, so they cannot collide.
   - **Not changed** (finding 9): the person x cell uniqueness check runs
     twice (wrapper and recursive call). It is a paste and duplicated() over
     the rows, well under a second for a million rows, small next to
     posterior_epred. Kept for the clearer error message of the wrapper.
   - **PCA message** (finding 10): the incomplete-rows message now counts
     "id x rater combinations" in the facet path, via a `unit_label`
     attribute on the internal unit column.
5c. Second code review (2026-11-16, `dev/code-review/code-review-2026-11-16.md`).
   All nine findings verified. Fixed 2, 3, 4, 5, 7, with tests. The full
   suite passes.
   - Grouping factors (2): `.group_factors()` expands `id/rater` to `id` and
     `id:rater` and reads the factor of `gr(rater, ...)`. Halo detection
     and the connectivity check now see these terms.
   - Fixed effects (4, 5): matched by exact names, built as brms does
     (`.brms_rename()`, a copy of the replacement table in brms:::rename(),
     brms 2.23). No prefix matching, so `rater_exp` next to `rater` is
     ignored. If k or k - 1 coefficients carry the prefix but match
     neither naming, the contrasts at fitting differed (options(contrasts))
     and the call stops. brms does not store option contrasts, so this
     cannot be recovered, only caught. Also a 0.4.0 bug, in NEWS.
   - Halo message (7): names the term as brms expands it (`halo_term`). The
     infit caveat comes from `.infit_halo_note(by)`: 2% of cells, 1-2% of
     facet levels, 1% of item levels (null run, step 5).
   - `item_restscore_post()` (3): facet columns are those between the item
     column and gamma_obs.
   - Open (user decision): 1, a minimum number of shared persons for Q3
     cell pairs, and 6, building the compact ESF weights without the full
     n_persons x k x C arrays. 8 and 9 not changed.
   - Finding 1 study (2026-11-16):
     `dev/q3_min_shared.R` + `_summary.R`. Known-parameter shortcut (no
     brms, replicated data from .cond_q3_draw() with the true parameters),
     8 raters x 4 items, 2 raters per person, rater-pair overlap 2-30
     persons, so 240 of 496 cell pairs share fewer than 10 persons. null
     (200 reps), halo SD 0.6 and 1.0 for R1 only (100 each), 400 draws.
     Minimum m in 0, 3, 5, 8, 10, 15, 20 for the max-T family. Outcomes:
     FWER, power for R1's same-rater pairs, per-pair calibration by
     shared-person bin, share of max-T maxima from small pairs.
     RUN 2026-11-16 (4 min). Results:
     - Per pair under null, one-sided mid-p < .05: .049 (3-4 shared),
       .047 (5-9), .048 (10-19), .048 (20-49), .050 (50+). With 2 shared
       persons .001: Q3 = +/-1 and P(rep = obs) is usually large, so these
       pairs almost never flag. 4% of them have undefined Q3. (The smoke
       note that mid-p >= .25 for these pairs was wrong: it can be small
       when the replicated Q3 is almost always -1.)
     - Max-T: FWER .05 for every minimum 0-15, .06 at 20 (MC SE .015).
       Power for R1's 6 same-rater pairs: halo SD 0.6 .03, SD 1.0 .14,
       the same for minimum 0-15, .17 at 20 (a smaller family, not the
       small pairs). 95th percentile of the max-T reference 3.40 with or
       without the small pairs. Pairs with < 10 shared persons give the
       maximum in 10% of draws, about their share of the family, so they
       do not dominate it.
     - Decision: no minimum. The reviewer's mechanism does not show. The
       low power of pair-level q3_p_adj for a one-rater halo is the known
       weakness that the pooled halo test (6b) addresses.
     - Limitation: known parameters, so no posterior variation in the
       observed Q3.
   - Finding 6 (fixed 2026-11-16, user approved). .cond_item_weights() is
     used by all conditional statistics (restscore, infit, Q3) for every
     model, so the aim was output identical to the old code. New
     `.cond_weight_groups()`: the hash is the all-missing value plus each
     observation's difference from it, summed per person in item order
     (rowsum), and each group's weights come from its first person's
     observation rows. No n_persons x k x C arrays. .cond_draw_setup()
     uses it and stops with an internal error if a group's items differ
     from the answered items. `.cond_item_weights()` is now a wrapper that
     expands to the full array, used only in tests.
     Check against the old file (scratchpad copy, sourced into an
     environment): identical() for the .cond_draw_setup() result and for
     seeded .cond_q3_draw(), .cond_infit_draw() and .cond_restscore_draw()
     on complete PCM, complete dichotomous, 50% planned missing, DIF,
     sparse many-facet, person-specific weights (halo-like) and a large
     case (1000 persons, 600 cells, 40 per person). All identical.
     Peak memory of one setup in the large case 270 -> 205 MB. The rest
     is the compact suffix tables and the result (59 MB). Time the same.
     Full test suite passes.
5d. Third code review (2026-11-16,
   `dev/code-review/code-review-2026-11-16-third-round.md`). All eight
   findings verified and fixed, with tests (6 after the user's yes). The
   full suite passes.
   - Session contrasts (1): with no contrasts stored on the factor,
     `.item_contrast_candidates()` tries R's default (treatment, or poly
     for ordered factors) and then the current `options(contrasts)` if it
     differs. One match is used. Two matches (labels 0..3 give the same
     names under treatment and contr.sum) stop with an error that asks to
     set the option used at fitting. Treatment-coded models read under a
     sum option work again, as in 0.4.0.
   - `cell` as item or facet name (2): rejected in `.resolve_facet()`,
     since the output has a `cell` column.
   - Silent zero item effects (3): `.facet_draws()` with thresholds by
     facet stops when no item effects are found but the item variable is
     in the formula (fixed or random). Zeros only when the model has no
     item term.
   - thres() regex (4): the facet name is escaped and must be followed by
     `,` or `)`, so `rater` does not match `thres(gr = rater.grp)`.
   - Restscore plot order (5): factor cells are plotted in level order,
     as in the table. Character labels are still sorted.
   - `.expected_rating()` (7): cumulative, sratio and cratio vectorised over
     draws, checked against `.compute_cat_probs()` for logit and probit.
   - Test-only helpers (8): `.cond_item_weights()` and
     `.factor_contrast_matrix()` moved to `tests/testthat/helper-internal.R`.
   - Q3 pairs (6): the facet branch of `q3_statistic()` built all
     k(k-1)/2 cell pairs for every draw before dropping unshared ones. Now
     `.q3_shared_pairs()` finds the pairs with shared persons from the
     data, and `.q3_conditional()` and `.q3_theta()` (the theta branch,
     moved out of `q3_statistic()` unchanged) get `shared_only = TRUE` and
     compute only those. `.cond_q3_draw()` gets a `keep` argument. The
     any-non-NA filter is kept. Output `identical()` to the old code for
     both methods on three mock designs (connected, unshared rater pairs,
     one rater per person), and for the standard call. The gain depends
     on the design: 1000 persons, 30 raters x 20 items, 2 raters each
     share 91% of the 179,700 pairs, 300 persons 53%, 200 persons with 50
     raters x 10 items 17%. With most pairs shared, output size is
     inherent, and `ndraws_use` is the way to bound it.

6. ~~Vignette~~ (done 2026-11-16, `vignettes/mfrm-rasch-analysis.Rmd.orig`,
   knitted in the container with the host fits, added to precompile.R and
   the pkgdown Articles menu). As built:
   - Structure: design and connectivity (crossprod of person x rater),
     models, contr_sum_zero vs contr.sum (prior-only fits: contr.sum gives
     SD 1.16, 1.15, 1.16, 2.05, contr_sum_zero about 1.0 each), severity +
     fair average + RMU caution, rater fit (by rater, by item, by cell,
     restscore), halo (Q3 by pair type with and without the halo term, the
     same-rater pairs plotted, PCA, loo), rater-specific scales, the
     dichotomous model, a truth reveal at the end, null-run summary.
   - The truth is revealed only at the end, so the diagnostics are read
     first, as a user would.
   - Q3 plot shows only the same-rater pairs. Sorted by absolute
     difference, the top pairs are mostly negative pairs that share only
     25 persons, which hides the halo.
   - Found while writing: the comment in fit-mfrm-vignette-models.R gave
     the contr.sum prior SD of the first three levels as s * sqrt(3/4).
     Correct is s (1.15). Fixed.
   - Halo message: infit_statistic() now adds that infit was conservative
     in simulations with the halo modelled (about 2% at .05). Same in the
     help page.
   - New bib entries: Linacre 1989, Eckes 2015, Myford & Wolfe 2003, Wang &
     Wilson 2005, Morris 2025. Details written from memory, check before
     release.
   - Citation check against the PDFs (2026-11-16, articles/MFRM/):
     Myford & Wolfe 2003 and Wang & Wilson 2005 bib details correct.
     Two passages were miscited and are reworded in .Rmd.orig and .Rmd:
     (a) Wang & Wilson model the person x rater effect (with a variance per
     rater) and call it severity variation, not halo. The halo definition
     now cites Myford & Wolfe 2004 (Part II, added to the bib).
     (b) Overfit as central tendency is in Part II, not Part I, and Part II
     shows central tendency can also give infit > 1. The unsupported "depend
     on the ratings of other raters" explanation was removed.
     Linacre 1989, Eckes 2015 and Morris 2025 are still unchecked (no PDFs).
6b. Pooled halo test (drafted 2026-11-16, NOT RUN).
   `dev/halo_test_calibration.R` + `_summary.R`. Motivation: in the null
   run a halo that was not modelled raised the mean same-rater Q3 clearly,
   but the pair-level max-T found it in only 4-64% of datasets.
   - Statistics, all from conditional `q3_statistic(facet_var =)` draws of
     a model without the halo term, one-sided mid-p (small = more
     same-rater dependence than replicated):
     `mean_same` (mean Q3 over same-rater pairs, obs vs rep per draw),
     `contrast` (same-rater mean minus mean of all other pairs; a halo
     pushes the other pairs down because residuals given the total sum to
     zero, so about twice the signal), per rater (mean over the rater's
     pairs, max-T p_adj over raters with `.ppp_maxt_upper()`), and the
     existing pair-level q3_p_adj for comparison.
   - Design: the null-run rating design (5 items, 8 raters, 2 per person),
     RM and RSM, n = 300 and 600, 200 reps, conditions null, misfit (R2
     discrimination 0.4, RSM also R5 thresholds x 1.8, no halo), halo03,
     halo06, halo_two (SD 0.6 for R3 and R7 only), halo06_model (halo term
     in the model). 4,800 fits, very roughly 15-20 h on 10 workers.
     **Quick version the same day (user: too extensive for now, keep the
     full one for later):** `dev/halo_test_quick.R`, n = 300 only,
     conditions null (200 reps), misfit, halo03 and halo_two (100 reps
     each). 1,000 fits, about 1.5-2 h. Results in
     dev/results/halo_test_quick*, summary with
     `RUN=halo_test_quick Rscript dev/halo_test_calibration_summary.R`.
     The full version, `dev/halo_test_calibration.R`, is unchanged. Answers calibration, specificity and pooled vs
     pair-level power. Leaves out the halo left after the halo term, and
     sample size.
   - Checked without fitting (container): the generator gives the intended
     rater behaviour (R5 middle categories 79% against 54-61%, R2
     correlation with the person mean .37 against .45-.50). On the vignette
     fits: model without halo term, mean_same and contrast p = 0, four of six
     raters flagged. Model with halo term, p = .31 and .32, none flagged.
   - Observation: the inconsistent vignette rater R6 has strongly negative
     same-rater Q3 (-0.17 against -0.09 replicated). An inconsistent rater
     may mask a halo in the pooled mean. The misfit condition checks the
     false-alarm side, not masking.
   - If the pooled test works, candidate package form: columns or a
     `halo` element in `q3_post()` output when `pair_type` is present, or a
     small `halo_statistic()`. Decide after the results.
   - **Quick run results (2026-11-16, 1,000 fits, 0 errors, Rhat < 1.05 in
     all; summary `RUN=halo_test_quick Rscript
     dev/halo_test_calibration_summary.R`, figures
     dev/results/halo_test_quick_ecdf.png and _power.png).** Rejection
     rates at .05, n = 300:

     | model | cond | mean_same | contrast | rater_any | pair_same |
     |---|---|---|---|---|---|
     | RM | null | .075 | .075 | .075 | .035 |
     | RM | misfit | .12 | .11 | .13 | .00 |
     | RM | halo03 | .30 | .29 | .20 | .02 |
     | RM | halo_two | .31 | .31 | .15 | .02 |
     | RSM | null | .04 | .04 | .025 | .015 |
     | RSM | misfit | .25 | .25 | .36 | .01 |
     | RSM | halo03 | .86 | .85 | .56 | .05 |
     | RSM | halo_two | .78 | .78 | .75 | .08 |

     - Null: RSM at or below .05. RM .075 (15 of 200, 95% CI .043-.12,
       P(X >= 15 | .05) = .08), not clearly above .05. The ECDF has no
       lower-tail excess, only a mid-range bump. Per-rater flag rates under
       the null are .000-.020 per rater.
     - mean_same and contrast are the same test: p correlate .998 and agree
       on p < .05 in 99-100% of replications. The contrast doubles the
       observed and the replicated difference alike. Keep only mean_same.
     - Power: pooling helps a lot. RSM .86 (SD 0.3, all raters) and .78
       (SD 0.6, two raters) against .05 and .08 for the pair-level max-T.
       RM .30 and .31 against .02. Dichotomous ratings carry little halo
       signal at n = 300.
     - halo_two, per rater: RSM halo raters R3 and R7 flagged in .46, their
       neighbours (R2, R4, R6, R8, which share persons with them) .06-.12,
       R1 and R5 .00-.01. RM .07-.08 for the halo raters.
     - **Specificity fails.** Rater misfit without a halo raises the pooled
       test (RM .12, RSM .25), and the per-rater flags often point at the
       wrong rater. The inconsistent rater R2 is not flagged (RM .02, RSM
       .00). Its neighbours R1 and R3 are (RM .05 and .04, RSM .14 and
       .13), with same-rater Q3 above replicated. A plausible reason: with
       R2's ratings noisy, the partner rater's ratings carry more of the
       person's location than the model expects, which the conditional
       residuals show as dependence within the partner (as for item fit,
       see memory "item fit is relative to LD"). The central-tendency rater
       R5 is flagged directly (.16 adjusted, p < .05 in .43), as Myford and
       Wolfe (2004) lead us to expect: a rater's compressed scale gives a
       person x rater residual pattern like a halo.
     - Conclusion: the pooled same-rater Q3 is a sensitive test of extra
       dependence within raters, but not of halo specifically. A flagged
       rater can have a halo, central tendency, or share persons with an
       inconsistent rater. Not yet decided whether to add it. Candidate
       follow-ups (not run, ask first): (a) does a model with
       thres(gr = rater) remove the R5 flag, (b) does rater infit
       (by = "facet") flag R2 in the same data, so that infit plus the
       pooled Q3 together separate the cases.
   - **Follow-up script written 2026-11-16 (user chose the fuller
     version):** `dev/halo_followup.R`, summary
     `dev/halo_followup_summary.R`. Same generator and seeds as the quick
     run, so the datasets are identical. Arm common: the quick run's model
     refitted on RM and RSM misfit and halo_two (400 fits), adds rater
     infit (two-sided max-T over raters) and checks that the Q3 results
     reproduce. Arm thres: RSM with `thres(gr = rater)` on null, misfit and
     halo_two (300 fits), with Q3, rater infit and the raters' threshold
     spread. 100 replications each, so the thres null has MC SE about
     .022 (can be extended, finished replications are skipped). Flags are
     tabulated by rater role: halo rater, inconsistent R2,
     central-tendency R5, and raters that share persons with each.
     Smoke test passed (14 fits). The thres fits were faster than the
     common RSM fits (0.75 times), so about 1.5 h on 10 workers. The
     common arm reproduced a fresh smoke run of the quick script exactly
     (the saved quick smoke results were older than the script's last
     edit, so they did not match). Run on the host 2026-11-16, 700 fits,
     no errors.
   - **Follow-up results (n = 300, 100 replications per condition, flag
     = p_adj < .05 per rater):**
     - Reproducibility: the common arm matched the quick run in all 800
       pooled p-values and all 3200 per-rater p_adj values.
     - Convergence: R-hat at most 1.05 in every fit. Common RSM fits had
       divergent transitions in 19-23% of replications, RM 5-6%, the
       thres fits none (and took 47 s against 79 s).
     - Calibration with `thres(gr = rater)` under the null: pooled
       mean_same rejects at .03, any Q3 rater flag .01, any infit rater
       flag .01.
     - Central tendency (R5, RSM misfit): with common thresholds, infit
       flags R5 as overfitting in .90 and Q3 flags it in .16 (all of
       these also infit-flagged). With rater thresholds both flags vanish
       (.00 and .01), and R5's threshold spread is estimated as 4.18
       (generating 4.32, others 2.2-2.3). So (a) holds: the thres model
       removes the R5 flag and describes the effect in the thresholds.
     - Inconsistency (R2): infit flags R2 as underfitting in .73 (RM) and
       1.00 (RSM, both models). Q3 never flags R2 (.00-.02) but flags its
       two partner raters at .13-.135 in RSM and .045 in RM, also in the
       thres model. So (b) holds, and the partner flags are Q3-only. R2's
       thresholds are narrowed in the thres model (1.68), as noise pushes
       ratings to the extreme categories.
     - Halo (R3, R7): Q3 flags each halo rater in .46 (RSM) and .075 (RM),
       infit nearly never (.00-.005). Raters sharing persons with a halo
       rater are Q3-flagged in .08-.09 (RSM). Pooled mean_same is .78
       (common) and .76 (thres). The halo raters' thresholds are slightly
       narrower (2.16).
     - Pooled mean_same in RSM misfit drops from .25 to .11 with the thres
       model. The rest is the R2 partner effect, so the pooled test is
       still not halo-specific.
     - Reading the flags together: infit overfit = central tendency (or
       fit rater thresholds), infit underfit = inconsistency, Q3 only =
       halo or a partner of an inconsistent rater. The two Q3-only cases
       are told apart by whether a partner rater is infit-underfit. Q3
       power for a single halo rater at n = 300 is moderate in RSM (.46)
       and low in RM.
     - Figure `dev/results/halo_followup_flags.png`.
   - **User documentation of the halo results (2026-11-16):** vignette
     `mfrm-rasch-analysis` (prose only, edited in both `.Rmd.orig` and the
     knitted `.Rmd`, no re-knit needed): a paragraph in the halo section
     that same-rater Q3 is not halo-specific, a sentence in the
     rater-specific scales section on calibration and central tendency,
     and a new subsection "Which rater effect a flag points to" at the end
     of "How well the checks work" with the design, a flag-rate table and
     the reading rules. Help pages: a paragraph in the many-facet section
     of `q3_statistic()` and of `infit_statistic()`. The pooled halo test
     is still not in the package, so the vignette describes the per-rater
     Q3 test as a simulation tool and points to averaging `q3_diff` by
     rater. Checked in the vignette fits (not in the vignette): with the
     halo term, R3 (central tendency, also R6's main partner) has the
     highest mean same-rater q3_diff (0.06), which drops to 0.02 with
     rater thresholds. R6's own is negative (-0.11).
   - **Fourth code review (2026-11-16):**
     `dev/code-review/code-review-2026-11-16-fourth-round.md`. Findings
     1-5 fixed: stop on unequal threshold counts with `thres(gr = )`
     (set `thres(K, gr = rater)`), `contr_sum_zero()` callable by name,
     hash sentinel for undefined weights, covariates named like the facet
     no longer trip the contrasts error, halo message for the theta
     method. 6-8 (efficiency, reuse, prefix matching in
     `posterior_to_prior()`/`item_parameters()`) left for later.
7. Later scope: facet map, interactions, agreement, calibration.

Each new or changed exported function: tests, NEWS (brief), `_pkgdown.yml`,
README check.

## Open questions

- Halo term: message wording, and whether the conditional results with the
  halo absorbed answer the question users ask. Null run (2026-11-15): with
  the halo modelled, infit is conservative (about .02 at .05), restscore
  and Q3 near nominal. The message should say so.
- Pooled halo statistic (see null run results under step 5).
- Facet-specific thresholds: verify the engine against iarm with a
  cell-as-item layout.
- Max-T adjustment with many cells.
- Rule for fixed vs random small facets: recommend fixed with
  `contr_sum_zero()` below about 5-6 elements? Decide from the vignette fits,
  not a fixed rule.
