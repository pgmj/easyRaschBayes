# Code review 2026-11-16, fourth round (/code-review high)

Scope: uncommitted diff on top of 0.4.0 (0.4.0.9001, MFRM work after the
third-round fixes). Eight findings, none fixed yet. Findings are as the
reviewer reported them. They have not been checked by hand except where
noted.

## Correctness

1. **Rater location from a different threshold set.**
   `R/facet_parameters.R:281`. With `thres(gr = rater)`, a rater's
   location is the mean of its own thresholds. If brms estimates fewer
   thresholds for a rater (a severe rater never uses the top category),
   `vapply(tau_f, rowMeans, ...)` averages K-1 thresholds instead of K.
   That rater looks more lenient, and the centring shifts every other
   rater. `.expected_rating()` also caps that rater's fair average one
   category lower without a warning.

2. **`contr_sum_zero()` fails when called by name as a contrasts function.**
   `R/facets.R:67`. It has no `contrasts`/`sparse` arguments, so
   `options(contrasts = c("contr_sum_zero", "contr.poly"))` or
   `contrasts.arg = list(f = "contr_sum_zero")` stops with "unused argument
   (contrasts = contrasts)". The reviewer reproduced this. The new error
   message in `.item_effect_draws()` points users to exactly this route.
   `contr.sum()` takes `(n, contrasts = TRUE, sparse = FALSE)`.

3. **NaN weights hash like a missing item.** `R/conditional_score.R:135`.
   An answered item with undefined weights (NaN in `lo`, for example
   `prob[, 1]` underflowing to 0) gets the same sentinel (1e4) as a missing
   item. Two persons with different answered sets can share a group. The
   internal check tests only the representative person. If that is the
   person who did not answer, the other person's response on that item is
   dropped from E/V, and `.cond_sample()` stops with "Conditional sampling
   failed". Follow-up to third-round Finding 6.

4. **Prefix count in the new fallback error.** `R/item_draws.R:133`. A
   covariate whose name starts with the facet name (`b_rater_exp`) can make
   `n_prefix` equal k or k-1. `facet_parameters()` on
   `rating | thres(gr = rater) ~ 1 + item + rater_exp + (1 | id)` with 2
   raters then stops with "Could not match the fixed effects of 'rater'"
   instead of returning NULL. This is the ambiguity the NEWS entry calls
   fixed.

5. **No halo message for `method = "theta"`.** `R/brms_localdep.R:248`.
   The halo message, and the halo caveat in `infit_statistic()`, appear
   only for the conditional method. `posterior_epred()` includes the
   `(1 | rater:id)` effects, so theta residuals are also computed after the
   halo is absorbed. The docs say only the conditional method accounts for
   the halo.

## Efficiency

6. **Full k x k correlation matrix per draw in many-facet Q3.**
   `R/conditional_score.R:594`, and `.q3_theta()` at
   `R/brms_localdep.R:355`. Each draw runs `cor(..., "pairwise.complete.obs")`
   over all cell pairs and subsets to `keep` afterwards. 200 cells give
   about 40k pairs, most of them all NA. Compute only the kept pairs.

## Reuse

7. **`.q3_theta()` duplicates draw selection.** `R/brms_localdep.R:276`.
   It re-implements the `ndraws_use` check, warning and sampling instead of
   calling `.select_draw_ids()`. The copies already differ (NULL vs
   `seq_len(total)` for all draws).

8. **Prefix matching of item effects remains elsewhere.**
   `R/post_helpers.R:950` (`posterior_to_prior()`) and
   `R/item_parameters.R:249` still use `grep("^b_", item_name)`. With
   `y ~ 0 + item + item_order + (1 | id)`, `posterior_to_prior()` writes an
   item prior for `b_item_order`. Root fix: route all detection through the
   exact-name `.item_effect_draws()`.

## Status

Checked 2026-11-16. Findings 1-5 fixed the same day (user approved, stop
chosen for 1). 6-8 left for later. Full test suite passes, `check_man`
clean.

Fixes:

1. `.facet_draws()` stops when the facet levels have different numbers of
   thresholds, naming the counts and suggesting `thres(K, gr = rater)`.
   Help of `facet_parameters()` and the vignette's rater-scales section
   say to set the number. The vignette fit has all four categories for
   every rater, so it is unaffected.
2. `contr_sum_zero(n, contrasts = TRUE, sparse = FALSE)`. `contrasts =
   FALSE` gives the identity, `sparse = TRUE` stops. Tests by name via
   `contrasts<-` and `options(contrasts)`.
3. Rows with undefined weights (`is.na(lo[, 1])`) get their own hash
   sentinel 2e4. Checked against the old code: the collision needs a
   person whose every row is NaN (a person whose only answered item has
   prob[, 1] = 0, against a person with no answers). The old code gave
   both key 1, the new code separates them. Very unlikely in practice.
4. `.item_effect_draws()` and `.item_tau_draws()` take `other_vars`
   (names(model$data)). Coefficients of a longer variable name starting
   with the factor name are left out of the prefix count.
   `facet_parameters()` and `.model_tau_draws()` pass it. Other callers of
   `.item_tau_draws()` are unchanged (outside the MFRM work, as Finding 8).
5. `.facet_halo_message()` takes `method`, with a theta wording. `q3` and
   `infit` with `facet_var` show it for both methods. The infit caveat
   from the null run stays conditional only. Help text of both updated.

| # | Check | Verdict |
|---|---|---|
| 1 | `make_standata()` with `thres(gr = rater)`, R2 never uses category 4: `nthres` = 3, 2 | confirmed |
| 2 | `contrasts(x) <- "contr_sum_zero"` stops with "unused argument (contrasts = contrasts)"; `contr.sum` has `(n, contrasts, sparse)` | confirmed |
| 3 | code read: NaN in `lo` gets the missing-item sentinel. Needs `prob[, 1]` to underflow to 0 (log-odds beyond about 700), and then the reference subtraction also gives Inf - Inf, so the exact failure path was not traced | plausible, rare |
| 4 | code read: `.item_tau_draws()` (item_draws.R:152) calls `.item_effect_draws()` for the facet. With 2 raters and one `b_rater_*` covariate, `n_prefix` = 1 = k - 1 | confirmed, narrow |
| 5 | code read: message only when `method == "conditional"`. `posterior_epred()` includes `(1 \| rater:id)` by default | confirmed (message and docs) |
| 6 | efficiency, small for the vignette's 24 cells | valid, optional |
| 7 | `.q3_theta()` is the old theta branch moved unchanged (third round Finding 6), so the duplicate predates this diff | valid, optional |
| 8 | `posterior_to_prior()` and `item_parameters()` code predates 0.4.0.9001 | valid, outside the MFRM work |
