# easyRaschBayes code review, 2026-11-13

Review of the uncommitted MFRM work in 0.4.0.9001 (`facet_var` in the fit
functions, `R/facets.R`, `R/facet_parameters.R`, contrast-aware
`.item_levels()`). Run as `/code-review high` in the background. The
findings below are the reviewer's. Status added 2026-11-13 after fixing;
design notes are under step 5b in `dev/mfrm-design-plan.md`.

| # | Status |
|---|--------|
| 1 | Fixed. ESF over answered items, rest-score distributions only at needed totals. Measured: 13 MB instead of ~1.3 GB per draw on the example. |
| 2 | Fixed. Split once in `q3_post()`, one tibble in the theta branch, unshared cell pairs dropped. |
| 3 | Fixed. Exact family names, hurdle stops. Test added. |
| 4 | Fixed. Contrasts built lazily via `.item_contrasts()`. Test with a 100-level ordered random item. |
| 5 | Fixed. Factors in brms order from both methods. Test added. |
| 6 | Fixed. Error on colliding cell labels, integer-coded PCA units. Test added. |
| 7 | Fixed. `pair_type` in `$hdi`, NEWS corrected. |
| 8 | Fixed. `contrasts`, `n_facet` dropped. `facet_levels` kept (now used for output order), `rank_deficiency` kept (computed for the warning anyway, read by tests). |
| 9 | Not changed. The check is cheap next to posterior_epred, and the wrapper gives the clearer message. |
| 10 | Fixed. Message counts person x facet combinations. Test added. |

## High (feasibility on realistic rater designs)

1. **Conditional path builds the ESF over every item x facet cell.**
   `R/conditional_score.R:267`. Each cell is treated as an item, including
   cells a person never had, so ESF tables are sized to all cells for every
   weight group. Example: 10 items x 20 raters, 5 categories, 2 raters per
   person gives k = 200 and M1 = 801. About 190 rater pairs means
   `suffix_all` is 190 x 201 x 801 doubles, about 245 MB per draw. With a
   `(1 | rater:id)` halo term every person is their own group (U = 1000),
   about 1.3 GB per draw, and `.cond_esf`/`.cond_expected_probs` loop over
   all 200 cells per person. `infit_statistic()`,
   `item_restscore_statistic()` and `q3_statistic()` with `facet_var` become
   infeasible. Suggested fix: build the ESF over each group's answered cells
   only.
2. **`q3_post()` is quadratic in the number of cell pairs.**
   `R/post_helpers.R:575`. `q3_statistic(facet_var =)` gives about
   cells^2/2 pairs, and `rep_mat`/`obs_mean` scan the whole long table once
   per pair. 40 cells give 780 pairs, with 4000 draws 3.1M rows and about
   5e9 comparisons. The theta branch also builds 780 tibbles. Most pairs are
   "neither" pairs that no person shares (all NA). Suggested fix: one
   `pivot_wider()` or a split by `item_pair`, and consider dropping
   unshared pairs.

## Medium (wrong or inconsistent output)

3. **`facet_parameters()` accepts hurdle_acat models.**
   `R/facet_parameters.R:155`. `grepl("acat|cumul|sratio|cratio")` also
   matches the package's own `hurdle_acat`, which some brms versions report
   by name (see `R/brms_hpcm.R:2358`). `.expected_rating()` then uses only
   the acat thresholds, so `fair_avg` silently ignores the hurdle
   probability of a zero. Suggested fix: explicit `stop()` for hurdle
   families.
4. **`.item_levels()` now builds contrasts for every caller.**
   `R/item_draws.R:34`. This affects item_parameters, plots, hpcm, cicc, RMU
   and the facet code, including random-item and `thres(gr = item)` models
   that never use contrasts. An ordered item factor with 96 or more levels
   in a `(1 | item)` model now errors in `contr.poly()`. A stored contrasts
   matrix whose row count differs from the levels stops in the new `nrow`
   check. Suggested fix: compute contrasts lazily, only in the k-1
   fixed-effect branch of `.item_effect_draws()`.
5. **Conditional and theta methods return differently typed and ordered
   `by = "item"`/`"facet"` output.** `R/facets.R:332`. Conditional gives
   character labels sorted alphabetically (rater 1, 10, 11, 12, 2, ...),
   theta gives the original column type in its own order. Binding or joining
   the two fails or misaligns, and `by = "item"` no longer follows the item
   order used elsewhere.
6. **Cell labels can collide.** `R/facets.R:207`.
   `paste(item, facet, sep = ":")` has no escaping, so item `a:b` with rater
   `c` and item `a` with rater `b:c` both give `a:b:c`. `grid_lab` then has
   duplicated levels and `factor()` errors. In `plot_residual_pca()` the
   same pattern for person x facet units can merge two units, which gives a
   spurious one-response-per-cell error or mixed residual rows.
7. **`q3_post()` adds `pair_type` to the summary table but not to `hdi`.**
   `R/post_helpers.R:585`. `infit_post()` adds the extra columns to both.
   NEWS says `q3_post()` keeps the item and facet columns, but it keeps
   neither. Filtering `$hdi` by `pair_type == "same rater"` fails.

## Low (cleanup)

8. **Unused fields in `.resolve_facet()`.** `R/facets.R:216`. `contrasts`,
   `facet_levels`, `n_facet` and `rank_deficiency` are returned but never
   read. Computing the contrasts calls `.factor_contrast_matrix()`, which
   can `stop()`, so a fit function with `facet_var` can fail on a value it
   never uses. Drop them, or compute them only in `facet_parameters()`.
9. **Duplicate validation.** `R/facets.R:276`. `.facet_setup()` runs
   `.check_unique_person_item()` on person x cell, then the recursive public
   call runs the same check again, along with response and variable
   validation. Costly on millions of rows. Suggested fix: an internal entry
   point that skips checks already done.
10. **Misleading message in the `plot_residual_pca()` facet path.**
    `R/brms_globalfit.R:233`. `fr` is computed only for the connectivity
    warning. The inner call's incomplete-rows message counts person x rater
    units as "persons" ("uses the 180 persons who answered every item"
    when 180 is person x rater rows).
