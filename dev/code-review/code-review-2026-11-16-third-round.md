# easyRaschBayes code review, 2026-11-16 (third round)

`/code-review high` of the uncommitted MFRM work in 0.4.0.9001, run after
the fixes logged in `dev/code-review/code-review-2026-11-16.md` (second round). Findings are the
reviewer's. The reviewer confirmed 1 and 3 in R. All were verified by
reading the code, and fixes are described in `dev/mfrm-design-plan.md`,
item 5d.

| # | Status |
|---|--------|
| 1 | Fixed: default coding first, then the session option, error if both match |
| 2 | Fixed: `cell` rejected as item or facet name |
| 3 | Fixed: error when the item term is in the formula but not found |
| 4 | Fixed: escaped name, followed by `,` or `)` |
| 5 | Fixed: factor level order in the plot |
| 6 | Fixed: only shared pairs computed, output unchanged |
| 7 | Fixed: vectorised over draws |
| 8 | Fixed: moved to `tests/testthat/helper-internal.R` |

## Can give wrong results

1. **Session contrasts used for models fitted with treatment coding.**
   `R/item_draws.R:99`. When the factor stores no contrasts, the extractor
   decodes with the reading session's `options(contrasts)`. A model fitted
   with the default treatment coding and read after
   `options(contrasts = c("contr.sum", "contr.poly"))` (or
   `afex::set_sum_contrasts()`) stops with a "refit" error. This worked in
   0.4.0. With item labels "0".."3" the names match under both codings, and
   the true effects 0, 1, 2, 3 come back as 1, 2, 3, -6 with no error
   (confirmed). Suggested fix: when no contrasts are stored, try treatment
   first and the session option second, use whichever names match, and stop
   only if neither or both match. NEWS currently says "Treatment-coded models
   are unchanged", which this contradicts.
2. **An item or facet variable named `cell`.** `R/facets.R:353`.
   `.label_cells()` always adds `cell`, so `item_var = cell` gives
   `cell...1`, `cell...2` (confirmed). `infit_post()` then drops the item
   and facet columns. Restscore and Q3 are affected too. Suggested fix:
   reject the name in `.resolve_facet()`, or use an internal column name.
3. **Silent zero item effects.** `R/item_draws.R:109`. The old count
   mismatch error is gone. When exact names fail and the prefix count is
   neither k nor k-1, `.item_effect_draws()` returns NULL silently, and the
   `thres(gr = rater)` branch of `.facet_draws()` (facet_parameters.R:272)
   replaces it with zeros. Fair averages then ignore item effects without a
   warning.
4. **Unescaped facet name in the thres() regex.** `R/facets.R:137`.
   `gr\s*=\s*rater\b` matches `thres(gr = rater.grp)`, because `.` is a
   word boundary, so `roles$thres` is TRUE and `facet_parameters()` stops
   with a misleading label error. Suggested fix: escape the name and
   require a following `,` or `)`.
5. **Restscore plot order differs from the table.** `R/post_helpers.R:338`.
   The plot sorts character cell labels alphabetically (crit10 before
   crit2), and the table uses the factor order.

## Efficiency

6. **Q3 builds all cell pairs before dropping unshared ones.**
   `R/brms_localdep.R:244`. 600 cells give 179,700 pairs, and with 4000
   draws that is about 719M rows before the filter. Suggested fix: find the
   shared pairs once from x_obs and compute only those.
7. **`.expected_rating()` loops over draws in R** for cumulative, sratio and
   cratio (`R/facet_parameters.R:345`). 4000 draws x 20 items x 30 raters
   makes 2.4M calls for the fair average. Suggested fix: vectorise over
   draws as the acat branch does.

## Cleanup

8. `.cond_item_weights()` (`R/conditional_score.R:163`) and
   `.factor_contrast_matrix()` (`R/facets.R:92`) are only called by tests.
   Move them to a test helper or test the real functions.
