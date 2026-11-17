# easyRaschBayes code review, 2026-11-16

Second review of the uncommitted MFRM work in 0.4.0.9001, after the fixes
from `dev/code-review/code-review-2026-11-13.md`. Covers `R/facets.R`,
`R/facet_parameters.R`, `R/conditional_score.R`, `R/item_draws.R`,
`R/brms_fitstats.R`, `R/brms_localdep.R`, `R/brms_globalfit.R` and
`R/post_helpers.R`. Run as `/code-review high` in the background. The
findings below are the reviewer's and have not been verified. The reviewer
checked three things in R: a data-column attribute survives assignment,
`terms()` keeps `(1 | id/rater)` and `(1 | gr(rater, ...))` as single term
labels, and selecting the same tibble column twice gives duplicate columns
without an error.

| # | Status |
|---|--------|
| 1 | Studied, not changed (2026-11-16, `dev/q3_min_shared.R`). Small-n pairs neither flag spuriously (2 shared persons: .001 at .05) nor inflate the max-T reference (95th percentile 3.40 with and without them). A minimum changes neither FWER (.05) nor power. No minimum needed. |
| 2 | Verified (`terms()` keeps `1 \| id/rater` and `1 \| gr(rater, by = g)`). Fixed: `.group_factors()` expands nesting and reads gr(). Test added. |
| 3 | Verified. Fixed: extra columns are those between the item column and gamma_obs. Test added. |
| 4 | Verified. Fixed with 5: fixed effects matched by exact brms names (`.brms_rename()`, copy of brms:::rename()). Test added. |
| 5 | Verified. Fixed: coefficients that match neither the k nor the k-1 naming of the stored contrasts stop with an error that points to options(contrasts). Cannot be recovered, brms does not store option contrasts. Test added. |
| 6 | Fixed (2026-11-16). New `.cond_weight_groups()` hashes and builds weights from the observations only. Output identical to the old code on 7 designs. Peak 270 -> 205 MB at n = 1000, 600 cells. |
| 7 | Verified. Fixed: message uses the term as found in the formula (`halo_term`), the infit caveat is built by the caller from `by` (cells 2%, facet 1-2%, item 1%). Test added. |
| 8 | Verified, not changed. Works and the message test covers it. An internal argument would change the exported signature. |
| 9 | Verified, not changed. Repeated checks cost little next to the draws. |

## Can give wrong results

1. **Q3 keeps cell pairs with very few shared persons.**
   `R/brms_localdep.R:244`. The filter drops only pairs that are NA in every
   draw. With 2 shared persons q3 and q3_rep are +/-1 in every draw. In
   `q3_post()` these pairs get extreme standardised values that inflate the
   max-T reference distribution, so `q3_p_adj` loses power for every other
   pair, and they show as spurious same-rater or same-item dependence.
   Suggested fix: require a minimum number of shared persons.
2. **`.formula_roles()` misreads nested and `gr()` syntax.**
   `R/facets.R:130`. Grouping factors come from the raw term label, so
   `(1 | id/rater)` gives group `id/rater` and `(1 | gr(rater, ...))` gives
   `gr(rater, ...)`. With `y ~ 1 + item + (1 | id/rater)` (brms expands it
   to `(1|id) + (1|id:rater)`) halo is FALSE and the halo message never
   appears. With `(1 | gr(rater, by = g))` `roles$random` is FALSE, the
   connectivity check is skipped and `facet_parameters()` stops with
   "'rater' is not in the model".
3. **Many-facet output detected by the first column name.**
   `R/post_helpers.R:318`. `item_restscore_post()` treats any result whose
   first column is `cell` as many-facet and takes columns 2:3 as item and
   facet. A standard model whose item variable is named `cell` gets
   `extra_cols = c("gamma_obs", "gamma_exp")`, so the summary has duplicated
   `gamma_obs`/`gamma_exp` columns and only one copy is rounded.
   `infit_post()` uses a safer test (columns before `draw`).
4. **Fixed effects found by name prefix.**
   `R/item_draws.R:73`. `.item_effect_draws()` matches `^b_<name>(.+)$`,
   which `facet_parameters()` now also applies to the facet. The
   `.formula_roles()` comment already notes the ambiguity ("b_r would match
   b_rater"). A fixed facet `rater` plus a covariate `rater_exp`, or an item
   variable that is a prefix of the facet name, adds extra columns. The count
   then takes the wrong branch (k instead of k-1), so the call errors or,
   with non-plain labels, maps coefficients to the wrong elements without a
   warning.
5. **Contrasts taken from the session at extraction time.**
   `R/facets.R:107`. When the factor stores no contrasts,
   `.contrast_from_spec()` uses `getOption("contrasts")` now, not the option
   in force at fitting. A model fitted under
   `options(contrasts = c("contr.sum", "contr.poly"))` and extracted in a
   fresh session is decoded as treatment coding. With labels such as
   `item 1`, `.check_label_map()` skips its check, so `item_parameters()` and
   `facet_parameters()` return wrong locations with no error.

## Claims and messages

6. **Peak memory still scales with all cells.**
   `R/conditional_score.R:289`. The compact ESF layout is built after
   `.cond_item_weights()` has allocated the full n_persons x k x C array `L`,
   its copy `Lh` and the full U x k x C array `W`. A halo model (U =
   n_persons) with 20 items x 30 raters = 600 cells still allocates three
   n_persons x 600 x C arrays per draw. NEWS says the methods "use far less
   memory" in many-facet designs. Suggested fix: build `W_c` directly from
   the answered rows.
7. **Halo message wording and calibration figure.**
   `R/facets.R:304`. `.facet_halo_message()` always prints the term as
   `(1 | <facet>:<person>)`, so a user with `(1 | id:rater)` is told the
   model has `(1 | rater:id)`. The extra text is chosen by comparing
   `fun_name` strings, and for `infit_statistic(by = "facet")` it gives the
   per-cell figure ("2% of fitting cells") for a per-rater statistic (pooled
   rater infit was 1-2% of raters). Renaming the function or calling the
   helper from a wrapper silently drops the caveat.

## Code quality

8. **Message choice passed through a data-column attribute.**
   `R/brms_globalfit.R:241`. The facet branch sets a `unit_label` attribute
   on a data column and the inner function reads `attr(persons, ...)`. Any
   later step that rebuilds `persons` (`factor()`, `unique()`,
   `as.character()`, dplyr) drops it and the wrong message prints.
   Suggested fix: an internal argument such as `.unit_label`.
9. **Checks repeated within and across calls.**
   `R/facets.R:288`. `.facet_setup()` runs `.check_unique_person_item()` on
   person x cell, then the recursive call (or `.infit_conditional()`) runs it
   again. `plot_residual_pca()` also runs it twice. Every function repeats
   `.resolve_facet()` with its eigen-based rank check. A workflow of infit,
   restscore, Q3 and PCA runs the uniqueness check about 8 times and the
   (I+F) x P table, eigen decomposition and connectivity warning 4 times.
   Suggested fix: check once and pass a flag to the inner call.
