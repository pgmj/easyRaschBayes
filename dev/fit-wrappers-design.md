# Fitting wrappers for supported Rasch models (design)

Status: design only, 2026-11-06. Nothing implemented. Option B from the
discussion of rstanarm-style precompiled models: keep brms, but give users
one function that builds the supported models correctly. Options A
(precompiled brms Stan code inserted into an empty brmsfit) and C (upstream
brms fix for custom ordinal families) stay open. B is the interface A would
later sit behind.

## Goals

- Users give wide data (one column per item, scored from 0, as in
  easyRasch2 and the `plot_*` data functions) and get a `brmsfit` that every
  function in the package accepts without further arguments.
- Only supported models can be built: PCM, rating scale model (RSM),
  dichotomous Rasch model (RM), hurdle PCM (hPCM).
- The usual mistakes are handled inside: wide-to-long reshaping, the +1
  coding for `acat`, item and person variable names, priors, the hurdle
  family's stanvars.
- The result is a plain `brmsfit`, so `loo()`, `pp_check()`, `update()` and
  all existing functions keep working. Nothing is lost compared with
  calling `brm()` directly.

Non-goals: precompiled models, new Stan code, MFRM (a later `fit_mfrm()` can
follow the same pattern, see dev/mfrm-design-plan.md).

## API

One exported function, plus one exported helper:

```r
fit_rasch(
  data,
  model     = c("auto", "pcm", "rsm", "rm", "hpcm"),
  items     = NULL,     # item columns; NULL = all except id and keep
  id        = NULL,     # person id column; NULL = row numbers
  keep      = NULL,     # person-level columns to carry along (e.g. "gender")
  max_score = NULL,     # per-item maximum score; NULL = observed maximum
  prior     = NULL,     # brmsprior; NULL = default_prior_rasch(model)
  chains = 4, iter = 2000, cores = getOption("mc.cores", 1L),
  ...                   # passed to brms::brm(): seed, file, backend, control, ...
)

default_prior_rasch(model = c("pcm", "rsm", "rm", "hpcm"))
```

`model = "auto"` picks `"rm"` when every item is 0/1 and `"pcm"` otherwise.
It never picks `"rsm"` or `"hpcm"`, which are modelling choices, not data
properties.

Why one function and not `fit_pcm()`, `fit_rsm()`, ...: a single entry point
with `model` is easier to document and to extend, and `"auto"` covers the
common case. Thin aliases can be added later if wanted. (Open question 1.)

### Models built

Long data columns are always `id`, `item`, `response`, which are the
package's defaults, so later calls need no `item_var` or `person_var`.

| `model` | Response coding | brms formula | family |
|---|---|---|---|
| `"pcm"` | x + 1 | `response \| thres(nthres, gr = item) ~ 1 + (1 \| id)` | `acat` |
| `"rsm"` | x + 1 | `response ~ 1 + item + (1 \| id)` | `acat` |
| `"rm"` | x | `response ~ 0 + item + (1 \| id)` | `bernoulli` |
| `"hpcm"` | x | `bf(response \| thres(nthres, gr = item) ~ 1 + (1 \|g\| id), hu ~ 0 + factor(item) + (1 \|g\| id))` | `hurdle_acat()` + stanvars |

Item effects are fixed effects in RSM and RM. The vignettes and some
examples use `(1 | item)`, which shrinks item difficulties toward each other
(the null calibration run showed a small bias in expected gamma from this
for `dich_ri`). Fixed items are the Rasch convention and are what CML
estimates. All extractors already support `0 + item` and treatment-coded
`item`. (Open question 2.)

`thres(nthres, gr = item)`: brms sets the number of thresholds per item
from the observed maximum unless told otherwise. An item whose top category
nobody used then silently gets fewer categories. The wrapper passes the
number of thresholds per item as a data column, from `max_score` or the
observed maximum. **Verify** that brms accepts a per-group `nthres` column in
`thres()` for `acat` and for the custom hurdle family, and that
`.item_tau_draws()` reads the result correctly.

### Default priors (`default_prior_rasch()`)

As in the vignettes, all on the logit scale:

- thresholds: `normal(0, 3)` (class `Intercept`)
- fixed item effects (RSM, RM): `normal(0, 3)` (class `b`)
- person SD: `normal(0, 3)`, half-normal through the lower bound
- hPCM: also `normal(0, 3)` on `b` and `sd` for `dpar = "hu"`, and
  `lkj(2)` on `cor`

`prior` replaces the defaults entirely, so users can start from
`default_prior_rasch(model)` and change single rows. No silent merging.

### Input checks

Stop with a clear message when:

- an item column is not integer-valued or has negative values
- fewer than 2 items
- an item has only one observed category
- `model = "rm"` and a response is above 1
- `model = "hpcm"` and an item has no responses above 0, or only one
  positive category (no severity thresholds)
- `max_score` is below an observed score

Warn when:

- an item has a category with zero responses below its maximum (the
  threshold then rests on the prior)
- persons have no responses (dropped, with their number)
- `keep` columns vary within a person (they must be person-level)

Missing responses are dropped in the long data, as brms would.

### hPCM and the brms version

CRAN brms 2.23 emits invalid Stan code for this family. Before sampling,
call `brms::make_stancode()` (no compilation) and, if Stan's parser or the
compile fails, stop with the install line for the patched branch. Keep
the check in one internal function so it can be removed when brms is
fixed.

### What the fit carries

The returned object is a `brmsfit` with one extra list element:

```r
fit$easyRaschBayes <- list(
  model     = "pcm",
  items     = c("q1", ...),      # in data order
  offset    = 1L,                # added to responses
  max_score = c(q1 = 3, ...),
  data_long = <long data incl. keep columns>,
  call      = match.call()
)
```

`data_long` matters because brms keeps only formula variables in
`fit$data`, which is why `dif_statistic()` and `plot_icc()` need the data
passed again. With `data_long` stored, phase 2 can let them find group
variables there: `dif_statistic(fit, group_var = gender)` would work
without `data =`. **Check** that an extra element survives `update()`
(brms may rebuild the object) and does not upset `brms::` methods.

`items` in data order also fixes a small annoyance: brms orders items as
`levels(factor(item))`, so output for `q1..q10` comes out as q1, q10, q2.
The wrapper makes `item` a factor with levels in data order, and the
extractors already use `levels(factor(...))`, so output follows the data.
**Check** that brms keeps factor level order for `thres(gr = item)` and
for `0 + item`.

### Exported helper for long data? (Open question 3)

Users need long data for `dif_statistic(data =)` and `plot_icc(dif_data =)`
today. Either export `rasch_long(data, items, id, keep, offset)` (the
reshaping step on its own), or rely on `fit$easyRaschBayes$data_long` once
phase 2 is in. Exporting it costs little and helps users who fit with
`brm()` directly.

## Phases

1. `fit_rasch()`, `default_prior_rasch()`, internal `.rasch_spec()` that
   returns formula, family, stanvars, prior and long data without sampling.
   Tests for `.rasch_spec()` (no Stan needed): coding, formulas, priors,
   nthres per item, checks and warnings, `"auto"` choice, factor level order.
   One sampled smoke test per model with `skip_on_cran()` and a tiny
   `iter`, checking that `infit_statistic()`, `item_parameters()` and
   `plot_targeting()` run on the result. hPCM test skipped unless the
   installed brms can build the family.
2. `dif_statistic()` and `plot_icc()` look up group variables in
   `fit$easyRaschBayes$data_long` when `data`/`dif_data` is `NULL`.
3. Vignettes switch to `fit_rasch()`, showing the equivalent `brm()` call
   once for transparency. `fit-vignette-models.R` uses it too. README quick
   example uses it.

Docs each phase: NEWS (brief), `_pkgdown.yml` (new "Fitting" section),
README function table.

## Effort

Phase 1 about one session, phase 2 a small one, phase 3 part of the next
vignette re-render. No host simulation needed. The smoke tests need small
fits, which the container can run.

## Open questions

1. One `fit_rasch(model = )` or separate `fit_pcm()`, `fit_rsm()`,
   `fit_rm()`, `fit_hpcm()`?
2. Fixed item effects for RSM and RM (Rasch convention, proposed), or
   random as in the current examples? Fixed changes the dichotomous example
   in the README and docs.
3. Export `rasch_long()`?
4. Default sampling settings: brms defaults (4 chains, 2000 iterations), or
   more for the hurdle model?
5. Name: `fit_rasch()` fits the package's snake_case style. Anything closer
   to easyRasch2's `RM*` names would clash with that style.
