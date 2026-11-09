# Conditional infit calibration: note for a possible future study

Status 2026-11-09: parked. Nothing decided. Discuss the design anew before
writing any code.

## Where things stand

- The conditional infit ppp (default in `infit_statistic()`) is conservative
  under the null. In `dev/null_restscore_calibration.R` the two-sided rejection
  at .05 was 0.5% to 2% of items. The theta method is worse.
- `dev/null_infit_adj_calibration.R` (1200 fits, 2026-11-09) tested a two-sided
  max-T `p_adj` for infit, as in `item_restscore_statistic()`. Pooled familywise
  error at .05 was .012 for conditional and .001 for theta. The unadjusted
  standardised per-item p (`p_marg`) gave .10 to .30. So the adjustment works,
  but it inherits the conservativeness. Not added to `infit_post()`.
- Cause, measured directly: across null datasets the posterior-mean observed
  infit varies less than the replicated draws say. The ratio SD across datasets
  / mean replicated SD is .65 to .96 for conditional (mean .82) and .39 to .69
  for theta. It is lowest for the easiest and hardest dichotomous items. Hence
  middle items are the most extreme item about 3 times as often as end items
  under the null (dichotomous argmax chi-square p < 1e-6; PCM p = .03 to .07).
  Results in `dev/results/null_infit_adj*.rds` and `null_infit_adj_ecdf.png`.
- The cause is double use of the data: the item parameters are fitted to the
  same responses.

## Candidate fixes (none chosen)

1. Hjort, Dahl and Steinbakk (2006, JASA 101, 1157-1174) calibration: simulate
   datasets from the fitted model, refit, recompute ppp, read the observed ppp
   against that distribution. Recording the per-item minimum in each refit gives
   a max-T reference in the same loop. Cost about 66 s per brms fit here, so
   B = 200 is about 20-25 min on 10 cores per analysis. A null check is a double
   loop (200 x 200 = 40,000 fits, about 3 days).
   Cheaper variants (my suggestions, untested):
   - refit with CML (psychotools) inside the loop
   - rescale the replicated reference by a per-item SD ratio from a small B
   - calibrate only the minimum, for the act-on-the-top-item workflow
2. CML parametric bootstrap reference, as in easyRasch2.
3. Item parameters not informed by the item's own responses.
4. Document only (current state).

The SD ratio above is a natural check for any fix: it should move to about 1
for every item.

## Suggested starting point (user, 2026-11-09)

Begin with a small detection power comparison between conditional infit and
item-restscore, before investing in calibration. Whether calibrating infit is
worth it depends on whether it detects anything restscore does not. Design to
be discussed when we return to this.

Outfit could be a third arm at little cost (`infit_statistic(outfit = TRUE)`
already computes it for both methods). Expectation, untested: the power
disadvantage from the user's frequentist infit vs outfit study (calibrated
under WY, never more powerful, worst off target) carries over, since it comes
from the statistic. Double use may make outfit even more conservative than
infit for off-target items.

## Related open questions

- Sequential error rate over several rounds of item removal.
- A family across statistics (infit, restscore, Q3) rather than within each.
