# cran-comments.md

## R CMD check results

0 errors | 0 warnings | 0 note

## Test environments

* local: macOS 26.7.2, R 4.6.1 and Ubuntu 24.04.5 LTS (aarch64), R 4.6.1
* rhub::rhub_check(): ubuntu-next, linux-r-devel
* devtools::check_win_devel()

## Notes

* Update version 0.3.0 -> 0.4.0.

* The item fit and local dependence functions (`item_restscore_statistic()`,
  `infit_statistic()`, `q3_statistic()`) now condition on each person's
  total score. Their results change, and the previous behaviour is kept as
  `method = "theta"` where it applies. NEWS.md lists these under
  "Breaking changes".

* Several bug fixes, including `dif_statistic()` for polytomous models and
  rating scale models.

* Two new datasets, `phq9` and `phq9_nhanes`.

## Downstream dependencies

There are currently no downstream dependencies for this package.
