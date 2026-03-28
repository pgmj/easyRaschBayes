# cran-comments.md

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: macOS 26.3.1(a), R 4.5.3
* devtools::check_win_devel()
* devtools::check_mac_release()

## Notes

* Update version 0.1.0 -> 0.2.0 

* Improved performance in 4 of the main functions and harmonized the output 
  to be similar across functions.
* Added 3 post-processing functions
* Added 4 new analysis functions
* Added 3 new pre-analysis plotting functions

* Notes since prior submission:

* The package depends on 'brms', which requires a C++ toolchain for
  Stan model compilation. All examples that include `brmsfit` models and 
  the vignette use `\donttest{}` or pre-computed results to avoid 
  long-running MCMC sampling during checks.

* The vignette is designed to display code without evaluation on CRAN
  (model fitting results are loaded from a pre-saved file that is not
  included in the package tarball). This is noted in the vignette text.

## Downstream dependencies

There are currently no downstream dependencies for this package.
