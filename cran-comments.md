# cran-comments.md

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: macOS 26.3, R 4.5.2
* R-hub: linux, R-devel
* R-hub: windows, R-devel
* R-hub: macOS-arm64, R-devel
* R-hub: m1-san, R-devel

## Notes

* This is a resubmission. In this version I have:

* Fixed unexecutable code in examples for `plot_residual_pca` 
  (missing comma) and `RMUreliability` (uncommented line).
* Replaced all `\dontrun{}` with `\donttest{}` in examples, 
  since they can be executed but take longer than 5 seconds 
  due to MCMC sampling. I have reduced the sampling in the 
  examples to reduce their processing time.

* Notes since prior submission:

* "DIF", "Rasch", and "infit" are established psychometric terms,
  not misspellings.

* The package depends on 'brms', which requires a C++ toolchain for
  Stan model compilation. All examples and the vignette use
  `\donttest{}` or pre-computed results to avoid long-running MCMC
  sampling during checks.

* The vignette is designed to display code without evaluation on CRAN
  (model fitting results are loaded from a pre-saved file that is not
  included in the package tarball). This is noted in the vignette text.

## Downstream dependencies

There are currently no downstream dependencies for this package.
