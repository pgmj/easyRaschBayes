# cran-comments.md

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: macOS 26.3.1(a), R 4.5.3 - success
* check_win_devel: 1 NOTE (resolved)
* R-hub: linux, R-devel - success
* R-hub: windows, R-devel - failed
* R-hub: macOS-arm64, R-devel - failed

## Notes

Can't make sense of the error message from R-hub for windows/macOS-arm64:

> ** installing vignettes
> ** testing if installed package can be loaded from temporary location
> ** testing if installed package can be loaded from final location
> ** testing if installed package keeps a record of temporary installation path
> * SHA256 sums
> * creating tarball
> Error in if (custom.bin) { : argument is of length zero

* Update version 0.1.0 -> 0.2.0 

* Improved performance in 4 of the main functions and harmonized the output to be more 
  similar across functions.
* Added 3 post-processing functions
* Added 4 new analysis functions
* Added 3 new pre-analysis plotting functions

* Notes since prior submission:

* The package depends on 'brms', which requires a C++ toolchain for
  Stan model compilation. All examples and the vignette use
  `\donttest{}` or pre-computed results to avoid long-running MCMC
  sampling during checks.

* The vignette is designed to display code without evaluation on CRAN
  (model fitting results are loaded from a pre-saved file that is not
  included in the package tarball). This is noted in the vignette text.

## Downstream dependencies

There are currently no downstream dependencies for this package.
