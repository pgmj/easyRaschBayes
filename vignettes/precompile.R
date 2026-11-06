# Pre-compile the vignettes locally.
# Run from the package root directory:
#   Rscript vignettes/precompile.R
# Both vignettes need the fits from vignettes/fit-vignette-models.R.

library(knitr)
#devtools::install(".", upgrade = "never")  # Install from local source
# Knit from within vignettes/ so fig.path = "figures/..." resolves correctly
withr::with_dir("vignettes", {
  knit(
    input  = "pcm-rasch-analysis.Rmd.orig",
    output = "pcm-rasch-analysis.Rmd"
  )
  if (file.exists("fits/fit_hpcm.rds")) {
    knit(
      input  = "hpcm-rasch-analysis.Rmd.orig",
      output = "hpcm-rasch-analysis.Rmd"
    )
  }
})
