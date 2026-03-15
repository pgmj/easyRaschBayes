# Pre-compile the vignette locally.
# Run from the package root directory:
#   Rscript vignettes/precompile.R

library(knitr)
devtools::install(".", upgrade = "never")  # Install from local source
# Knit from within vignettes/ so fig.path = "figures/pcm-" resolves correctly
withr::with_dir("vignettes", {
  knit(
    input  = "pcm-rasch-analysis.Rmd.orig",
    output = "pcm-rasch-analysis.Rmd"
  )
})
