# Pre-compile the vignette locally.
# Run from the package root directory:
#   Rscript vignettes/precompile.R

library(knitr)
devtools::install(".", upgrade = "never")  # Install from local source
knit(
  input  = "vignettes/pcm-rasch-analysis.Rmd.orig",
  output = "vignettes/pcm-rasch-analysis.Rmd"
)
