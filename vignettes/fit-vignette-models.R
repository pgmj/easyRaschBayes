# Fit the models used in the vignettes and save them to vignettes/fits/:
#   fit_pcm_phq9.rds    pcm-rasch-analysis.Rmd.orig
#   dif_gender_phq9.rds pcm-rasch-analysis.Rmd.orig (dif_statistic() result)
#   fit_pcm_nhanes.rds  hpcm-rasch-analysis.Rmd.orig
#   fit_hpcm.rds        hpcm-rasch-analysis.Rmd.orig
# Run on the host from the package root:
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript vignettes/fit-vignette-models.R
# Set FIT=phq9 (or FIT=nhanes) to fit only one vignette's models, or
# FIT=dif to refit only the DIF model from the saved fit_pcm_phq9.rds.
# Requires the patched brms branch for hurdle_acat():
#   remotes::install_github("rpsychologist/brms@fix-custom-ordinal-grouped-thres")
# Then knit the vignette with vignettes/precompile.R.

library(brms)
library(dplyr)
library(tidyr)
library(tibble)
devtools::load_all(".")

which_fits <- Sys.getenv("FIT", "all")
dir.create("vignettes/fits", showWarnings = FALSE)
fits <- list()

if (which_fits %in% c("all", "phq9")) {
  df_pcm <- phq9 %>%
    rownames_to_column("id") %>%
    select(id, gender, starts_with("q")) %>%
    pivot_longer(starts_with("q"), names_to = "item",
                 values_to = "response") %>%
    mutate(response = response + 1)
  t0 <- Sys.time()
  fits$pcm_phq9 <- brm(
    response | thres(gr = item) ~ 1 + (1 | id),
    data   = df_pcm,
    family = acat,
    prior  = prior("normal(0, 3)", class = "Intercept") +
      prior("normal(0, 3)", class = "sd", group = "id"),
    chains = 4, cores = 4, iter = 2000, seed = 2024
  )
  saveRDS(fits$pcm_phq9, "vignettes/fits/fit_pcm_phq9.rds")
  message("PCM (phq9) done in ", format(Sys.time() - t0))
}

if (which_fits %in% c("all", "phq9", "dif")) {
  df_dif <- phq9 %>%
    rownames_to_column("id") %>%
    select(id, gender, starts_with("q")) %>%
    pivot_longer(starts_with("q"), names_to = "item",
                 values_to = "response") %>%
    mutate(response = response + 1) %>%
    filter(!is.na(gender))
  fit_pcm_phq9 <- if (!is.null(fits$pcm_phq9)) fits$pcm_phq9 else
    readRDS("vignettes/fits/fit_pcm_phq9.rds")
  t0 <- Sys.time()
  dif <- dif_statistic(fit_pcm_phq9, group_var = gender, data = df_dif,
                       rope = 0, cores = 4, refresh = 0, seed = 2024)
  saveRDS(dif, "vignettes/fits/dif_gender_phq9.rds")
  fits$dif_phq9 <- dif$dif_model
  message("DIF model (phq9) done in ", format(Sys.time() - t0))
}

if (which_fits %in% c("all", "nhanes")) {
df_hpcm <- phq9_nhanes[, 1:9] %>%
  rownames_to_column("id") %>%
  pivot_longer(!id, names_to = "item", values_to = "response") %>%
  mutate(response_pcm = response + 1)

t0 <- Sys.time()
fit_pcm <- brm(
  response_pcm | thres(gr = item) ~ 1 + (1 | id),
  data   = df_hpcm,
  family = acat,
  prior  = prior("normal(0, 3)", class = "Intercept") +
    prior("normal(0, 3)", class = "sd", group = "id"),
  chains = 4, cores = 4, iter = 2000, seed = 2024
)
saveRDS(fit_pcm, "vignettes/fits/fit_pcm_nhanes.rds")
message("PCM done in ", format(Sys.time() - t0))

prior_hpcm <- prior("normal(0, 3)", class = "Intercept") +
  prior("normal(0, 3)", class = "b", dpar = "hu") +
  prior("normal(0, 3)", class = "sd", group = "id") +
  prior("normal(0, 3)", class = "sd", group = "id", dpar = "hu") +
  prior("lkj(2)", class = "cor")

t0 <- Sys.time()
fit_hpcm <- brm(
  bf(
    response | thres(gr = item) ~ 1 + (1 |g| id),
    hu ~ 0 + factor(item) + (1 |g| id)
  ),
  data     = df_hpcm,
  family   = hurdle_acat(),
  stanvars = hurdle_acat_stanvars(),
  prior    = prior_hpcm,
  chains   = 4, cores = 4, iter = 2000, seed = 2024
)
saveRDS(fit_hpcm, "vignettes/fits/fit_hpcm.rds")
message("hPCM done in ", format(Sys.time() - t0))
fits$pcm_nhanes <- fit_pcm
fits$hpcm <- fit_hpcm
}

# Convergence summary
for (f in fits) {
  print(summary(f)$fixed[, c("Rhat", "Bulk_ESS")])
  message("max Rhat: ", round(max(rhat(f), na.rm = TRUE), 3),
          ", divergences: ", sum(nuts_params(f)$Value[nuts_params(f)$Parameter == "divergent__"]))
}
