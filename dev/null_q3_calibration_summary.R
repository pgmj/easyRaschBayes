# Summary of dev/null_q3_calibration.R
#
# Run from the package root:
#   Rscript dev/null_q3_calibration_summary.R
# With SMOKE=1 it reads the smoke-test results instead.
#
# Under the null, a calibrated p-value is uniform. Reported per condition
# and method:
#   mean, sd     uniform has 0.5 and 0.289
#   rej05_one    per-pair rejection, ppp > .95 (local dependence)
#   rej05_two    per-pair two-sided rejection at .05, 2 * min(p, 1 - p)
#   fwer05_one   proportion of replications with any pair ppp > .95
# Methods: theta_ppp (q3_statistic()), cond_ppp (conditional prototype),
# cond_ppp_adj (one-sided max-T over pairs, stored as 1 - adjusted p;
# fwer05_one should be near .05, mean/sd/rej05_two are not reported).
# Pairs: 45 (dichotomous, 10 items) or 28 (PCM, 8 items).
# MC standard error at 200 replications: about 0.015 for an FWER of .05,
# and much smaller for per-pair rates (pairs are pooled, but correlated).

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
smoke <- Sys.getenv("SMOKE") == "1"
name <- if (smoke) "null_q3_smoke" else "null_q3"
res <- readRDS(file.path("dev", "results", paste0(name, ".rds")))

cat("Replications per condition:\n")
print(res |> distinct(cond, rep) |> count(cond))
cat("\nConvergence and timing (median seconds):\n")
print(res |> distinct(cond, rep, rhat_max, n_div, fit_sec, theta_sec, cond_sec) |>
        group_by(cond) |>
        summarise(rhat_gt_101 = mean(rhat_max > 1.01), any_div = mean(n_div > 0),
                  fit_sec = median(fit_sec), theta_sec = median(theta_sec),
                  cond_sec = median(cond_sec), .groups = "drop"))

cat("\nMean Q3 over pairs (observed and replicated) and SD across draws:\n")
print(res |> group_by(cond) |>
        summarise(theta_obs = mean(theta_q3_obs), theta_rep = mean(theta_q3_rep),
                  cond_obs = mean(cond_q3_obs), cond_rep = mean(cond_q3_rep),
                  theta_sd_obs = mean(theta_sd_obs), theta_sd_rep = mean(theta_sd_rep),
                  cond_sd_obs = mean(cond_sd_obs), cond_sd_rep = mean(cond_sd_rep),
                  .groups = "drop") |>
        as.data.frame(), digits = 3)

long <- res |>
  select(cond, model_type, n, rep, item_pair, theta_ppp, cond_ppp, cond_ppp_adj) |>
  pivot_longer(c(theta_ppp, cond_ppp, cond_ppp_adj),
               names_to = "method", values_to = "p") |>
  filter(!is.na(p))

tab <- long |>
  group_by(cond, method) |>
  summarise(mean = mean(p), sd = sd(p),
            rej05_one = mean(p > .95),
            rej05_two = mean(pmin(1, 2 * pmin(p, 1 - p)) < .05),
            .groups = "drop") |>
  # the max-T column is 1 - adjusted p, not uniform under the null: only
  # rej05_one and fwer05_one are meaningful for it
  mutate(across(c(mean, sd, rej05_two),
                ~ if_else(method == "cond_ppp_adj", NA_real_, .x))) |>
  left_join(long |> group_by(cond, method, rep) |>
              summarise(any05 = any(p > .95), .groups = "drop") |>
              group_by(cond, method) |>
              summarise(fwer05_one = mean(any05), .groups = "drop"),
            by = c("cond", "method"))
cat("\nCalibration:\n")
print(as.data.frame(tab), digits = 3)
saveRDS(tab, file.path("dev", "results", paste0(name, "_table.rds")))

# ECDF minus the uniform CDF for the per-pair ppp. A calibrated p-value
# stays near 0. A conservative one bends below 0 at small p and above 0 at
# large p.
grid <- seq(0, 1, by = 0.01)
ecdf_df <- long |>
  filter(method != "cond_ppp_adj") |>
  group_by(cond, method) |>
  reframe(p_grid = grid, diff = ecdf(p)(grid) - grid) |>
  mutate(method = recode(method,
                         theta_ppp = "Theta-based (q3_statistic)",
                         cond_ppp  = "Conditional on total score"))

p <- ggplot(ecdf_df, aes(p_grid, diff, colour = method, linetype = method)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ cond, ncol = 2) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = "ppp", y = "ECDF minus uniform CDF", colour = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path("dev", "results", paste0(name, "_ecdf.png")), p,
       width = 8, height = 7, dpi = 150)
cat("\nSaved", paste0(name, "_ecdf.png"), "\n")
