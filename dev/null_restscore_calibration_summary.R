# Summary of dev/null_restscore_calibration.R (restscore and infit)
#
# Run from the package root:
#   Rscript dev/null_restscore_calibration_summary.R
# With SMOKE=1 it reads the smoke-test results instead.
#
# Under the null, a calibrated p-value is uniform. Reported per condition
# and method:
#   mean, sd         uniform has 0.5 and 0.289
#   rej05, rej10     per-item two-sided rejection rate at .05 and .10.
#                    For all methods except ppp_adj the two-sided p is
#                    2 * min(p, 1 - p).
#   fwer05           proportion of replications with any item rejected at .05
# Methods: ppp (restscore, conditional), ppp_old (restscore, 0.3.0
# theta-based), ppp_adj (restscore, max-T adjusted), infit_ppp_cond and
# infit_ppp_theta (infit, P(observed > replicated)).
# ppp_adj is already two-sided and adjusted, so its rej05 is per item and its
# fwer05 should be close to .05. The unadjusted methods are expected to have
# fwer05 well above .05 with 8 to 10 items.
# MC standard error at 200 replications: about 0.015 at a rate of .05.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
smoke <- Sys.getenv("SMOKE") == "1"
name <- if (smoke) "null_restscore_smoke" else "null_restscore"
res <- readRDS(file.path("dev", "results", paste0(name, ".rds")))

cat("Replications per condition:\n")
print(res |> distinct(cond, rep) |> count(cond))
cat("\nConvergence (max Rhat > 1.01, any divergences) and timing:\n")
print(res |> distinct(cond, rep, rhat_max, n_div, fit_sec, new_sec) |>
        group_by(cond) |>
        summarise(rhat_gt_101 = mean(rhat_max > 1.01), any_div = mean(n_div > 0),
                  fit_sec = median(fit_sec), restscore_sec = median(new_sec),
                  .groups = "drop"))

long <- res |>
  transmute(cond, model_type, n, rep, item,
            ppp = ppp, ppp_old = ppp_old, ppp_adj = ppp_adj,
            infit_ppp_cond = infit_ppp_cond,
            infit_ppp_theta = infit_ppp_theta) |>
  pivot_longer(c(ppp, ppp_old, ppp_adj, infit_ppp_cond, infit_ppp_theta),
               names_to = "method", values_to = "p") |>
  filter(!is.na(p)) |>
  mutate(p_two = if_else(method == "ppp_adj", p, pmin(1, 2 * pmin(p, 1 - p))))

tab <- long |>
  group_by(cond, method) |>
  summarise(mean = mean(p), sd = sd(p),
            rej05 = mean(p_two < .05), rej10 = mean(p_two < .10),
            .groups = "drop") |>
  left_join(long |> group_by(cond, method, rep) |>
              summarise(any05 = any(p_two < .05), .groups = "drop") |>
              group_by(cond, method) |>
              summarise(fwer05 = mean(any05), .groups = "drop"),
            by = c("cond", "method"))
cat("\nCalibration:\n")
print(as.data.frame(tab), digits = 3)
saveRDS(tab, file.path("dev", "results", paste0(name, "_table.rds")))

# Mean gamma_obs - gamma_exp, which should be close to 0 under the null.
cat("\nMean gamma_obs - gamma_exp by condition:\n")
print(res |> group_by(cond) |>
        summarise(mean_diff = mean(gamma_obs - gamma_exp),
                  se = sd(gamma_obs - gamma_exp) / sqrt(n()), .groups = "drop"))

# ECDF minus the uniform CDF for the one-sided per-item p-values. A
# calibrated p-value stays near 0. A conservative one bends below 0 at small
# p and above 0 at large p.
grid <- seq(0, 1, by = 0.01)
ecdf_df <- long |>
  filter(method != "ppp_adj") |>
  group_by(cond, method) |>
  reframe(p_grid = grid, diff = ecdf(p)(grid) - grid) |>
  mutate(statistic = if_else(grepl("^infit", method), "Infit", "Restscore"),
         method = if_else(method %in% c("ppp", "infit_ppp_cond"),
                          "Conditional (0.3.0.9000)", "Theta-based (0.3.0)"),
         panel = paste(statistic, cond))

okabe_ito <- c("#E69F00", "#0072B2")
p <- ggplot(ecdf_df, aes(p_grid, diff, colour = method, linetype = method,
                         shape = method)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.7) +
  geom_point(data = filter(ecdf_df, p_grid %in% seq(0, 1, by = 0.1)), size = 1.8) +
  facet_wrap(~ panel, ncol = 3) +
  scale_colour_manual(values = okabe_ito) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Posterior predictive p-value (one-sided)",
       y = "ECDF minus uniform CDF",
       colour = NULL, linetype = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path("dev", "results", paste0(name, "_ecdf.png")), p,
       width = 9, height = 3 * ceiling(length(unique(ecdf_df$panel)) / 3) + 1,
       dpi = 150)
cat("\nPlot saved to", file.path("dev", "results", paste0(name, "_ecdf.png")), "\n")
