# Summary of dev/null_infit_adj_calibration.R
#
# Run from the package root:
#   Rscript dev/null_infit_adj_calibration_summary.R
# With SMOKE=1 it reads the smoke-test results instead.
#
# Per condition and method (conditional, theta):
#   Per item, with uniform targets in brackets
#     ppp_sd       SD of ppp (0.289)
#     ppp_rej05    two-sided rejection of ppp at .05, 2 * min(p, 1 - p) (.05)
#     marg_sd      SD of p_marg (0.289)
#     marg_rej05   p_marg < .05 (.05)
#   Per replication (familywise error, target .05, MC SE about 0.015)
#     fwer_ppp     any item with two-sided ppp < .05 (no adjustment)
#     fwer_marg    any item with p_marg < .05 (no adjustment)
#     fwer_adj     any item with p_adj < .05
#     fwer_adj_draw  any item with p_adj_draw < .05
#   argmax_chisq   chi-square p for whether each item is equally often the
#                  most extreme (smallest p_adj), pooled over sample sizes
# Plot: ECDF minus uniform for ppp (two-sided) and p_marg.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
smoke <- Sys.getenv("SMOKE") == "1"
name <- if (smoke) "null_infit_adj_smoke" else "null_infit_adj"
res <- readRDS(file.path("dev", "results", paste0(name, ".rds")))

cat("Replications per condition:\n")
print(res |> distinct(cond, rep) |> count(cond))
cat("\nConvergence, and largest difference between p_adj and the package",
    "helper .ppp_maxt() (should be 0):\n")
print(res |> group_by(cond) |>
        summarise(rhat_gt_101 = mean(rhat_max > 1.01), any_div = mean(n_div > 0),
                  check = max(abs(p_adj - p_adj_pkg), na.rm = TRUE),
                  .groups = "drop"))

cat("\nMean infit and SD across draws (observed, replicated):\n")
print(res |> group_by(cond, method) |>
        summarise(obs = mean(infit_obs), rep = mean(infit_rep),
                  sd_obs = mean(sd_obs), sd_rep = mean(sd_rep), .groups = "drop") |>
        as.data.frame(), digits = 3)

res <- res |> mutate(p_two = pmin(1, 2 * pmin(ppp, 1 - ppp)))

per_item <- res |>
  group_by(cond, method) |>
  summarise(ppp_sd = sd(ppp), ppp_rej05 = mean(p_two < .05),
            marg_sd = sd(p_marg), marg_rej05 = mean(p_marg < .05),
            .groups = "drop")

fwer <- res |>
  group_by(cond, method, rep) |>
  summarise(a_ppp = any(p_two < .05), a_marg = any(p_marg < .05),
            a_adj = any(p_adj < .05), a_draw = any(p_adj_draw < .05),
            .groups = "drop") |>
  group_by(cond, method) |>
  summarise(fwer_ppp = mean(a_ppp), fwer_marg = mean(a_marg),
            fwer_adj = mean(a_adj), fwer_adj_draw = mean(a_draw),
            .groups = "drop")

set.seed(1)
argmax <- res |>
  group_by(model_type, cond, method, rep) |>
  mutate(u = runif(n())) |>
  arrange(p_adj, u, .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  count(model_type, method, item) |>
  right_join(res |> distinct(model_type, method, item),
             by = c("model_type", "method", "item")) |>
  mutate(n = coalesce(n, 0L)) |>
  group_by(model_type, method) |>
  summarise(top_share = max(n) / sum(n), uniform = 1 / n(),
            argmax_chisq = suppressWarnings(chisq.test(n)$p.value),
            .groups = "drop")

tab <- left_join(per_item, fwer, by = c("cond", "method")) |>
  arrange(method, cond)
cat("\nCalibration:\n")
print(as.data.frame(tab), digits = 3)
cat("\nWhich item is the most extreme (smallest p_adj):\n")
print(as.data.frame(argmax), digits = 3)
cat("\nPooled over conditions:\n")
print(as.data.frame(res |> group_by(method, cond, rep) |>
        summarise(a_adj = any(p_adj < .05), a_draw = any(p_adj_draw < .05),
                  .groups = "drop") |>
        group_by(method) |>
        summarise(fwer_adj = mean(a_adj), fwer_adj_draw = mean(a_draw))),
      digits = 3)
saveRDS(list(table = tab, argmax = argmax),
        file.path("dev", "results", paste0(name, "_table.rds")))

grid <- seq(0, 1, by = 0.01)
ecdf_df <- res |>
  select(cond, method, p_two, p_marg) |>
  pivot_longer(c(p_two, p_marg), names_to = "pvalue", values_to = "p") |>
  group_by(cond, method, pvalue) |>
  reframe(p_grid = grid, diff = ecdf(p)(grid) - grid) |>
  mutate(pvalue = recode(pvalue, p_two = "ppp, two-sided",
                         p_marg = "Standardised, per item"))

p <- ggplot(ecdf_df, aes(p_grid, diff, colour = pvalue, linetype = pvalue)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.7) +
  facet_grid(method ~ cond) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = "p", y = "ECDF minus uniform CDF", colour = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path("dev", "results", paste0(name, "_ecdf.png")), p,
       width = 11, height = 6, dpi = 150)
cat("\nSaved", paste0(name, "_ecdf.png"), "\n")
