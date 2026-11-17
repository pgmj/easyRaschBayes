# Summary of dev/halo_test_calibration.R (pooled halo test). DRAFT.
# Works for the full design (dev/halo_test_calibration.R) and the quick
# one (dev/halo_test_quick.R).
#
# Run from the package root:
#   Rscript dev/halo_test_calibration_summary.R                      full run
#   RUN=halo_test_quick Rscript dev/halo_test_calibration_summary.R  quick run
# With SMOKE=1 it reads the smoke-test results of that run instead.
#
# All p-values are one-sided: small = more same-rater dependence than
# replicated. Under the null a calibrated p is uniform, so p < .05 in 5% of
# replications. MC SE at a rate of .05: about .015 with 200 replications
# (null), .022 with 100 (misfit and the halo conditions).
#
# Tables:
#   1. Convergence and timing per condition.
#   2. Rejection rate at .05 per condition and statistic:
#      mean_same, contrast         p < .05
#      rater_any                   any rater with p_adj < .05
#      pair_same                   any same-rater pair with q3_p_adj < .05
#      Under null and misfit these are false-alarm rates (should be <= .05
#      for the adjusted ones and about .05 for the pooled ones). Under the
#      halo conditions they are power.
#   3. halo_two: per-rater flags (p_adj < .05) for the two halo raters and
#      the six others.
#   4. misfit: per-rater p for the inconsistent and central-tendency
#      raters (do they look like halo?).
#   5. Mean observed minus replicated for the pooled statistics.
# Figures:
#   halo_test_ecdf.png    ECDF minus uniform for the pooled p under null,
#                         misfit and halo06_model.
#   halo_test_power.png   rejection rate by condition, pooled vs pair-level.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
smoke <- Sys.getenv("SMOKE") == "1"
run <- Sys.getenv("RUN", "halo_test")
name <- if (smoke) paste0(run, "_smoke") else run
res <- readRDS(file.path("dev", "results", paste0(name, ".rds")))
# Files written before 2026-11-16 (first smoke test) named the means obs and
# rep, which clashed with the replication number `rep` (saved as rep.1).
fix_names <- function(df) {
  if ("rep.1" %in% names(df)) {
    names(df)[names(df) == "obs"] <- "obs_mean"
    names(df)[names(df) == "rep.1"] <- "rep_mean"
  }
  df
}
res$pooled <- fix_names(res$pooled)
res$per_rater <- fix_names(res$per_rater)
cond_levels <- c("null", "misfit", "halo06_model", "halo03", "halo06",
                 "halo_two")

# --- 1. Convergence and timing ---
cat("Replications, convergence and timing per condition:\n")
print(res$info |>
        group_by(id) |>
        summarise(reps = n(), rhat_med = median(rhat_max),
                  rhat_gt105 = mean(rhat_max > 1.05),
                  any_div = mean(n_div > 0),
                  fit_sec = median(fit_sec), stat_sec = median(stat_sec),
                  .groups = "drop"),
      n = Inf)

# --- 2. Rejection rates ---
rates <- bind_rows(
  res$pooled |>
    group_by(model_type, cond, n, statistic) |>
    summarise(rej05 = mean(p < .05), mean_p = mean(p), sd_p = sd(p),
              .groups = "drop"),
  res$per_rater |>
    group_by(model_type, cond, n, rep) |>
    summarise(any = any(p_adj < .05), .groups = "drop") |>
    group_by(model_type, cond, n) |>
    summarise(rej05 = mean(any), .groups = "drop") |>
    mutate(statistic = "rater_any"),
  res$pair |>
    group_by(model_type, cond, n) |>
    summarise(rej05 = mean(min_p_adj_same < .05), .groups = "drop") |>
    mutate(statistic = "pair_same")
) |>
  mutate(cond = factor(cond, levels = cond_levels),
         statistic = factor(statistic, levels = c("mean_same", "contrast",
                                                  "rater_any", "pair_same"))) |>
  arrange(model_type, cond, n, statistic)
cat("\nRejection rate at .05 (false alarms under null, misfit and",
    "halo06_model, power otherwise):\n")
print(as.data.frame(rates |>
                      select(model_type, cond, n, statistic, rej05) |>
                      pivot_wider(names_from = statistic,
                                  values_from = rej05)),
      digits = 3)
cat("\nPooled p: mean and SD (uniform: 0.5 and 0.289):\n")
print(as.data.frame(rates |>
                      filter(statistic %in% c("mean_same", "contrast"),
                             cond %in% c("null", "misfit", "halo06_model")) |>
                      select(model_type, cond, n, statistic, mean_p, sd_p)),
      digits = 3)

# --- 3. halo_two: which raters are flagged ---
cat("\nhalo_two: per-rater flag rate (p_adj < .05) for halo and other",
    "raters:\n")
print(res$per_rater |>
        filter(cond == "halo_two") |>
        mutate(rater_type = if_else(halo_sd_true > 0, "halo rater",
                                    "other rater")) |>
        group_by(model_type, n, rater_type) |>
        summarise(flag_rate = mean(p_adj < .05), mean_p = mean(p),
                  .groups = "drop"),
      n = Inf)

# --- 4. misfit: do non-halo rater problems look like halo? ---
cat("\nmisfit: per-rater mean Q3 difference and flag rate by rater type:\n")
print(res$per_rater |>
        filter(cond == "misfit") |>
        mutate(rater_type = case_when(
          disc_true < 1 ~ "inconsistent (R2)",
          stretch_true > 1 ~ "central tendency (R5)",
          TRUE ~ "other")) |>
        group_by(model_type, n, rater_type) |>
        summarise(q3_diff = mean(obs_mean - rep_mean), mean_p = mean(p),
                  p05 = mean(p < .05), flag_adj = mean(p_adj < .05),
                  .groups = "drop"),
      n = Inf)

# --- 5. Size of the pooled signal ---
cat("\nPooled statistics: mean observed minus replicated:\n")
print(res$pooled |>
        mutate(cond = factor(cond, levels = cond_levels)) |>
        group_by(model_type, cond, n, statistic) |>
        summarise(diff = mean(obs_mean - rep_mean), .groups = "drop") |>
        pivot_wider(names_from = statistic, values_from = diff) |>
        arrange(model_type, cond, n),
      n = Inf)
saveRDS(rates, file.path("dev", "results", paste0(name, "_table.rds")))

# --- Figures (Okabe-Ito + line type + shape, theme_bw) ---
okabe_ito <- c("#E69F00", "#0072B2", "#009E73", "#CC79A7")

grid <- seq(0, 1, by = 0.01)
ecdf_df <- res$pooled |>
  filter(cond %in% c("null", "misfit", "halo06_model")) |>
  group_by(model_type, cond, n, statistic) |>
  reframe(p_grid = grid, diff = ecdf(p)(grid) - grid) |>
  mutate(panel = paste0(toupper(model_type), ", ", cond),
         series = paste0(statistic, ", n = ", n))
p1 <- ggplot(ecdf_df, aes(p_grid, diff, colour = series, linetype = series,
                          shape = series)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.7) +
  geom_point(data = filter(ecdf_df, p_grid %in% seq(0, 1, by = 0.1)),
             size = 1.8) +
  facet_wrap(~ panel, ncol = 2) +
  scale_colour_manual(values = okabe_ito) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(x = "Halo test p-value (one-sided)", y = "ECDF minus uniform CDF",
       colour = NULL, linetype = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
n_panels <- length(unique(ecdf_df$panel))
ggsave(file.path("dev", "results", paste0(name, "_ecdf.png")), p1,
       width = 8, height = 2.6 * ceiling(n_panels / 2) + 1, dpi = 150)

power_df <- rates |>
  filter(cond %in% c("null", "halo03", "halo06", "halo_two")) |>
  mutate(n = paste("n =", n))
p2 <- ggplot(power_df, aes(cond, rej05, colour = statistic,
                           linetype = statistic, shape = statistic,
                           group = statistic)) +
  geom_hline(yintercept = .05, colour = "grey50", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_grid(model_type ~ n, labeller = labeller(model_type = toupper)) +
  scale_colour_manual(values = okabe_ito) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Rejection rate at .05", colour = NULL,
       linetype = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path("dev", "results", paste0(name, "_power.png")), p2,
       width = 8, height = 6, dpi = 150)
cat("\nFigures saved to dev/results/", name, "_ecdf.png and _power.png\n",
    sep = "")
