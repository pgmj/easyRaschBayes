# Summary of dev/halo_followup.R. DRAFT.
#
# Run from the package root:
#   Rscript dev/halo_followup_summary.R             full run
#   SMOKE=1 Rscript dev/halo_followup_summary.R     smoke run
#
# Flags: Q3 per rater, one-sided max-T p_adj < .05 (more same-rater
# dependence than replicated). Infit per rater, two-sided max-T
# infit_p_adj < .05, split into underfit (observed infit above replicated)
# and overfit. MC SE at a rate of .05 with 100 replications: about .022.
#
# Tables:
#   1. Convergence and timing per condition.
#   2. Reproducibility: arm common against dev/results/halo_test_quick.rds
#      (same datasets, fits and Q3 seeds), pooled p and per-rater p_adj.
#   3. Pooled Q3 rejection rate at .05 per arm and condition (thres null:
#      calibration with rater thresholds in the model).
#   4. Per-rater flag rates by rater role: Q3, infit underfit, infit
#      overfit.
#   5. Joint flags per rater role: share of replications with Q3 only,
#      infit only, both, neither.
#   6. thres arm: threshold spread (third minus first threshold) per rater
#      role, which should single out the central-tendency rater.
# Figure:
#   halo_followup_flags.png  per-rater flag rates by role, Q3 vs infit.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
smoke <- Sys.getenv("SMOKE") == "1"
name <- if (smoke) "halo_followup_smoke" else "halo_followup"
res <- readRDS(file.path("dev", "results", paste0(name, ".rds")))

# Rater roles. Each person is rated by two adjacent raters (cyclic), so a
# rater shares persons with the raters on either side.
neighbours <- function(r) {
  i <- as.integer(sub("R", "", r))
  paste0("R", c((i - 2L) %% 8L + 1L, i %% 8L + 1L))
}
role_of <- function(rater, cond, model_type) {
  out <- rep("other", length(rater))
  if (cond == "misfit") {
    out[rater %in% neighbours("R2")] <- "shares persons with R2"
    if (model_type == "rsm") {
      out[rater %in% neighbours("R5")] <- "shares persons with R5"
      out[rater == "R5"] <- "central tendency (R5)"
    }
    out[rater == "R2"] <- "inconsistent (R2)"
  }
  if (cond == "halo_two") {
    out[rater %in% c(neighbours("R3"), neighbours("R7"))] <-
      "shares persons with a halo rater"
    out[rater %in% c("R3", "R7")] <- "halo (R3, R7)"
  }
  out
}
role_levels <- c("halo (R3, R7)", "shares persons with a halo rater",
                 "inconsistent (R2)", "shares persons with R2",
                 "central tendency (R5)", "shares persons with R5", "other")
pr <- res$per_rater |>
  rowwise() |>
  mutate(role = role_of(rater, cond, model_type)) |>
  ungroup() |>
  mutate(role = factor(role, levels = role_levels),
         q3_flag = p_adj < .05,
         infit_under = infit_p_adj < .05 & infit_obs > infit_rep,
         infit_over = infit_p_adj < .05 & infit_obs < infit_rep,
         infit_flag = infit_p_adj < .05)
cond_levels <- c("null", "misfit", "halo_two")

# --- 1. Convergence and timing ---
cat("Replications, convergence and timing per condition:\n")
print(res$info |>
        group_by(id) |>
        summarise(reps = n(), rhat_med = median(rhat_max),
                  rhat_gt105 = mean(rhat_max > 1.05),
                  any_div = mean(n_div > 0), fit_sec = median(fit_sec),
                  q3_sec = median(q3_sec), infit_sec = median(infit_sec),
                  .groups = "drop"),
      n = Inf)

# --- 2. Reproducibility of the quick run (arm common) ---
quick_file <- file.path("dev", "results",
                        if (smoke) "halo_test_quick_smoke.rds" else
                          "halo_test_quick.rds")
if (file.exists(quick_file)) {
  quick <- readRDS(quick_file)
  cmp <- res$pooled |>
    filter(arm == "common") |>
    inner_join(quick$pooled, by = c("model_type", "cond", "n", "rep",
                                    "statistic"),
               suffix = c("", "_quick"))
  cmp_r <- res$per_rater |>
    filter(arm == "common") |>
    inner_join(quick$per_rater, by = c("model_type", "cond", "n", "rep",
                                       "rater"),
               suffix = c("", "_quick"))
  cat("\nReproducibility against the quick run (arm common):\n")
  print(as.data.frame(bind_rows(
    cmp |> summarise(what = "pooled p", matched = n(),
                     identical = mean(abs(p - p_quick) < 1e-12),
                     max_abs_diff = max(abs(p - p_quick)),
                     cor = cor(p, p_quick)),
    cmp_r |> summarise(what = "per-rater p_adj", matched = n(),
                       identical = mean(abs(p_adj - p_adj_quick) < 1e-12),
                       max_abs_diff = max(abs(p_adj - p_adj_quick)),
                       cor = cor(p_adj, p_adj_quick))
  )), digits = 3)
} else {
  cat("\nNo quick-run results found for the reproducibility check.\n")
}

# --- 3. Pooled Q3 ---
cat("\nPooled Q3, rejection rate at .05, and smallest pair-level q3_p_adj",
    "over same-rater pairs < .05:\n")
pooled_tbl <- res$pooled |>
  filter(statistic == "mean_same") |>
  group_by(arm, model_type, cond) |>
  summarise(reps = n(), mean_same = mean(p < .05), mean_p = mean(p),
            .groups = "drop") |>
  left_join(res$pair |>
              group_by(arm, model_type, cond) |>
              summarise(pair_same = mean(min_p_adj_same < .05),
                        .groups = "drop"),
            by = c("arm", "model_type", "cond")) |>
  left_join(pr |>
              group_by(arm, model_type, cond, rep) |>
              summarise(q3_any = any(q3_flag), infit_any = any(infit_flag),
                        .groups = "drop") |>
              group_by(arm, model_type, cond) |>
              summarise(q3_rater_any = mean(q3_any),
                        infit_rater_any = mean(infit_any),
                        .groups = "drop"),
            by = c("arm", "model_type", "cond")) |>
  mutate(cond = factor(cond, levels = cond_levels)) |>
  arrange(arm, model_type, cond)
print(as.data.frame(pooled_tbl), digits = 3)

# --- 4. Per-rater flag rates by role ---
flag_tbl <- pr |>
  group_by(arm, model_type, cond, role) |>
  summarise(raters = n_distinct(rater),
            q3_diff = mean(obs_mean - rep_mean),
            q3_flag = mean(q3_flag),
            infit_obs = mean(infit_obs), infit_rep = mean(infit_rep),
            infit_under = mean(infit_under), infit_over = mean(infit_over),
            .groups = "drop") |>
  mutate(cond = factor(cond, levels = cond_levels)) |>
  arrange(arm, model_type, cond, role)
cat("\nPer-rater flag rates by role (rate per rater, p_adj < .05):\n")
print(as.data.frame(flag_tbl), digits = 3)

# --- 5. Joint flags ---
cat("\nJoint flags per rater role (share of rater x replication):\n")
print(as.data.frame(pr |>
                      mutate(pattern = case_when(
                        q3_flag & infit_flag ~ "both",
                        q3_flag ~ "q3_only",
                        infit_flag ~ "infit_only",
                        TRUE ~ "neither")) |>
                      count(arm, model_type, cond, role, pattern) |>
                      group_by(arm, model_type, cond, role) |>
                      mutate(share = n / sum(n)) |>
                      select(-n) |>
                      pivot_wider(names_from = pattern, values_from = share,
                                  values_fill = 0) |>
                      ungroup() |>
                      mutate(cond = factor(cond, levels = cond_levels)) |>
                      arrange(arm, model_type, cond, role)),
      digits = 3)

# --- 6. Threshold spread (thres arm) ---
if (!is.null(res$thresholds) && nrow(res$thresholds)) {
  spread <- res$thresholds |>
    group_by(cond, rep, rater) |>
    summarise(spread = max(location[threshold == max(threshold)]) -
                min(location[threshold == 1]), .groups = "drop") |>
    rowwise() |>
    mutate(role = role_of(rater, cond, "rsm")) |>
    ungroup() |>
    mutate(role = factor(role, levels = role_levels),
           cond = factor(cond, levels = cond_levels))
  cat("\nthres arm: threshold spread (third minus first), generating",
      "value 2.4 (4.32 for R5 in misfit):\n")
  print(as.data.frame(spread |>
                        group_by(cond, role) |>
                        summarise(mean = mean(spread), sd = sd(spread),
                                  .groups = "drop") |>
                        arrange(cond, role)),
        digits = 3)
}
saveRDS(list(pooled = pooled_tbl, flags = flag_tbl),
        file.path("dev", "results", paste0(name, "_table.rds")))

# --- Figure: per-rater flag rates by role (theme_bw, colour + line type +
# shape) ---
okabe_ito <- c("#E69F00", "#0072B2", "#009E73")
fig_df <- flag_tbl |>
  filter(cond != "null") |>
  select(arm, model_type, cond, role, q3_flag, infit_under, infit_over) |>
  pivot_longer(c(q3_flag, infit_under, infit_over),
               names_to = "statistic", values_to = "rate") |>
  mutate(statistic = factor(statistic,
                            levels = c("q3_flag", "infit_under",
                                       "infit_over"),
                            labels = c("Q3 same-rater", "infit underfit",
                                       "infit overfit")),
         panel = paste0(toupper(model_type), ", ", arm, ", ", cond),
         role = factor(role, levels = rev(role_levels)))
p <- ggplot(fig_df, aes(rate, role, colour = statistic, shape = statistic,
                        linetype = statistic, group = statistic)) +
  geom_vline(xintercept = .05, colour = "grey50", linetype = "dashed") +
  geom_path(linewidth = 0.6) +
  geom_point(size = 2.2) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = okabe_ito) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Flag rate per rater (p_adj < .05)", y = NULL, colour = NULL,
       shape = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
n_panels <- length(unique(fig_df$panel))
ggsave(file.path("dev", "results", paste0(name, "_flags.png")), p,
       width = 9, height = 2.4 * ceiling(n_panels / 2) + 1, dpi = 150)
cat("\nFigure saved to dev/results/", name, "_flags.png\n", sep = "")
