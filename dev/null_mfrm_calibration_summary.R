# Summary of dev/null_mfrm_calibration.R (many-facet fit statistics)
#
# Run from the package root:
#   Rscript dev/null_mfrm_calibration_summary.R
# With SMOKE=1 it reads the smoke-test results instead.
#
# Under the null a calibrated p-value is uniform. Halo settings "none" and
# "modelled" are null conditions. "unmodelled" (halo in the data, not in the
# model) is a power check and is reported separately.
#
# Reported per condition and statistic:
#   mean, sd     uniform has 0.5 and 0.289
#   rej05        per-unit two-sided rejection rate at .05, with
#                p_two = 2 * min(p, 1 - p) for the ppp columns. The p_adj
#                columns are already two-sided (restscore, infit) or one-sided
#                for dependence (Q3), and are used as they are.
#   fwer05       proportion of replications with any unit rejected at .05
#                (only meaningful for the adjusted p-values, should be <= .05)
# Units: cells (40, restscore and infit), raters (8, pooled infit), items
# (5, pooled infit), cell pairs by pair type (Q3).
# MC standard error at 200 replications: about 0.015 at a rate of .05 for
# fwer05. Per-unit rates pool 40 cells x 200 reps, so their SE is smaller,
# but cells within a replication are not independent.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
smoke <- Sys.getenv("SMOKE") == "1"
name <- if (smoke) "null_mfrm_smoke" else "null_mfrm"
res <- readRDS(file.path("dev", "results", paste0(name, ".rds")))

two_sided <- function(p) pmin(1, 2 * pmin(p, 1 - p))
rates <- function(df, p_col, adjusted = FALSE, group_cols = c("cond", "halo")) {
  df <- df |> filter(!is.na(.data[[p_col]]))
  df$p <- df[[p_col]]
  df$p_two <- if (adjusted) df$p else two_sided(df$p)
  per_unit <- df |>
    group_by(across(all_of(group_cols))) |>
    summarise(units = n(), mean = mean(p), sd = sd(p),
              rej05 = mean(p_two < .05), .groups = "drop")
  fwer <- df |>
    group_by(across(all_of(c(group_cols, "rep")))) |>
    summarise(any05 = any(p_two < .05), .groups = "drop") |>
    group_by(across(all_of(group_cols))) |>
    summarise(fwer05 = mean(any05), .groups = "drop")
  left_join(per_unit, fwer, by = group_cols) |>
    mutate(statistic = p_col, .before = 1)
}

# --- Replications, convergence, timing ---
cat("Replications, convergence and timing per condition:\n")
print(res$info |>
        group_by(cond) |>
        summarise(reps = n(), rhat_gt_101 = mean(rhat_max > 1.01),
                  any_div = mean(n_div > 0),
                  fit_sec = median(fit_sec), stat_sec = median(stat_sec),
                  .groups = "drop"),
      n = Inf)

# --- Calibration tables ---
cell_tab <- bind_rows(
  rates(res$cells, "rs_ppp"),
  rates(res$cells, "rs_p_adj", adjusted = TRUE),
  rates(res$cells, "infit_ppp")
)
group_tab <- bind_rows(
  rates(res$groups, "infit_ppp", group_cols = c("cond", "halo", "level")),
  rates(res$groups, "infit_p_adj", adjusted = TRUE,
        group_cols = c("cond", "halo", "level"))
)
q3_tab <- bind_rows(
  rates(res$q3, "q3_ppp", group_cols = c("cond", "halo", "pair_type")),
  rates(res$q3, "q3_p_adj", adjusted = TRUE,
        group_cols = c("cond", "halo", "pair_type"))
)

for (h in c("null", "power")) {
  keep <- if (h == "null") c("none", "modelled") else "unmodelled"
  cat("\n==== ", if (h == "null") "NULL CONDITIONS" else
        "POWER CHECK: halo in data, not in model", " ====\n", sep = "")
  cat("\nPer cell (restscore, infit):\n")
  print(as.data.frame(filter(cell_tab, halo %in% keep)), digits = 3)
  cat("\nInfit pooled by rater and by item:\n")
  print(as.data.frame(filter(group_tab, halo %in% keep)), digits = 3)
  cat("\nQ3 by pair type (q3_ppp: P(obs > rep), q3_p_adj one-sided):\n")
  print(as.data.frame(filter(q3_tab, halo %in% keep)), digits = 3)
}
saveRDS(list(cells = cell_tab, groups = group_tab, q3 = q3_tab),
        file.path("dev", "results", paste0(name, "_table.rds")))

# --- Expected vs observed gamma, Q3 centre ---
cat("\nMean gamma_obs - gamma_exp per condition (about 0 under the null):\n")
print(res$cells |> group_by(cond) |>
        summarise(mean_diff = mean(gamma_obs - gamma_exp),
                  se = sd(gamma_obs - gamma_exp) / sqrt(n()),
                  .groups = "drop"), n = Inf)
cat("\nMean observed and replicated Q3 by pair type:\n")
print(res$q3 |> group_by(cond, pair_type) |>
        summarise(q3_obs = mean(q3_obs), q3_rep = mean(q3_rep),
                  .groups = "drop"), n = Inf)

# --- Facet parameters: coverage and bias for rater severity ---
cat("\nRater severity: 95% HDCI coverage, bias and RMSE (centred truth):\n")
print(res$facet |>
        group_by(cond) |>
        summarise(coverage = mean(truth >= hdci_lower & truth <= hdci_upper),
                  bias = mean(location - truth),
                  rmse = sqrt(mean((location - truth)^2)),
                  .groups = "drop"), n = Inf)
cat("\nFair minus observed average, by rater (rotation design: each rater",
    "sees a different person mix, so these differ by sampling only):\n")
print(res$facet |>
        group_by(cond, rater) |>
        summarise(diff = mean(fair_avg - observed_avg), .groups = "drop") |>
        group_by(cond) |>
        summarise(max_abs_diff = max(abs(diff)), .groups = "drop"), n = Inf)

# --- ECDF minus uniform CDF, null conditions, one-sided ppp ---
# A calibrated p-value stays near 0. A conservative one bends below 0 at
# small p and above 0 at large p.
grid <- seq(0, 1, by = 0.01)
ecdf_long <- bind_rows(
  res$cells |> filter(halo != "unmodelled") |>
    transmute(cond, halo, statistic = "Restscore, cell", p = rs_ppp),
  res$cells |> filter(halo != "unmodelled") |>
    transmute(cond, halo, statistic = "Infit, cell", p = infit_ppp),
  res$groups |> filter(halo != "unmodelled", level == "rater") |>
    transmute(cond, halo, statistic = "Infit, rater", p = infit_ppp)
) |>
  filter(!is.na(p)) |>
  mutate(model = sub("_.*", "", cond),
         n = sub(".*_n", "n = ", cond),
         halo = factor(halo, levels = c("none", "modelled"),
                       labels = c("No halo", "Halo modelled")))
ecdf_df <- ecdf_long |>
  group_by(model, n, halo, statistic) |>
  reframe(p_grid = grid, diff = ecdf(p)(grid) - grid) |>
  mutate(panel = paste(toupper(model), n, halo, sep = ", "))

okabe_ito <- c("#E69F00", "#0072B2", "#009E73")
p <- ggplot(ecdf_df, aes(p_grid, diff, colour = statistic,
                         linetype = statistic, shape = statistic)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.7) +
  geom_point(data = filter(ecdf_df, p_grid %in% seq(0, 1, by = 0.1)),
             size = 1.8) +
  facet_wrap(~ panel, ncol = 2) +
  scale_colour_manual(values = okabe_ito) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = "Posterior predictive p-value (one-sided)",
       y = "ECDF minus uniform CDF",
       colour = NULL, linetype = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
n_panels <- length(unique(ecdf_df$panel))
ggsave(file.path("dev", "results", paste0(name, "_ecdf.png")), p,
       width = 8, height = 2.6 * ceiling(n_panels / 2) + 1, dpi = 150)
cat("\nPlot saved to", file.path("dev", "results", paste0(name, "_ecdf.png")),
    "\n")
