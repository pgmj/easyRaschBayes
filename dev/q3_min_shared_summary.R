# Summary of dev/q3_min_shared.R (minimum shared persons for Q3 cell
# pairs). Run from the package root:
#   Rscript dev/q3_min_shared_summary.R        (SMOKE=1 for the smoke run)
# MC SE at a rate of .05: about .015 with 200 null replicates, .022 with
# 100 per halo condition.
#
# Tables:
#   1. Per-pair calibration under null by shared-person bin: rate of mid-p
#      < .05 (should be about .05) and of q3_ppp > .95 (as q3_post()).
#   2. Max-T by minimum m: familywise rate (null), power for R1's
#      same-rater pairs (halo06, halo10), size of the family, 95th
#      percentile of the max-T reference, and the share of draws whose
#      maximum comes from a pair with fewer than 10 shared persons.
# Figure: q3_min_shared_power.png, power and FWER against m.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
smoke <- Sys.getenv("SMOKE") == "1"
name <- if (smoke) "q3_min_shared_smoke" else "q3_min_shared"
res <- readRDS(file.path("dev", "results", paste0(name, ".rds")))

cat("Replicates per condition and seconds per replicate:\n")
print(res$info |> group_by(cond) |>
        summarise(reps = n(), sec = median(sec), .groups = "drop"))

bins <- c(0, 2, 4, 9, 19, 49, Inf)
cat("\n1. Null, per pair, by number of shared persons:\n")
print(as.data.frame(
  res$pair |>
    filter(cond == "null") |>
    mutate(shared = cut(n_shared, bins, labels = c("2", "3-4", "5-9",
                                                   "10-19", "20-49",
                                                   "50+"))) |>
    group_by(shared) |>
    # Q3 is undefined (NA) when a pair's residuals are constant, about 4%
    # of the pairs with 2 shared persons
    summarise(pairs = n_distinct(pair), na = mean(is.na(q3_obs)),
              midp05 = mean(midp < .05, na.rm = TRUE),
              ppp95 = mean(q3_ppp > .95, na.rm = TRUE),
              mean_q3_obs = mean(q3_obs, na.rm = TRUE),
              mean_q3_rep = mean(q3_rep, na.rm = TRUE), .groups = "drop")
), digits = 3)

fam <- res$family |>
  group_by(cond, min_shared) |>
  summarise(n_family = n_family[1],
            rate = if (cond[1] == "null") mean(any_flag) else
              mean(target_flag),
            tmax_q95 = mean(tmax_q95), tmax_small = mean(tmax_small),
            .groups = "drop")
cat("\n2. Max-T by minimum shared persons (rate = FWER under null,",
    "power for R1's same-rater pairs under halo):\n")
print(as.data.frame(fam |> select(cond, min_shared, n_family, rate) |>
                      pivot_wider(names_from = cond, values_from = rate) |>
                      left_join(fam |> filter(cond == "null") |>
                                  select(min_shared, tmax_q95, tmax_small),
                                by = "min_shared")),
      digits = 3)
saveRDS(fam, file.path("dev", "results", paste0(name, "_table.rds")))

okabe_ito <- c(null = "#E69F00", halo06 = "#0072B2", halo10 = "#009E73")
p <- ggplot(fam, aes(min_shared, rate, colour = cond, linetype = cond,
                     shape = cond)) +
  geom_hline(yintercept = .05, colour = "grey50", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = okabe_ito) +
  scale_linetype_manual(values = c(null = "solid", halo06 = "dashed",
                                   halo10 = "dotted")) +
  scale_shape_manual(values = c(null = 16, halo06 = 17, halo10 = 15)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Minimum number of shared persons in the max-T family",
       y = "Rate of q3_p_adj < .05\n(null: any pair, halo: R1's pairs)",
       colour = NULL, linetype = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path("dev", "results", paste0(name, "_power.png")), p,
       width = 6, height = 4.5, dpi = 150)
cat("\nFigure saved to dev/results/", name, "_power.png\n", sep = "")
