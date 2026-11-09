# Summary of dev/null_q3_global_calibration.R
#
# Run from the package root:
#   Rscript dev/null_q3_global_calibration_summary.R
# With SMOKE=1 it reads the smoke-test results instead.
#
# Per condition, q3_statistic() method and max-T variant:
#   fwer05       proportion of replications with any pair adjusted p < .05
#                (target .05, MC SE about 0.015 at 200 replications)
#   fwer10       the same at .10
#   flag_min, flag_max
#                lowest and highest per-pair rate of adjusted p < .05 over
#                replications. Even spread means each pair is near
#                fwer05 / number of pairs.
#   share_max    share of the replications' false flags that fall on the
#                most-flagged pair
#   cor_sd       correlation between a pair's flag rate and its mean
#                replicated SD: positive means the variant flags the most
#                variable pairs more often
# Plot: per-pair flag rate against the pair's mean replicated SD.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
smoke <- Sys.getenv("SMOKE") == "1"
name <- if (smoke) "null_q3_global_smoke" else "null_q3_global"
res <- readRDS(file.path("dev", "results", paste0(name, ".rds")))

cat("Replications per condition:\n")
print(res |> distinct(cond, rep) |> count(cond))
cat("\nConvergence, and largest difference between p_pair_std and",
    "q3_post()'s q3_p_adj (should be <= 0.0005, rounding):\n")
print(res |> distinct(cond, rep, rhat_max, n_div, check_cond, check_theta) |>
        group_by(cond) |>
        summarise(rhat_gt_101 = mean(rhat_max > 1.01), any_div = mean(n_div > 0),
                  check_cond = max(check_cond), check_theta = max(check_theta),
                  .groups = "drop"))

long <- res |>
  select(cond, method, rep, item_pair, sd_rep, starts_with("p_")) |>
  pivot_longer(starts_with("p_"), names_to = "variant", values_to = "p",
               names_prefix = "p_")

fwer <- long |>
  group_by(cond, method, variant, rep) |>
  summarise(any05 = any(p < .05), any10 = any(p < .10), .groups = "drop") |>
  group_by(cond, method, variant) |>
  summarise(fwer05 = mean(any05), fwer10 = mean(any10), .groups = "drop")

per_pair <- long |>
  group_by(cond, method, variant, item_pair) |>
  summarise(flag = mean(p < .05), sd_rep = mean(sd_rep), .groups = "drop")

balance <- per_pair |>
  group_by(cond, method, variant) |>
  summarise(flag_min = min(flag), flag_max = max(flag),
            share_max = if (sum(flag) > 0) max(flag) / sum(flag) else NA_real_,
            cor_sd = suppressWarnings(cor(flag, sd_rep)),
            .groups = "drop")

tab <- left_join(fwer, balance, by = c("cond", "method", "variant")) |>
  arrange(method, cond, variant)
cat("\nFamilywise error and spread of false flags over pairs:\n")
print(as.data.frame(tab), digits = 3)
saveRDS(tab, file.path("dev", "results", paste0(name, "_table.rds")))

p <- ggplot(per_pair, aes(sd_rep, flag, colour = variant, shape = variant)) +
  geom_point(alpha = 0.7) +
  facet_grid(method ~ cond, scales = "free_x") +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  scale_shape_manual(values = c(16, 17, 15, 3)) +
  labs(x = "Mean SD of replicated Q3 for the pair",
       y = "Proportion of replications with adjusted p < .05",
       colour = "Variant", shape = "Variant") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path("dev", "results", paste0(name, "_balance.png")), p,
       width = 11, height = 6, dpi = 150)
cat("\nSaved", paste0(name, "_balance.png"), "\n")
