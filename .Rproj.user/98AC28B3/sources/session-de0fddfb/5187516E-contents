library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggh4x)

# ============================================================
# Configuration
# ============================================================
n_values <- c(20, 50, 100)

# ============================================================
# Load Results
# ============================================================
all_data <- list()

for (n in n_values) {
  filename <- sprintf("rdata/gmm2_2d_multiple_Y_n%d.RData", n)
  
  if (file.exists(filename)) {
    load(filename)
    all_data[[paste0("n", n)]] <- list(n = n, results = gmm_2d_multiple_Y_results)
  }
}

# ============================================================
# Combine Data
# ============================================================
df_combined <- bind_rows(
  lapply(names(all_data), function(key) {
    res <- all_data[[key]]$results
    n <- all_data[[key]]$n
    
    bind_rows(
      tibble(method = "THAMES", diff = res$thames_diff, time = res$time_thms, n = n),
      tibble(method = "MixTHAMES", diff = res$mixturethames_diff, time = res$time_mixthms, n = n),
      tibble(method = "ECMLE", diff = res$ecmle_diff, time = res$time_ecmle, n = n),
      tibble(method = "ePWK", diff = res$epwk_diff, time = res$time_pwk, n = n)
    )
  })
) %>%
  filter(is.finite(diff)) %>%
  mutate(
    ratio = exp(diff),
    method = factor(method, levels = c("THAMES", "MixTHAMES", "ECMLE", "ePWK")),
    n = factor(paste0("n = ", n), levels = paste0("n = ", n_values))
  )

# ============================================================
# Summary Statistics
# ============================================================
summary_stats <- df_combined %>%
  group_by(n, method) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    median_ratio = median(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    rmse_ratio = sqrt(mean((ratio - 1)^2, na.rm = TRUE)),
    mean_time = mean(time, na.rm = TRUE),
    n_valid = sum(is.finite(ratio)),
    .groups = "drop"
  )

print(summary_stats, n = Inf)

# ============================================================
# Mean Time Summary
# ============================================================
mean_time_by_method <- df_combined %>%
  filter(is.finite(time)) %>%
  group_by(n, method) %>%
  summarise(mean_time = mean(time, na.rm = TRUE), .groups = "drop")

print(mean_time_by_method, n = Inf)

# ============================================================
# Y-axis Limits
# ============================================================
y_scales <- list(
  scale_y_continuous(limits = c(0.95, 1.8)),
  scale_y_continuous(limits = c(0.85, 1.10)),
  scale_y_continuous(limits = c(0.95, 1.10))
)

# ============================================================
# Boxplot
# ============================================================
p_box <- ggplot(df_combined, aes(x = method, y = ratio, fill = method)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, outlier.alpha = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  facet_grid2(. ~ n, scales = "free_y", independent = "y", axes = "y") +
  facetted_pos_scales(y = y_scales) +
  scale_fill_manual(values = c(
    "THAMES" = "orchid", "MixTHAMES" = "skyblue",
    "ECMLE" = "salmon", "ePWK" = "lightgreen"
  )) +
  labs(y = "Estimated / Exact", x = NULL) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

print(p_box)