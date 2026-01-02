library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggh4x)

# ============================================================
# Configuration
# ============================================================
file_mapping <- list(
  list(d = 2,  n = 20,  file = "rdata/banana_d2_multiple_Y_n20_results.RData"),
  list(d = 2,  n = 50,  file = "rdata/banana_d2_multiple_Y_n50_results.RData"),
  list(d = 2,  n = 100, file = "rdata/banana_d2_multiple_Y_n100_results.RData"),
  list(d = 5,  n = 20,  file = "rdata/banana_d5_multiple_Y_n20_results.RData"),
  list(d = 5,  n = 50,  file = "rdata/banana_d5_multiple_Y_n50_results.RData"),
  list(d = 5,  n = 100, file = "rdata/banana_d5_multiple_Y_n100_results.RData"),
  list(d = 10, n = 20,  file = "rdata/banana_d10_multiple_Y_n20_results.RData"),
  list(d = 10, n = 50,  file = "rdata/banana_d10_multiple_Y_n50_results.RData"),
  list(d = 10, n = 100, file = "rdata/banana_d10_multiple_Y_n100_results.RData")
)

d_values <- c(2, 5, 10)
n_values <- c(20, 50, 100)

# ============================================================
# Load Results
# ============================================================
all_data <- list()

for (item in file_mapping) {
  if (file.exists(item$file)) {
    temp_env <- new.env()
    load(item$file, envir = temp_env)
    result <- get(ls(temp_env)[1], envir = temp_env)
    all_data[[paste0("d", item$d, "_n", item$n)]] <- list(
      d = item$d, n = item$n, results = result
    )
  }
}

# ============================================================
# Plot 1: All 4 methods with wide y-axis limits
# ============================================================
df_combined <- bind_rows(
  lapply(names(all_data), function(key) {
    res <- all_data[[key]]$results
    d_val <- all_data[[key]]$d
    n_val <- all_data[[key]]$n
    
    bind_rows(
      tibble(method = "THAMES", evidence = exp(as.numeric(res$thames)), 
             time = res$time_thms, d = d_val, n = n_val),
      tibble(method = "MixTHAMES", evidence = exp(res$mixturethames), 
             time = res$time_mixthms, d = d_val, n = n_val),
      tibble(method = "ECMLE", evidence = exp(res$ecmle), 
             time = res$time_ecmle, d = d_val, n = n_val),
      tibble(method = "ePWK", evidence = exp(as.numeric(res$epwk)), 
             time = res$time_pwk, d = d_val, n = n_val)
    )
  })
) %>%
  filter(is.finite(evidence)) %>%
  mutate(
    method = factor(method, levels = c("THAMES", "MixTHAMES", "ECMLE", "ePWK")),
    d = factor(paste0("d = ", d), levels = paste0("d = ", d_values)),
    n = factor(paste0("n = ", n), levels = paste0("n = ", n_values))
  )

y_scales <- list(
  scale_y_continuous(limits = c(0.95, 1.15)),
  scale_y_continuous(limits = c(0.97, 1.05)),
  scale_y_continuous(limits = c(0.98, 1.02)),
  scale_y_continuous(limits = c(0.9, 1e7)),
  scale_y_continuous(limits = c(0.95, 5e4)),
  scale_y_continuous(limits = c(0.97, 5e4)),
  scale_y_continuous(limits = c(0.95, 4e28)),
  scale_y_continuous(limits = c(0.95, 2e255)),
  scale_y_continuous(limits = c(0.97, 1.5e251))
)

summary_stats <- df_combined %>%
  group_by(d, n, method) %>%
  summarise(
    mean_evidence = mean(evidence, na.rm = TRUE),
    median_evidence = median(evidence, na.rm = TRUE),
    sd_evidence = sd(evidence, na.rm = TRUE),
    rmse = sqrt(mean((evidence - 1)^2, na.rm = TRUE)),
    mean_time = mean(time, na.rm = TRUE),
    n_valid = sum(is.finite(evidence)),
    .groups = "drop"
  )

print(summary_stats, n = 50)

p_box <- ggplot(df_combined, aes(x = method, y = evidence, fill = method)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, outlier.alpha = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  facet_grid2(d ~ n, scales = "free_y", independent = "y", axes = "y") +
  facetted_pos_scales(y = y_scales) +
  scale_fill_manual(values = c(
    "THAMES" = "orchid", "MixTHAMES" = "skyblue",
    "ECMLE" = "salmon", "ePWK" = "lightgreen"
  )) +
  labs(y = "Evidence", x = NULL) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

print(p_box)


# ============================================================
# Plot 2: All 4 methods with moderate y-axis limits
# ============================================================
df_combined <- bind_rows(
  lapply(names(all_data), function(key) {
    res <- all_data[[key]]$results
    d_val <- all_data[[key]]$d
    n_val <- all_data[[key]]$n
    
    bind_rows(
      tibble(method = "THAMES", evidence = exp(res$thames), 
             time = res$time_thms, d = d_val, n = n_val),
      tibble(method = "MixTHAMES", evidence = exp(res$mixturethames), 
             time = res$time_mixthms, d = d_val, n = n_val),
      tibble(method = "ECMLE", evidence = exp(res$ecmle), 
             time = res$time_ecmle, d = d_val, n = n_val),
      tibble(method = "ePWK", evidence = exp(res$epwk), 
             time = res$time_pwk, d = d_val, n = n_val)
    )
  })
) %>%
  filter(is.finite(evidence)) %>%
  mutate(
    method = factor(method, levels = c("THAMES", "MixTHAMES", "ECMLE", "ePWK")),
    d = factor(paste0("d = ", d), levels = paste0("d = ", d_values)),
    n = factor(paste0("n = ", n), levels = paste0("n = ", n_values))
  )

y_scales <- list(
  scale_y_continuous(limits = c(0.95, 1.15)),
  scale_y_continuous(limits = c(0.97, 1.05)),
  scale_y_continuous(limits = c(0.98, 1.01)),
  scale_y_continuous(limits = c(0.6, 200)),
  scale_y_continuous(limits = c(0.7, 100)),
  scale_y_continuous(limits = c(0.7, 50)),
  scale_y_continuous(limits = c(0, 150000)),
  scale_y_continuous(limits = c(0, 3000)),
  scale_y_continuous(limits = c(0, 2000))
)

summary_stats <- df_combined %>%
  group_by(d, n, method) %>%
  summarise(
    mean_evidence = mean(evidence, na.rm = TRUE),
    median_evidence = median(evidence, na.rm = TRUE),
    sd_evidence = sd(evidence, na.rm = TRUE),
    rmse = sqrt(mean((evidence - 1)^2, na.rm = TRUE)),
    mean_time = mean(time, na.rm = TRUE),
    n_valid = sum(is.finite(evidence)),
    .groups = "drop"
  )

print(summary_stats, n = 50)

p_box <- ggplot(df_combined, aes(x = method, y = evidence, fill = method)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, outlier.alpha = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  facet_grid2(d ~ n, scales = "free_y", independent = "y", axes = "y") +
  facetted_pos_scales(y = y_scales) +
  scale_fill_manual(values = c(
    "THAMES" = "orchid", "MixTHAMES" = "skyblue",
    "ECMLE" = "salmon", "ePWK" = "lightgreen"
  )) +
  labs(y = "Evidence", x = NULL) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

print(p_box)


# ============================================================
# Plot 3: Only MixTHAMES and ECMLE with tight y-axis limits
# ============================================================
df_combined <- bind_rows(
  lapply(names(all_data), function(key) {
    res <- all_data[[key]]$results
    d_val <- all_data[[key]]$d
    n_val <- all_data[[key]]$n
    
    bind_rows(
      tibble(method = "MixTHAMES", evidence = exp(res$mixturethames), 
             time = res$time_mixthms, d = d_val, n = n_val),
      tibble(method = "ECMLE", evidence = exp(res$ecmle), 
             time = res$time_ecmle, d = d_val, n = n_val)
    )
  })
) %>%
  filter(is.finite(evidence)) %>%
  mutate(
    method = factor(method, levels = c("MixTHAMES", "ECMLE")),
    d = factor(paste0("d = ", d), levels = paste0("d = ", d_values)),
    n = factor(paste0("n = ", n), levels = paste0("n = ", n_values))
  )

y_scales <- list(
  scale_y_continuous(limits = c(0.97, 1.03)),
  scale_y_continuous(limits = c(0.97, 1.03)),
  scale_y_continuous(limits = c(0.98, 1.02)),
  scale_y_continuous(limits = c(0.6, 1.4)),
  scale_y_continuous(limits = c(0.7, 1.3)),
  scale_y_continuous(limits = c(0.7, 1.3)),
  scale_y_continuous(limits = c(0, 3)),
  scale_y_continuous(limits = c(0, 3)),
  scale_y_continuous(limits = c(0, 3))
)

summary_stats <- df_combined %>%
  group_by(d, n, method) %>%
  summarise(
    mean_evidence = mean(evidence, na.rm = TRUE),
    median_evidence = median(evidence, na.rm = TRUE),
    sd_evidence = sd(evidence, na.rm = TRUE),
    rmse = sqrt(mean((evidence - 1)^2, na.rm = TRUE)),
    mean_time = mean(time, na.rm = TRUE),
    n_valid = sum(is.finite(evidence)),
    .groups = "drop"
  )

print(summary_stats, n = 50)

p_box <- ggplot(df_combined, aes(x = method, y = evidence, fill = method)) +
  geom_boxplot(width = 0.65, outlier.size = 0.6, outlier.alpha = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  facet_grid2(d ~ n, scales = "free_y", independent = "y", axes = "y") +
  facetted_pos_scales(y = y_scales) +
  scale_fill_manual(values = c("MixTHAMES" = "skyblue", "ECMLE" = "salmon")) +
  labs(y = "Evidence", x = NULL) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

print(p_box)