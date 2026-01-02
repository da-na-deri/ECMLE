load("rdata/ex1_lvl.Rdata")
load("rdata/ex1_result.Rdata")


library(ggplot2)
library(dplyr)
library(tidyr)

df <- as.data.frame(ex1_lvl) %>%
  pivot_longer(
    cols = everything(),
    names_to = "HPD",
    values_to = "log_evidence"
  ) %>%
  mutate(HPD = as.numeric(HPD))

true_log_evidence <- ex1_result$log_exact_marginal_likelihood

plot1 <- ggplot(df, aes(x = factor(HPD), y = log_evidence, fill = factor(HPD))) +
  geom_boxplot(width = 0.6, outlier.size = 1) +
  geom_hline(
    yintercept = true_log_evidence,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    x = "HPD Levels",
    y = "Log Evidence",
    fill = "HPD"
  ) +
  scale_fill_brewer(palette = "Set3") +  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(plot1)

#========================== EXP2 =======================================
load("rdata/ex2_lvl.Rdata")
load("rdata/ex2_result.Rdata")


df <- as.data.frame(ex2_lvl) %>%
  pivot_longer(
    cols = everything(),
    names_to = "HPD",
    values_to = "log_evidence"
  ) %>%
  mutate(HPD = as.numeric(HPD))

true_log_evidence <- ex2_result$log_exact_marginal_likelihood

plot2 <- ggplot(df, aes(x = factor(HPD), y = log_evidence, fill = factor(HPD))) +
  geom_boxplot(width = 0.6, outlier.size = 1) +
  geom_hline(
    yintercept = true_log_evidence,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    x = "HPD Levels",
    y = "Log Evidence",
    fill = "HPD"
  ) +
  scale_fill_brewer(palette = "Set3") +   
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(plot2)