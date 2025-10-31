rm(list=ls())
library(data.table)
library(tidyverse)
library(plyr)
library(dplyr)
library(survival)
library(Matrix)
library(glmnet)
library(GGally)
library(ggplot2)
library(ggsci)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_1022_TTFD10")

files <- list.files(pattern = "^cox_td_OS_.*\\.rds$")

results <- data.frame(
  model_id = character(),
  exp_coef = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (f in files) {
  model <- readRDS(f)
  sum_model <- summary(model)
  
  if ("has_TTFD10" %in% rownames(sum_model$coefficients)) {
    exp_coef <- sum_model$conf.int["has_TTFD10", "exp(coef)"]
    lower <- sum_model$conf.int["has_TTFD10", "lower .95"]
    upper <- sum_model$conf.int["has_TTFD10", "upper .95"]
    pval <- sum_model$coefficients["has_TTFD10", "Pr(>|z|)"]
    
    results <- rbind(results, data.frame(
      model_id = gsub(".rds", "", f),
      exp_coef = exp_coef,
      lower_95 = lower,
      upper_95 = upper,
      p_value = pval
    ))
  }
}

results$model_id <- factor(results$model_id, levels = results$model_id)

results <- results %>%
  mutate(
    suffix = gsub(".*_OS_(.*?)_\\d{8}", "\\1", model_id)
  ) %>%
  arrange(p_value)  

results$suffix <- factor(results$suffix, levels = results$suffix)

results <- results %>%
  mutate(p_label = formatC(p_value, format = "e", digits = 2))

ggplot(results, aes(x = suffix, y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(
    aes(label = p_label, x = suffix, y = upper_95 + 0.015),  # Slightly shift labels to the right to avoid overlapping with error bars
    hjust = 0,
    size = 3.5
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.2))) +  # Upper y-limit 0.5 with extra padding
  labs(
    x = NULL,
    y = "Hazard Ratio of TTFD10",
    title = "TTFD10 in OS Models"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),  # Increase right margin
    plot.title = element_text(hjust = 0.5)
  )

# Construct display string: HR (95% CI), keep 4 decimals
results <- results %>%
  mutate(
    suffix = gsub(".*_OS_(.*?)_\\d{8}", "\\1", model_id),
    hr_ci_label = sprintf("%.4f (%.4f, %.4f)", exp_coef, lower_95, upper_95),
    p_label = formatC(p_value, format = "e", digits = 2)
  ) %>%
  arrange(p_value)

results$suffix <- factor(results$suffix, levels = results$suffix)

# Plotting
ggplot(results, aes(x = suffix, y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  
  # Add HR (95% CI) labels, slightly to the right, right-aligned
  geom_text(
    aes(label = hr_ci_label, y = upper_95 + 0.02),
    hjust = 0,
    size = 3.5,
    family = "mono"  # Monospaced font for neater alignment
  )+
  
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.4))) +
  labs(
    x = NULL,
    y = "Hazard Ratio of TTFD10",
    title = "TTFD10 in OS Models"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(5.5, 100, 5.5, 5.5, "pt"),
    plot.title = element_text(hjust = 0.5)
  )


####PFS
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_1022_TTFD10")

files <- list.files(pattern = "^cox_td_PFS_.*\\.rds$")

results <- data.frame(
  model_id = character(),
  exp_coef = numeric(),
  lower_95 = numeric(),
  upper_95 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (f in files) {
  model <- readRDS(f)
  sum_model <- summary(model)
  
  if ("has_TTFD10" %in% rownames(sum_model$coefficients)) {
    exp_coef <- sum_model$conf.int["has_TTFD10", "exp(coef)"]
    lower <- sum_model$conf.int["has_TTFD10", "lower .95"]
    upper <- sum_model$conf.int["has_TTFD10", "upper .95"]
    pval <- sum_model$coefficients["has_TTFD10", "Pr(>|z|)"]
    
    results <- rbind(results, data.frame(
      model_id = gsub(".rds", "", f),
      exp_coef = exp_coef,
      lower_95 = lower,
      upper_95 = upper,
      p_value = pval
    ))
  }
}

results$model_id <- factor(results$model_id, levels = results$model_id)

results <- results %>%
  mutate(
    suffix = gsub(".*_PFS_(.*?)_\\d{8}", "\\1", model_id)
  ) %>%
  arrange(p_value)  

results$suffix <- factor(results$suffix, levels = results$suffix)

results <- results %>%
  mutate(p_label = formatC(p_value, format = "e", digits = 2))

ggplot(results, aes(x = suffix, y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(
    aes(label = p_label, x = suffix, y = upper_95 + 0.015),  # Slightly shift labels to the right to avoid overlapping with error bars
    hjust = 0,
    size = 3.5
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.2))) +  # Upper y-limit 0.5 with extra padding
  labs(
    x = NULL,
    y = "Hazard Ratio of TTFD10",
    title = "TTFD10 in PFS Models"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),  # ✅ Increase right margin
    plot.title = element_text(hjust = 0.5)
  )

ggplot(results, aes(x = suffix, y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(aes(label = p_label), hjust = -0.1, size = 3.5, color = "black") +  # ✅ Add p-value labels
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +  # Allow p-value labels to extend outside the panel
  labs(
    x = NULL,
    y = "Hazard Ratio of TTFD10",
    title = "Hazard Ratios and P-values of has_TTFD10 in PFS Models"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(5.5, 30, 5.5, 5.5, "pt")  # Add right margin to show p-values
  )

ggplot(results, aes(x = suffix, y = exp_coef)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Reference line at HR = 1
  coord_flip() +
  labs(
    x = " ",
    y = "Hazard Ratio of TTFD10",
    title = " "
  ) +
  theme_minimal(base_size = 14)

# Construct display string: HR (95% CI), keep 4 decimals
results <- results %>%
  mutate(
    suffix = gsub(".*_PFS_(.*?)_\\d{8}", "\\1", model_id),
    hr_ci_label = sprintf("%.4f (%.4f, %.4f)", exp_coef, lower_95, upper_95),
    p_label = formatC(p_value, format = "e", digits = 2)
  ) %>%
  arrange(p_value)

results$suffix <- factor(results$suffix, levels = results$suffix)

# Plotting
ggplot(results, aes(x = suffix, y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  
  # Add HR (95% CI) labels, slightly to the right, right-aligned
  geom_text(
    aes(label = hr_ci_label, y = upper_95 + 0.02),
    hjust = 0,
    size = 3.5,
    family = "mono"  # Monospaced font for neater alignment
  ) +
  
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.4))) +
  labs(
    x = NULL,
    y = "Hazard Ratio of TTFD10",
    title = "TTFD10 in PFS Models"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.margin = margin(5.5, 100, 5.5, 5.5, "pt"),
    plot.title = element_text(hjust = 0.5)
  )