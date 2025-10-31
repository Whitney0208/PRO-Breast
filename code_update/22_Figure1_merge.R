### updated version for 20_11_Cox_Time_dependent_HR_plot.R

### TTFD10

### OS
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
library(patchwork)
library(purrr)
library(tibble)

# setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10")
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_1001_TTFD10")

files <- list.files(pattern = "^cox_td_OS_.*\\.rds$")

bold_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(5.5, 100, 5.5, 5.5, "pt")
  )

# Collect model results -------------------------------------------------------
results <- purrr::map_dfr(files, function(f){
  m <- readRDS(f)
  sm <- summary(m)
  if ("has_TTFD10" %in% rownames(sm$coefficients)) {
    tibble(
      model_id = sub("\\.rds$", "", f),
      exp_coef = sm$conf.int["has_TTFD10", "exp(coef)"],
      lower_95 = sm$conf.int["has_TTFD10", "lower .95"],
      upper_95 = sm$conf.int["has_TTFD10", "upper .95"],
      p_value  = sm$coefficients["has_TTFD10", "Pr(>|z|)"]
    )
  } else {
    tibble()  # return empty tibble (safer than NULL for map_dfr)
  }
})

# Extract suffix and map to full symptom names --------------------------------
suffix_map <- c(
  FA="Fatigue", NV="Nausea/Vomit", PA="Pain", DY="Dyspnea", SL="Insomnia",
  AP="Appetite loss", CO="Constipation", DI="Diarrhea", FI="Financial",
  PF2="Physical Functioning", RF2="Role Functioning", EF="Emotional Functioning",
  CF="Cognitive Functioning", SF="Social Functioning", GHS="Global Health Status / QOL"
)

results <- results %>%
  mutate(
    suffix = sub(".*_OS_(.*?)_\\d{8}", "\\1", model_id),
    suffix = as.character(suffix),
    suffix_full = suffix_map[suffix],
    p_label = formatC(p_value, format = "e", digits = 2),
    hr_ci_label = sprintf("%.4f (%.4f, %.4f)", exp_coef, lower_95, upper_95)
  )

## =========================
## A) Order by p-value (ascending): smaller p at the top
## =========================
# Compute adaptive y-limits to avoid clipping labels
yrng <- range(c(results$lower_95, results$upper_95), na.rm = TRUE)
ypad <- if (diff(yrng) == 0) 0.1 * max(yrng, na.rm = TRUE) else 0.15 * diff(yrng)
ylims <- c(max(0, yrng[1] - ypad), yrng[2] + ypad)

p_plot <- ggplot(results,
                 aes(x = reorder(suffix_full, -p_value), y = exp_coef)) +
  geom_point(size = 2.5, color = "#4DBBD5") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = p_label,
                y = pmin(upper_95 + 0.05 * diff(ylims), ylims[2])),
            hjust = 0, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Hazard Ratio", title = "Overall Survival") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(p_plot)

## =========================
## B) Order by HR (descending): larger HR at the top
## =========================
yrng <- range(c(results$lower_95, results$upper_95), na.rm = TRUE)
ypad <- if (diff(yrng) == 0) 0.1 * max(yrng, na.rm = TRUE) else 0.15 * diff(yrng)
ylims <- c(max(0, yrng[1] - ypad), yrng[2] + ypad)

hr_plot <- ggplot(results,
                  aes(x = reorder(suffix_full, exp_coef), y = exp_coef)) +
  geom_point(size = 2.5, color = "#4DBBD5") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = hr_ci_label,
                y = pmin(upper_95 + 0.05 * diff(ylims), ylims[2])),
            hjust = 0, size = 3.5, family = "mono") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Hazard Ratio", title = "Overall Survival") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 100, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(hr_plot)

hr_plot_OS <- hr_plot + bold_theme

### PFS
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_1022_TTFD10")
files <- list.files(pattern = "^cox_td_PFS_.*\\.rds$")

results <- purrr::map_dfr(files, function(f){
  model <- readRDS(f)
  sm <- summary(model)
  if ("has_TTFD10" %in% rownames(sm$coefficients)) {
    tibble(
      model_id = sub("\\.rds$", "", f),
      exp_coef = sm$conf.int["has_TTFD10", "exp(coef)"],
      lower_95 = sm$conf.int["has_TTFD10", "lower .95"],
      upper_95 = sm$conf.int["has_TTFD10", "upper .95"],
      p_value  = sm$coefficients["has_TTFD10", "Pr(>|z|)"]
    )
  } else {
    tibble()
  }
})

# Map suffix to full symptom names -------------------------------------------
suffix_map <- c(
  FA = "Fatigue", NV = "Nausea/Vomit", PA = "Pain", DY = "Dyspnea", SL = "Insomnia",
  AP = "Appetite loss", CO = "Constipation", DI = "Diarrhea", FI = "Financial",
  PF2 = "Physical Functioning", RF2 = "Role Functioning", EF = "Emotional Functioning",
  CF = "Cognitive Functioning", SF = "Social Functioning", GHS = "Global Health Status / QOL"
)

results <- results %>%
  mutate(
    suffix = sub(".*_PFS_(.*?)_\\d{8}", "\\1", model_id),
    suffix = as.character(suffix),
    suffix_full = suffix_map[suffix],
    p_label = formatC(p_value, format = "e", digits = 2),
    hr_ci_label = sprintf("%.4f (%.4f, %.4f)", exp_coef, lower_95, upper_95)
  )

## =========================
## A) Order by p-value (ascending)
## =========================
# Adaptive y-limits (same approach as OS, avoids clipping)
yrng <- range(c(results$lower_95, results$upper_95), na.rm = TRUE)
ypad <- if (diff(yrng) == 0) 0.1 * max(yrng, na.rm = TRUE) else 0.15 * diff(yrng)
ylims <- c(max(0, yrng[1] - ypad), yrng[2] + ypad)

p_plot <- ggplot(results,
                 aes(x = reorder(suffix_full, -p_value), y = exp_coef)) +
  geom_point(size = 2.5, color = "#E64B35") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = p_label,
                y = pmin(upper_95 + 0.05 * diff(ylims), ylims[2])),
            hjust = 0, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Hazard Ratio",
       title = "Progression Free Survival") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(p_plot)

## =========================
## B) Order by HR (descending)
## =========================
yrng <- range(c(results$lower_95, results$upper_95), na.rm = TRUE)
ypad <- if (diff(yrng) == 0) 0.1 * max(yrng, na.rm = TRUE) else 0.15 * diff(yrng)
ylims <- c(max(0, yrng[1] - ypad), yrng[2] + ypad)

hr_plot <- ggplot(results,
                  aes(x = reorder(suffix_full, exp_coef), y = exp_coef)) +
  geom_point(size = 2.5, color = "#E64B35") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = hr_ci_label,
                y = pmin(upper_95 + 0.05 * diff(ylims), ylims[2])),
            hjust = 0, size = 3.5, family = "mono") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Hazard Ratio",
       title = "Progression Free Survival") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 100, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(hr_plot)

hr_plot_PFS <- hr_plot + bold_theme
hr_plot_PFS

# Stack OS over PFS with panel tags a/b
final_plot <- hr_plot_OS / hr_plot_PFS + plot_annotation(tag_levels = "a")

# Show final composed figure
print(final_plot)