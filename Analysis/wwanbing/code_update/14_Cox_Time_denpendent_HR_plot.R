#Updated for R script 15-16

### TTFD10

###OS
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

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10")

files <- list.files(pattern = "^cox_td_OS_.*\\.rds$")

# Collect results --------------------------------------------------------------
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
    NULL
  }
})

# Extract suffix + map full names ----------------------------------------------------
suffix_map <- c(
  FA="Fatigue", NV="Nausea/Vomit", PA="Pain", DY="Dyspnea", SL="Insomnia",
  AP="Appetite loss", CO="Constipation", DI="Diarrhea", FI="Financial",
  PF2="Physical Functioning", RF2="Role Functioning", EF="Emotional Functioning",
  CF="Cognitive Functioning", SF="Social Functioning", GHS="Global Health Status / QOL"
)

results <- results %>%
  mutate(
    suffix = sub(".*_OS_(.*?)_\\d{8}", "\\1", model_id),
    suffix = as.character(suffix),                     # Key: convert to character before mapping
    suffix_full = suffix_map[suffix],
    p_label = formatC(p_value, format = "e", digits = 2),
    hr_ci_label = sprintf("%.4f (%.4f, %.4f)", exp_coef, lower_95, upper_95)
  )

## =========================
## A) Sort by p-value (small -> large), smallest at the top
## =========================
p_plot <- ggplot(results,
                 aes(x = reorder(suffix_full, -p_value),  # minus ensures smallest p on top
                     y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = p_label, y = upper_95 + 0.015),
            hjust = 0, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Hazard Ratio of TTFD10",
       title = "TTFD10 in OS Models (sorted by p-value)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(p_plot)

## =========================
## B) Sort by HR (large -> small), largest at the top
## =========================
hr_plot <- ggplot(results,
                  aes(x = reorder(suffix_full, exp_coef),  # minus = descending
                      y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = hr_ci_label, y = upper_95 + 0.02),
            hjust = 0, size = 3.5, family = "mono") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5),
                     expand = expansion(mult = c(0, 0.4))) +
  labs(x = NULL, y = "Hazard Ratio of TTFD10",
       title = "TTFD10 in OS Models (sorted by HR)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 100, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(hr_plot)

### PFS
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10")

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
  }
})

# Map symptom full names ----------------------------------------------------------
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
## A) Sort by p-value (small -> large)
## =========================
p_plot <- ggplot(results,
                 aes(x = reorder(suffix_full, -p_value), y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = p_label, y = upper_95 + 0.015),
            hjust = 0, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Hazard Ratio of TTFD10",
       title = "TTFD10 in PFS Models (sorted by p-value)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(p_plot)

## =========================
## B) Sort by HR (large -> small)
## =========================
hr_plot <- ggplot(results,
                  aes(x = reorder(suffix_full, exp_coef), y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = hr_ci_label, y = upper_95 + 0.02),
            hjust = 0, size = 3.5, family = "mono") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5),
                     expand = expansion(mult = c(0, 0.4))) +
  labs(x = NULL, y = "Hazard Ratio of TTFD10",
       title = "TTFD10 in PFS Models (sorted by HR)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 100, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(hr_plot)


### TTFD10P
# ============== Common: mapping table ==============
suffix_map <- c(
  FA="Fatigue", NV="Nausea/Vomit", PA="Pain", DY="Dyspnea", SL="Insomnia",
  AP="Appetite loss", CO="Constipation", DI="Diarrhea", FI="Financial",
  PF2="Physical Functioning", RF2="Role Functioning", EF="Emotional Functioning",
  CF="Cognitive Functioning", SF="Social Functioning", GHS="Global Health Status / QOL"
)

# ============== Common: read & tidy function ==============
read_td <- function(dir_path, pattern, term){  # term: "has_TTFD10P"
  setwd(dir_path)
  files <- list.files(pattern = pattern)
  
  purrr::map_dfr(files, function(f){
    m  <- readRDS(f)
    sm <- summary(m)
    if (term %in% rownames(sm$coefficients)) {
      tibble(
        model_id = sub("\\.rds$", "", f),
        exp_coef = sm$conf.int[term, "exp(coef)"],
        lower_95 = sm$conf.int[term, "lower .95"],
        upper_95 = sm$conf.int[term, "upper .95"],
        p_value  = sm$coefficients[term, "Pr(>|z|)"]
      )
    }
  })
}

# ============== 1) OS models ==============
os_dir <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10P"
os_res <- read_td(os_dir, "^cox_td_OS_.*\\.rds$", "has_TTFD10P") %>%
  mutate(
    suffix = sub(".*_OS_(.*?)_\\d{8}", "\\1", model_id),
    suffix = as.character(suffix),                     # Key: convert to character before mapping
    suffix_full = suffix_map[suffix],
    p_label = formatC(p_value, format = "e", digits = 2),
    hr_ci_label = sprintf("%.4f (%.4f, %.4f)", exp_coef, lower_95, upper_95)
  )

# ---- OS: sort by p-value (small→large), most significant at the top ----
os_p_plot <- ggplot(os_res, aes(x = reorder(suffix_full, -p_value), y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(aes(label = p_label, y = upper_95 + 0.015), hjust = 0, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.3), expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Hazard Ratio of TTFD10P", title = "TTFD10P in OS Models (sorted by p-value)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))
print(os_p_plot)

# ---- OS: sort by HR (large→small), largest at the top ----
os_hr_plot <- ggplot(os_res, aes(x = reorder(suffix_full, exp_coef), y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(aes(label = hr_ci_label, y = upper_95 + 0.02), hjust = 0, size = 3.5, family = "mono") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.3), expand = expansion(mult = c(0, 0.4))) +
  labs(x = NULL, y = "Hazard Ratio of TTFD10P", title = "TTFD10P in OS Models (sorted by HR)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 100, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))
print(os_hr_plot)

# ============== 2) PFS models ==============
pfs_dir <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10P"
pfs_res <- read_td(pfs_dir, "^cox_td_PFS_.*\\.rds$", "has_TTFD10P") %>%
  mutate(
    suffix = sub(".*_PFS_(.*?)_\\d{8}", "\\1", model_id),
    suffix = as.character(suffix),
    suffix_full = suffix_map[suffix],
    p_label = formatC(p_value, format = "e", digits = 2),
    hr_ci_label = sprintf("%.4f (%.4f, %.4f)", exp_coef, lower_95, upper_95)
  )

# ---- PFS: sort by p-value (small→large)----
pfs_p_plot <- ggplot(pfs_res, aes(x = reorder(suffix_full, -p_value), y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(aes(label = p_label, y = upper_95 + 0.015), hjust = 0, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.3), expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Hazard Ratio of TTFD10P", title = "TTFD10P in PFS Models (sorted by p-value)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))
print(pfs_p_plot)

# ---- PFS: sort by HR (large→small)----
pfs_hr_plot <- ggplot(pfs_res, aes(x = reorder(suffix_full, exp_coef), y = exp_coef)) +
  geom_point(size = 2.5, color = "steelblue") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(aes(label = hr_ci_label, y = upper_95 + 0.02), hjust = 0, size = 3.5, family = "mono") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.3), expand = expansion(mult = c(0, 0.4))) +
  labs(x = NULL, y = "Hazard Ratio of TTFD10P", title = "TTFD10P in PFS Models (sorted by HR)") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 100, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))
print(pfs_hr_plot)