### updated version for 20_6_Cox_Time_denpendent_HR_plot.R

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
## A) Sort by p-value (small -> large), smallest at top
## =========================
# Add these lines before plotting
yrng <- range(c(results$lower_95, results$upper_95), na.rm = TRUE)
ypad <- if (diff(yrng) == 0) 0.1 * max(yrng, na.rm = TRUE) else 0.15 * diff(yrng)
ylims <- c(max(0, yrng[1] - ypad), yrng[2] + ypad)

# Plot (everything else unchanged; only adjust two places: geom_text y and scale_y_continuous limits)
p_plot <- ggplot(results,
                 aes(x = reorder(suffix_full, -p_value), y = exp_coef)) +
  geom_point(size = 2.5, color = "#4DBBD5") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = p_label,
                y = pmin(upper_95 + 0.05 * diff(ylims), ylims[2])),  # Adaptive position to avoid overflow
            hjust = 0, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims, expand = expansion(mult = c(0, 0.05))) +  # ★ Adaptive range
  labs(x = NULL, y = "Hazard Ratio",
       title = "Overall Survival") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(p_plot)

## =========================
## B) Sort by HR (large -> small), largest at top
## =========================
# Compute range before plotting
yrng <- range(c(results$lower_95, results$upper_95), na.rm = TRUE)
ypad <- if (diff(yrng) == 0) 0.1 * max(yrng, na.rm = TRUE) else 0.15 * diff(yrng)
ylims <- c(max(0, yrng[1] - ypad), yrng[2] + ypad)

# Plot
hr_plot <- ggplot(results,
                  aes(x = reorder(suffix_full, exp_coef),  # from large to small
                      y = exp_coef)) +
  geom_point(size = 2.5, color = "#4DBBD5") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = hr_ci_label,
                y = pmin(upper_95 + 0.05 * diff(ylims), ylims[2])), # Adaptive position
            hjust = 0, size = 3.5, family = "mono") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims,
                     expand = expansion(mult = c(0, 0.05))) +       # ★ Adaptive range
  labs(x = NULL, y = "Hazard Ratio",
       title = "Overall Survival") +
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
  geom_point(size = 2.5, color = "#E64B35") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = p_label, y = upper_95 + 0.015),
            hjust = 0, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Hazard Ratio",
       title = "Progression Free Survival") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))

print(p_plot)

## =========================
## B) Sort by HR (large -> small)
## =========================
hr_plot <- ggplot(results,
                  aes(x = reorder(suffix_full, exp_coef), y = exp_coef)) +
  geom_point(size = 2.5, color = "#E64B35") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                width = 0.2, color = "gray40") +
  geom_text(aes(label = hr_ci_label, y = upper_95 + 0.02),
            hjust = 0, size = 3.5, family = "mono") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 0.5),
                     expand = expansion(mult = c(0, 0.4))) +
  labs(x = NULL, y = "Hazard Ratio",
       title = "Progression Free Survival") +
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

# Compute adaptive range once before the two plots
yrng  <- range(c(os_res$lower_95, os_res$upper_95), na.rm = TRUE)
ypad  <- if (diff(yrng) == 0) 0.1 * max(yrng, na.rm = TRUE) else 0.15 * diff(yrng)
ylims <- c(max(0, yrng[1] - ypad), yrng[2] + ypad)

# ---- OS: sort by p-value (small->large), most significant on top ----
os_p_plot <- ggplot(os_res, aes(x = reorder(suffix_full, -p_value), y = exp_coef)) +
  geom_point(size = 2.5, color = "#4DBBD5") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(
    aes(label = p_label, y = pmin(upper_95 + 0.05 * diff(ylims), ylims[2])),
    hjust = 0, size = 3.5
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Hazard Ratio", title = "Overall Survival") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))
print(os_p_plot)

# ---- OS: sort by HR (large->small), largest on top ----
os_hr_plot <- ggplot(os_res, aes(x = reorder(suffix_full, exp_coef), y = exp_coef)) +
  geom_point(size = 2.5, color = "#4DBBD5") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(
    aes(label = hr_ci_label, y = pmin(upper_95 + 0.05 * diff(ylims), ylims[2])),
    hjust = 0, size = 3.5, family = "mono"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Hazard Ratio", title = "Overall Survival") +
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

# —— Place before the two PFS plots ——
yrng_pfs  <- range(c(pfs_res$lower_95, pfs_res$upper_95), na.rm = TRUE)
ypad_pfs  <- if (diff(yrng_pfs) == 0) 0.1 * max(yrng_pfs, na.rm = TRUE) else 0.15 * diff(yrng_pfs)
ylims_pfs <- c(max(0, yrng_pfs[1] - ypad_pfs), yrng_pfs[2] + ypad_pfs)

# ---- PFS: sort by p-value (small->large) ----
pfs_p_plot <- ggplot(pfs_res, aes(x = reorder(suffix_full, -p_value), y = exp_coef)) +
  geom_point(size = 2.5, color = "#E64B35") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(
    aes(label = p_label, y = pmin(upper_95 + 0.05 * diff(ylims_pfs), ylims_pfs[2])),
    hjust = 0, size = 3.5
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims_pfs, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Hazard Ratio", title = "Progression Free Survival") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 60, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))
print(pfs_p_plot)

# ---- PFS: sort by HR (large->small) ----
pfs_hr_plot <- ggplot(pfs_res, aes(x = reorder(suffix_full, exp_coef), y = exp_coef)) +
  geom_point(size = 2.5, color = "#E64B35") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2, color = "gray40") +
  geom_text(
    aes(label = hr_ci_label, y = pmin(upper_95 + 0.05 * diff(ylims_pfs), ylims_pfs[2])),
    hjust = 0, size = 3.5, family = "mono"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = ylims_pfs, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Hazard Ratio", title = "Progression Free Survival") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 100, 5.5, 5.5, "pt"),
        plot.title = element_text(hjust = 0.5))
print(pfs_hr_plot)