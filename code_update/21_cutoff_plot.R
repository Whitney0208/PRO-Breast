## =========================
## Tumor size vs. symptom deterioration (EN)
## =========================

rm(list=ls())

library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(haven)
library(tidyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggsci)
library(reshape2)
library(survival)
library(survminer)
library(RColorBrewer)
library(purrr)    # used for map_dfr
library(tibble)   # used for tibble()
library(scales)   # used for scale_y_log10 labels

## Does tumor size correlate with symptom deterioration?
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/data")

tumor_all <- read_csv("PRO_TTFD_tumor05AUG2025.csv")

Monolix_2001 <- read_csv("Breast_Monolix_2001_tumor.csv")
Monolix_2004 <- read_csv("Breast_Monolix_2004_tumor.csv")

Monolix_2001$Year <- "2001"
Monolix_2004$Year <- "2004"

Monolix_2001$UID <- as.character(Monolix_2001$UID)
Monolix_2004$UID <- as.character(Monolix_2004$UID)

Monolix <- bind_rows(Monolix_2001, Monolix_2004)

PRO_2001 <- read_csv("PRO_Size_2001.csv")
PRO_2004 <- read_csv("PRO_Size_2004.csv")

PRO_2001$Year <- "2001"
PRO_2004$Year <- "2004"

PRO_2001$UID <- as.character(PRO_2001$UID)
PRO_2004$UID <- as.character(PRO_2004$UID)

PRO <- bind_rows(PRO_2001, PRO_2004)

## Left-join censored variables from tumor_all into PRO

# length(unique(PRO$UID))        # 356
# length(unique(tumor_all$UID))  # 436

censor_cols <- grep("^censored", colnames(tumor_all), value = TRUE)
censor_df <- tumor_all %>%
  select(UID, all_of(censor_cols)) %>%
  distinct()

PRO <- PRO %>%
  left_join(censor_df, by = "UID")

Monolix$UID <- as.character(Monolix$UID)

vars_to_merge <- c("UID", "AGE", "MENOS", "CHEMO", "HORMON", "ERS", "PGRS", "ECOG")

monolix_clinical <- Monolix %>%
  select(all_of(vars_to_merge)) %>%
  distinct()

PRO <- PRO %>%
  left_join(monolix_clinical, by = "UID")

# Column structure reference (truncated)
# [1] "UID" "ID" "TID" "ORGAN2" "TTFD10_AP" ... "ECOG"

## 1) Preprocessing: map censored (0=event, 1=censor) to Surv event indicator (1=event)
dat_prep <- PRO %>%
  mutate(
    ratio_val = nadir_ratio_TS_TTFD10_AP,  # change to your chosen metric if needed
    time_val  = TTFD10_AP,
    event_val = 1 - censored_AP            # key: 0->1 (event), 1->0 (censored)
  ) %>%
  drop_na(ratio_val, time_val, event_val)
# If adding covariates, also drop NAs for AGE/MENOS/... etc.

## Quick cross-check (should be 0→1 and 1→0):
with(dat_prep, table(censored_AP, event_val))
# expected:
# censored_AP  0  1
# event_val    1  0

## 2) Cox models across cutoffs (to draw HR/p-value vs. cutoff)
covs <- c("AGE","MENOS","CHEMO","HORMON","ERS","PGRS","ECOG")
cuts <- seq(1.2, 2.5, by = 0.05)

fit_one_cut <- function(cut){
  dd <- dat_prep %>%
    mutate(Ratio_cat = factor(
      ifelse(ratio_val >= cut, paste0("≥", cut), paste0("<", cut)),
      levels = c(paste0("<", cut), paste0("≥", cut))
    ))
  
  # Unadjusted model (as in your current code). To adjust, uncomment covariates line.
  # fml <- as.formula(paste0("Surv(time_val, event_val) ~ Ratio_cat + ",
  #                          paste(covs, collapse = " + ")))
  fml <- as.formula("Surv(time_val, event_val) ~ Ratio_cat")
  fit <- coxph(fml, data = dd)
  term <- paste0("Ratio_cat", "≥", cut)
  
  tibble(
    cutoff   = cut,
    HR       = exp(coef(fit)[term]),
    conf.low = exp(stats::confint(fit)[term, 1]),
    conf.high= exp(stats::confint(fit)[term, 2]),
    pval     = summary(fit)$coefficients[term, "Pr(>|z|)"],
    n_low    = sum(dd$Ratio_cat == paste0("<",cut)),
    n_high   = sum(dd$Ratio_cat == paste0("≥",cut))
  )
}
res <- map_dfr(cuts, fit_one_cut)

## 3) Plots (no patchwork)
p_hr <- ggplot(res, aes(cutoff, HR)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line() + geom_point() +
  labs(x = "Tumor nadir ratio cutoff", y = "Hazard Ratio (≥cut vs <cut)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 13, face = "bold"),
    axis.text.y  = element_text(size = 13, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )

p_p <- ggplot(res, aes(cutoff, pval)) +
  geom_hline(yintercept = 0.05, linetype = 2) +
  geom_line() + geom_point() +
  scale_y_log10(
    breaks = c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(x = "Tumor nadir ratio cutoff", y = "p-value") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 13, face = "bold"),
    axis.text.y  = element_text(size = 13, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )

print(p_hr)
print(p_p)

## 4) KM plot at cutoff = 1.85
##    Helper that saves the full figure (curve + risk table)
plot_km_cut <- function(
    cut,
    df = dat_prep,
    outdir = "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/Figure_summary/cutoff_plot"
){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  dd <- df %>%
    mutate(group = factor(ifelse(ratio_val >= cut, paste0("≥", cut), paste0("<", cut)),
                          levels = c(paste0("<", cut), paste0("≥", cut))))
  
  fit <- survfit(Surv(time_val, event_val) ~ group, data = dd)
  
  p <- ggsurvplot(
    fit, data = dd,
    pval        = TRUE,
    conf.int    = TRUE,
    risk.table  = TRUE,              # show risk table
    risk.table.height = 0.25,
    risk.table.y.text = TRUE,
    risk.table.title  = "No. at risk",
    tables.theme = theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold")
      ),
    legend.title = "Group",
    legend.labs  = c(paste0("<", cut), paste0("≥", cut)),
    xlab = "Time to appetite deterioration (days)",
    ylab = "Probability",
    palette = c("#E69F00", "#56B4E9"),
    ggtheme = theme_bw(base_size = 14),
    xlim = c(0, 300),
    break.time.by = 100
  )
  
  ## Bold axis ticks/titles on the KM panel
  p$plot <- p$plot +
    theme(
      axis.text.x  = element_text(size = 13, face = "bold"),
      axis.text.y  = element_text(size = 13, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold")
    )
  
  ## Save the combined figure (curve + risk table) using ggexport
  outfile <- file.path(outdir, paste0("KM_cut_", cut, ".tiff"))
  ggexport(
    plotlist = list(p$plot, p$table + theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))),
    ncol = 1, nrow = 2,
    filename = outfile,
    width = 2000, height = 1600, res = 300
  )
  message("Saved: ", outfile)
  
  ## Also print to viewer
  print(p)
  invisible(p)
}

plot_km_cut(1.85)