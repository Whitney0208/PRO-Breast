## =========================
## 19_tumor_size_analysis.R (EN version)
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
library(forcats)

## Do tumor size and its change correlate with symptom deterioration?
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/data")

tumor_all <- read_csv("PRO_TTFD_tumor_20251015.csv")

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

# Left-join censored variables from tumor_all into PRO
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

# Column structure check (reference)
# [1] "UID" "ID" "TID" "ORGAN2" "TTFD10_AP" "TTFD10P_AP" "TTFD10_CO" "TTFD10P_CO"
# ... (omitted for brevity; unchanged from your listing)

## =========================
## Hazard ratio plot (single-run version)
## =========================
symptoms <- c("AP", "CO", "DI", "DY", "FA", "FI", "GHS", "NV",
              "PA", "SL", "PF2", "RF2", "EF", "CF", "SF")

symptom_names <- c(
  "FA"  = "Fatigue",
  "NV"  = "Nausea/Vomit",
  "PA"  = "Pain",
  "DY"  = "Dyspnea",
  "SL"  = "Insomnia",
  "AP"  = "Appetite loss",
  "CO"  = "Constipation",
  "DI"  = "Diarrhea",
  "FI"  = "Financial",
  "PF2" = "Physical",
  "RF2" = "Role",
  "EF"  = "Emotional",
  "CF"  = "Cognitive",
  "SF"  = "Social",
  "GHS" = "GHS/QOL"
)

results <- data.frame()

for (sym in symptoms) {
  ttfd_col <- paste0("TTFD10_", sym)
  cens_col <- paste0("censored_", sym)
  # If using nadir ratio instead, swap the next line:
  # ratio_col <- paste0("nadir_ratio_TS_TTFD10_", sym)
  ratio_col <- paste0("TS_TTFD10_", sym)
  
  dat <- PRO %>%
    filter(!is.na(.data[[ttfd_col]]),
           !is.na(.data[[ratio_col]]),
           !is.na(AGE),
           !is.na(ECOG),
           .data[[ttfd_col]] != -999,
           .data[[cens_col]] != -999,
           .data[[ratio_col]] != -999,
           AGE != -999,
           !MENOS %in% "UNK",
           !CHEMO %in% "UNK",
           !HORMON %in% "UNK",
           !ERS %in% "UNK",
           !PGRS %in% "UNK",
           !ECOG %in% "UNK") %>%
    mutate(event = ifelse(.data[[cens_col]] == 0, 1, 0))
  
  if (nrow(dat) >= 30) {
    formula_str <- as.formula(paste0(
      "Surv(", ttfd_col, ", event) ~ ", ratio_col,
      " + AGE + MENOS + CHEMO + HORMON + ERS + PGRS + ECOG"
    ))
    
    cox_model <- coxph(formula_str, data = dat)
    sum_model <- summary(cox_model)
    
    hr    <- sum_model$coefficients[1, "exp(coef)"]
    lower <- sum_model$conf.int[1, "lower .95"]
    upper <- sum_model$conf.int[1, "upper .95"]
    pval  <- sum_model$coefficients[1, "Pr(>|z|)"]
    
    results <- rbind(results, data.frame(
      Symptom  = sym,
      HR       = hr,
      Lower_CI = lower,
      Upper_CI = upper,
      P_value  = pval
    ))
  }
}

results$Symptom <- symptom_names[results$Symptom]

results <- results %>%
  arrange(HR) %>%
  mutate(Symptom = fct_inorder(Symptom),
         Significance = ifelse(P_value < 0.05, "p < 0.05", "ns"),
         Label = sprintf("%.2f (%.2f–%.2f)", HR, Lower_CI, Upper_CI))

ggplot(results, aes(x = Symptom, y = HR, ymin = Lower_CI, ymax = Upper_CI)) +
  geom_pointrange(aes(color = Significance), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  geom_text(aes(label = Label), hjust = -0.1, size = 4) +   # Right-side HR (95% CI)
  scale_color_manual(values = c("p < 0.05" = "red", "ns" = "gray40")) +
  labs(
    title = "Hazard Ratio of Tumor Growth for Symptom Deterioration",
    y = "Hazard Ratio (HR)", x = "Symptom",
    color = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 80, 5.5, 5.5))

## =========================
## Updated forest plot (cleaned labels & stars)
## =========================

# Abbrev -> Full name (duplicate kept for clarity in this block)
symptom_names <- c(
  "FA"="Fatigue","NV"="Nausea/Vomit","PA"="Pain","DY"="Dyspnea",
  "SL"="Insomnia","AP"="Appetite loss","CO"="Constipation","DI"="Diarrhea",
  "FI"="Financial","PF2"="Physical","RF2"="Role","EF"="Emotional",
  "CF"="Cognitive","SF"="Social","GHS"="GHS/QOL"
)

# Key change: as.character() to avoid factor misalignment after mapping
results_viz <- results %>%
  mutate(
    Symptom_chr = as.character(Symptom),
    Symptom_full = dplyr::coalesce(               # If already full name, keep as-is
      unname(symptom_names[Symptom_chr]),         # Map by name index
      Symptom_chr
    ),
    SigFlag = ifelse(P_value < 0.05, "p < 0.05", "ns"),
    Stars   = case_when(P_value < 0.001 ~ "***",
                        P_value < 0.01  ~ "**",
                        P_value < 0.05  ~ "*",
                        TRUE ~ ""),
    Label   = sprintf("%.2f (%.2f–%.2f)%s",
                      HR, Lower_CI, Upper_CI,
                      ifelse(Stars=="","", paste0("  ", Stars)))
  ) %>%
  select(-Symptom, -Symptom_chr) %>%
  rename(Symptom = Symptom_full) %>%
  mutate(Symptom = fct_reorder(Symptom, HR, .desc = FALSE))  # Larger HRs at the top

# Colors
col_sig <- "#D81B60"; col_ns <- "#6B7280"

ggplot(results_viz, aes(y = Symptom, x = HR, color = SigFlag)) +
  geom_vline(xintercept = 1, linetype = "11", linewidth = 0.6, color = "#9aa0a6") +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI),
                 height = 0, linewidth = 1.1, alpha = 0.95) +
  geom_point(size = 3.5) +
  geom_text(aes(label = Label, x = pmax(Upper_CI, HR)),
            hjust = -0.1, size = 4.2, color = "black") +
  scale_color_manual(values = c("p < 0.05" = col_sig, "ns" = col_ns),
                     name = "Significance") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(title = "Hazard Ratio of Tumor Growth for Symptom Deterioration",
       x = "Hazard Ratio (HR)", y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    plot.title.position = "plot",
    axis.line.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    panel.grid.major.x = element_line(color = "#eeeeee"),
    plot.margin = margin(10, 120, 10, 10)
  ) +
  coord_cartesian(clip = "off")

## =========================
## Extended analysis: add organ involvement flags into the Cox models
## =========================

# Organ mapping
organ_map <- function(x){
  dplyr::case_when(
    x %in% c("Liver") ~ "Liver",
    x %in% c("Lymph node","Lymph nodes") ~ "Lymph",
    x %in% c("Lung") ~ "Lung",
    x %in% c("Bone") ~ "Bone",
    x %in% c("Breast") ~ "Breast",
    x %in% c("Pleura","Pleural effusion") ~ "Pleura",
    TRUE ~ "Other"
  )
}

# For each UID, create 0/1 flags for organ involvement (a UID can have multiple organs)
organ_flags <- tumor_all %>%
  dplyr::mutate(Organ_category = organ_map(ORGAN)) %>%
  dplyr::select(UID, Organ_category) %>%
  dplyr::distinct() %>%
  dplyr::mutate(flag = 1L) %>%
  tidyr::pivot_wider(
    names_from  = Organ_category,
    values_from = flag,
    values_fill = 0L
  )

# Include "Other" as well
major_organs <- c("Liver","Lymph","Lung","Bone","Breast","Pleura","Other")

# Merge organ flags into PRO (create missing columns as 0)
for (oc in major_organs) {
  if (!oc %in% names(organ_flags)) organ_flags[[oc]] <- 0L
}
PRO <- PRO %>%
  dplyr::left_join(organ_flags %>% dplyr::select(UID, dplyr::all_of(major_organs)), by = "UID")

## =========================
## Define symptoms and display names
## =========================
symptoms <- c("AP","CO","DI","DY","FA","FI","GHS","NV","PA","SL","PF2","RF2","EF","CF","SF")
symptom_names <- c(
  FA="Fatigue", NV="Nausea/Vomit", PA="Pain", DY="Dyspnea", SL="Insomnia",
  AP="Appetite loss", CO="Constipation", DI="Diarrhea", FI="Financial",
  PF2="Physical Functioning", RF2="Role Functioning", EF="Emotional Functioning",
  CF="Cognitive Functioning", SF="Social Functioning", GHS="Global Health Status / QOL"
)

## =========================
## Cox loop: include organ flags (0/1) with clinical covariates in the model
## =========================
results_symptom <- data.frame()  # A) HR of the tumor-size ratio per symptom
results_organ   <- data.frame()  # B) HR of organ flags per symptom

for (sym in symptoms) {
  ttfd_col  <- paste0("TTFD10_", sym)
  cens_col  <- paste0("censored_", sym)
  ratio_col <- paste0("Ratio_TS_TTFD10_", sym)
  
  dat <- PRO %>%
    dplyr::filter(
      !is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
      !is.na(.data[[cens_col]]), .data[[cens_col]] != -999,
      !is.na(.data[[ratio_col]]), .data[[ratio_col]] != -999,
      !is.na(AGE), AGE != -999,
      !MENOS  %in% "UNK",
      !CHEMO  %in% "UNK",
      !HORMON %in% "UNK",
      !ERS    %in% "UNK",
      !PGRS   %in% "UNK",
      !ECOG   %in% "UNK"
    ) %>%
    dplyr::mutate(event = ifelse(.data[[cens_col]] == 0, 1, 0))
  
  # Fill missing organ columns with 0 (including "Other")
  for (oc in major_organs) {
    if (!oc %in% names(dat)) dat[[oc]] <- 0L
    dat[[oc]][is.na(dat[[oc]])] <- 0L
  }
  
  if (nrow(dat) < 30) next
  
  # Model formula: ratio + 7 organ flags (include clinical covariates here if desired)
  rhs <- paste(
    ratio_col,
    # "AGE", "MENOS", "CHEMO", "HORMON", "ERS", "PGRS", "ECOG",
    paste(major_organs, collapse = " + "),
    sep = " + "
  )
  fml <- as.formula(paste0("Surv(", ttfd_col, ", event) ~ ", rhs))
  
  cox_model <- survival::coxph(fml, data = dat)
  sum_model <- summary(cox_model)
  
  ## ------ A) HR for the tumor-size ratio (safer by name lookup) ------ ##
  if (ratio_col %in% rownames(sum_model$conf.int)) {
    results_symptom <- rbind(
      results_symptom,
      data.frame(
        Symptom  = sym,
        HR       = sum_model$conf.int[ratio_col, "exp(coef)"],
        Lower_CI = sum_model$conf.int[ratio_col, "lower .95"],
        Upper_CI = sum_model$conf.int[ratio_col, "upper .95"],
        P_value  = sum_model$coefficients[ratio_col, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    )
  }
  
  ## ------ B) HRs for organ flags (including "Other") ------ ##
  organ_rows <- intersect(major_organs, rownames(sum_model$conf.int))
  if (length(organ_rows) > 0) {
    organ_df <- data.frame(
      Symptom = sym,
      Organ   = organ_rows,
      HR      = sum_model$conf.int[organ_rows, "exp(coef)"],
      Lower_CI= sum_model$conf.int[organ_rows, "lower .95"],
      Upper_CI= sum_model$conf.int[organ_rows, "upper .95"],
      P_value = sum_model$coefficients[organ_rows, "Pr(>|z|)"],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    results_organ <- rbind(results_organ, organ_df)
  }
}

## =========================
## Forest plot for symptom-level HRs
## =========================
if (!"forcats" %in% .packages()) library(forcats)

results_symptom$Symptom <- symptom_names[results_symptom$Symptom]
results_symptom <- results_symptom %>%
  dplyr::arrange(HR) %>%
  dplyr::mutate(
    Symptom = forcats::fct_inorder(Symptom),
    SigFlag = ifelse(P_value < 0.05, "p < 0.05", "ns"),
    Stars   = dplyr::case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Label   = sprintf("%.2f (%.2f–%.2f)%s", HR, Lower_CI, Upper_CI,
                      ifelse(Stars=="","", paste0("  ", Stars)))
  )

col_sig <- "#D81B60"  # stronger magenta
col_ns  <- "#374151"  # dark gray

p_forest <- ggplot2::ggplot(results_symptom, aes(y = Symptom, x = HR, color = SigFlag)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "11", linewidth = 0.7, color = "#9aa0a6") +
  ggplot2::geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0, linewidth = 1.1) +
  ggplot2::geom_point(size = 3.6) +
  ggplot2::geom_text(aes(label = Label, x = pmax(Upper_CI, HR)),
                     hjust = -0.10, size = 4.2, color = "black") +
  ggplot2::scale_color_manual(values = c("p < 0.05" = col_sig, "ns" = col_ns), name = "Significance") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.22))) +
  ggplot2::labs(title = " ",
                x = "Hazard Ratio", y = NULL) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title.position = "plot",
    axis.line.y   = ggplot2::element_blank(),
    axis.ticks.y  = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(color = "#eeeeee"),
    plot.margin = ggplot2::margin(10, 120, 10, 10),
    axis.text.x = ggplot2::element_text(size = 13, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 13, face = "bold")
  ) +
  ggplot2::coord_cartesian(clip = "off")

print(p_forest)