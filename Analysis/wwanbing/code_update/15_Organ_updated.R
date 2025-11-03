### Organ: Liver, Lymph nodes, Lung, Bone, Breast, Pleura ：Pleura + Pleural effusion
rm(list=ls())
library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(haven)
library(tidyr)
library(readr)
library(ggplot2)
library(GGally)
library(ggsci)
library(ggpubr) 
library(reshape2)
library(survival)
library(survminer)
library(forcats)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_Celgene_2001_107/wwanbing/output")
tumor_2001 <- read_csv("Breast_Celgene_2001_107_tumor.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_SanofiU_2004_135/wwanbing/output")
tumor_2004 <- read_csv("Breast_SanofiU_2004_135_tumor_18MAR2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
Merged_TTFD <- read_csv("Merged_PRO_TTFD_28JUN2025.csv")


#Merge 2001 and 2004
tumor_2001$Study <- "2001"
tumor_2004$Study <- "2004"

all(colnames(tumor_2001) == colnames(tumor_2004))

tumor_2001 <- tumor_2001 %>%
  rename(UID = LUID)

tumor_2001 <- tumor_2001 %>%
  mutate(
    UID = as.character(as.numeric(UID)),
    ID = as.character(as.numeric(ID))
  )

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = as.numeric(SIZE))

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = ifelse(is.na(SIZE), -999, SIZE))


#Merge TTFD
# Select all variables starting with TTFD10, TTFD10P, censored, censored10P
vars_to_merge <- colnames(Merged_TTFD)[grepl("^(TTFD10|TTFD10P|censored|censored10P)_", colnames(Merged_TTFD))]

vars_to_merge <- c("UID", vars_to_merge)

symptoms_TTFD <- Merged_TTFD[, vars_to_merge]

# Merge into tumor_ALL
tumor_2001 <- left_join(tumor_2001, symptoms_TTFD, by = "UID")
tumor_2004 <- left_join(tumor_2004, symptoms_TTFD, by = "UID")

tumor_all <- bind_rows(tumor_2001, tumor_2004)

###Do target or non-target lesions more strongly correlate with symptom deterioration? 
tumor_all$LESTYPE <- ifelse(tumor_all$LESTYPE == 1, 1, 0)
table(tumor_all$LESTYPE)

# Step 1: Organ categorization (new variable Organ_clean)
tumor_all <- tumor_all %>%
  mutate(
    Organ_category = case_when(
      ORGAN %in% c("Liver") ~ "Liver",
      ORGAN %in% c("Lymph node", "Lymph nodes") ~ "Lymph",
      ORGAN %in% c("Lung") ~ "Lung",
      ORGAN %in% c("Bone") ~ "Bone",
      ORGAN %in% c("Breast") ~ "Breast",
      ORGAN %in% c("Pleura", "Pleural effusion") ~ "Pleura",
      TRUE ~ "Other"
    )
  )

#> table(tumor_all$Organ_category)
#Bone Breast  Liver   Lung  Lymph  Other Pleura 
#907    415   2214   1734   1564   2970    363 

# Step 2: Extract ever-metastasis status by ID (ever)
recode_organ <- function(x) {
  dplyr::case_when(
    x %in% c("Liver") ~ "Liver",
    x %in% c("Lymph node","Lymph nodes") ~ "Lymph",
    x %in% c("Lung") ~ "Lung",
    x %in% c("Bone","Bone marrow","Skull") ~ "Bone",
    x %in% c("Breast") ~ "Breast",
    x %in% c("Pleura","Pleural effusion") ~ "Pleura",
    TRUE ~ "Other"
  )
}

# Ever-by-ID organ flag matrix (seven columns)
organ_flags <- tumor_all %>%
  filter(!is.na(ORGAN)) %>%
  mutate(Organ_clean = recode_organ(ORGAN)) %>%
  distinct(ID, Organ_clean) %>%
  mutate(flag = 1L) %>%
  tidyr::pivot_wider(names_from = Organ_clean, values_from = flag, values_fill = 0L)

# Vector of organ columns to be used
organs <- c("Liver","Lymph","Lung","Bone","Breast","Pleura","Other")

# Step 3: Merge organ flags back to TTFD data
ttfd_cols <- grep("^TTFD10_", names(tumor_all), value = TRUE)
censor_cols <- grep("^censored_", names(tumor_all), value = TRUE)

# Use AP as example
symptom <- "AP"
ttfd_col <- paste0("TTFD10_", symptom)
censored_col <- paste0("censored_", symptom)

df <- tumor_all %>%
  dplyr::select(ID, UID, dplyr::all_of(ttfd_col), dplyr::all_of(censored_col)) %>%
  dplyr::left_join(organ_flags, by = "ID") %>%
  dplyr::filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
                !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)

# Safety check: ensure seven organ columns exist
for (org in c("Liver","Lymph","Lung","Bone","Breast","Pleura","Other")) {
  if (!org %in% names(df)) df[[org]] <- 0L
  df[[org]] <- ifelse(is.na(df[[org]]), 0L, as.integer(df[[org]]))
}

surv_obj <- Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
fit <- coxph(surv_obj ~ Liver + Lymph + Lung + Bone + Breast + Pleura + Other, data = df)
sm  <- summary(fit)

forest_df <- data.frame(
  Term   = rownames(sm$coefficients),
  HR     = sm$coefficients[, "exp(coef)"],
  LCL    = sm$conf.int[, "lower .95"],
  UCL    = sm$conf.int[, "upper .95"],
  p_value= sm$coefficients[, "Pr(>|z|)"],
  row.names = NULL
) %>%
  dplyr::mutate(
    Sig   = ifelse(p_value < 0.05, "Significant", "Not Significant"),
    Term  = gsub("^Other$", "Other organ", Term),
    OR_CI = sprintf("%.2f (%.2f, %.2f)", HR, LCL, UCL),
    Label = forcats::fct_reorder(Term, HR)   # One-time ordering by HR
  )

# Dynamically place right-side text + dynamic x-range
label_x_pos <- max(forest_df$UCL, na.rm = TRUE) * 1.05
x_max <- max(label_x_pos * 1.15, max(forest_df$UCL, na.rm = TRUE) * 1.1)

ggplot(forest_df, aes(x = HR, y = Label, color = Sig)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2) +
  geom_text(aes(x = label_x_pos, label = OR_CI), hjust = 0, size = 3.5, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Significant" = "#d73027", "Not Significant" = "gray40")) +
  scale_x_log10(limits = c(0.6, x_max)) +
  labs(
    title = paste0("Forest Plot: Organ Involvement and ", symptom, " Deterioration"),
    x = "Hazard Ratio (HR, log scale)", y = NULL, color = "p < 0.05"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    plot.margin = margin(t = 10, r = 140, b = 10, l = 10) # wider right margin
  ) +
  coord_cartesian(clip = "off")   # allow right-side text to overflow the panel

### 15 symptoms

plot_forest_for_symptom <- function(symptom, tumor_all, organ_flags,
                                    organs = c("Liver","Lymph","Lung","Bone","Breast","Pleura","Other")) {
  ttfd_col     <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  # Prepare data
  df <- tumor_all %>%
    dplyr::select(ID, UID, dplyr::all_of(ttfd_col), dplyr::all_of(censored_col)) %>%
    dplyr::left_join(organ_flags, by = "ID") %>%
    dplyr::filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
                  !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  # Safety check: ensure all seven organ columns exist
  for (org in organs) {
    if (!org %in% names(df)) df[[org]] <- 0L
    df[[org]] <- ifelse(is.na(df[[org]]), 0L, as.integer(df[[org]]))
  }
  
  # Cox
  surv_obj <- survival::Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
  fml <- as.formula(paste("surv_obj ~", paste(organs, collapse = " + ")))
  fit <- survival::coxph(fml, data = df)
  sm  <- summary(fit)
  
  # Summary
  forest_df <- data.frame(
    Term    = rownames(sm$coefficients),
    HR      = sm$coefficients[, "exp(coef)"],
    LCL     = sm$conf.int[,   "lower .95"],
    UCL     = sm$conf.int[,   "upper .95"],
    p_value = sm$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  ) %>%
    dplyr::mutate(
      Term  = gsub("^Other$", "Other organ", Term),
      Sig   = ifelse(p_value < 0.05, "Significant", "Not Significant"),
      OR_CI = sprintf("%.2f (%.2f, %.2f)", HR, LCL, UCL),
      Label = forcats::fct_reorder(Term, HR)   # Order by HR
    )
  
  # Dynamic right-side text position & x-axis range
  label_x_pos <- max(forest_df$UCL, na.rm = TRUE) * 1.05
  x_max <- max(label_x_pos * 1.15, max(forest_df$UCL, na.rm = TRUE) * 1.1)
  
  ggplot2::ggplot(forest_df, ggplot2::aes(x = HR, y = Label, color = Sig)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = LCL, xmax = UCL), height = 0.2) +
    ggplot2::geom_text(ggplot2::aes(x = label_x_pos, label = OR_CI),
                       hjust = 0, size = 3.5, color = "black") +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    ggplot2::scale_color_manual(values = c("Significant" = "#d73027", "Not Significant" = "gray40")) +
    ggplot2::scale_x_log10(limits = c(0.6, x_max)) +
    ggplot2::labs(
      title = paste0("Forest Plot: Organ Involvement and ", symptom, " Deterioration"),
      x = "Hazard Ratio (HR, log scale)", y = NULL, color = "p < 0.05"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 10),
      plot.margin  = ggplot2::margin(t = 10, r = 140, b = 10, l = 10)  # wider right margin
    ) +
    ggplot2::coord_cartesian(clip = "off")   # allow right-side text to overflow the panel
}

symptoms <- c("AP","CO","DI","DY","FA","FI","GHS","NV","PA","SL","PF2","RF2","EF","CF","SF")

symptom_names <- c(
  "FA"  = "Fatigue",
  "NV"  = "Nausea_Vomit",
  "PA"  = "Pain",
  "DY"  = "Dyspnea",
  "SL"  = "Insomnia",
  "AP"  = "Appetite_loss",
  "CO"  = "Constipation",
  "DI"  = "Diarrhea",
  "FI"  = "Financial",
  "PF2" = "Physical_Functioning",
  "RF2" = "Role_Functioning",
  "EF"  = "Emotional_Functioning",
  "CF"  = "Cognitive_Functioning",
  "SF"  = "Social_Functioning",
  "GHS" = "Global_Health_Status_QOL"
)

out_dir <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/Figure_summary/Organ/Forest_plot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (symptom in names(symptom_names)) {
  message("Plotting & saving: ", symptom)
  p <- plot_forest_for_symptom(symptom, tumor_all, organ_flags)
  
  # Build file name: abbreviation + full name
  file_name <- paste0(symptom, "_", symptom_names[[symptom]], ".png")
  
  ggsave(
    filename = file.path(out_dir, file_name),
    plot     = p,
    width    = 7, height = 5, dpi = 300
  )
}

############### Hazard Ratio heatmap

# 7 organs, 15 symptoms, and full-name mapping
organs <- c("Liver","Lymph","Lung","Bone","Breast","Pleura","Other")
symptoms <- c("AP","CO","DI","DY","FA","FI","GHS","NV","PA","SL","PF2","RF2","EF","CF","SF")
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
  "PF2" = "Physical Functioning",
  "RF2" = "Role Functioning",
  "EF"  = "Emotional Functioning",
  "CF"  = "Cognitive Functioning",
  "SF"  = "Social Functioning",
  "GHS" = "Global Health Status / QOL"
)

#==== 1) Single symptom: extract HR/p for 7 organs (auto-skip organs with no variation) ====
get_hr_for_symptom <- function(symptom, tumor_all, organ_flags, organs_vec = organs) {
  ttfd_col     <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  df <- tumor_all %>%
    select(ID, UID, all_of(ttfd_col), all_of(censored_col)) %>%
    left_join(organ_flags, by = "ID") %>%
    filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
           !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  # Ensure organ columns exist
  for (o in organs_vec) {
    if (!o %in% names(df)) df[[o]] <- 0L
    df[[o]] <- ifelse(is.na(df[[o]]), 0L, as.integer(df[[o]]))
  }
  
  # Keep only organs with variation
  active_orgs <- organs_vec[vapply(organs_vec, function(o) {
    s <- df[[o]]; sum(s, na.rm = TRUE) > 0 && sum(s, na.rm = TRUE) < nrow(df)
  }, logical(1))]
  
  # If none vary, return all NA
  if (length(active_orgs) == 0) {
    return(tibble(symptom = symptom, organ = organs_vec, HR = NA_real_, p = NA_real_))
  }
  
  surv_obj <- Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
  fml <- as.formula(paste("surv_obj ~", paste(active_orgs, collapse = " + ")))
  fit <- coxph(fml, data = df)
  sm  <- summary(fit)
  
  res <- tibble(
    organ = rownames(sm$coefficients),
    HR    = sm$coefficients[, "exp(coef)"],
    p     = sm$coefficients[, "Pr(>|z|)"],
  )
  
  # Fill missing organs as NA and order by organs_vec
  res_full <- tibble(organ = organs_vec) %>%
    left_join(res, by = "organ") %>%
    mutate(symptom = symptom, .before = 1)
  
  res_full
}

#==== 2) Aggregate 15×7 HR and p ====
hr_long <- map_dfr(symptoms, ~get_hr_for_symptom(.x, tumor_all, organ_flags, organs_vec = organs))

# Add symptom full names and plotting order (column order)
hr_long <- hr_long %>%
  mutate(symptom_full = symptom_names[symptom],
         symptom_full = factor(symptom_full, levels = symptom_names[symptoms]),
         organ = factor(organ, levels = organs))

#==== 3) Plot heatmap ====
# Color scheme: neutral at 0 (HR=1), <0 blue, >0 red; adjust as needed
heat <- ggplot(hr_long, aes(x = symptom_full, y = organ, fill = HR)) +
  geom_tile(color = "white") +
  # Outline cells with significance
  geom_tile(data = subset(hr_long, !is.na(p) & p < 0.05),
            fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient2(
    name = "Hazard Ratio",
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 1, na.value = "grey90"
  ) +
  labs(
    title = "Hazard Ratio Heatmap by Symptom and Organ",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(heat)

#==== 4) Save PNG / PDF ====
out_dir <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/Figure_summary/Organ/HR_heatmap"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "HR_heatmap_log2HR.png"),
       heat, width = 12, height = 5.5, dpi = 300)

ggsave(file.path(out_dir, "HR_heatmap_log2HR.pdf"),
       heat, width = 12, height = 5.5)


### Multiple-organ distribution plots
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_Celgene_2001_107/wwanbing/output")
tumor_2001 <- read_csv("Breast_Celgene_2001_107_tumor.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_SanofiU_2004_135/wwanbing/output")
tumor_2004 <- read_csv("Breast_SanofiU_2004_135_tumor_18MAR2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
Merged_TTFD <- read_csv("Merged_PRO_TTFD_28JUN2025.csv")


#Merge 2001 and 2004
tumor_2001$Study <- "2001"
tumor_2004$Study <- "2004"

all(colnames(tumor_2001) == colnames(tumor_2004))

tumor_2001 <- tumor_2001 %>%
  rename(UID = LUID)

tumor_2001 <- tumor_2001 %>%
  mutate(
    UID = as.character(as.numeric(UID)),
    ID = as.character(as.numeric(ID))
  )

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = as.numeric(SIZE))

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = ifelse(is.na(SIZE), -999, SIZE))


#Merge TTFD
# Select all variables starting with TTFD10, TTFD10P, censored, censored10P
vars_to_merge <- colnames(Merged_TTFD)[grepl("^(TTFD10|TTFD10P|censored|censored10P)_", colnames(Merged_TTFD))]

vars_to_merge <- c("UID", vars_to_merge)

symptoms_TTFD <- Merged_TTFD[, vars_to_merge]

# Merge into tumor_ALL
tumor_2001 <- left_join(tumor_2001, symptoms_TTFD, by = "UID")
tumor_2004 <- left_join(tumor_2004, symptoms_TTFD, by = "UID")

tumor_all <- bind_rows(tumor_2001, tumor_2004)

###Do target or non-target lesions more strongly correlate with symptom deterioration? 
tumor_all$LESTYPE <- ifelse(tumor_all$LESTYPE == 1, 1, 0)
table(tumor_all$LESTYPE)

# Step 1: Organ categorization (new variable Organ_clean)
tumor_all <- tumor_all %>%
  mutate(
    Organ_category = case_when(
      ORGAN %in% c("Liver") ~ "Liver",
      ORGAN %in% c("Lymph node", "Lymph nodes") ~ "Lymph",
      ORGAN %in% c("Lung") ~ "Lung",
      ORGAN %in% c("Bone") ~ "Bone",
      ORGAN %in% c("Breast") ~ "Breast",
      ORGAN %in% c("Pleura", "Pleural effusion") ~ "Pleura",
      TRUE ~ "Other"
    )
  )

#> table(tumor_all$Organ_category)
#Bone Breast  Liver   Lung  Lymph  Other Pleura 
#907    415   2214   1734   1564   2970    363 

# Step 2: Extract ever-metastasis status by ID (ever)
recode_organ <- function(x) {
  dplyr::case_when(
    x %in% c("Liver") ~ "Liver",
    x %in% c("Lymph node","Lymph nodes") ~ "Lymph",
    x %in% c("Lung") ~ "Lung",
    x %in% c("Bone","Bone marrow","Skull") ~ "Bone",
    x %in% c("Breast") ~ "Breast",
    x %in% c("Pleura","Pleural effusion") ~ "Pleura",
    TRUE ~ "Other"
  )
}

# Ever-by-ID organ flag matrix (seven columns)
organ_flags <- tumor_all %>%
  filter(!is.na(ORGAN)) %>%
  mutate(Organ_clean = recode_organ(ORGAN)) %>%
  distinct(ID, Organ_clean) %>%
  mutate(flag = 1L) %>%
  tidyr::pivot_wider(names_from = Organ_clean, values_from = flag, values_fill = 0L)

# Merge organ flags back to TTFD data
ttfd_cols <- grep("^TTFD10_", names(tumor_all), value = TRUE)
censor_cols <- grep("^censored_", names(tumor_all), value = TRUE)

### 15 symptoms
organs <- c("Liver","Lymph","Lung","Bone","Breast","Pleura","Other")

# Symptom abbreviation -> full name (for filenames)
symptom_names <- c(
  "FA" = "Fatigue", "NV" = "Nausea/Vomit", "PA" = "Pain", "DY" = "Dyspnea",
  "SL" = "Insomnia", "AP" = "Appetite loss", "CO" = "Constipation", 
  "DI" = "Diarrhea", "FI" = "Financial", "PF2" = "Physical", "RF2" = "Role",
  "EF" = "Emotional", "CF" = "Cognitive", "SF" = "Social", "GHS" = "GHS/QOL"
)

# Output directory
out_dir_km <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/Figure_summary/Organ/KM_curve"
dir.create(out_dir_km, showWarnings = FALSE, recursive = TRUE)

# ---------- Function: produce a 3×3 KM panel for a single symptom ----------
plot_km_panel_for_symptom <- function(symptom, tumor_all, organ_flags,
                                      organs_vec = organs, ncol = 3, nrow = 3) {
  ttfd_col     <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  df <- tumor_all %>%
    dplyr::select(ID, UID, dplyr::all_of(ttfd_col), dplyr::all_of(censored_col)) %>%
    dplyr::left_join(organ_flags, by = "ID") %>%
    dplyr::filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
                  !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  # Safety check: ensure each organ column exists
  for (org in organs_vec) {
    if (!org %in% names(df)) df[[org]] <- 0L
    df[[org]] <- ifelse(is.na(df[[org]]), 0L, as.integer(df[[org]]))
  }
  
  plot_list <- list()
  for (org in organs_vec) {
    tmp <- data.frame(
      time  = df[[ttfd_col]],
      event = 1 - df[[censored_col]],
      group = factor(df[[org]], levels = c(0,1), labels = c("No","Yes"))
    )
    
    if (length(unique(na.omit(tmp$group))) < 2) {
      g <- ggplot() +
        ggtitle(paste0(org, " (insufficient groups)")) +
        theme_void()
      plot_list[[org]] <- g
      next
    }
    
    surv_diff <- survdiff(Surv(time, event) ~ group, data = tmp)
    pval <- 1 - pchisq(surv_diff$chisq, df = 1)
    
    fit <- survfit(Surv(time, event) ~ group, data = tmp)
    
    g <- ggsurvplot(
      fit, data = tmp,
      legend.title = paste0(org, " involvement"),
      legend.labs  = levels(tmp$group),
      palette      = c("#4DBBD5FF", "#E64B35FF"),
      xlab = paste0("Time to ", symptom, " deterioration"),
      ylab = "Survival Probability",
      title = paste0(org, " (p = ", format.pval(pval, digits = 3), ")"),
      ggtheme = theme_minimal(base_size = 12),
      pval = FALSE
    )
    plot_list[[org]] <- g$plot
  }
  
  panel <- ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow)
  panel <- annotate_figure(
    panel,
    top = text_grob(
      paste0("KM Curves for ", symptom, " by Organ Involvement"),
      face = "bold", size = 14
    )
  )
  return(panel)
}

# ---------- Batch plotting & saving ----------
symptoms <- c("AP","CO","DI","DY","FA","FI","GHS","NV","PA","SL","PF2","RF2","EF","CF","SF")

for (symptom in symptoms) {
  message("Plotting KM panel for: ", symptom)
  panel_plot <- plot_km_panel_for_symptom(symptom, tumor_all, organ_flags,
                                          organs_vec = organs, ncol = 3, nrow = 3)
  
  file_name <- paste0(symptom, "_", symptom_names[[symptom]], "_KMpanel.png")
  ggsave(
    filename = file.path(out_dir_km, file_name),
    plot     = panel_plot,
    width    = 11, height = 9, dpi = 300
  )
}