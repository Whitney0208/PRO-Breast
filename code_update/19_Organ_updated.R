### Organ: Liver, Lymph nodes, Lung, Bone, Breast, Pleura ：Pleura + Pleural effusion

### update code for organ HR heatmap

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
Merged_TTFD <- read_csv("Merged_PRO_TTFD_24SEP2025.csv")


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

# Step 1: Organ classification
tumor_all <- tumor_all %>%
  mutate(
    Organ_category = case_when(
      ORGAN %in% c("Liver") ~ "Liver",
      ORGAN %in% c("Lymph node", "Lymph nodes") ~ "Lymph nodes",
      ORGAN %in% c("Lung") ~ "Lung",
      ORGAN %in% c("Bone") ~ "Bone",
      ORGAN %in% c("Breast") ~ "Breast",
      ORGAN %in% c("Pleura", "Pleural effusion") ~ "Pleura",
      TRUE ~ "Other"
    )
  )

#> table(tumor_all$Organ_category)
#Bone Breast  Liver   Lung  Lymph nodes  Other Pleura 
#907    415   2214   1734     1564       2970    363 

# Step 2: Extract whether there has ever been an organ metastasis (ever) by ID
recode_organ <- function(x) {
  dplyr::case_when(
    x %in% c("Liver") ~ "Liver",
    x %in% c("Lymph node","Lymph nodes") ~ "Lymph nodes",
    x %in% c("Lung") ~ "Lung",
    x %in% c("Bone","Bone marrow","Skull") ~ "Bone",
    x %in% c("Breast") ~ "Breast",
    x %in% c("Pleura","Pleural effusion") ~ "Pleura",
    TRUE ~ "Other"
  )
}

# The organ marker matrix of Ever-by-ID (seven columns)
organ_flags <- tumor_all %>%
  filter(!is.na(ORGAN)) %>%
  mutate(Organ_clean = recode_organ(ORGAN)) %>%
  distinct(ID, Organ_clean) %>%
  mutate(flag = 1L) %>%
  tidyr::pivot_wider(names_from = Organ_clean, values_from = flag, values_fill = 0L)

# The column vectors of the organs needed
organs <- c("Liver","Lymph nodes","Lung","Bone","Breast","Pleura","Other")

# Step 3: Merge the TTFD data of organ markers
ttfd_cols <- grep("^TTFD10_", names(tumor_all), value = TRUE)
censor_cols <- grep("^censored_", names(tumor_all), value = TRUE)

### 15 symptoms

plot_forest_for_symptom <- function(symptom, tumor_all, organ_flags,
                                    organs = c("Liver","Lymph nodes","Lung","Bone","Breast","Pleura","Other")) {
  ttfd_col     <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  df <- tumor_all %>%
    dplyr::select(ID, UID, dplyr::all_of(ttfd_col), dplyr::all_of(censored_col)) %>%
    dplyr::left_join(organ_flags, by = "ID") %>%
    dplyr::filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
                  !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  for (org in organs) {
    if (!org %in% names(df)) df[[org]] <- 0L
    df[[org]] <- ifelse(is.na(df[[org]]), 0L, as.integer(df[[org]]))
  }
  
  # Cox
  surv_obj <- survival::Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
  
  # Put backticks around column names to avoid formula parsing failures caused by Spaces or special characters
  org_terms <- paste0("`", organs, "`")
  fml <- as.formula(paste("surv_obj ~", paste(org_terms, collapse = " + ")))
  
  fit <- survival::coxph(fml, data = df)
  sm  <- summary(fit)
  
  forest_df <- data.frame(
    Term    = rownames(sm$coefficients),
    HR      = sm$coefficients[, "exp(coef)"],
    LCL     = sm$conf.int[,   "lower .95"],
    UCL     = sm$conf.int[,   "upper .95"],
    p_value = sm$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  ) %>%
    dplyr::mutate(
      Term  = gsub("`", "", Term),                 # Remove the backticks
      Term  = gsub("^Other$", "Other organ", Term),
      Sig   = ifelse(p_value < 0.05, "Significant", "Not Significant"),
      OR_CI = sprintf("%.2f (%.2f, %.2f)", HR, LCL, UCL),
      Label = forcats::fct_reorder(Term, HR)
    )
  
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
      plot.margin  = ggplot2::margin(t = 10, r = 140, b = 10, l = 10) 
    ) +
    ggplot2::coord_cartesian(clip = "off")   
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
  "PF2" = "Physical Functioning",
  "RF2" = "Role Functioning",
  "EF"  = "Emotional Functioning",
  "CF"  = "Cognitive Functioning",
  "SF"  = "Social Functioning",
  "GHS" = "Global Health Status / QOL"
)

out_dir <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/Figure_summary/Organ/Forest_plot"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (symptom in names(symptom_names)) {
  message("Plotting & saving: ", symptom)
  p <- plot_forest_for_symptom(symptom, tumor_all, organ_flags)
  
  file_name <- paste0(symptom, "_", symptom_names[[symptom]], ".png")
  
  ggsave(
    filename = file.path(out_dir, file_name),
    plot     = p,
    width    = 7, height = 5, dpi = 300
  )
}

############### Hazard Ratio heatmap
organs <- c("Liver","Lymph nodes","Lung","Bone","Breast","Pleura","Other")
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

#==== 1) Single symptom: Extract HR/p values from 7 organs (automatically ignore organs without variations) ====
get_hr_for_symptom <- function(symptom, tumor_all, organ_flags, organs_vec = organs) {
  ttfd_col     <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  df <- tumor_all %>%
    select(ID, UID, all_of(ttfd_col), all_of(censored_col)) %>%
    left_join(organ_flags, by = "ID") %>%
    filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
           !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  for (o in organs_vec) {
    if (!o %in% names(df)) df[[o]] <- 0L
    df[[o]] <- ifelse(is.na(df[[o]]), 0L, as.integer(df[[o]]))
  }
  
  active_orgs <- organs_vec[vapply(organs_vec, function(o) {
    s <- df[[o]]; sum(s, na.rm = TRUE) > 0 && sum(s, na.rm = TRUE) < nrow(df)
  }, logical(1))]
  
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
  
  # Fill in the missing organs as NA and sort them by organs_vec
  res_full <- tibble(organ = organs_vec) %>%
    left_join(res, by = "organ") %>%
    mutate(symptom = symptom, .before = 1)
  
  res_full
}

#==== 2) Summarize the HR and p-value====
hr_long <- map_dfr(symptoms, ~get_hr_for_symptom(.x, tumor_all, organ_flags, organs_vec = organs))

# Add the full name of the symptom and the order of drawing (column order)
hr_long <- hr_long %>%
  mutate(symptom_full = symptom_names[symptom],
         symptom_full = factor(symptom_full, levels = symptom_names[symptoms]),
         organ = factor(organ, levels = organs))

#==== 3) heatmap ====
# Neutral is 0 (HR=1), <0 is blue, >0 is red
heat <- ggplot(hr_long, aes(x = symptom_full, y = organ, fill = HR)) +
  geom_tile(color = "white") +
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

#==== 4) save PNG / PDF ====
out_dir <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/Figure_summary/Organ/HR_heatmap"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "HR_heatmap_log2HR.png"),
       heat, width = 12, height = 5.5, dpi = 300)

ggsave(file.path(out_dir, "HR_heatmap_log2HR.pdf"),
       heat, width = 12, height = 5.5)

library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(ggplot2)

organs <- c("Liver","Lymph nodes","Lung","Bone","Breast","Pleura")
symptoms <- c("AP","CO","DI","DY","FA","FI","GHS","NV",
              "PA","SL","PF2","RF2","EF","CF","SF")
symptom_names <- c(
  FA="Fatigue", NV="Nausea/Vomit", PA="Pain", DY="Dyspnea", SL="Insomnia",
  AP="Appetite loss", CO="Constipation", DI="Diarrhea", FI="Financial",
  PF2="Physical Functioning", RF2="Role Functioning", EF="Emotional Functioning",
  CF="Cognitive Functioning", SF="Social Functioning", GHS="Global Health Status / QOL"
)

# ===== 1)ID- Organ Marker matrix =====
recode_organ <- function(x){
  case_when(
    x %in% c("Liver") ~ "Liver",
    x %in% c("Lymph node","Lymph nodes") ~ "Lymph nodes",
    x %in% c("Lung") ~ "Lung",
    x %in% c("Bone","Bone marrow","Skull") ~ "Bone",
    x %in% c("Breast") ~ "Breast",
    x %in% c("Pleura","Pleural effusion") ~ "Pleura",
    TRUE ~ "Other"
  )
}

organ_flags_raw <- tumor_all %>%
  filter(!is.na(ORGAN)) %>%
  mutate(Organ_clean = recode_organ(ORGAN)) %>%
  distinct(ID, Organ_clean) %>%
  mutate(flag = 1L) %>%
  pivot_wider(names_from = Organ_clean, values_from = flag, values_fill = 0L)

organ_flags <- organ_flags_raw
for (o in organs) if (!o %in% names(organ_flags)) organ_flags[[o]] <- 0L
organ_flags <- organ_flags[, c("ID", organs)]

# ===== 2) For individual symptoms, Cox regression was used to obtain HR =====
get_hr_for_symptom <- function(symptom, tumor_all, organ_flags, organs_vec = organs){
  ttfd_col     <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  df <- tumor_all %>%
    select(ID, all_of(ttfd_col), all_of(censored_col)) %>%
    left_join(organ_flags, by = "ID") %>%
    filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
           !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  for(o in organs_vec){
    if(!o %in% names(df)) df[[o]] <- 0L
    df[[o]] <- ifelse(is.na(df[[o]]), 0L, as.integer(df[[o]]))
  }
  
  active_orgs <- organs_vec[vapply(organs_vec, function(o){
    s <- df[[o]]; sum(s, na.rm=TRUE) > 0 && sum(s, na.rm=TRUE) < nrow(df)
  }, logical(1))]
  
  if(length(active_orgs) == 0){
    return(tibble(symptom = symptom, organ = organs_vec, HR = NA_real_, p = NA_real_))
  }
  
  surv_obj <- Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
  fit <- coxph(as.formula(paste("surv_obj ~", paste(active_orgs, collapse = " + "))), data = df)
  sm  <- summary(fit)
  
  res <- tibble(
    organ = rownames(sm$coefficients),
    HR    = sm$coefficients[,"exp(coef)"],
    p     = sm$coefficients[,"Pr(>|z|)"]
  )
  
  tibble(organ = organs_vec) %>%
    left_join(res, by = "organ") %>%
    mutate(symptom = symptom, .before = 1)
}

hr_long <- map_dfr(symptoms, ~get_hr_for_symptom(.x, tumor_all, organ_flags, organs_vec = organs)) %>%
  mutate(
    symptom_full = factor(symptom_names[symptom], levels = symptom_names[symptoms]),
    organ = factor(organ, levels = organs)
  )

ggplot(hr_long, aes(x = symptom_full, y = organ, fill = HR)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    name = "Hazard Ratio",
    colours = c("#f7fbff", "#deebf7", "white", "#d73027","#a50026"), 
    values = scales::rescale(c(0.5,0.8,1,1.3,2.5)), 
    limits = c(0.5, 2.5),
    oob = scales::squish
  ) +
  labs(
    title = "Hazard Ratio Heatmap by Symptom and Organ",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, 
                               size = 14, face = "bold"), 
    axis.text.y = element_text(size = 14, face = "bold"),  
    panel.grid = element_blank(),
    legend.position = "right"
  )

