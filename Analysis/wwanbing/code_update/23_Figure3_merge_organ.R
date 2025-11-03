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
Merged_TTFD <- read_csv("Merged_PRO_TTFD_01OCT2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output")
Merge_EQLQ_24APR2025 <- read_csv("Merge_EQLQ_24APR2025.csv")

Merged_TTFD <- Merged_TTFD %>%
  left_join(
    Merge_EQLQ_24APR2025 %>% distinct(UID, STDID),
    by = "UID"
  )

# Merge 2001 and 2004
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


# Merge TTFD
# Select all variables starting with TTFD10, TTFD10P, censored, censored10P
vars_to_merge <- colnames(Merged_TTFD)[grepl("^(TTFD10|TTFD10P|censored|censored10P)_", colnames(Merged_TTFD))]

vars_to_merge <- c("UID", vars_to_merge)

symptoms_TTFD <- Merged_TTFD[, vars_to_merge]

# Merge into tumor_ALL
tumor_2001 <- left_join(tumor_2001, symptoms_TTFD, by = "UID")
tumor_2004 <- left_join(tumor_2004, symptoms_TTFD, by = "UID")

tumor_all <- bind_rows(tumor_2001, tumor_2004)

# UIDs to remove, where all TTFD = NA
uids_remove <- c("006089-000-999-034",
                 "006089-000-999-039",
                 "006089-000-999-167",
                 "006089-000-999-177")

# Remove rows for the specified UIDs
tumor_all <- tumor_all %>%
  filter(!UID %in% uids_remove)

# Check removal result
length(unique(tumor_all$UID))

###Do target or non-target lesions more strongly correlate with symptom deterioration? 
tumor_all$LESTYPE <- ifelse(tumor_all$LESTYPE == 1, 1, 0)
table(tumor_all$LESTYPE)

# Step 1: Organ categorization (create a new variable Organ_clean)
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
# Bone Breast Liver Lung Lymph_nodes Other Pleura 
# 907   415    2214  1734 1564        2970  363 

# Step 2: Extract whether an organ was ever involved (ever) at the ID level
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

# Ever-by-ID organ flag matrix (seven columns)
organ_flags <- tumor_all %>%
  filter(!is.na(ORGAN)) %>%
  mutate(Organ_clean = recode_organ(ORGAN)) %>%
  distinct(ID, Organ_clean) %>%
  mutate(flag = 1L) %>%
  tidyr::pivot_wider(names_from = Organ_clean, values_from = flag, values_fill = 0L)

# Vector of organ columns needed
organs <- c("Liver","Lymph nodes","Lung","Bone","Breast","Pleura","Other")

# Step 3: Merge organ flags back to TTFD data
ttfd_cols <- grep("^TTFD10_", names(tumor_all), value = TRUE)
censor_cols <- grep("^censored_", names(tumor_all), value = TRUE)

### 15 symptoms

plot_forest_for_symptom <- function(symptom, tumor_all, organ_flags,
                                    organs = c("Liver","Lymph nodes","Lung","Bone","Breast","Pleura","Other")) {
  ttfd_col     <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  # Prepare data
  df <- tumor_all %>%
    dplyr::select(ID, UID, dplyr::all_of(ttfd_col), dplyr::all_of(censored_col)) %>%
    dplyr::left_join(organ_flags, by = "ID") %>%
    dplyr::filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
                  !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  # Safety: ensure all seven organ columns are present
  for (org in organs) {
    if (!org %in% names(df)) df[[org]] <- 0L
    df[[org]] <- ifelse(is.na(df[[org]]), 0L, as.integer(df[[org]]))
  }
  
  # Cox
  surv_obj <- survival::Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
  
  # Add backticks to column names to avoid formula parsing errors from spaces/special chars
  org_terms <- paste0("`", organs, "`")
  fml <- as.formula(paste("surv_obj ~", paste(org_terms, collapse = " + ")))
  fit <- survival::coxph(fml, data = df)
  sm  <- summary(fit)
  
  # Summarize
  forest_df <- data.frame(
    Term    = rownames(sm$coefficients),
    HR      = sm$coefficients[, "exp(coef)"],
    LCL     = sm$conf.int[,   "lower .95"],
    UCL     = sm$conf.int[,   "upper .95"],
    p_value = sm$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  ) %>%
    dplyr::mutate(
      Term  = gsub("`", "", Term),                 # Remove backticks
      Term  = gsub("^Other$", "Other organ", Term),
      Sig   = ifelse(p_value < 0.05, "Significant", "Not Significant"),
      OR_CI = sprintf("%.2f (%.2f, %.2f)", HR, LCL, UCL),
      Label = forcats::fct_reorder(Term, HR)
    )
  
  # Dynamic right-side label position & x-axis range
  # Dynamic right-side label position
  label_x_pos <- max(forest_df$UCL, na.rm = TRUE) * 1.05
  x_max <- max(label_x_pos * 1.15, max(forest_df$UCL, na.rm = TRUE) * 1.1)
  
  ggplot2::ggplot(forest_df, ggplot2::aes(x = HR, y = Label, color = Sig)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = LCL, xmax = UCL), height = 0.2) +
    ggplot2::geom_text(ggplot2::aes(x = label_x_pos, label = OR_CI),
                       hjust = 0, size = 4, color = "black") +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    ggplot2::scale_color_manual(values = c("Significant" = "#d73027", "Not Significant" = "gray40")) +
    ggplot2::scale_x_log10(limits = c(0.6, x_max)) +
    ggplot2::labs(
      title = paste0(symptom_names[symptom]),
      x = "Hazard Ratio", y = NULL, color = "p < 0.05"
    ) +
    ggplot2::theme_minimal(base_size = 15) +   # Increase base font size
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 16, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 16, face = "bold"),
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      plot.margin  = ggplot2::margin(t = 10, r = 140, b = 10, l = 10)
    ) +
    ggplot2::coord_cartesian(clip = "off")
}

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

out_dir <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/Figure_summary/Organ/Forest_plot_update"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (symptom in names(symptom_names)) {
  message("Plotting & saving: ", symptom)
  p <- plot_forest_for_symptom(symptom, tumor_all, organ_flags)
  
  # Build filename: abbreviation + full name
  file_name <- paste0(symptom, "_", symptom_names[[symptom]], ".png")
  
  ggsave(
    filename = file.path(out_dir, file_name),
    plot     = p,
    width    = 7, height = 5, dpi = 300
  )
}


############### Hazard Ratio heatmap

# Mapping for 7 organs and 15 symptoms with full names
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

#==== 1) Single symptom: extract HR/p for 7 organs (automatically ignore organs with no variation) ====
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
  
  # Fill back missing organs as NA and order by organs_vec
  res_full <- tibble(organ = organs_vec) %>%
    left_join(res, by = "organ") %>%
    mutate(symptom = symptom, .before = 1)
  
  res_full
}

####################
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(ggplot2)

# ===== Configuration (excluding "Other") =====
organs <- c("Liver","Lymph nodes","Lung","Bone","Breast","Pleura")
symptoms <- c("AP","CO","DI","DY","FA","FI","GHS","NV",
              "PA","SL","PF2","RF2","EF","CF","SF")
symptom_names <- c(
  FA="Fatigue", NV="Nausea/Vomit", PA="Pain", DY="Dyspnea", SL="Insomnia",
  AP="Appetite loss", CO="Constipation", DI="Diarrhea", FI="Financial",
  PF2="Physical Functioning", RF2="Role Functioning", EF="Emotional Functioning",
  CF="Cognitive Functioning", SF="Social Functioning", GHS="Global Health Status / QOL"
)

# ===== 1) ID–organ flag matrix =====
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

# ===== 2) Cox per single symptom to extract HR =====
get_hr_for_symptom <- function(symptom, tumor_all, organ_flags, organs_vec = organs){
  ttfd_col     <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  df <- tumor_all %>%
    select(ID, all_of(ttfd_col), all_of(censored_col)) %>%
    left_join(organ_flags, by = "ID") %>%
    filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
           !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  # Ensure organ columns exist
  for(o in organs_vec){
    if(!o %in% names(df)) df[[o]] <- 0L
    df[[o]] <- ifelse(is.na(df[[o]]), 0L, as.integer(df[[o]]))
  }
  
  # Keep only organs with variation
  active_orgs <- organs_vec[vapply(organs_vec, function(o){
    s <- df[[o]]; sum(s, na.rm=TRUE) > 0 && sum(s, na.rm=TRUE) < nrow(df)
  }, logical(1))]
  
  if(length(active_orgs) == 0){
    return(tibble(symptom = symptom, organ = organs_vec, HR = NA_real_, p = NA_real_))
  }
  
  surv_obj <- Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
  
  # Key: add backticks to active_orgs to avoid errors with "Lymph nodes"
  org_terms <- paste0("`", active_orgs, "`")
  fit <- coxph(as.formula(paste("surv_obj ~", paste(org_terms, collapse = " + "))), data = df)
  sm  <- summary(fit)
  
  res <- tibble(
    organ = gsub("`", "", rownames(sm$coefficients)),  # Remove backticks
    HR    = sm$coefficients[, "exp(coef)"],
    p     = sm$coefficients[, "Pr(>|z|)"]
  )
  
  tibble(organ = organs_vec) %>%
    left_join(res, by = "organ") %>%
    mutate(symptom = symptom, .before = 1)
}

# ===== 3) Aggregate 15×6 & draw heatmap =====
hr_long <- map_dfr(symptoms, ~get_hr_for_symptom(.x, tumor_all, organ_flags, organs_vec = organs)) %>%
  mutate(
    symptom_full = factor(symptom_names[symptom], levels = symptom_names[symptoms]),
    organ = factor(organ, levels = organs)
  )

ggplot(hr_long, aes(x = symptom_full, y = organ, fill = HR)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    name = "Hazard Ratio",
    colours = c("#f7fbff", "#deebf7", "white", "#d73027","#a50026"), # dark blue–light blue–white–light red–dark red
    values = scales::rescale(c(0.5,0.8,1,1.3,2.5)), # map by HR range
    limits = c(0.5, 2.5),
    oob = scales::squish
  ) +
  labs(title = "Hazard Ratio Heatmap by Symptom and Organ",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = 13, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"))

