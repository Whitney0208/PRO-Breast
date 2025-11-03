###Do target or non-target lesions more strongly correlate with symptom deterioration? 
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
library(reshape2)
library(survival)
library(survminer)

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

# For each ID, take the earliest TIME record (baseline)
baseline_lesions <- tumor_all %>%
  filter(TIME != -999) %>%
  group_by(ID) %>%
  filter(TIME == min(TIME, na.rm = TRUE)) %>%
  ungroup()

baseline_lesions$LESTYPE<-as.integer(baseline_lesions$LESTYPE)

# baseline target tumor
lesion_summary <- baseline_lesions %>%
  distinct(UID, TIME, LESIONID, LESTYPE) %>%                 # Remove duplicates to avoid double counting
  summarise(
    target_n = n_distinct(LESIONID[LESTYPE == 1]),    # Count only target lesions
    .by = UID
  ) %>%
  mutate(target_n = replace_na(target_n, 0))          # Set to 0 if no target lesions


vars_to_merge <- colnames(Merged_TTFD)[grepl("^(TTFD10|censored)_", colnames(Merged_TTFD))]
vars_to_merge <- c("UID", vars_to_merge)
symptoms_TTFD <- Merged_TTFD[, vars_to_merge]

# Merge into BCR_1997 and BCR_2000
lesion_summary <- left_join(lesion_summary, symptoms_TTFD, by = "UID")

#Cox model + KM curve
ttfd_cols <- names(lesion_summary)[grepl("^TTFD", names(lesion_summary)) & !grepl("censored", names(lesion_summary))]
model_list <- list()

df<-lesion_summary

for (ttfd in ttfd_cols) {
  # suffix
  suffix <- sub("TTFD10_", "", ttfd)
  censored_col <- paste0("censored_", suffix)
  
  # clean data
  clean_df <- df[!is.na(df[[ttfd]]) & df[[ttfd]] != -999 &
                   !is.na(df[[censored_col]]) & df[[censored_col]] != -999, ]
  
  clean_df <- clean_df %>%
    mutate(target_cat = case_when(
      target_n == 0 ~ "0",
      target_n == 1 ~ "1",
      target_n == 2 ~ "2",
      target_n >= 3 ~ "3+"
    )) %>%
    mutate(target_cat = factor(target_cat, levels = c("0", "1", "2", "3+")))
  
  surv_obj <- Surv(time = clean_df[[ttfd]], event = 1 - clean_df[[censored_col]])
  model <- coxph(surv_obj ~ target_cat, data = clean_df)
  
  model_list[[paste0("model_", suffix)]] <- summary(model)
  
  cat("\n=== Cox Model for suffix:", suffix, "===\n")
  print(summary(model))
}


model_list <- list()

for (ttfd in ttfd_cols) {
  # suffix and censor columns
  suffix <- sub("^TTFD10_", "", ttfd)
  censored_col <- paste0("censored_", suffix)
  
  # clean data
  clean_df <- df[!is.na(df[[ttfd]]) & df[[ttfd]] != -999 &
                   !is.na(df[[censored_col]]) & df[[censored_col]] != -999, ]
  
  if (nrow(clean_df) < 10 || length(unique(clean_df[[ttfd]])) < 2) next
  
  # survival object
  surv_obj <- Surv(time = clean_df[[ttfd]], event = 1 - clean_df[[censored_col]])
  
  # -------- Grouping variable target_cat (0,1,2,3+) --------
  clean_df <- clean_df %>%
    mutate(target_cat = case_when(
      target_n == 0 ~ "0",
      target_n == 1 ~ "1",
      target_n == 2 ~ "2",
      target_n >= 3 ~ "3+"
    )) %>%
    mutate(target_cat = factor(target_cat, levels = c("0", "1", "2", "3+")))
  
  # Cox model (using grouped variable)
  model <- coxph(surv_obj ~ target_cat, data = clean_df)
  sm <- summary(model)
  model_list[[paste0("model_", suffix)]] <- sm
  
  cat("\n=== Cox Model for", suffix, "===\n")
  print(sm)
  
  # -------- KM survival curves (grouped by target_cat) --------
  fit_km <- survfit(surv_obj ~ target_cat, data = clean_df)
  
  g <- ggsurvplot(
    fit_km, data = clean_df,
    pval = TRUE,
    risk.table = TRUE,
    break.time.by = 50,  # <-- New addition: Display a time scale every 50 
    legend.title = "Target lesion count",
    legend.labs = levels(clean_df$target_cat),
    xlab = paste0("Time to symptom worsening (", ttfd, ")"),
    ylab = "Survival probability",
    title = paste0("KM: ", ttfd, " by target lesion count"),
    ggtheme = theme_minimal(),
    risk.table.height = 0.22
  )
  
  print(g)
}