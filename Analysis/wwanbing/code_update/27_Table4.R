rm(list=ls())
library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(haven)
library(tidyr)
library(purrr)
library(readr)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

# Merged_TTFD <- read_csv("Merged_PRO_TTFD_01OCT2025.csv")

Merge_EQLQ_24APR2025 <- read_csv("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/Merge_EQLQ_24APR2025.csv")

# Merged_TTFD <- Merged_TTFD %>%
#   left_join(
#     Merge_EQLQ_24APR2025 %>% distinct(UID, STDID),
#     by = "UID"
#   )
# 
# df0 <- Merged_TTFD
# 
# colnames(Merge_EQLQ_24APR2025)

# ---- 1 Data cleaning: ensure -999/-999.0 in DV are treated as missing ----
eqlq <- Merge_EQLQ_24APR2025 %>%
  mutate(DV = ifelse(DV %in% c(-999, -999.0), NA, DV)) %>%
  select(UID, TIME, FLAG, DV, STDID)

# ---- 2 Summarize total DV count, missing count, and missing percentage for each FLAG ----
dv_summary <- eqlq %>%
  group_by(FLAG) %>%
  summarise(
    Total_DV = n(),                            # total number of rows
    Missing_DV = sum(is.na(DV)),               # number of missing values
    Missing_Fraction = round(100 * Missing_DV / Total_DV, 2)  # missing percentage (%)
  ) %>%
  arrange(as.numeric(FLAG))                    # arrange by FLAG order

# ---- 3 View results ----
print(dv_summary, n = 30)





