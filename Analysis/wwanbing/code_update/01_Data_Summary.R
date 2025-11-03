library(data.table)
library(dplyr)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output")

eqlq <- fread("Merge_EQLQ_04JUN2025.csv")
bcr <- fread("Merged_BCR_21APR2025.csv")
tumor <- fread("Merged_tumor_21APR2025.csv")

# 1. how many patients in study, number and proportion

study_uid_summary <- eqlq %>%
  distinct(STDID, UID) %>%
  group_by(STDID) %>%
  summarise(n_patients = n(), .groups = "drop") %>%
  mutate(proportion = round(n_patients / sum(n_patients), 4))

print(study_uid_summary)

# 2. How many subjects do not have EROTC-QLQ C30 data by study (DV = -999 at all time points), 
#    number and proportion 

# Step 1: find the ID with all DV = -999
no_response_uid <- eqlq %>%
  group_by(STDID, UID) %>%
  summarise(all_missing = all(DV == -999), .groups = "drop") %>%
  filter(all_missing)

# Step 2: Total UID in each study
total_uid_by_study <- eqlq %>%
  distinct(STDID, UID) %>%
  group_by(STDID) %>%
  summarise(total_UID = n(), .groups = "drop")

# Step 3: the number of ID with all DV = -999 in each study
missing_uid_by_study <- no_response_uid %>%
  group_by(STDID) %>%
  summarise(missing_UID = n(), .groups = "drop")

study_missing_summary <- left_join(total_uid_by_study, missing_uid_by_study, by = "STDID") %>%
  mutate(missing_UID = ifelse(is.na(missing_UID), 0, missing_UID),
         proportion_missing = round(missing_UID / total_UID, 4))

print(study_missing_summary)

# 3. How many data points are missing in each study (DV = -999) 

# Step 1: Count missing DV points in each study
missing_data_by_study <- eqlq %>%
  filter(DV == -999) %>%
  group_by(STDID) %>%
  summarise(missing_count = n(), .groups = "drop")

# Step 2: Count total data points in each study
total_data_by_study <- eqlq %>%
  group_by(STDID) %>%
  summarise(total_count = n(), .groups = "drop")

# Step 3: Join and calculate proportion of missing
missing_summary <- left_join(total_data_by_study, missing_data_by_study, by = "STDID") %>%
  mutate(
    missing_count = ifelse(is.na(missing_count), 0, missing_count),
    missing_proportion = round(missing_count / total_count, 4)
  )

print(missing_summary)

# 4. How many data points are missing by question (30 questions). 

# Step 1: Count missing responses per question (DV == -999)
missing_by_question <- eqlq %>%
  filter(DV == -999) %>%
  group_by(FLAG) %>%
  summarise(missing_count = n(), .groups = "drop")

# Step 2: Count total responses per question
total_by_question <- eqlq %>%
  group_by(FLAG) %>%
  summarise(total_count = n(), .groups = "drop")

# Step 3: Join and calculate proportion
question_missing_summary <- left_join(total_by_question, missing_by_question, by = "FLAG") %>%
  mutate(
    missing_count = ifelse(is.na(missing_count), 0, missing_count),
    missing_proportion = round(missing_count / total_count, 4)
  ) %>%
  arrange(FLAG)

print(question_missing_summary)




