library(data.table)
library(dplyr)

eqlq <- fread("Merge_EQLQ_04JUN2025.csv") ###update PFS
bcr <- fread("Merged_BCR_21APR2025.csv")
tumor <- fread("Merged_tumor_21APR2025.csv")

###FA dataset(10 12 18)

#exclude TIME=-999 & DV =-999

# Step 1: Total number of rows before exclusion
total_rows <- nrow(eqlq)

# Step 2: Rows to be excluded (TIME = -999 or DV = -999)
excluded_rows <- eqlq %>%
  filter(TIME == -999 | DV == -999)

excluded_count <- nrow(excluded_rows)
excluded_proportion <- round(excluded_count / total_rows, 4)

# cat("Total excluded rows:", excluded_count, "\n")
# cat("Proportion of excluded rows:", excluded_proportion, "\n")

# Step 3: Excluded rows per study
excluded_by_study <- excluded_rows %>%
  group_by(STDID) %>%
  summarise(excluded_count = n(), .groups = "drop")

# Step 4: Total rows per study
total_by_study <- eqlq %>%
  group_by(STDID) %>%
  summarise(total_count = n(), .groups = "drop")

# Step 5: Merge and calculate proportion per study
excluded_summary <- left_join(total_by_study, excluded_by_study, by = "STDID") %>%
  mutate(
    excluded_count = ifelse(is.na(excluded_count), 0, excluded_count),
    excluded_proportion = round(excluded_count / total_count, 4)
  )

print(excluded_summary)

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999)

# Subset to only FLAG 10, 12, 18
fatigue_raw <- eqlq2 %>%
  filter(FLAG %in% c(10, 12, 18))

Fatigue <- fatigue_raw %>%
  filter(FLAG %in% c(10, 12, 18)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 3) %>%  # ensure all 3 flags present
  ungroup()

Fatigue2 <- Fatigue %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 3) * 100, 2)
  ) %>%
  ungroup()

Fatigue2 <- Fatigue2 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> FA
Fatigue2 <- Fatigue2 %>%
  mutate(FLAG = "FA")   

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Fatigue_", toupper(current_date), ".csv")

write.csv(Fatigue2, new_filename, quote=FALSE,row.names=FALSE)






