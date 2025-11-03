library(data.table)
library(dplyr)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_24APR2025.csv")
# bcr <- fread("Merged_BCR_21APR2025.csv")
# tumor <- fread("Merged_tumor_21APR2025.csv")

###GHS dataset(QL2 29 30)

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

# Subset to only FLAG 10, 12, 18
QL2_1 <- eqlq2 %>%
  filter(FLAG %in% c(29,30))

QL2_2 <- QL2_1 %>%
  filter(FLAG %in% c(29,30)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 2) %>%  # ensure all 3 flags present
  ungroup()

QL2_3 <- QL2_2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 6) * 100, 2)
  ) %>%
  ungroup()

QL2_3 <- QL2_3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> QL2
QL2_4 <- QL2_3 %>%
  mutate(FLAG = "FA")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_GHS_", toupper(current_date), ".csv")

write.csv(QL2_4, new_filename, quote=FALSE,row.names=FALSE)



