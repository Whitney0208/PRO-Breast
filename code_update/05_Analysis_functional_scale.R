library(data.table)
library(dplyr)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_24APR2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

###PF2

# Subset to only FLAG 1, 2, 3, 4 ,5
PF2_1 <- eqlq2 %>%
  filter(FLAG %in% c(1, 2, 3, 4, 5))

PF2_2 <- PF2_1 %>%
  filter(FLAG %in% c(1, 2, 3, 4, 5)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 5) %>%  # ensure all 5 flags present
  ungroup()

PF2_3 <- PF2_2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round((1 - ((mean(as.numeric(DV)) - 1) / 3)) * 100, 2)
  ) %>%
  ungroup()

PF2_3 <- PF2_3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> NV
PF2_4 <- PF2_3 %>%
  mutate(FLAG = "PF2")

PF2_5 <- PF2_4 %>% dplyr::select(-QSFLAG, -DV) %>% distinct()

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Physical_functioning(revised)_", toupper(current_date), ".csv")

write.csv(PF2_5, new_filename, quote=FALSE,row.names=FALSE)


###RF2 6 7
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

# Subset to only FLAG 6,7
RF2_1 <- eqlq2 %>%
  filter(FLAG %in% c(6,7))

RF2_2 <- RF2_1 %>%
  filter(FLAG %in% c(6,7)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 2) %>%  # ensure all 2 flags present
  ungroup()

RF2_3 <- RF2_2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round((1 - ((mean(as.numeric(DV)) - 1) / 3)) * 100, 2)
  ) %>%
  ungroup()

RF2_3 <- RF2_3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> NV
RF2_4 <- RF2_3 %>%
  mutate(FLAG = "RF2")

RF2_5 <- RF2_4 %>% dplyr::select(-QSFLAG, -DV) %>% distinct()

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Role_functioning(revised)_", toupper(current_date), ".csv")

write.csv(RF2_5, new_filename, quote=FALSE,row.names=FALSE)

###EF
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

# Subset to only FLAG 21 to 24
EF1 <- eqlq2 %>%
  filter(FLAG %in% c(21,22,23,24))

EF2 <- EF1 %>%
  filter(FLAG %in% c(21,22,23,24)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 4) %>%  # ensure all 2 flags present
  ungroup()

EF3 <- EF2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round((1 - ((mean(as.numeric(DV)) - 1) / 3)) * 100, 2)
  ) %>%
  ungroup()

EF3 <- EF3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> EF
EF4 <- EF3 %>%
  mutate(FLAG = "EF")

EF5 <- EF4 %>% dplyr::select(-QSFLAG, -DV) %>% distinct()

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Emotional_functioning_", toupper(current_date), ".csv")

write.csv(EF5, new_filename, quote=FALSE,row.names=FALSE)


###CF
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

# Subset to only FLAG 20, 25
CF1 <- eqlq2 %>%
  filter(FLAG %in% c(20,25))

CF2 <- CF1 %>%
  filter(FLAG %in% c(20,25)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 2) %>%  # ensure all 2 flags present
  ungroup()

CF3 <- CF2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round((1 - ((mean(as.numeric(DV)) - 1) / 3)) * 100, 2)
  ) %>%
  ungroup()

CF3 <- CF3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> NV
CF4 <- CF3 %>%
  mutate(FLAG = "CF")

CF5 <- CF4 %>% dplyr::select(-QSFLAG, -DV) %>% distinct()

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Cognitive_functioning_", toupper(current_date), ".csv")

write.csv(CF5, new_filename, quote=FALSE,row.names=FALSE)

###SF
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

# Subset to only FLAG 26,27
SF1 <- eqlq2 %>%
  filter(FLAG %in% c(26,27))

SF2 <- SF1 %>%
  filter(FLAG %in% c(26,27)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 2) %>%  # ensure all 2 flags present
  ungroup()

SF3 <- SF2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round((1 - ((mean(as.numeric(DV)) - 1) / 3)) * 100, 2)
  ) %>%
  ungroup()

SF3 <- SF3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> NV
SF4 <- SF3 %>%
  mutate(FLAG = "SF")

SF5 <- SF4 %>% dplyr::select(-QSFLAG, -DV) %>% distinct()

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Social_functioning_", toupper(current_date), ".csv")

write.csv(SF5, new_filename, quote=FALSE,row.names=FALSE)
