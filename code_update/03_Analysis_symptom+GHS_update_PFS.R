library(data.table)
library(dplyr)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

###Nausea and vomiting (NV 14 15)

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999)

# Subset to only FLAG 14, 15
NV <- eqlq2 %>%
  filter(FLAG %in% c(14, 15))

NV2 <- NV %>%
  filter(FLAG %in% c(14,15)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 2) %>%  # ensure all 2 flags present
  ungroup()

NV3 <- NV2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 3) * 100, 2)
  ) %>%
  ungroup()

NV3 <- NV3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> NV
NV4 <- NV3 %>%
  mutate(FLAG = "NV")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Nausea_and_vomiting_", toupper(current_date), ".csv")

write.csv(NV4, new_filename, quote=FALSE,row.names=FALSE)


### Pain (PA 9 19)
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

PA <- eqlq2 %>%
  filter(FLAG %in% c(9, 19))

PA2 <- PA %>%
  filter(FLAG %in% c(9,19)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 2) %>%  # ensure all 2 flags present
  ungroup()

PA3 <- PA2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 3) * 100, 2)
  ) %>%
  ungroup()

PA3 <- PA3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> PA
PA4 <- PA3 %>%
  mutate(FLAG = "PA") 

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Pain_", toupper(current_date), ".csv")

write.csv(PA4, new_filename, quote=FALSE,row.names=FALSE)

### Dyspnoea (DY 8)
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

DY <- eqlq2 %>%
  filter(FLAG %in% c(8))

DY2 <- DY %>%
  filter(FLAG %in% c(8)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 1) %>%  # ensure all 2 flags present
  ungroup()

DY3 <- DY2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 3) * 100, 2)
  ) %>%
  ungroup()

DY3 <- DY3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> DY
DY4 <- DY3 %>%
  mutate(FLAG = "DY")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Dyspnoea_", toupper(current_date), ".csv")

write.csv(DY4, new_filename, quote=FALSE,row.names=FALSE)

### Insomnia (SL 11)
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

SL <- eqlq2 %>%
  filter(FLAG %in% c(11))

SL2 <- SL %>%
  filter(FLAG %in% c(11)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 1) %>%  # ensure all 2 flags present
  ungroup()

SL3 <- SL2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 3) * 100, 2)
  ) %>%
  ungroup()

SL3 <- SL3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> SL
SL4 <- SL3 %>%
  mutate(FLAG = "SL")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Insomnia_", toupper(current_date), ".csv")

write.csv(SL4, new_filename, quote=FALSE,row.names=FALSE)

### Appetite loss (AP 13)
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

AP <- eqlq2 %>%
  filter(FLAG %in% c(13))

AP2 <- AP %>%
  filter(FLAG %in% c(13)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 1) %>%  # ensure all 2 flags present
  ungroup()

AP3 <- AP2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 3) * 100, 2)
  ) %>%
  ungroup()

AP3 <- AP3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> AP
AP4 <- AP3 %>%
  mutate(FLAG = "AP") 

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Appetite_loss_", toupper(current_date), ".csv")

write.csv(AP4, new_filename, quote=FALSE,row.names=FALSE)

### Constipation (CO 16)
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

CO <- eqlq2 %>%
  filter(FLAG %in% c(16))

CO2 <- CO %>%
  filter(FLAG %in% c(16)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 1) %>%  # ensure all 2 flags present
  ungroup()

CO3 <- CO2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 3) * 100, 2)
  ) %>%
  ungroup()

CO3 <- CO3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> CO
CO4 <- CO3 %>%
  mutate(FLAG = "CO") 

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Constipation_", toupper(current_date), ".csv")

write.csv(CO4, new_filename, quote=FALSE,row.names=FALSE)

### Diarrhoea (DI 17)
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

DI <- eqlq2 %>%
  filter(FLAG %in% c(17))

DI2 <- DI %>%
  filter(FLAG %in% c(17)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 1) %>%  # ensure all 2 flags present
  ungroup()

DI3 <- DI2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 3) * 100, 2)
  ) %>%
  ungroup()

DI3 <- DI3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> DI
DI4 <- DI3 %>%
  mutate(FLAG = "DI") 

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Diarrhoea_", toupper(current_date), ".csv")

write.csv(DI4, new_filename, quote=FALSE,row.names=FALSE)

### Financial difficulties (FI 28)
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

FI <- eqlq2 %>%
  filter(FLAG %in% c(28))

FI2 <- FI %>%
  filter(FLAG %in% c(28)) %>%
  group_by(UID, TIME, FLAG) %>%
  slice(1) %>%   # keep only the first entry per UID-TIME-FLAG
  ungroup() %>%
  group_by(UID, TIME) %>%
  filter(n_distinct(FLAG) == 1) %>%  # ensure all 2 flags present
  ungroup()

FI3 <- FI2 %>%
  group_by(UID, TIME) %>%
  mutate(
    DV2 = round(((mean(as.numeric(DV)) - 1) / 3) * 100, 2)
  ) %>%
  ungroup()

FI3 <- FI3 %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "UNK", .)))

###delete DV, rename FLAG -> FI
FI4 <- FI3 %>%
  mutate(FLAG = "FI")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Financial_difficulties_", toupper(current_date), ".csv")

write.csv(FI4, new_filename, quote=FALSE,row.names=FALSE)


###GHS dataset(QL2 29 30)
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/")
eqlq <- fread("Merge_EQLQ_04JUN2025.csv")

# Exclude rows where TIME == -999 or DV == -999
eqlq2 <- eqlq %>%
  filter(TIME != -999 & DV != -999 & DV != "U")

# Subset to only FLAG 29 30
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
  mutate(FLAG = "GHS")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_GHS_", toupper(current_date), ".csv")

write.csv(QL2_4, new_filename, quote=FALSE,row.names=FALSE)



###update PFS dataset(04_Analysis_symptom.R)
#delete DV, distinct dataset
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS_old")

NV <- fread("Breast_PRO_Nausea_and_vomiting_18JUN2025.csv")
PA <- fread("Breast_PRO_Pain_18JUN2025.csv")
DY <- fread("Breast_PRO_Dyspnoea_18JUN2025.csv")
SL <- fread("Breast_PRO_Insomnia_18JUN2025.csv")
AP <- fread("Breast_PRO_Appetite_loss_18JUN2025.csv")  
CO <- fread("Breast_PRO_Constipation_18JUN2025.csv")
DI <- fread("Breast_PRO_Diarrhoea_18JUN2025.csv")
FI <- fread("Breast_PRO_Financial_difficulties_18JUN2025.csv")
GHS <- fread("Breast_PRO_GHS_18JUN2025.csv")

NV2 <- NV[, !c("QSFLAG", "DV")] %>% distinct()
PA2 <- PA[, !c("QSFLAG", "DV")] %>% distinct()
DY2 <- DY[, !c("QSFLAG", "DV")] %>% distinct()
SL2 <- SL[, !c("QSFLAG", "DV")] %>% distinct()
AP2 <- AP[, !c("QSFLAG", "DV")] %>% distinct()
CO2 <- CO[, !c("QSFLAG", "DV")] %>% distinct()
DI2 <- DI[, !c("QSFLAG", "DV")] %>% distinct()
FI2 <- FI[, !c("QSFLAG", "DV")] %>% distinct()
GHS2 <- GHS[, !c("QSFLAG", "DV")] %>% distinct()

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")

current_date <- format(Sys.Date(), "%d%b%Y")

write.csv(NV2, paste0("Breast_PRO_Nausea_and_vomiting_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(PA2, paste0("Breast_PRO_Pain_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(DY2, paste0("Breast_PRO_Dyspnoea_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(SL2, paste0("Breast_PRO_Insomnia_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(AP2, paste0("Breast_PRO_Appetite_loss_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(CO2, paste0("Breast_PRO_Constipation_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(DI2, paste0("Breast_PRO_Diarrhoea_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(FI2, paste0("Breast_PRO_Financial_difficulties_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(GHS2, paste0("Breast_PRO_GHS_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)


