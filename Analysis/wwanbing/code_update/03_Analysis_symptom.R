library(data.table)
library(dplyr)

eqlq <- fread("Merge_EQLQ_24APR2025.csv")

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

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Nausea_and_vomiting_", toupper(current_date), ".csv")

write.csv(NV4, new_filename, quote=FALSE,row.names=FALSE)


### Pain (PA 9 19)
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

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Pain_", toupper(current_date), ".csv")

write.csv(PA4, new_filename, quote=FALSE,row.names=FALSE)

### Dyspnoea (DY 8)
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

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Dyspnoea_", toupper(current_date), ".csv")

write.csv(DY4, new_filename, quote=FALSE,row.names=FALSE)

### Insomnia (SL 11)
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

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Insomnia_", toupper(current_date), ".csv")

write.csv(SL4, new_filename, quote=FALSE,row.names=FALSE)

### Appeite loss (AP 13)
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

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Apprtite_loss_", toupper(current_date), ".csv")

write.csv(AP4, new_filename, quote=FALSE,row.names=FALSE)

### Constipation (CO 16)
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

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Constipation_", toupper(current_date), ".csv")

write.csv(CO4, new_filename, quote=FALSE,row.names=FALSE)

### Diarrhoea (DI 17)
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

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Diarrhoea_", toupper(current_date), ".csv")

write.csv(DI4, new_filename, quote=FALSE,row.names=FALSE)

### Financial difficulties (FI 28)
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
  mutate(FLAG = "FA")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_PRO_Financial_difficulties_", toupper(current_date), ".csv")

write.csv(FI4, new_filename, quote=FALSE,row.names=FALSE)




