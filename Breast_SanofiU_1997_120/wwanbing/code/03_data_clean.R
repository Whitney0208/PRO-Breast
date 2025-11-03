rm(list=ls())

library(plyr)
library(dplyr)
library(tidyr)
library(haven)
library(stringr)

adae <- read_sas(data_file='adae.sas7bdat') # AE
adcm <- read_sas(data_file='adcm.sas7bdat') # hormone therapy
addm <- read_sas(data_file='addm.sas7bdat')
adds <- read_sas(data_file='adds.sas7bdat') #disposition
addv <- read_sas(data_file='addv.sas7bdat') #not useful
adef <- read_sas(data_file='adef.sas7bdat') #week death
adeg <- read_sas(data_file='adeg.sas7bdat') #ECG
ades <- read_sas(data_file='ades.sas7bdat') # dose intensity
adex <- read_sas(data_file='adex.sas7bdat') # dose
adfa <- read_sas(data_file='adfa.sas7bdat') # tumor size, receptor, metastasis
adlb <- read_sas(data_file='adlb.sas7bdat') # lab test, not included in this dataset
admh <- read_sas(data_file='admh.sas7bdat') #anxiety or depression
adoc <- read_sas(data_file='adoc.sas7bdat') #prior treatment
adpe <- read_sas(data_file='adpe.sas7bdat') #overall assessment
adpea <- read_sas(data_file='adpea.sas7bdat') #disease free survival
adpr <- read_sas(data_file='adpr.sas7bdat') #breast cancer surgery
adqs <- read_sas(data_file='adqs.sas7bdat') 
adsc <- read_sas(data_file='adsc.sas7bdat') #clinical characteristics
adsl <- read_sas(data_file='adsl.sas7bdat') #not useful
adsv <- read_sas(data_file='adsv.sas7bdat') #plan day
advs <- read_sas(data_file='advs.sas7bdat')

addm2 <- addm %>%
  transmute(
    UID     = RUSUBJID,
    ID      = RSUBJID,
    AGE     = AGE,
    SEX     = SEX,
    RACE    = RACE,
    ETHNIC  = -999,               
    REGION  = REGION,
    ARM     = 0,                  
    EMPLOY  = -999,            
    MENOS   = PRE_POST,
    LESION1 = NBLYNOI
  ) %>%
  mutate(
    SEX = case_when(
      SEX == "M"   ~ 0,
      SEX == "F" ~ 1,
      TRUE            ~ -999
    ),
    RACE = case_when(
      RACE == "CAUCASIAN" ~ 1,
      RACE == "BLACK"     ~ 2,
      RACE == "" | is.na(RACE) ~ -999,
      TRUE ~ 3
    )
  )


adsc_wide <- adsc %>%
  filter(SCTESTCD %in% c("DIAGC", "PSWHO0", "CHEMONY", "SURGNY", "RADIONY", "HORMONY", "HT", "WT0")) %>%
  select(RUSUBJID, SCTESTCD, SCORRES) %>%
  pivot_wider(names_from = SCTESTCD, values_from = SCORRES)

adsc2 <- adsc_wide %>%
  transmute(
    UID    = RUSUBJID,
    
    STAGE = case_when(
      DIAGC %in% c("pT1,pN1,M1", "pT2,pN2,M1", "pT2,pN2,M0", "pT3,pN1,M0") ~ "III",
      DIAGC %in% c("pT1,pN2,M0", "pT2,pN1,M0", "pT2,pN1,M1") ~ "II",
      DIAGC %in% c("pT1,pN1,M0") ~ "IV",
      TRUE ~ "-999"
    ),
    
    ECOG = case_when(
      PSWHO0 == "0" ~ 0,
      PSWHO0 == "1" ~ 1,
      PSWHO0 %in% c("2", "3") ~ 2,
      TRUE ~ -999
    ),
    
    CHEMO   = case_when(CHEMONY == "YES" ~ 1, CHEMONY == "NO" ~ 2, TRUE ~ -999),
    SURGERY = case_when(SURGNY == "YES" ~ 1, SURGNY == "NO" ~ 2, TRUE ~ -999),
    RADIO   = case_when(RADIONY == "YES" ~ 1, RADIONY == "NO" ~ 2, TRUE ~ -999),
    HORMON  = case_when(HORMONY == "YES" ~ 1, HORMONY == "NO" ~ 2, TRUE ~ -999),
    
    HEIGHT = ifelse(is.na(as.numeric(HT)), -999, as.numeric(HT)),
    WEIGHT = ifelse(is.na(as.numeric(WT0)), -999, as.numeric(WT0))
  )

adfa_wide <- adfa %>%
  filter(FATEST %in% c("Estrogen Receptors Immunohist", 
                       "Progesterone Receptors Immunohist")) %>%
  select(RUSUBJID, FATEST, FAORRES) %>%
  pivot_wider(names_from = FATEST, values_from = FAORRES)

adfa2 <- adfa_wide %>%
  transmute(
    UID = RUSUBJID,
    ERS = case_when(
      `Estrogen Receptors Immunohist` == "POSITIVE" ~ 1,
      `Estrogen Receptors Immunohist` == "NEGATIVE" ~ 2,
      `Estrogen Receptors Immunohist` == "NOT ASSESSABLE" ~ 3,
      TRUE ~ -999
    ),
    PGRS = case_when(
      `Progesterone Receptors Immunohist` == "POSITIVE" ~ 1,
      `Progesterone Receptors Immunohist` == "NEGATIVE" ~ 2,
      `Progesterone Receptors Immunohist` == "NOT ASSESSABLE" ~ 3,
      TRUE ~ -999
    )
  )

# HER2 adpea, 0 = 'Negative', 1 = 'Positive', -1 = 'Unknown'; Initial her2 status, 1= Positive, 2= Negative, 3= Unknown
adpea3 <- adpea %>%
  transmute(
    UID = RUSUBJID,
    HER2 = case_when(
      `HER2NEU` == "1" ~ 1,
      `HER2NEU` == "0" ~ 2,
      `HER2NEU` == "-1" ~ 3,
      TRUE ~ -999
    )
  )

admh2 <- admh %>%
  select(RUSUBJID, MHICD9CD) %>%
  group_by(RUSUBJID) %>%
  summarise(
    MHDEPRESSION = if_else(any(MHICD9CD == "311"), 1, 0),
    MDANXIETY = if_else(any(MHICD9CD %in% c("300.00", "300.01")), 1, 0)
  ) %>%
  ungroup() %>%
  rename(UID = RUSUBJID)

adpea2 <- adpea %>%
  transmute(
    UID = RUSUBJID,
    
    DTHDY = ifelse(is.na(DROVER), -999, DROVER * 30),
    
    DTH = case_when(
      is.na(DTHDY) ~ -999,
      !is.na(DTHDY) ~ as.numeric(EV_OVER)  # make sure EV_OVER== 0 or 1
    ),
    
    PFSDY = ifelse(is.na(drsurfd), -999, drsurfd * 30),
    
    PFS = case_when(
      is.na(BCRYN) ~ -999,
      TRUE ~ as.numeric(BCRYN)
    )
  )

# EQLQ C30 Score
adqs2 <- adqs %>%
  filter(QSSCAT == "C30 Score")
         
eqlq_c30 <- adqs2 %>%
  select(RUSUBJID, QSTESTCD, QSSTRESC, QSDY) %>%
  mutate(
    QSFLAG = QSTESTCD,
    FLAG = as.numeric(gsub("[^0-9]", "", QSTESTCD)),
    DV = case_when(
      QSSTRESC %in% c("No", "Not at All")      ~ 1,
      QSSTRESC == "A Little"                   ~ 2,
      QSSTRESC == "Quite a Bit"               ~ 3,
      QSSTRESC %in% c("Yes", "Very Much")      ~ 4,
      str_detect(QSSTRESC, "^[1-7]$")           ~ as.numeric(QSSTRESC),
      QSSTRESC == "" | is.na(QSSTRESC)          ~ -999,
      TRUE                                     ~ -999
    ),
    TIME = QSDY
  ) %>%
  select(RUSUBJID, QSFLAG, FLAG, DV, TIME) %>%
  distinct()

# Long format
eqlq_long <- eqlq_c30 %>%
  mutate(
    DV = ifelse(is.na(DV), -999, DV),
    TIME
  ) %>%
  select(RUSUBJID, QSFLAG, FLAG, TIME, DV) %>%
  distinct()

# 120 duplicate ID, use slice(1)
eqlq_long2 <- eqlq_long %>%
  group_by(RUSUBJID, FLAG, TIME) %>%
  arrange(
    (DV == -999 & any(DV != -999)),  # Put "-999 of groups mixed with -999" at the back
    row_number()                     # The rest of the order does not move
  ) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(RUSUBJID, TIME, FLAG)

#wide format
eqlq_wide <- eqlq_long2 %>%
  arrange(RUSUBJID, TIME, FLAG) %>%  
  group_by(RUSUBJID) %>%
  mutate(
    all_na = all(is.na(TIME)),
    TIME_valid = ifelse(is.na(TIME), Inf, TIME),
    baseline_time = ifelse(all_na, NA, min(TIME_valid))
  ) %>%
  # Retain the entire group of 30 items corresponding to the baseline time or retain the first 30 items when all nas are present
  mutate(row_idx = row_number()) %>%
  filter(
    (!all_na & TIME == baseline_time) |
      (all_na & row_idx <= 30)
  ) %>%
  ungroup() %>%
  select(RUSUBJID, FLAG, DV) %>%
  mutate(FLAG = paste0("B", sprintf("%02d", FLAG))) %>%
  tidyr::pivot_wider(names_from = FLAG, values_from = DV) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), -999, .x)))

# Merge long-format and wide-format datasets
eqlq_final <- eqlq_long %>%
  left_join(eqlq_wide, by = "RUSUBJID") %>%
  rename(UID = RUSUBJID)  # Rename ID variable

# Ensure consistent missing value handling
eqlq_final <- eqlq_final %>%
  mutate(across(everything(), ~ ifelse(.x == "", NA, .x))) %>%  # Convert empty strings to NA
  mutate(
    across(where(is.numeric), ~ replace_na(.x, -999)),  
    across(where(is.character), ~ replace_na(.x, "-999"))  
  )

final_data <- addm2 %>%
  left_join(adsc2,   by = "UID") %>%
  left_join(adfa2, by = "UID") %>%
  left_join(admh2,   by = "UID") %>%
  left_join(adpea2,    by = "UID") %>% 
  left_join(adpea3,    by = "UID") %>%
  distinct()

final_data2 <- final_data %>%
  full_join(eqlq_final, by = "UID")   # Perform a full join to keep all IDs

dataset2 <- final_data2 %>%
  mutate(across(everything(), ~ ifelse(.x == "", NA, .x))) %>%  # Convert empty strings to NA
  mutate(
    across(where(is.numeric), ~ replace_na(.x, -999)),  
    across(where(is.character), ~ replace_na(.x, "-999"))  
  )%>%
  arrange(UID, TIME, FLAG)

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_SanofiU_1997_120_", toupper(current_date), ".csv")

write.csv(dataset2, new_filename, quote=FALSE,row.names=FALSE)










