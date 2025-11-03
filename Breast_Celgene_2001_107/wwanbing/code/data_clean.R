library(plyr)
library(dplyr)
library(haven)
library(tidyr)
library(purrr)

# check dataset
adex <- read_sas(data_file='adex.sas7bdat') #adverse effect
base <- read_sas(data_file='base.sas7bdat') #baseline target lesion number and tumor size
brst <- read_sas(data_file='brst.sas7bdat') #baseline individual tumor location and size
cen_lab2 <- read_sas(data_file='cen_lab2.sas7bdat') # lab tests 2
cen_lab3 <- read_sas(data_file='cen_lab3.sas7bdat') # lab tests 3
cen_labs <- read_sas(data_file='cen_labs.sas7bdat') #lab test
cmed <- read_sas(data_file='cmed.sas7bdat') # comed
cpro <- read_sas(data_file='cpro.sas7bdat') # procedure
demo <- read_sas(data_file='demo.sas7bdat') # demographics
dose <- read_sas(data_file='dose.sas7bdat') #dose
echo <- read_sas(data_file='echo.sas7bdat') #echo
ekgr <- read_sas(data_file='ekgr.sas7bdat') #EKG
elig <- read_sas(data_file='elig.sas7bdat') #eligibility
eosr <- read_sas(data_file='eosr.sas7bdat') #treatment cycle
eosr <- read_sas(data_file='eosr.sas7bdat') #EOSR signature
eqlq <- read_sas(data_file='eqlq.sas7bdat') #EORTC questionnaire
fact <- read_sas(data_file='fact.sas7bdat') #fact questionnaire, use it
lesn <- read_sas(data_file='lesn.sas7bdat') #lesion size and location
loc_labs <- read_sas(data_file='loc_labs.sas7bdat') #lab test
mehx <- read_sas(data_file='mehx.sas7bdat') #medical history
muga <- read_sas(data_file='muga.sas7bdat') #LVEF
ntle <- read_sas(data_file='ntle.sas7bdat') #non-target lesion
phex <- read_sas(data_file='phex.sas7bdat')
phon <- read_sas(data_file='phon.sas7bdat') #phone contact
prth <- read_sas(data_file='prth.sas7bdat') #drug and response
prtx <- read_sas(data_file='prtx.sas7bdat')
pspn <- read_sas(data_file='pspn.sas7bdat') #physician assessment
ptss <- read_sas(data_file='ptss.sas7bdat') #other disease
resp <- read_sas(data_file='resp.sas7bdat') #lesion response
rlps <- read_sas(data_file='rlps.sas7bdat') #relpase and ER/HER2 status
scan <- read_sas(data_file='scan.sas7bdat') # organ lesion
surv <- read_sas(data_file='surv.sas7bdat') #survival
targ <- read_sas(data_file='targ.sas7bdat') #target lesion size
toxy <- read_sas(data_file='toxy.sas7bdat') #Toxicity
vitl <- read_sas(data_file='vitl.sas7bdat') #vitals
wclesion <- read_sas(data_file='wclesion.sas7bdat') #CT measured lesion size
wcresp <- read_sas(data_file='wcresp.sas7bdat') #CT measured total lesion size

# all dataset in list
# ds_list <- list(
#   adex=adex, base=base, brst=brst, cen_lab2=cen_lab2, cen_lab3=cen_lab3,
#   cen_labs=cen_labs, cmed=cmed, cpro=cpro, demo=demo, dose=dose,
#   echo=echo, ekgr=ekgr, elig=elig, eosr=eosr, eqlq=eqlq, fact=fact,
#   lesn=lesn, loc_labs=loc_labs, mehx=mehx, muga=muga, ntle=ntle,
#   phex=phex, phon=phon, prth=prth, pspn=pspn, ptss=ptss, resp=resp,
#   rlps=rlps, scan=scan, surv=surv, targ=targ, toxy=toxy, vitl=vitl,
#   wclesion=wclesion, wcresp=wcresp
# )

# dim() for each dataset
# dataset_dims <- sapply(ds_list, dim)
# print(dataset_dims)
# 
# # str() for each dataset
# for (name in names(ds_list)) {
#   cat("\n====================\n")
#   cat("Structure of dataset:", name, "\n")
#   str(ds_list[[name]])
# }
# 
# # Match variable names between data_spec.xlsx and all datasets
# 
# ###import Data_spec.xlsx
# library(readxl)
# data_spec <- read_excel("Data_spec.xlsx", sheet = 1)
# 
# # Creates/initializes the Source column, all NA
# data_spec$Source <- NA_character_
# 
# library(purrr)
# # Extract the column names label of all data sets and summarize them into ds_meta
# # map2_dfr() All the data sets are grouped into one large table (data frame).
# # For columns without labels, set to the empty string "".
# 
# ds_meta <- map2_dfr(ds_list, names(ds_list), function(df, dsname) {
#   tibble(
#     dataset = dsname,
#     colname = names(df),
#     label = map_chr(df, ~ attr(.x, "label") %||% "")
#   )
# })


# demo2: Extract UID, ID, AGE, RACE, HEIGHT, MENOS, ARM, and ETHNIC, retaining RUSUBJ_ID for merging.
demo2 <- demo %>%
  transmute(
    UID       = RSUBJ_ID,
    ID        = RSUBJ_ID,
    AGE,
    SEX       = 1,         # FEMALE (all subjects are female)
    RACE      = RACE_GEN,
    HEIGHT    = DEMO_004,
    MENOS     = DEMO_005,
    ARM       = ELIG_007,
    ETHNIC    = -999,       # Assign -999 for ETHNIC uniformly
    WBC   = -999,
    NEUT  = -999,
    CREAT = -999,
    ALT   = -999,
    BILI  = -999,
    ALP   = -999,
    AST   = -999,
    ALB   = -999
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, -999)))

# base2: Extract LESION1 and retain RUSUBJ_ID.
base2 <- base %>%
  transmute(
    ID = RSUBJ_ID,
    LESION1   = BASE_001
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, -999)))

# brst2: Extract CHEMO, SURGERY, RADIO, HORMON, ERS, and PGRS while retaining RUSUBJ_ID.
brst2 <- brst %>%
  transmute(
    ID = RSUBJ_ID,
    CHEMO     = BRST_009,  # Baseline Prior chemotherapy (1=Yes, 2=No)
    SURGERY   = BRST_006,  # Baseline Surgery history (1=Yes, 2=No)
    RADIO     = BRST_007,  # Baseline radiotherapy (1=Yes, 2=No)
    HORMON    = BRST_008,  # Initial hormonal therapy (1=Yes, 2=No)
    ERS       = BRST_010,  # Initial ER status (1=Positive, 2=Negative, 3=Unknown)
    PGRS      = BRST_011   # Initial PgR status (1=Positive, 2=Negative, 3=Unknown)
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, -999)))

#prtx2: ANTHRAC
prtx2 <- prtx %>%
  transmute(
    ID = RSUBJ_ID,
    ANTHRAC = ifelse(PRTX_020 == "U", -999, PRTX_020), 
  )%>%
  mutate(ANTHRAC = ifelse(ANTHRAC %in% c(1, -999), ANTHRAC, 0))

# ptss2: From ptss, derive MHDEPRESSION and MHANXIETY based on PTSS_001, and retain RUSUBJ_ID.
ptss2 <- ptss %>%
  mutate(
    MHDEPRESSION = if_else(PTSS_001 == "DEPRESSION", 1, 0),
    MHANXIETY    = if_else(PTSS_001 == "ANXIETY",    1, 0)
  ) %>%
  transmute(
    ID = RSUBJ_ID,
    MHDEPRESSION,
    MHANXIETY
  )

ptss3 <- ptss2 %>%
  group_by(ID) %>%
  summarise(
    MHDEPRESSION = max(MHDEPRESSION),  # keep 1
    MHANXIETY    = max(MHANXIETY)      # keep 1
  ) %>%
  ungroup()

# Find IDs that are in prtx2 but not in ptss3
missing_ids <- anti_join(prtx2 %>% select(ID), ptss3 %>% select(ID), by = "ID")

# Create a new dataframe with missing IDs and set MHDEPRESSION and MHANXIETY to -999
missing_data <- missing_ids %>%
  mutate(MHDEPRESSION = -999, MHANXIETY = -999)

# Append the missing IDs to ptss3
ptss4 <- ptss3 %>%
  bind_rows(missing_data)

# pspn2: Extract ECOG and retain RUSUBJ_ID.
pspn2 <- pspn %>%
  transmute(
    ID = RSUBJ_ID,
    VISIT_ID = VISIT_ID,
    ECOG = PSPN_001
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, -999)))%>%
  filter(VISIT_ID == 1)

# Find IDs that are in prtx2 but not in pspn2
missing_ids <- anti_join(prtx2 %>% select(ID), pspn2 %>% select(ID), by = "ID")

# Create a new dataframe with missing IDs and set ECOG to -999
missing_data <- missing_ids %>%
  mutate(ECOG = -999, VISIT_ID = 1)%>%
  mutate(ECOG = as.character(ECOG))

# Append the missing IDs to ptss3
pspn3 <- pspn2 %>%
  bind_rows(missing_data)%>%
  transmute(
    ID = ID,
    ECOG = ECOG
  )

# phex2: Extract WEIGHT and retain RUSUBJ_ID.
phex2 <- phex %>%
  transmute(
    ID = RSUBJ_ID,
    VISIT_ID = VISIT_ID,
    WEIGHT = PHEX_002
  ) %>%
  mutate(across(everything(), ~ replace_na(.x, -999))) %>%
  filter(VISIT_ID == 1)

###delete VISIT_ID
phex3 <- phex2 %>%
  transmute(
    ID = ID,
    WEIGHT = WEIGHT
  ) 

# Finally, merge all datasets by ID
final_data <- demo2 %>%
  left_join(base2, by = "ID") %>%
  left_join(brst2, by = "ID") %>%
  left_join(ptss4, by = "ID") %>%
  left_join(phex3, by = "ID") %>%
  left_join(pspn3, by = "ID") %>%
  left_join(prtx2, by = "ID")

# add survival: DTHDY DTH PFSDY PFS
survive <- surv %>% mutate(ID = RSUBJ_ID, DTHDY = SURVDAY_002, PFSDY = SURVDAY_004) %>%
  mutate(DTH = ifelse(SURV_001 == 3, 1, 0),
         PFS = ifelse(SURV_003 == 1, 1, 0),
         PFSDY = ifelse(is.na(PFSDY), -999, PFSDY)) %>%
  dplyr::select(ID, DTHDY, DTH, PFSDY, PFS)

# Merge final_data with survive, ensuring all IDs are retained
final_data <- final_data %>%
  full_join(survive, by = "ID") %>%
  mutate(across(everything(), ~ replace_na(.x, -999)))  # Replace all NA with -999

# eqlq: QSFLAG FLAG DV TIME B01--B30
# Convert the dataset to long format and replace NA with -999
eqlq_long <- eqlq %>%
  #filter(VISIT_ID == 1) %>%  # Keep only rows where VISIT_ID equals 1
  pivot_longer(cols = EQLQ_001:EQLQ_030,  # Convert only EQLQ_001 to EQLQ_030 into long format
               names_to = "QSFLAG", 
               values_to = "DV") %>%
  mutate(
    FLAG = as.numeric(substr(QSFLAG, 6, 8)),  # Extract numerical question number
    TIME = EQLQDAY_000,  # Assign TIME from EQLQDAY_000
    DV = replace_na(DV, -999)  # Replace NA values in DV with -999
  ) %>%
  select(RSUBJ_ID, QSFLAG, FLAG, TIME, DV)%>% # Keep necessary columns
  distinct()  # Remove duplicated observations

# Convert EQLQ_001 ~ EQLQ_030 to B01 ~ B30 and join with eqlq_long by RSUBJ_ID
eqlq_wide <- eqlq %>%
  filter(VISIT_ID == 1) %>%
  select(RSUBJ_ID, EQLQ_001:EQLQ_030) %>%  # Select only EQLQ_001 ~ EQLQ_030
  rename_with(~ paste0("B", sprintf("%02d", as.numeric(substr(.x, 6, 8)))), starts_with("EQLQ_"))  # Replace NA values with -999

eqlq_wide <- eqlq_wide %>%
  mutate(across(everything(), ~ ifelse(.x == "", NA, .x))) %>%  # Convert empty strings to NA
  mutate(
    across(where(is.numeric), ~ replace_na(.x, -999)),  
    across(where(is.character), ~ replace_na(.x, "-999"))  
  )

# Merge the long format data with the wide format transformed variables
eqlq_final <- eqlq_long %>%
  left_join(eqlq_wide, by = "RSUBJ_ID")%>%
  rename(ID = RSUBJ_ID)

# Replace NA values with -999
eqlq_final <- eqlq_final %>%
  mutate(across(everything(), ~ ifelse(.x == "", NA, .x))) %>%  # Convert empty strings to NA
  mutate(
    across(where(is.numeric), ~ replace_na(.x, -999)),  
    across(where(is.character), ~ replace_na(.x, "-999"))  
  )

###merge data
dataset1 <- final_data %>%
  full_join(eqlq_final, by = "ID") %>%  # Perform a full join to keep all IDs
  mutate(across(everything(), ~ replace_na(.x, -999)))  # Replace all NA values with -999

write.csv(dataset1, file='Breast_Celgene_2001_107_updated.csv', quote=FALSE,row.names=FALSE)





