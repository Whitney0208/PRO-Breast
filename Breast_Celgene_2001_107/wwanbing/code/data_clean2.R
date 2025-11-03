###process surv and eqlq datset
library(plyr)
library(dplyr)
library(haven)
library(tidyr)

eqlq <- read_sas(data_file='eqlq.sas7bdat') #EORTC questionnaire
surv <- read_sas(data_file='surv.sas7bdat') #survival

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
  filter(VISIT_ID == 1) %>%  # Keep only rows where VISIT_ID equals 1
  pivot_longer(cols = EQLQ_001:EQLQ_030,  # Convert only EQLQ_001 to EQLQ_030 into long format
               names_to = "QSFLAG", 
               values_to = "DV") %>%
  mutate(
    FLAG = as.numeric(substr(QSFLAG, 6, 8)),  # Extract numerical question number
    TIME = EQLQDAY_000,  # Assign TIME from EQLQDAY_000
    DV = replace_na(DV, -999)  # Replace NA values in DV with -999
  ) %>%
  select(RSUBJ_ID, QSFLAG, FLAG, TIME, DV)  # Keep necessary columns

# Convert EQLQ_001 ~ EQLQ_030 to B01 ~ B30 and join with eqlq_long by RSUBJ_ID
eqlq_wide <- eqlq %>%
  filter(VISIT_ID == 1) %>%
  select(RSUBJ_ID, EQLQ_001:EQLQ_030) %>%  # Select only EQLQ_001 ~ EQLQ_030
  rename_with(~ paste0("B", sprintf("%02d", as.numeric(substr(.x, 6, 8)))), starts_with("EQLQ_"))  # Replace NA values with -999

# Merge the long format data with the wide format transformed variables
eqlq_final <- eqlq_long %>%
  left_join(eqlq_wide, by = "RSUBJ_ID")%>%
  rename(ID = RSUBJ_ID)

# Replace NA values with -999
eqlq_final <- eqlq_final %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, -999))) %>%  # Replace NAs in numeric columns
  mutate(across(where(is.character), ~ replace_na(.x, "-999")))  # Replace NAs in character columns











































