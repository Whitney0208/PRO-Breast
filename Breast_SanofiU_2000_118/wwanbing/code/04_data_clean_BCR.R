rm(list=ls())
library(plyr)
library(dplyr)
library(haven)

bcrelap <- read_sas(data_file='bcrelap.sas7bdat') # relapse day, regional or distant, which organ
demog <- read_sas(data_file='demog.sas7bdat')

# 1. create BCRYN, filter BCR1
bcryn_df <- bcrelap %>%
  filter(EVENT_ID == "BCR1") %>%
  distinct(RSUBJID) %>%
  mutate(BCRYN = 1)

# 2. combine two datasets
demog_with_bcryn <- demog %>%
  select(RUSUBJID, RSUBJID) %>%
  left_join(bcryn_df, by = "RSUBJID") %>%
  mutate(BCRYN = ifelse(is.na(BCRYN), 0, BCRYN))%>%
  distinct()

# 3. Coding Yes=1 No=0
final_df <- demog_with_bcryn %>%
  left_join(bcrelap %>% 
              filter(EVENT_ID == "BCR1") %>%
              select(RSUBJID, 
                     DISTANT = DSRELAP, 
                     REGIONAL = RGRELAP, 
                     LOCAL = LCRELAP),
            by = "RSUBJID") %>%
  mutate(across(c(DISTANT, REGIONAL, LOCAL), ~ case_when(
    . == "Yes" ~ 1,
    . == "No" ~ 0,
    is.na(.) | . == "" | . == "." ~ -999,
    TRUE ~ -999
  )))%>%
  distinct()%>%
  rename(
    UID = RUSUBJID,
    ID = RSUBJID)%>%
  mutate(across(c(DISTANT, REGIONAL, LOCAL), ~ ifelse(. == -999, 0, .)))

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_SanofiU_2000_118_BCR_", toupper(current_date), ".csv")

write.csv(final_df, new_filename, quote=FALSE,row.names=FALSE)