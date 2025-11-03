###TTFD10: Delete some variables in OS/PFS (OS:stage)

# Significant covariates for PFS: disease stage, age, arm, menopause (clinically important), ECOG, ERS, HER2
# Significant covariates for OS: disease stage, arm, employment status, lesion, ECOG, HER2.
# 
# Run cox proportional hazard model with TTFD (time to event) + selected significant covariates

library(data.table)
library(tidyverse)
library(plyr)
library(dplyr)
library(survival)
library(Matrix)
library(glmnet)
library(GGally)


select = dplyr::select
filter = dplyr::filter
rename = dplyr::rename
mutate = dplyr::mutate
relocate = dplyr::relocate
summarise = dplyr::summarise

###AP
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

####============================================================================####
#### LASSO for OS

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT,TTFD10_AP,censored_AP,TTFD10P_AP,censored10P_AP) %>%
  distinct()%>%
  mutate(
    censored_AP = as.character(censored_AP),
    censored10P_AP = as.character(censored10P_AP)
  )

#OS[OS == -999] <- NA

# rescale continuous variables
set.seed(1234)
OS_scaled <- OS %>%
  # Convert all character to factor
  mutate(across(where(is.character), as.factor)) %>%
  # Scale numeric columns except DTHDY and DTH
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  # Remove rows with NA
  na.omit()

# Significant covariates for OS: disease stage, arm, employment status, lesion, ECOG, HER2.
dat_td <- OS_scaled %>% 
  select(UID, DTHDY, DTH, TTFD10_AP, censored_AP, TTFD10P_AP, censored10P_AP,
         STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)

# Step 1: create tstart and tstop
base <- dat_td %>%
  filter(DTHDY >= 0) %>%     # filter DTHDY >= 0
  mutate(tstart = 0, tstop = DTHDY)

# Step 2: Use tmerge to split time intervals and add time-dependent covariates
td_data <- tmerge(
  data1 = base,
  data2 = base,
  id = UID,
  death = event(tstop, DTH)  # DTH must be numeric
)

td_data <- tmerge(
  td_data,
  base,
  id = UID,
  has_TTFD10 = tdc(TTFD10_AP),
  has_TTFD10P = tdc(TTFD10P_AP)
)

# Step 3: construct Time-dependent Cox model
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 +
                     ARM + EMPLOY + LESION1 + ECOG + HER2,
                   data = td_data)

# Step 4: summarize results
summary(cox_td_OS)

####============================================================================####
#### LASSO for PFS

# Extract necessary variables for PFS analysis
PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_AP,censored_AP,TTFD10P_AP,censored10P_AP) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_AP),
    censored10P = as.character(censored10P_AP)
  )

# Optional: Replace invalid values (-999) with NA if needed
# PFS[PFS == -999] <- NA

# Normalize continuous covariates and convert categorical variables to factor
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()  # Remove rows with missing values

# Select variables for time-dependent modeling
# Significant covariates for PFS: disease stage, age, arm, menopause (clinically important), ECOG, ERS, HER2
dat_td <- PFS_scaled %>% 
  select(UID, PFSDY, PFS, TTFD10_AP, censored_AP, TTFD10P_AP, censored10P_AP,
         STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)

# Step 1: Create tstart and tstop columns
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)

# Step 2: Use tmerge to define survival intervals and time-dependent covariates
td_data <- tmerge(
  data1 = base,
  data2 = base,
  id = UID,
  progression = event(tstop, PFS)  # PFS must be numeric: 1 = event, 0 = censored
)

td_data <- tmerge(
  td_data,
  base,
  id = UID,
  has_TTFD10 = tdc(TTFD10_AP),
  has_TTFD10P = tdc(TTFD10P_AP)
)

# Step 3: Build the time-dependent Cox model
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 +
                      STAGE + AGE+ ARM + MENOS + ECOG + ERS + HER2,
                    data = td_data)

# Step 4: Summarize results
summary(cox_td_PFS)

summary(cox_td_OS)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_AP_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_AP_", date_str, ".rds")))

###CO
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

####============================================================================####
#### LASSO for OS
OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_CO, censored_CO, TTFD10P_CO, censored10P_CO) %>%
  distinct() %>%
  mutate(
    censored_CO = as.character(censored_CO),
    censored10P_CO = as.character(censored10P_CO)
  )

set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()

dat_td <- OS_scaled %>% 
  select(UID, DTHDY, DTH, TTFD10_CO, censored_CO, TTFD10P_CO, censored10P_CO,
         STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)

base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)

td_data <- tmerge(
  data1 = base,
  data2 = base,
  id = UID,
  death = event(tstop, DTH)
)
td_data <- tmerge(
  td_data,
  base,
  id = UID,
  has_TTFD10 = tdc(TTFD10_CO),
  has_TTFD10P = tdc(TTFD10P_CO)
)

cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 +
                     ARM + EMPLOY + LESION1 + ECOG + HER2,
                   data = td_data)
summary(cox_td_OS)

####============================================================================####
#### LASSO for PFS

# Extract necessary variables for PFS analysis
PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_CO, censored_CO, TTFD10P_CO, censored10P_CO) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_CO),
    censored10P = as.character(censored10P_CO)
  )

set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()

dat_td <- PFS_scaled %>% 
  select(UID, PFSDY, PFS, TTFD10_CO, censored_CO, TTFD10P_CO, censored10P_CO,
         STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)

base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)

td_data <- tmerge(
  data1 = base,
  data2 = base,
  id = UID,
  progression = event(tstop, PFS)
)
td_data <- tmerge(
  td_data,
  base,
  id = UID,
  has_TTFD10 = tdc(TTFD10_CO),
  has_TTFD10P = tdc(TTFD10P_CO)
)

cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 +
                      STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2,
                    data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_CO_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_CO_", date_str, ".rds")))

###DI
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_DI, censored_DI, TTFD10P_DI, censored10P_DI) %>%
  distinct() %>%
  mutate(
    censored_DI = as.character(censored_DI),
    censored10P_DI = as.character(censored10P_DI)
  )

set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()

dat_td <- OS_scaled %>% 
  select(UID, DTHDY, DTH, TTFD10_DI, censored_DI, TTFD10P_DI, censored10P_DI,
         STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)

base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)

td_data <- tmerge(
  data1 = base,
  data2 = base,
  id = UID,
  death = event(tstop, DTH)
)
td_data <- tmerge(
  td_data,
  base,
  id = UID,
  has_TTFD10 = tdc(TTFD10_DI),
  has_TTFD10P = tdc(TTFD10P_DI)
)

cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 +
                     ARM + EMPLOY + LESION1 + ECOG + HER2,
                   data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_DI, censored_DI, TTFD10P_DI, censored10P_DI) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_DI),
    censored10P = as.character(censored10P_DI)
  )

set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()

dat_td <- PFS_scaled %>% 
  select(UID, PFSDY, PFS, TTFD10_DI, censored_DI, TTFD10P_DI, censored10P_DI,
         STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)

base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)

td_data <- tmerge(
  data1 = base,
  data2 = base,
  id = UID,
  progression = event(tstop, PFS)
)
td_data <- tmerge(
  td_data,
  base,
  id = UID,
  has_TTFD10 = tdc(TTFD10_DI),
  has_TTFD10P = tdc(TTFD10P_DI)
)

cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 +
                      STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2,
                    data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_DI_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_DI_", date_str, ".rds")))

###DY
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv") 

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_DY, censored_DY, TTFD10P_DY, censored10P_DY) %>%
  distinct() %>%
  mutate(
    censored_DY = as.character(censored_DY),
    censored10P_DY = as.character(censored10P_DY)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_DY, censored_DY, TTFD10P_DY, censored10P_DY, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_DY), has_TTFD10P = tdc(TTFD10P_DY))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_DY, censored_DY, TTFD10P_DY, censored10P_DY) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_DY),
    censored10P = as.character(censored10P_DY)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_DY, censored_DY, TTFD10P_DY, censored10P_DY, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_DY), has_TTFD10P = tdc(TTFD10P_DY))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_DY_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_DY_", date_str, ".rds")))

###FA
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_FA, censored_FA, TTFD10P_FA, censored10P_FA) %>%
  distinct() %>%
  mutate(
    censored_FA = as.character(censored_FA),
    censored10P_FA = as.character(censored10P_FA)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_FA, censored_FA, TTFD10P_FA, censored10P_FA, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_FA), has_TTFD10P = tdc(TTFD10P_FA))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_FA, censored_FA, TTFD10P_FA, censored10P_FA) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_FA),
    censored10P = as.character(censored10P_FA)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_FA, censored_FA, TTFD10P_FA, censored10P_FA, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_FA), has_TTFD10P = tdc(TTFD10P_FA))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_FA_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_FA_", date_str, ".rds")))

##FI
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_FI, censored_FI, TTFD10P_FI, censored10P_FI) %>%
  distinct() %>%
  mutate(
    censored_FI = as.character(censored_FI),
    censored10P_FI = as.character(censored10P_FI)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_FI, censored_FI, TTFD10P_FI, censored10P_FI, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_FI), has_TTFD10P = tdc(TTFD10P_FI))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_FI, censored_FI, TTFD10P_FI, censored10P_FI) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_FI),
    censored10P = as.character(censored10P_FI)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_FI, censored_FI, TTFD10P_FI, censored10P_FI, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_FI), has_TTFD10P = tdc(TTFD10P_FI))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_FI_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_FI_", date_str, ".rds")))


###GHS
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_GHS, censored_GHS, TTFD10P_GHS, censored10P_GHS) %>%
  distinct() %>%
  mutate(
    censored_GHS = as.character(censored_GHS),
    censored10P_GHS = as.character(censored10P_GHS)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_GHS, censored_GHS, TTFD10P_GHS, censored10P_GHS, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_GHS), has_TTFD10P = tdc(TTFD10P_GHS))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_GHS, censored_GHS, TTFD10P_GHS, censored10P_GHS) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_GHS),
    censored10P = as.character(censored10P_GHS)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_GHS, censored_GHS, TTFD10P_GHS, censored10P_GHS, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_GHS), has_TTFD10P = tdc(TTFD10P_GHS))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_GHS_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_GHS_", date_str, ".rds")))

###NV
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_NV, censored_NV, TTFD10P_NV, censored10P_NV) %>%
  distinct() %>%
  mutate(
    censored_NV = as.character(censored_NV),
    censored10P_NV = as.character(censored10P_NV)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_NV, censored_NV, TTFD10P_NV, censored10P_NV, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_NV), has_TTFD10P = tdc(TTFD10P_NV))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_NV, censored_NV, TTFD10P_NV, censored10P_NV) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_NV),
    censored10P = as.character(censored10P_NV)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_NV, censored_NV, TTFD10P_NV, censored10P_NV, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_NV), has_TTFD10P = tdc(TTFD10P_NV))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_NV_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_NV_", date_str, ".rds")))

###PA
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_PA, censored_PA, TTFD10P_PA, censored10P_PA) %>%
  distinct() %>%
  mutate(
    censored_PA = as.character(censored_PA),
    censored10P_PA = as.character(censored10P_PA)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_PA, censored_PA, TTFD10P_PA, censored10P_PA, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_PA), has_TTFD10P = tdc(TTFD10P_PA))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_PA, censored_PA, TTFD10P_PA, censored10P_PA) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_PA),
    censored10P = as.character(censored10P_PA)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_PA, censored_PA, TTFD10P_PA, censored10P_PA, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_PA), has_TTFD10P = tdc(TTFD10P_PA))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_PA_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_PA_", date_str, ".rds")))

###SL
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_SL, censored_SL, TTFD10P_SL, censored10P_SL) %>%
  distinct() %>%
  mutate(
    censored_SL = as.character(censored_SL),
    censored10P_SL = as.character(censored10P_SL)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_SL, censored_SL, TTFD10P_SL, censored10P_SL, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_SL), has_TTFD10P = tdc(TTFD10P_SL))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_SL, censored_SL, TTFD10P_SL, censored10P_SL) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_SL),
    censored10P = as.character(censored10P_SL)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_SL, censored_SL, TTFD10P_SL, censored10P_SL, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_SL), has_TTFD10P = tdc(TTFD10P_SL))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_SL_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_SL_", date_str, ".rds")))

###PF2
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_PF2, censored_PF2, TTFD10P_PF2, censored10P_PF2) %>%
  distinct() %>%
  mutate(
    censored_PF2 = as.character(censored_PF2),
    censored10P_PF2 = as.character(censored10P_PF2)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_PF2, censored_PF2, TTFD10P_PF2, censored10P_PF2, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_PF2), has_TTFD10P = tdc(TTFD10P_PF2))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_PF2, censored_PF2, TTFD10P_PF2, censored10P_PF2) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_PF2),
    censored10P = as.character(censored10P_PF2)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_PF2, censored_PF2, TTFD10P_PF2, censored10P_PF2, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_PF2), has_TTFD10P = tdc(TTFD10P_PF2))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_PF2_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_PF2_", date_str, ".rds")))

###RF2
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_RF2, censored_RF2, TTFD10P_RF2, censored10P_RF2) %>%
  distinct() %>%
  mutate(
    censored_RF2 = as.character(censored_RF2),
    censored10P_RF2 = as.character(censored10P_RF2)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_RF2, censored_RF2, TTFD10P_RF2, censored10P_RF2, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_RF2), has_TTFD10P = tdc(TTFD10P_RF2))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_RF2, censored_RF2, TTFD10P_RF2, censored10P_RF2) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_RF2),
    censored10P = as.character(censored10P_RF2)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_RF2, censored_RF2, TTFD10P_RF2, censored10P_RF2, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_RF2), has_TTFD10P = tdc(TTFD10P_RF2))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_RF2_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_RF2_", date_str, ".rds")))

###EF
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_EF, censored_EF, TTFD10P_EF, censored10P_EF) %>%
  distinct() %>%
  mutate(
    censored_EF = as.character(censored_EF),
    censored10P_EF = as.character(censored10P_EF)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_EF, censored_EF, TTFD10P_EF, censored10P_EF, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_EF), has_TTFD10P = tdc(TTFD10P_EF))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_EF, censored_EF, TTFD10P_EF, censored10P_EF) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_EF),
    censored10P = as.character(censored10P_EF)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_EF, censored_EF, TTFD10P_EF, censored10P_EF, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_EF), has_TTFD10P = tdc(TTFD10P_EF))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_EF_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_EF_", date_str, ".rds")))

###CF
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_CF, censored_CF, TTFD10P_CF, censored10P_CF) %>%
  distinct() %>%
  mutate(
    censored_CF = as.character(censored_CF),
    censored10P_CF = as.character(censored10P_CF)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_CF, censored_CF, TTFD10P_CF, censored10P_CF, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_CF), has_TTFD10P = tdc(TTFD10P_CF))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_CF, censored_CF, TTFD10P_CF, censored10P_CF) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_CF),
    censored10P = as.character(censored10P_CF)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_CF, censored_CF, TTFD10P_CF, censored10P_CF, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_CF), has_TTFD10P = tdc(TTFD10P_CF))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_CF_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_CF_", date_str, ".rds")))

###SF
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

OS <- TTFD %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT, TTFD10_SF, censored_SF, TTFD10P_SF, censored10P_SF) %>%
  distinct() %>%
  mutate(
    censored_SF = as.character(censored_SF),
    censored10P_SF = as.character(censored10P_SF)
  )
set.seed(1234)
OS_scaled <- OS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  na.omit()
dat_td <- OS_scaled %>%
  select(UID, DTHDY, DTH, TTFD10_SF, censored_SF, TTFD10P_SF, censored10P_SF, STAGE, ARM, EMPLOY, LESION1, ECOG, HER2)
base <- dat_td %>%
  filter(DTHDY >= 0) %>%
  mutate(tstart = 0, tstop = DTHDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, death = event(tstop, DTH))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_SF), has_TTFD10P = tdc(TTFD10P_SF))
cox_td_OS <- coxph(Surv(tstart, tstop, DTH) ~ has_TTFD10 + ARM + EMPLOY + LESION1 + ECOG + HER2, data = td_data)
summary(cox_td_OS)

PFS <- TTFD %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT, TTFD10_SF, censored_SF, TTFD10P_SF, censored10P_SF) %>%
  distinct() %>%
  mutate(
    censored = as.character(censored_SF),
    censored10P = as.character(censored10P_SF)
  )
set.seed(1234)
PFS_scaled <- PFS %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  na.omit()
dat_td <- PFS_scaled %>%
  select(UID, PFSDY, PFS, TTFD10_SF, censored_SF, TTFD10P_SF, censored10P_SF, STAGE, AGE, ARM, MENOS, ECOG, ERS, HER2)
base <- dat_td %>%
  filter(PFSDY >= 0) %>%
  mutate(tstart = 0, tstop = PFSDY)
td_data <- tmerge(data1 = base, data2 = base, id = UID, progression = event(tstop, PFS))
td_data <- tmerge(td_data, base, id = UID, has_TTFD10 = tdc(TTFD10_SF), has_TTFD10P = tdc(TTFD10P_SF))
cox_td_PFS <- coxph(Surv(tstart, tstop, PFS) ~ has_TTFD10 + STAGE + AGE + ARM + MENOS + ECOG + ERS + HER2, data = td_data)
summary(cox_td_PFS)

save_path <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Time_dependent_Cox_Model_0703_TTFD10"
date_str <- format(Sys.Date(), "%Y%m%d")

saveRDS(cox_td_OS, file = file.path(save_path, paste0("cox_td_OS_SF_", date_str, ".rds")))
saveRDS(cox_td_PFS, file = file.path(save_path, paste0("cox_td_PFS_SF_", date_str, ".rds")))