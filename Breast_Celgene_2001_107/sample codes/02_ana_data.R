## Project: explore Merck dataset NCT00409188
## Objective: Generate Analysis dataset for PRO outcomes
## Author: Jiawei Zhou
## Date: Sept 15, 2024
rm(list=ls())
setwd('C:\\JZ\\UNC_faculty\\Research\\Patient Centric\\PRO Oncology\\Merck_PDS\\NCT00409188\\Data')
library(plyr)
library(dplyr)
library(haven)
library(tidyr)

# Load covariate data
load('dat7.RData')

# Load analysis data
qs <- read_sas(data_file='qs.sas7bdat')
qs2 <- qs %>% select(USUBJID, QSDY, QSBLFL, QSTEST, QSCAT, QSSCAT, QSSTRESN) %>%
  filter(is.na(QSSTRESN) == F) %>%
  mutate(USUBJID = as.numeric(as.character(USUBJID)),
  QSFLAG = ifelse(QSTEST == 'Mobility', 'EQ5D1',
           ifelse(QSTEST == 'Self-Care', 'EQ5D2',
           ifelse(QSTEST == 'Usual Activities', 'EQ5D3',
           ifelse(QSTEST == 'Pain/Discomfort', 'EQ5D4',
           ifelse(QSTEST == 'Anxiety/Depression', 'EQ5D5',
           ifelse(QSTEST == 'Best Imaginable Health State', 'EQ5D6',
           ifelse(QSTEST == 'Loss of Appetite (O)','LCSSO1',
           ifelse(QSTEST == 'Fatigue (O)', 'LCSSO2',
           ifelse(QSTEST == 'Cough (O)', 'LCSSO3',
           ifelse(QSTEST == 'Dyspnea (O)', 'LCSSO4',
           ifelse(QSTEST == 'Hemoptysis (O)', 'LCSSO5',
           ifelse(QSTEST == 'Pain (O)', 'LCSSO6',
           ifelse(QSTEST == 'Loss of Appetite (S)', 'LCSSS1',
           ifelse(QSTEST == 'Fatigue (S)', 'LCSSS2',
           ifelse(QSTEST == 'Cough (S)', 'LCSSS3',
           ifelse(QSTEST == 'Dyspnea (S)', 'LCSSS4',
           ifelse(QSTEST == 'Hemoptysis (S)', 'LCSSS5',
           ifelse(QSTEST == 'Pain (S)', 'LCSSS6',
           ifelse(QSTEST == 'Lung Cancer Symptoms (S)', 'LCSSS7',
           ifelse(QSTEST == 'Ability to Carry out Activities (S)', 'LCSSS8',
           ifelse(QSTEST == 'Quality of Life Today (S)', 'LCSSS9', 'Unknown'))))))))))))))))))))))

bsl.qs <- qs2 %>% filter(QSBLFL == 'Y') %>%
  select(USUBJID, QSFLAG, QSSTRESN) %>%
  group_by(USUBJID) %>%
  pivot_wider(names_from = QSFLAG, names_glue = "B_{QSFLAG}",values_from = QSSTRESN)

dat8 <- merge(dat7, bsl.qs, by='USUBJID', all=T)
dat8 <- dat8 %>% replace(is.na(.), -999)
save(dat8, file='dat8.RData')

# generate analysis data
qs3 <- qs2 %>% select(USUBJID, QSFLAG, QSDY, QSSTRESN) %>%
  mutate(TIME = QSDY, DV = QSSTRESN,
         FLAG = as.numeric(factor(QSFLAG))) %>%
  select(USUBJID, QSFLAG, FLAG, TIME, DV)

# merge to final dataset
dat9 <- merge(dat8,qs3, by='USUBJID', all=T )
dat9 <- dat9 %>% replace(is.na(.), -999)

# Add C column
dat9$C <- '.'
dat9[dat9$TIME == -999,]$C <- 'C'
dat9$MDV <- 0
dat9[dat9$DV == -999,]$MDV <- 1
dat9 <- dat9 %>% relocate(C,.before = USUBJID)
save(dat9, file='dat9.RData')
write.csv(dat9, file='Merck_Lung_9188_Data_15SEP2024.csv', quote=FALSE, row.names=FALSE)

