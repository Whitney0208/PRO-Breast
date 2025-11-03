## Project: Sanofi_2007_133 NSCLC dataset, get tumor, PRO, survival, appetite, weight
## Objective: Generate Covariate dataset for PRO outcomes
## Author: Jiawei Zhou
## Date: Jan 21, 2025
rm(list=ls())
setwd(' ')
library(plyr)
library(dplyr)
library(haven)
library(tidyr)

# 1. load demographics information
dm <- read_sas(data_file='dm.sas7bdat') # Important, demographics
length(unique(dm$RSUBJID))#N=455
dm2 <- dm %>% mutate(
  SUBJID = as.numeric(as.character(RSUBJID)),
  AGE = as.numeric(as.character(AGE)),
  SEX = factor(SEX, levels=c('M','F'), labels=c(0, 1)),
  RACE = factor(RACE, levels=c('CAUCASIAN/WHITE','OTHER'), labels=c(1,5)),
  REGION = factor(REGION, levels=c('Asia','Eastern Europe','North America','Other countries','Western Europe'), 
                  labels=c(1,2,3,4,5)),
  ARM = factor(ARMCD, levels=c('PLACEBO'), labels=c(0))
) %>% select(RUSUBJID, SUBJID, AGE, SEX, RACE, REGION, ARM) %>% distinct()

# 2. add disease stage
cd <- read_sas(data_file='cd.sas7bdat') # disease stage
cd2 <- cd %>% filter(CDTESTCD == 'STAGE') %>% select(RUSUBJID, CDSTRESC) %>% 
      mutate(STAGE = factor(CDSTRESC, levels=c('STAGE I','STAGE II',
                                               'STAGE III A','STAGE III B','STAGE IV',
                                               'UNKNOWN','UNKNOWN - M MISSING'),
                            labels = c(1, 2, 3, 3, 4, -999, -999))) %>% select(RUSUBJID, STAGE) %>% distinct()
#merge to demographics
dat1 <- merge(dm2, cd2, by='RUSUBJID')

# 3. add survival data
ds <- read_sas(data_file='ds.sas7bdat') 
# death day
dth <- ds %>% group_by(RUSUBJID) %>%
  summarise(STATE = last(DSSCAT), DTHWK = last(DSSTHWK), ALIVEDY = last(DSSTDY)) %>%
  mutate(DTH = ifelse(STATE == 'DEATH', 1,
                ifelse(STATE == 'LAST CONTACT', 0, -999))) %>%
  mutate(DTHDY = ifelse(DTH == 1, DTHWK*7, 
                 ifelse(DTH == 0, ALIVEDY, -999))) %>%
  select(RUSUBJID, DTH, DTHDY) %>% distinct()

#merge data
dat2 <- merge(dat1, dth, by='RUSUBJID')

# select disease progression time
pfs <- ds %>% filter(DSDECOD == 'DISEASE PROGRESSION') %>%
  group_by(RUSUBJID) %>%
  summarise(STATE = first(DSDECOD), STATE2 = first(DSSCAT), PFST = first(DSSTHWK)*7)
trtend.time <- dm %>% select(RUSUBJID, RFENDY)
pfs2 <- merge(pfs, trtend.time, by='RUSUBJID', all='T') %>%
  mutate(PFS = ifelse(is.na(STATE) == T, 0, 1)) %>%
  mutate(PFSDY = ifelse(PFS == 0, RFENDY,
                 ifelse(PFS == 1 & STATE2 == 'END OF TREATMENT', RFENDY, PFST))) %>%
  select(RUSUBJID, PFS, PFSDY)

## merge data
dat3 <- merge(dat2, pfs2, by='RUSUBJID')

#  4.liver and kidney lab values
lb <- read_sas(data_file='lb.sas7bdat') # lab test ***
lb2 <- lb %>% filter(LBBLFL == 'Y')
#WBC
WBC <- lb2 %>% filter(LBTESTCD == 'WBC' & is.na(LBSTRESN) == FALSE) %>% select(RUSUBJID, LBSTRESN) %>% dplyr::rename(WBC = LBSTRESN)
dat4 <- merge(dat3, WBC, by='RUSUBJID', all=T)
#NEUT
NEUT <- lb2 %>% filter(LBTESTCD == 'NEUT' & is.na(LBSTRESN) == FALSE) %>% select(RUSUBJID, LBSTRESN) %>% dplyr::rename(NEUT = LBSTRESN)
dat4 <- merge(dat4, NEUT, by='RUSUBJID', all=T)
# CREAT
CREAT <- lb2 %>% filter(LBTESTCD == 'CREAT' & is.na(LBSTRESN) == FALSE & LBCAT =='BIOCHEMISTRY') %>% select(RUSUBJID, LBSTRESN) %>% dplyr::rename(CREAT = LBSTRESN)
dat4 <- merge(dat4, CREAT, by='RUSUBJID', all=T)
dat4[is.na(dat4$CREAT) ==T,]$CREAT <- -999
# ALT
ALT <- lb2 %>% filter(LBTESTCD == 'ALT' & is.na(LBSTRESN) == FALSE) %>% select(RUSUBJID, LBSTRESN) %>% dplyr::rename(ALT = LBSTRESN)
dat4 <- merge(dat4, ALT, by='RUSUBJID', all=T)
dat4[is.na(dat4$ALT) ==T,]$ALT <- -999
# BILI
BILI <- lb2 %>% filter(LBTESTCD == 'BILI' & is.na(LBSTRESN) == FALSE) %>% select(RUSUBJID, LBSTRESN) %>% dplyr::rename(BILI = LBSTRESN)
dat4 <- merge(dat4, BILI, by='RUSUBJID', all=T)
dat4[is.na(dat4$BILI) ==T,]$BILI <- -999
# ALP
ALP <- lb2 %>% filter(LBTESTCD == 'ALP' & is.na(LBSTRESN) == FALSE) %>% select(RUSUBJID, LBSTRESN) %>% dplyr::rename(ALP = LBSTRESN)
dat4 <- merge(dat4, ALP, by='RUSUBJID', all=T)
dat4[is.na(dat4$ALP) ==T,]$ALP <- -999
# AST
AST <- lb2 %>% filter(LBTESTCD == 'AST' & is.na(LBSTRESN) == FALSE) %>% select(RUSUBJID, LBSTRESN) %>% dplyr::rename(AST = LBSTRESN)
dat4 <- merge(dat4, AST, by='RUSUBJID', all=T)
dat4[is.na(dat4$AST) ==T,]$AST <- -999
# ALB
ALB <- lb2 %>% filter(LBTESTCD == 'ALB' & is.na(LBSTRESN) == FALSE& LBCAT =='BIOCHEMISTRY') %>% select(RUSUBJID, LBSTRESN) %>% dplyr::rename(ALB = LBSTRESN) %>% distinct()
dat4 <- merge(dat4, ALB, by='RUSUBJID', all=T)
dat4[is.na(dat4$ALB) ==T,]$ALB <- -999

# 5. anxiety and depression
dat5 <- dat4
mh <- read_sas(data_file='mh.sas7bdat') # medical history
depres <- mh %>% filter(MHDECOD %in% c('Depression')) %>% select(RUSUBJID) %>% distinct()
dat5$MHDEPRESSION <- 0
dat5[dat5$RUSUBJID %in% depres$RUSUBJID,]$MHDEPRESSION <- 1

anxiety <- mh %>% filter(MHDECOD %in% c('Anxiety','Anxiety disorder')) %>% select(RUSUBJID) %>% distinct()
dat5$MHANXIETY <- 0
dat5[dat5$RUSUBJID %in% anxiety$RUSUBJID,]$MHANXIETY <- 1

# smoking
su <- read_sas(data_file='su.sas7bdat')
su2 <- su %>% select(RUSUBJID, SUSMKST) %>% distinct() %>% 
  mutate(SUOCCUR =  ifelse(SUSMKST == 'NEVER', 0,
                    ifelse(SUSMKST == 'CURRENT', 1, 2))) %>% select(RUSUBJID, SUOCCUR)
dat6 <- merge(dat5, su2, by='RUSUBJID', all=T)

# height
vs <- read_sas(data_file='vs.sas7bdat')
height <- vs %>% filter(VSTESTCD == 'HEIGHT') %>% select(RUSUBJID, VSSTRESN) %>% dplyr::rename(HEIGHT = VSSTRESN) %>% distinct()
dat7 <- merge(dat6, height, by='RUSUBJID', all=T)

# weight
weight <- vs %>% filter(VSTESTCD == 'WEIGHT' & VSBLFL == 'Y') %>% select(RUSUBJID, VSSTRESN) %>% dplyr::rename(WEIGHT = VSSTRESN) %>% distinct()
dat8 <- merge(dat7, weight, by='RUSUBJID', all=T)

#ECOG
ECOG <- vs %>% filter(VSTESTCD == 'ECOG' & VSBLFL == 'Y')%>% select(RUSUBJID, VSSTRESN) %>% dplyr::rename(ECOG = VSSTRESN) %>% distinct()
dat9 <- merge(dat8, ECOG, by='RUSUBJID', all=T)

write.csv(dat9, file='Cov_133data_21JAN2025.csv', quote=FALSE, row.names=FALSE)
