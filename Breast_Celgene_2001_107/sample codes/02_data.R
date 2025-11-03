# Generate preliminary data for breast cancer proposal
# study: NCT00046527	Breast_Celgene_2001_107
rm(list=ls())
setwd(' ')
library(plyr)
library(dplyr)
library(haven)
library(tidyr)

eqlq <- read_sas(data_file='eqlq.sas7bdat') #EORTC questionnaire
wclesion <- read_sas(data_file='wclesion.sas7bdat') #CT measured lesion size
surv <- read_sas(data_file='surv.sas7bdat') #survival

# Generate total EQLQ score
eqlq2 <- eqlq %>% filter(EQLQ_031 != 1)
# transform EQLD to numberic
i <- c(6:35)  
eqlq2[, i] <- apply(eqlq2[, i], 2, function(x) as.numeric(as.character(x)))
eqlq2$EQLQ_total <- rowSums(eqlq2[, 6:35])

eqlq3 <- eqlq2 %>% filter(is.na(EQLQ_total) == FALSE) %>% 
  mutate(ID = RSUBJ_ID, QDAY = EQLQDAY_000) %>%
  dplyr::select(ID, EQLQ_total, QDAY) %>% distinct()

# exclude subjects with only one measurement
freq <- eqlq3 %>% group_by(ID) %>%
  summarise(freq = n()) #N=225
exclude.id <- freq %>% filter(freq == 1) #N= 14

# data
eqlq4 <- eqlq3 %>% filter(!c(ID %in% exclude.id$ID))


# add survival
survive <- surv %>% mutate(ID = RSUBJ_ID, DTHDY = SURVDAY_002, PFSDY = SURVDAY_004) %>%
  mutate(DTH = ifelse(SURV_001 == 3, 1, 0),
         PFS = ifelse(SURV_003 == 1, 1, 0)) %>%
  dplyr::select(ID, DTHDY, DTH, PFSDY, PFS)

eqlq5 <- merge(eqlq4, survive, by='ID', all=T) %>% filter(is.na(QDAY) == FALSE) %>% filter(QDAY >-200)%>% arrange(ID, QDAY)

eqlq6 <- eqlq5 %>% group_by(ID) %>%
  mutate(DTHDY2 = ifelse(is.na(DTHDY) == T, last(QDAY), DTHDY),
         DTH2 = ifelse(is.na(DTH) == T, 0, DTH),
         PFSDY2 = ifelse(is.na(PFSDY) == T, last(QDAY), PFSDY),
         PFS2 = ifelse(is.na(PFS) == T, 0, PFS))


baseline <- eqlq6 %>% group_by(ID) %>%
  dplyr::summarise(BSL = first(EQLQ_total), FIRSTTIME = first(QDAY))
eqlq7 <- merge(eqlq6, baseline, by='ID', all=TRUE) %>% mutate(CFB = EQLQ_total - BSL) 
write.csv(eqlq7, file='PROdata.csv', row.names=FALSE, quote=FALSE)


# tumor
lesion <- wclesion %>% mutate(ID = RSUBJ_ID) %>%
  dplyr::select(ID, EXAMDAY, LESTYPE, TGTSIZE, LESLOC, LESNUM) %>% distinct()
write.csv(lesion, file='lesion.csv', row.names=FALSE, quote=FALSE)



