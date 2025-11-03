# Generate preliminary data for breast cancer proposal
rm(list=ls())
setwd(' ') #your path
library(plyr)
library(dplyr)
library(haven)

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






