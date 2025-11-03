rm(list=ls())

library(plyr)
library(dplyr)
library(haven)

adpea <- read_sas(data_file='adpea.sas7bdat') #disease free survival

bcr <- adpea %>%
  transmute(
    UID     = RUSUBJID,
    ID      = RSUBJID,
    BCRYN   = BCRYN,
    DISTANT = DISTANT,
    REGIONAL = REGIONAL,
    LOCAL   = LOCAL
  )

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_SanofiU_1997_120_BCR_", toupper(current_date), ".csv")

write.csv(bcr, new_filename, quote=FALSE,row.names=FALSE)