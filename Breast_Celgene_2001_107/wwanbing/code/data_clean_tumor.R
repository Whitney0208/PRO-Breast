library(plyr)
library(dplyr)
library(haven)
library(tidyr)

wclesion <- read_sas(data_file='wclesion.sas7bdat') #CT measured lesion size

wclesion2 <- wclesion %>%
  select(
    UID     = RSUBJ_ID,   # unique subject identifier
    ID       = RSUBJ_ID,   # subject identifier
    TIME     = EXAMDAY,    # time in days
    LESIONID = LESNUM,     # ID for each tumor
    LESTYPE  = LESTYPE,     # Lesion type
    SIZE     = TGTSIZE,   # Lesion longest diameter
    ORGAN    = LESLOC       # Lesion organ
  ) %>%
  mutate(
    # hot coding
    LESTYPE = if_else(LESTYPE == "Target", 1, 0)
  )

tumor <- wclesion2 %>%
  mutate(across(everything(), ~ replace_na(.x, -999)))  # Replace all NA values with -999

write.csv(tumor, file='Breast_Celgene_2001_107_tumor.csv', quote=FALSE,row.names=FALSE)

