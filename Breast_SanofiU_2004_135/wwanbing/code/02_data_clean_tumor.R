library(plyr)
library(dplyr)
library(haven)
library(tidyr)

tumasses <- read_sas(data_file='tumasses.sas7bdat')

tumasses2 <- tumasses %>%
  select(
    UID      = RUSUBJID,   # unique subject identifier
    ID       = RSUBJID,   # subject identifier
    TIME     = TADY,    # time in days
    LESIONID = TANUM,     # ID for each tumor
    LESTYPE  = TAREMES,     # Lesion type
    SIZE     = TAMEAS,   # Lesion longest diameter
    ORGAN    = TASITECD       # Lesion organ
  ) %>%
  mutate(
    # hot coding
    LESTYPE = if_else(LESTYPE == "Measurable (target)", 1, 2)
  )

#tumor <- tumasses2 %>%
#  mutate(across(everything(), ~ replace_na(.x, -999)))  # Replace all NA values with -999

tumor <- tumasses2 %>%
  mutate(across(everything(), ~ replace_na(.x, -999))) %>%  # Replace all NA values with -999
  mutate(across(everything(), ~ ifelse(.x == "." | .x == "", -999, .x)))  # Replace "." and " "

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_SanofiU_2004_135_tumor_", toupper(current_date), ".csv")

write.csv(tumor, new_filename, quote=FALSE,row.names=FALSE)

