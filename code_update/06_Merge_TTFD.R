library(data.table)
library(dplyr)
library(purrr)

###Merge TTFD
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates/")
Fatigue<-fread("TTFD+Covariates_FA_04JUN2025.csv")

FA <- Fatigue %>%
  select(-FLAG, -TIME, -STDID, -NCT, -DV2) %>%
  distinct()

NV<-fread("TTFD+Covariates_NV_18JUN2025.csv")
PA<-fread("TTFD+Covariates_PA_18JUN2025.csv")
DY<-fread("TTFD+Covariates_DY_18JUN2025.csv")
SL<-fread("TTFD+Covariates_SL_18JUN2025.csv")
AP<-fread("TTFD+Covariates_AP_18JUN2025.csv")
CO<-fread("TTFD+Covariates_CO_18JUN2025.csv")
DI<-fread("TTFD+Covariates_DI_18JUN2025.csv")
FI<-fread("TTFD+Covariates_FI_18JUN2025.csv")
GHS<-fread("TTFD+Covariates_GHS_18JUN2025.csv")

PF2<-fread("TTFD+Covariates_PF2_22OCT2025.csv")
RF2<-fread("TTFD+Covariates_RF2_22OCT2025.csv")
EF<-fread("TTFD+Covariates_EF_22OCT2025.csv")
CF<-fread("TTFD+Covariates_CF_22OCT2025.csv")
SF<-fread("TTFD+Covariates_SF_22OCT2025.csv")

clean_ttfd_data <- function(df) {
  df %>%
    select(-FLAG, -TIME, -STDID, -NCT, -DV2) %>%
    distinct()
}

NV <- clean_ttfd_data(NV)
PA <- clean_ttfd_data(PA)
DY <- clean_ttfd_data(DY)
SL <- clean_ttfd_data(SL)
AP <- clean_ttfd_data(AP)
CO <- clean_ttfd_data(CO)
DI <- clean_ttfd_data(DI)
FI <- clean_ttfd_data(FI)
GHS <- clean_ttfd_data(GHS)

PF2 <- clean_ttfd_data(PF2)
RF2 <- clean_ttfd_data(RF2)
EF <- clean_ttfd_data(EF)
CF <- clean_ttfd_data(CF)
SF <- clean_ttfd_data(SF)

FA <- dplyr::rename(FA,
                    TTFD10_FA = TTFD10,
                    censored_FA = censored,
                    TTFD10P_FA = TTFD10P,
                    censored10P_FA = censored10P
)

rename_ttfds <- function(df, suffix) {
  names(df)[names(df) == "TTFD10"] <- paste0("TTFD10_", suffix)
  names(df)[names(df) == "censored"] <- paste0("censored_", suffix)
  names(df)[names(df) == "TTFD10P"] <- paste0("TTFD10P_", suffix)
  names(df)[names(df) == "censored10P"] <- paste0("censored10P_", suffix)
  return(df)
}

NV <- rename_ttfds(NV, "NV")
PA <- rename_ttfds(PA, "PA")
DY <- rename_ttfds(DY, "DY")
SL <- rename_ttfds(SL, "SL")
AP <- rename_ttfds(AP, "AP")
CO <- rename_ttfds(CO, "CO")
DI <- rename_ttfds(DI, "DI")
FI <- rename_ttfds(FI, "FI")
GHS <- rename_ttfds(GHS, "GHS")

PF2 <- rename_ttfds(PF2, "PF2")
RF2 <- rename_ttfds(RF2, "RF2")
EF <- rename_ttfds(EF, "EF")
CF <- rename_ttfds(CF, "CF")
SF <- rename_ttfds(SF, "SF")

#Take the last four columns
join_last4 <- function(main_df, df_to_join) {
  suffix <- sub(".*_", "", names(df_to_join)[ncol(df_to_join)])  
  last4_cols <- c("UID", tail(colnames(df_to_join), 4))
  main_df <- left_join(main_df, df_to_join[, last4_cols, with = FALSE], by = "UID")
  return(main_df)
}

AP_df <- as.data.frame(AP)

merged_final <- AP_df %>%
  join_last4(CO) %>%
  join_last4(DI) %>%
  join_last4(DY) %>%
  join_last4(FA) %>%
  join_last4(FI) %>%
  join_last4(GHS) %>%
  join_last4(NV) %>%
  join_last4(PA) %>%
  join_last4(SL) %>%
  join_last4(PF2) %>%
  join_last4(RF2) %>%
  join_last4(EF) %>%
  join_last4(CF) %>%
  join_last4(SF) 

merged_final[is.na(merged_final)] <- -999

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Merged_PRO_TTFD_", toupper(current_date), ".csv")

write.csv(merged_final, new_filename, quote = FALSE, row.names = FALSE)

### 10/01/2025 add these codes to change DTHDY/DTH/PFSPY/PFS (we keep original data for them as Merge_EQLQ_24APR.CSV) 
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output")
eqlq<-fread("Merge_EQLQ_24APR2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
Merge_TTFD<-fread("Merged_PRO_TTFD_22OCT2025.csv")

eqlq_4 <- eqlq %>%
  select(UID, DTHDY, DTH, PFSDY, PFS) %>%
  distinct(UID, .keep_all = TRUE)

Merge_TTFD_1 <- Merge_TTFD %>%
  left_join(eqlq_4, by = "UID", suffix = c("", ".eqlq")) %>%
  mutate(
    DTHDY = DTHDY.eqlq,
    DTH   = DTH.eqlq,
    PFSDY = PFSDY.eqlq,
    PFS   = PFS.eqlq
  ) %>%
  select(-ends_with(".eqlq"))

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Merged_PRO_TTFD_", toupper(current_date), ".csv")

write.csv(Merge_TTFD_1, new_filename, quote = FALSE, row.names = FALSE)


