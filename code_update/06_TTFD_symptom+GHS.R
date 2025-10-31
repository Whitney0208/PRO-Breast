###FA
library(data.table)
library(dplyr)
library(purrr)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
Fatigue<-fread("Breast_PRO_Fatigue_updated_04JUN2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- Fatigue %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  summarise(
    # Define baseline DV2 score as the first time point
    baseline = first(DV2),
    
    # --- TTFD10: First time DV2 >= baseline + 10 ---
    TTFD10 = {
      idx <- which(DV2 >= baseline + 10)
      if (length(idx) > 0) TIME[idx[1]] else max(TIME)
    },
    
    # Censoring indicator for TTFD10 (1 = censored, 0 = not censored)
    censored = ifelse(any(DV2 >= baseline + 10), 0, 1),
    
    # --- TTFD10P: First time the difference (DV2 - baseline) increases by >= 10 between two consecutive visits ---
    TTFD10P = {
      delta <- DV2 - baseline            # Difference from baseline for all time points
      ddiff <- diff(delta)              # Difference between consecutive delta values
      if (any(ddiff >= 10)) {
        i <- which(ddiff >= 10)[1]      # Get the index of the first jump >= 10
        TIME[i + 1]                     # The corresponding TIME (i + 1 because diff shifts index)
      } else {
        max(TIME)                       # If no such jump, use last available TIME
      }
    },
    
    # Censoring indicator for TTFD10P
    censored10P = ifelse(any(diff(DV2 - baseline) >= 10), 0, 1)
  )


###export TTFD result
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_FA_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
Fatigue2 <- Fatigue %>%
  left_join(
    ttfd_result %>% 
      select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_FA_", toupper(current_date), ".csv")

write.csv(Fatigue2, new_filename, quote=FALSE,row.names=FALSE)

###NV
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
NV<-fread("Breast_PRO_Nausea_and_vomiting_updated_18JUN2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- NV %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  group_modify(~{
    baseline <- first(.x$DV2)
    
    # TTFD10
    idx <- which(.x$DV2 >= baseline + 10)
    TTFD10 <- if (length(idx) > 0) .x$TIME[idx[1]] else max(.x$TIME)
    censored <- ifelse(length(idx) > 0, 0, 1)
    
    # TTFD10P
    delta <- .x$DV2 - baseline
    ddiff <- diff(delta)
    if (any(ddiff >= 10)) {
      i <- which(ddiff >= 10)[1]
      TTFD10P <- .x$TIME[i + 1]
      censored10P <- 0
    } else {
      TTFD10P <- max(.x$TIME)
      censored10P <- 1
    }
    
    tibble(
      baseline = baseline,
      TTFD10 = TTFD10,
      censored = censored,
      TTFD10P = TTFD10P,
      censored10P = censored10P
    )
  }) %>% ungroup()

###export TTFD result
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD")
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_NV_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
NV2 <- NV %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_NV_", toupper(current_date), ".csv")

write.csv(NV2, new_filename, quote=FALSE,row.names=FALSE)

###PA
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
PA<-fread("Breast_PRO_Pain_updated_18JUN2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- PA %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  group_modify(~{
    baseline <- first(.x$DV2)
    
    # TTFD10
    idx <- which(.x$DV2 >= baseline + 10)
    TTFD10 <- if (length(idx) > 0) .x$TIME[idx[1]] else max(.x$TIME)
    censored <- ifelse(length(idx) > 0, 0, 1)
    
    # TTFD10P
    delta <- .x$DV2 - baseline
    ddiff <- diff(delta)
    if (any(ddiff >= 10)) {
      i <- which(ddiff >= 10)[1]
      TTFD10P <- .x$TIME[i + 1]
      censored10P <- 0
    } else {
      TTFD10P <- max(.x$TIME)
      censored10P <- 1
    }
    
    tibble(
      baseline = baseline,
      TTFD10 = TTFD10,
      censored = censored,
      TTFD10P = TTFD10P,
      censored10P = censored10P
    )
  }) %>% ungroup()

###export TTFD result
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD")
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_PA_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
PA2 <- PA %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_PA_", toupper(current_date), ".csv")

write.csv(PA2, new_filename, quote=FALSE,row.names=FALSE)


###DY
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
DY<-fread("Breast_PRO_Dyspnoea_updated_18JUN2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- DY %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  group_modify(~{
    baseline <- first(.x$DV2)
    
    # TTFD10
    idx <- which(.x$DV2 >= baseline + 10)
    TTFD10 <- if (length(idx) > 0) .x$TIME[idx[1]] else max(.x$TIME)
    censored <- ifelse(length(idx) > 0, 0, 1)
    
    # TTFD10P
    delta <- .x$DV2 - baseline
    ddiff <- diff(delta)
    if (any(ddiff >= 10)) {
      i <- which(ddiff >= 10)[1]
      TTFD10P <- .x$TIME[i + 1]
      censored10P <- 0
    } else {
      TTFD10P <- max(.x$TIME)
      censored10P <- 1
    }
    
    tibble(
      baseline = baseline,
      TTFD10 = TTFD10,
      censored = censored,
      TTFD10P = TTFD10P,
      censored10P = censored10P
    )
  }) %>% ungroup()

###export TTFD result
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD")
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_DY_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
DY2 <- DY %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_DY_", toupper(current_date), ".csv")

write.csv(DY2, new_filename, quote=FALSE,row.names=FALSE)

###SL
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
SL<-fread("Breast_PRO_Insomnia_updated_18JUN2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- SL %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  group_modify(~{
    baseline <- first(.x$DV2)
    
    # TTFD10
    idx <- which(.x$DV2 >= baseline + 10)
    TTFD10 <- if (length(idx) > 0) .x$TIME[idx[1]] else max(.x$TIME)
    censored <- ifelse(length(idx) > 0, 0, 1)
    
    # TTFD10P
    delta <- .x$DV2 - baseline
    ddiff <- diff(delta)
    if (any(ddiff >= 10)) {
      i <- which(ddiff >= 10)[1]
      TTFD10P <- .x$TIME[i + 1]
      censored10P <- 0
    } else {
      TTFD10P <- max(.x$TIME)
      censored10P <- 1
    }
    
    tibble(
      baseline = baseline,
      TTFD10 = TTFD10,
      censored = censored,
      TTFD10P = TTFD10P,
      censored10P = censored10P
    )
  }) %>% ungroup()

###export TTFD result
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD")
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_SL_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
SL2 <- SL %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_SL_", toupper(current_date), ".csv")

write.csv(SL2, new_filename, quote=FALSE,row.names=FALSE)


###AP
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
AP<-fread("Breast_PRO_Appetite_loss_updated_18JUN2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- AP %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  group_modify(~{
    baseline <- first(.x$DV2)
    
    # TTFD10
    idx <- which(.x$DV2 >= baseline + 10)
    TTFD10 <- if (length(idx) > 0) .x$TIME[idx[1]] else max(.x$TIME)
    censored <- ifelse(length(idx) > 0, 0, 1)
    
    # TTFD10P
    delta <- .x$DV2 - baseline
    ddiff <- diff(delta)
    if (any(ddiff >= 10)) {
      i <- which(ddiff >= 10)[1]
      TTFD10P <- .x$TIME[i + 1]
      censored10P <- 0
    } else {
      TTFD10P <- max(.x$TIME)
      censored10P <- 1
    }
    
    tibble(
      baseline = baseline,
      TTFD10 = TTFD10,
      censored = censored,
      TTFD10P = TTFD10P,
      censored10P = censored10P
    )
  }) %>% ungroup()

###export TTFD result
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD")
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_AP_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
AP2 <- AP %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_AP_", toupper(current_date), ".csv")

write.csv(AP2, new_filename, quote=FALSE,row.names=FALSE)


###CO
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
CO<-fread("Breast_PRO_Constipation_updated_18JUN2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- CO %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  group_modify(~{
    baseline <- first(.x$DV2)
    
    # TTFD10
    idx <- which(.x$DV2 >= baseline + 10)
    TTFD10 <- if (length(idx) > 0) .x$TIME[idx[1]] else max(.x$TIME)
    censored <- ifelse(length(idx) > 0, 0, 1)
    
    # TTFD10P
    delta <- .x$DV2 - baseline
    ddiff <- diff(delta)
    if (any(ddiff >= 10)) {
      i <- which(ddiff >= 10)[1]
      TTFD10P <- .x$TIME[i + 1]
      censored10P <- 0
    } else {
      TTFD10P <- max(.x$TIME)
      censored10P <- 1
    }
    
    tibble(
      baseline = baseline,
      TTFD10 = TTFD10,
      censored = censored,
      TTFD10P = TTFD10P,
      censored10P = censored10P
    )
  }) %>% ungroup()

###export TTFD result
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD")
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_CO_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
CO2 <- CO %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_CO_", toupper(current_date), ".csv")

write.csv(CO2, new_filename, quote=FALSE,row.names=FALSE)


###DI
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
DI<-fread("Breast_PRO_Diarrhoea_updated_18JUN2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- DI %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  group_modify(~{
    baseline <- first(.x$DV2)
    
    # TTFD10
    idx <- which(.x$DV2 >= baseline + 10)
    TTFD10 <- if (length(idx) > 0) .x$TIME[idx[1]] else max(.x$TIME)
    censored <- ifelse(length(idx) > 0, 0, 1)
    
    # TTFD10P
    delta <- .x$DV2 - baseline
    ddiff <- diff(delta)
    if (any(ddiff >= 10)) {
      i <- which(ddiff >= 10)[1]
      TTFD10P <- .x$TIME[i + 1]
      censored10P <- 0
    } else {
      TTFD10P <- max(.x$TIME)
      censored10P <- 1
    }
    
    tibble(
      baseline = baseline,
      TTFD10 = TTFD10,
      censored = censored,
      TTFD10P = TTFD10P,
      censored10P = censored10P
    )
  }) %>% ungroup()

###export TTFD result
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD")
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_DI_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
DI2 <- DI %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_DI_", toupper(current_date), ".csv")

write.csv(DI2, new_filename, quote=FALSE,row.names=FALSE)

###FI
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
FI<-fread("Breast_PRO_Financial_difficulties_updated_18JUN2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- FI %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  group_modify(~{
    baseline <- first(.x$DV2)
    
    # TTFD10
    idx <- which(.x$DV2 >= baseline + 10)
    TTFD10 <- if (length(idx) > 0) .x$TIME[idx[1]] else max(.x$TIME)
    censored <- ifelse(length(idx) > 0, 0, 1)
    
    # TTFD10P
    delta <- .x$DV2 - baseline
    ddiff <- diff(delta)
    if (any(ddiff >= 10)) {
      i <- which(ddiff >= 10)[1]
      TTFD10P <- .x$TIME[i + 1]
      censored10P <- 0
    } else {
      TTFD10P <- max(.x$TIME)
      censored10P <- 1
    }
    
    tibble(
      baseline = baseline,
      TTFD10 = TTFD10,
      censored = censored,
      TTFD10P = TTFD10P,
      censored10P = censored10P
    )
  }) %>% ungroup()

###export TTFD result
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD")
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_FI_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
FI2 <- FI %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_FI_", toupper(current_date), ".csv")

write.csv(FI2, new_filename, quote=FALSE,row.names=FALSE)


###GHS
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/updated_PFS")
GHS<-fread("Breast_PRO_GHS_updated_18JUN2025.csv")

ttfd_result <- GHS %>%
  group_by(UID) %>%
  arrange(TIME) %>%
  group_modify(~{
    baseline <- first(.x$DV2)
    
    # --- TTFD10: First time DV2 <= baseline - 10 ---
    idx <- which(.x$DV2 <= baseline - 10)
    TTFD10 <- if (length(idx) > 0) .x$TIME[idx[1]] else max(.x$TIME)
    censored <- ifelse(length(idx) > 0, 0, 1)
    
    # --- TTFD10P: Drop of >=10 between consecutive visits ---
    delta <- .x$DV2 - baseline
    ddiff <- diff(delta)
    if (any(ddiff <= -10)) {
      i <- which(ddiff <= -10)[1]
      TTFD10P <- .x$TIME[i + 1]
      censored10P <- 0
    } else {
      TTFD10P <- max(.x$TIME)
      censored10P <- 1
    }
    
    tibble(
      baseline = baseline,
      TTFD10 = TTFD10,
      censored = censored,
      TTFD10P = TTFD10P,
      censored10P = censored10P
    )
  }) %>% ungroup()

###export TTFD result
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD")
current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD_GHS_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
GHS2 <- GHS %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_GHS_", toupper(current_date), ".csv")

write.csv(GHS2, new_filename, quote=FALSE,row.names=FALSE)


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
  join_last4(SL)

merged_final[is.na(merged_final)] <- -999

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Merged_PRO_TTFD_", toupper(current_date), ".csv")

write.csv(merged_final, new_filename, quote = FALSE, row.names = FALSE)
