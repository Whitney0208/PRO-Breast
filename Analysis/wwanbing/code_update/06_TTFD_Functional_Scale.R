library(data.table)
library(dplyr)
library(purrr)

###PF2
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")
PF2<-fread("Breast_PRO_Physical_functioning(revised)_22OCT2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- PF2 %>%
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
new_filename <- paste0("TTFD_PF2_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
PF2_2 <- PF2 %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_PF2_", toupper(current_date), ".csv")

write.csv(PF2_2, new_filename, quote=FALSE,row.names=FALSE)


###RF2
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")
RF2<-fread("Breast_PRO_Role_functioning(revised)_22OCT2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- RF2 %>%
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
new_filename <- paste0("TTFD_RF2_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
RF2_2 <- RF2 %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_RF2_", toupper(current_date), ".csv")

write.csv(RF2_2, new_filename, quote=FALSE,row.names=FALSE)


###EF
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")
EF<-fread("Breast_PRO_Emotional_functioning_22OCT2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- EF %>%
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
new_filename <- paste0("TTFD_EF_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
EF_2 <- EF %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_EF_", toupper(current_date), ".csv")

write.csv(EF_2, new_filename, quote=FALSE,row.names=FALSE)


###CF
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")
CF<-fread("Breast_PRO_Cognitive_functioning_22OCT2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- CF %>%
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
new_filename <- paste0("TTFD_CF_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
CF_2 <- CF %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_CF_", toupper(current_date), ".csv")

write.csv(CF_2, new_filename, quote=FALSE,row.names=FALSE)

###SF
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/Functional_Scales")
SF<-fread("Breast_PRO_Social_functioning_22OCT2025.csv")

# Calculate TTFD10 and TTFD10P for each subject (UID)
ttfd_result <- SF %>%
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
new_filename <- paste0("TTFD_SF_", toupper(current_date), ".csv")

write.csv(ttfd_result, new_filename, quote=FALSE,row.names=FALSE)

###export TTFD+Covariate
SF_2 <- SF %>%
  left_join(
    ttfd_result %>% 
      dplyr::select(UID, TTFD10, censored, TTFD10P, censored10P),
    by = "UID"
  )

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("TTFD+Covariates_SF_", toupper(current_date), ".csv")

write.csv(SF_2, new_filename, quote=FALSE,row.names=FALSE)

















