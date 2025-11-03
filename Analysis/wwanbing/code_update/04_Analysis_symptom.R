#delete DV, distinct dataset

library(data.table)
library(dplyr)

FA <- fread("Breast_PRO_Fatigue_19MAY2025.csv")
NV <- fread("Breast_PRO_Nausea_and_vomiting_19MAY2025.csv")
PA <- fread("Breast_PRO_Pain_19MAY2025.csv")
DY <- fread("Breast_PRO_Dyspnoea_19MAY2025.csv")
SL <- fread("Breast_PRO_Insomnia_19MAY2025.csv")
AP <- fread("Breast_PRO_Appetite_loss_19MAY2025.csv")  
CO <- fread("Breast_PRO_Constipation_19MAY2025.csv")
DI <- fread("Breast_PRO_Diarrhoea_19MAY2025.csv")
FI <- fread("Breast_PRO_Financial_difficulties_19MAY2025.csv")
GHS <- fread("Breast_PRO_GHS_19MAY2025.csv")

FA2 <- FA[, !c("QSFLAG", "DV")] %>% distinct()
NV2 <- NV[, !c("QSFLAG", "DV")] %>% distinct()
PA2 <- PA[, !c("QSFLAG", "DV")] %>% distinct()
DY2 <- DY[, !c("QSFLAG", "DV")] %>% distinct()
SL2 <- SL[, !c("QSFLAG", "DV")] %>% distinct()
AP2 <- AP[, !c("QSFLAG", "DV")] %>% distinct()
CO2 <- CO[, !c("QSFLAG", "DV")] %>% distinct()
DI2 <- DI[, !c("QSFLAG", "DV")] %>% distinct()
FI2 <- FI[, !c("QSFLAG", "DV")] %>% distinct()
GHS2 <- GHS[, !c("QSFLAG", "DV")] %>% distinct()

current_date <- format(Sys.Date(), "%d%b%Y")

new_filename <- paste0("Breast_PRO_Fatigue_updated_", toupper(current_date), ".csv")
write.csv(FA2, new_filename, quote=FALSE,row.names=FALSE)

write.csv(NV2, paste0("Breast_PRO_Nausea_and_vomiting_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(PA2, paste0("Breast_PRO_Pain_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(DY2, paste0("Breast_PRO_Dyspnoea_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(SL2, paste0("Breast_PRO_Insomnia_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(AP2, paste0("Breast_PRO_Appetite_loss_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(CO2, paste0("Breast_PRO_Constipation_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(DI2, paste0("Breast_PRO_Diarrhoea_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(FI2, paste0("Breast_PRO_Financial_difficulties_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)
write.csv(GHS2, paste0("Breast_PRO_GHS_updated_", toupper(current_date), ".csv"), quote = FALSE, row.names = FALSE)

###update PFS
FA <- Fatigue2
FA2 <- FA %>%
  select(-QSFLAG, -DV) %>%
  distinct()

current_date <- format(Sys.Date(), "%d%b%Y")

new_filename <- paste0("Breast_PRO_Fatigue_updated_", toupper(current_date), ".csv")
write.csv(FA2, new_filename, quote=FALSE,row.names=FALSE)
