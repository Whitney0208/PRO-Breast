# Heatmap: Symptoms(censor variable) × Relapse Type
rm(list=ls())
library(data.table)
library(plyr)
library(dplyr)
library(haven)
library(tidyr)
library(readr)
library(ggplot2)
library(GGally)
library(ggsci)
library(reshape2)
library(survival)
library(dplyr)
library(stringr)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_SanofiU_1997_120/wwanbing/output")
BCR_1997 <- read_csv("Breast_SanofiU_1997_120_BCR_10APR2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_SanofiU_2000_118/wwanbing/output")
BCR_2000 <- read_csv("Breast_SanofiU_2000_118_BCR_21APR2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
Merged_TTFD <- read_csv("Merged_PRO_TTFD_28JUN2025.csv")

#Merge TTFD
# Select all variables starting with TTFD10, TTFD10P, censored, censored10P
vars_to_merge <- colnames(Merged_TTFD)[grepl("^(TTFD10|TTFD10P|censored|censored10P)_", colnames(Merged_TTFD))]

vars_to_merge <- c("UID", vars_to_merge)

symptoms_TTFD <- Merged_TTFD[, vars_to_merge]

# Merge into BCR_1997 and BCR_2000
df_1997 <- left_join(BCR_1997, symptoms_TTFD, by = "UID")
df_2000 <- left_join(BCR_2000, symptoms_TTFD, by = "UID")

#suffix
all_vars <- colnames(df_1997)
symptom_vars <- all_vars[grepl("^TTFD10_", all_vars)]
symptom_suffixes <- gsub("^TTFD10_", "", symptom_vars)
unique_symptoms <- unique(symptom_suffixes)

#Q2 Heatmap: Symptoms(censor variable) × Relapse Type. 
#   Color: % of patients with symptom deterioration, Rows: Symptoms, Columns: Relapse types 
#   (all, distant/local/regional) 

df_1997$Year <- "1997"
df_2000$Year <- "2000"
df <- bind_rows(df_1997, df_2000)

all_vars <- colnames(df_1997)
symptom_vars <- all_vars[grepl("^censored_", all_vars)]
symptom_suffixes <- gsub("^censored_", "", symptom_vars)


###new heatmap
# Step 1: Construct long-format data
df_all<- df

df_long <- df_all %>%
  select(UID, all_of(symptom_vars), DISTANT, REGIONAL, LOCAL) %>%
  pivot_longer(cols = all_of(symptom_vars), names_to = "Symptom", values_to = "Deteriorated") %>%
  mutate(
    Symptom = gsub("censored_", "", Symptom),
    RelapseType = case_when(
      DISTANT == 1 & REGIONAL == 1 & LOCAL == 1 ~ "Multiple",
      DISTANT == 1 & REGIONAL == 1 ~ "Multiple",
      DISTANT == 1 & LOCAL == 1 ~ "Multiple",
      REGIONAL == 1 & LOCAL == 1 ~ "Multiple",
      DISTANT == 1 ~ "Distant",
      REGIONAL == 1 ~ "Regional",
      LOCAL == 1 ~ "Local",
      TRUE ~ "NoRelapse"
    )
  )

# Add the tag "AllRelapse" (any relapse counts)
df_long <- df_long %>%
  mutate(AllRelapse = ifelse(RelapseType %in% c("Distant", "Regional", "Local", "Multiple"), "AllRelapse", "NoRelapse"))


# Step 2: Keep only Distant / Regional / Local for the heatmap columns (exclude Multiple)
df_long_filtered <- df_long %>%
  filter(RelapseType %in% c("Distant", "Regional", "Local"))

# Step 3: Compute the deterioration percentage for each RelapseType
relapse_specific <- df_long_filtered %>%
  group_by(Symptom, RelapseType) %>%
  summarise(
    total = n(),
    deteriorated = sum(Deteriorated == 1, na.rm = TRUE),
    percent = round(100 * deteriorated / total, 1),
    .groups = "drop"
  )

# Step 4: Compute the AllRelapse percentage separately
all_relapse <- df_long %>%
  filter(AllRelapse == "AllRelapse") %>%
  group_by(Symptom) %>%
  summarise(
    RelapseType = "AllRelapse",
    total = n(),
    deteriorated = sum(Deteriorated == 1, na.rm = TRUE),
    percent = round(100 * deteriorated / total, 1),
    .groups = "drop"
  )

# Step 5: Combine data
heatmap_data <- bind_rows(relapse_specific, all_relapse)

#suffix -> full name, rank by order
symptom_full_names <- c(
  "FA"  = "Fatigue",
  "NV"  = "Nausea/Vomit",
  "PA"  = "Pain",
  "DY"  = "Dyspnea",
  "SL"  = "Insomnia",
  "AP"  = "Appetite Loss",
  "CO"  = "Constipation",
  "DI"  = "Diarrhea",
  "FI"  = "Financial",
  "PF2" = "Physical Functioning",
  "RF2" = "Role Functioning",
  "EF"  = "Emotional Functioning",
  "CF"  = "Cognitive Functioning",
  "SF"  = "Social Functioning",
  "GHS" = "Global Health Status / QOL"
)

heatmap_data$Symptom <- symptom_full_names[as.character(heatmap_data$Symptom)]

# Sort the symptom (in descending order of AllRelapse)
symptom_order <- heatmap_data %>%
  filter(RelapseType == "AllRelapse") %>%
  arrange(percent) %>%
  pull(Symptom)

heatmap_data$Symptom <- factor(heatmap_data$Symptom, levels = symptom_order)

ggplot(heatmap_data, aes(x = RelapseType, y = Symptom, fill = percent)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#4575b4",    # Blue (low values)
    mid = "white",      # Mid value
    high = "#d73027",   # Orange-red (high values)
    midpoint = 50,
    limit = c(0, 100),
    name = "% Deteriorated"
  ) +
  geom_text(aes(label = paste0(percent, "%")), size = 3.5) +
  theme_minimal() +
  labs(
    title = "Symptom Deterioration by Relapse Type",
    x = "Relapse Type",
    y = "Symptom"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

heatmap_data$RelapseType <- recode(heatmap_data$RelapseType,
                                   "AllRelapse" = "All sites")

ggplot(heatmap_data, aes(x = RelapseType, y = Symptom, fill = percent)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#4575b4",    # Blue (low values)
    mid = "white",      # Mid value
    high = "#d73027",   # Orange-red (high values)
    midpoint = 50,
    limit = c(0, 100),
    name = "% Deteriorated"
  ) +
  geom_text(aes(label = paste0(percent, "%")), size = 3.5) +
  theme_minimal() +
  labs(
    title = "Symptom Deterioration by Relapse Type",
    x = "Relapse Type",
    y = "Symptom"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold")
  )