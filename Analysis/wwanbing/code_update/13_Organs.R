rm(list=ls())
library(data.table)
library(stringr)
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
library(survminer)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_Celgene_2001_107/wwanbing/output")
tumor_2001 <- read_csv("Breast_Celgene_2001_107_tumor.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_SanofiU_2004_135/wwanbing/output")
tumor_2004 <- read_csv("Breast_SanofiU_2004_135_tumor_18MAR2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
Merged_TTFD <- read_csv("Merged_PRO_TTFD_28JUN2025.csv")


#Merge 2001 and 2004
tumor_2001$Study <- "2001"
tumor_2004$Study <- "2004"

all(colnames(tumor_2001) == colnames(tumor_2004))

tumor_2001 <- tumor_2001 %>%
  rename(UID = LUID)

tumor_2001 <- tumor_2001 %>%
  mutate(
    UID = as.character(as.numeric(UID)),
    ID = as.character(as.numeric(ID))
  )

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = as.numeric(SIZE))

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = ifelse(is.na(SIZE), -999, SIZE))


#Merge TTFD
# Select all variables starting with TTFD10, TTFD10P, censored, censored10P
vars_to_merge <- colnames(Merged_TTFD)[grepl("^(TTFD10|TTFD10P|censored|censored10P)_", colnames(Merged_TTFD))]

vars_to_merge <- c("UID", vars_to_merge)

symptoms_TTFD <- Merged_TTFD[, vars_to_merge]

# Merge into tumor_ALL
tumor_2001 <- left_join(tumor_2001, symptoms_TTFD, by = "UID")
tumor_2004 <- left_join(tumor_2004, symptoms_TTFD, by = "UID")

tumor_all <- bind_rows(tumor_2001, tumor_2004)

###Do target or non-target lesions more strongly correlate with symptom deterioration? 
tumor_all$LESTYPE <- ifelse(tumor_all$LESTYPE == 1, 1, 0)
table(tumor_all$LESTYPE)

# Step 1: Organ classification (new variable Organ_clean)
tumor_all$Organ_category <- ifelse(tumor_all$ORGAN %in% c("Liver"), "Liver",
                                   ifelse(tumor_all$ORGAN %in% c("Lung"), "Lung",
                                          ifelse(tumor_all$ORGAN %in% c("Lymph node", "Lymph nodes"), "Lymph", "Other")))


# Step 2: Extract whether there was ever organ metastasis (by ID)
organ_flags <- tumor_all %>%
  filter(!is.na(ORGAN)) %>%
  mutate(Organ_clean = case_when(
    ORGAN %in% "Liver" ~ "Liver",
    ORGAN %in% "Lung"  ~ "Lung",
    ORGAN %in% c("Lymph node","Lymph nodes") ~ "Lymph",
    TRUE ~ "Other"            # Breast is also merged into Other
  )) %>%
  distinct(ID, Organ_clean) %>%
  mutate(flag = 1) %>%
  pivot_wider(names_from = Organ_clean, values_from = flag, values_fill = 0)


# Step 3: Merge organ flags back to TTFD data
ttfd_cols <- grep("^TTFD10_", names(tumor_all), value = TRUE)
censor_cols <- grep("^censored_", names(tumor_all), value = TRUE)

# Use AP as an example
symptom <- "AP"
ttfd_col <- paste0("TTFD10_", symptom)
censored_col <- paste0("censored_", symptom)

# Prepare analysis dataset
df <- tumor_all %>%
  select(ID, UID, all_of(ttfd_col), all_of(censored_col)) %>%
  left_join(organ_flags, by = "ID")

# Data cleaning
df <- df %>%
  filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
         !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)

# Fallback: ensure all four columns exist
for (org in c("Liver", "Lung", "Lymph", "Other")) {
  if (!org %in% names(df)) df[[org]] <- 0L
  df[[org]] <- ifelse(is.na(df[[org]]), 0L, as.integer(df[[org]]))
}

# Step 4: Build Cox model
surv_obj <- Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
fit <- coxph(surv_obj ~ Liver + Lung + Lymph + Other, data = df)
sm <- summary(fit)

# Step 5: Extract HR and related results
forest_df <- data.frame(
  Term = rownames(sm$coefficients),
  HR = sm$coefficients[, "exp(coef)"],
  LCL = sm$conf.int[, "lower .95"],
  UCL = sm$conf.int[, "upper .95"],
  p_value = sm$coefficients[, "Pr(>|z|)"]
)

# Add significance labels
forest_df <- forest_df %>%
  mutate(
    Sig = ifelse(p_value < 0.05, "Significant", "Not Significant"),
    Term = gsub("Other", "Other organ", Term),  # Beautify display
    Label = Term
  ) %>%
  arrange(desc(HR)) %>%
  mutate(Label = factor(Label, levels = unique(Label)))

forest_df <- forest_df %>%
  mutate(
    OR_CI = sprintf("%.2f (%.2f, %.2f)", HR, LCL, UCL),
    Label = factor(Label, levels = Label[order(HR)])  # Sort
  )

# Forest plot with OR + CI annotations
label_x_pos <- 1.8

ggplot(forest_df, aes(x = HR, y = Label, color = Sig)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2) +
  
  # ✅ Fix text to display at a uniform position
  geom_text(aes(x = label_x_pos, label = OR_CI), hjust = 0, size = 3.5, color = "black") +
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Significant" = "#d73027", "Not Significant" = "gray40")) +
  scale_x_log10(limits = c(0.6, 2.5)) +  # Slightly extend right side to ensure label_x_pos is visible
  labs(
    title = paste0("Forest Plot: Organ Involvement and ", symptom, " Deterioration"),
    x = "Hazard Ratio (HR, log scale)",
    y = NULL,
    color = "p < 0.05"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    plot.margin = margin(t = 10, r = 80, b = 10, l = 10)
  )


###Function for 15 symptoms
plot_forest_for_symptom <- function(symptom, tumor_all, organ_flags) {
  ttfd_col <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  df <- tumor_all %>%
    select(ID, UID, all_of(ttfd_col), all_of(censored_col)) %>%
    left_join(organ_flags, by = "ID") %>%
    filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
           !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  for (org in c("Liver", "Lung", "Lymph", "Other")) {
    if (!org %in% names(df)) df[[org]] <- 0L
    df[[org]] <- ifelse(is.na(df[[org]]), 0L, as.integer(df[[org]]))
  }
  
  # Cox model
  surv_obj <- Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
  fit <- coxph(surv_obj ~ Liver + Lung + Lymph + Other, data = df)
  sm <- summary(fit)
  
  # Construct forest_df
  forest_df <- data.frame(
    Term = rownames(sm$coefficients),
    HR = sm$coefficients[, "exp(coef)"],
    LCL = sm$conf.int[, "lower .95"],
    UCL = sm$conf.int[, "upper .95"],
    p_value = sm$coefficients[, "Pr(>|z|)"]
  ) %>%
    mutate(
      Sig = ifelse(p_value < 0.05, "Significant", "Not Significant"),
      Term = gsub("Other", "Other organ", Term),
      Label = Term,
      OR_CI = sprintf("%.2f (%.2f, %.2f)", HR, LCL, UCL),
      Label = factor(Label, levels = Label[order(HR)])
    )
  
  label_x_pos <- 1.8
  
  # Plot
  ggplot(forest_df, aes(x = HR, y = Label, color = Sig)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2) +
    geom_text(aes(x = label_x_pos, label = OR_CI), hjust = 0, size = 3.5, color = "black") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c("Significant" = "#d73027", "Not Significant" = "gray40")) +
    scale_x_continuous(limits = c(0.6, 2.5)) +
    labs(
      title = paste0("Forest Plot: Organ Involvement and ", symptom, " Deterioration"),
      x = "Hazard Ratio (HR)",
      y = NULL,
      color = "p < 0.05"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 10),
      plot.margin = margin(t = 10, r = 80, b = 10, l = 10)
    )
}


symptoms <- c("AP", "CO", "DI", "DY", "FA", "FI", "GHS", "NV", "PA", "SL", "PF2", "RF2", "EF", "CF", "SF")

for (symptom in symptoms) {
  print(paste("Plotting:", symptom))
  print(plot_forest_for_symptom(symptom, tumor_all, organ_flags))
}

### Distribution plot of multiple organs 
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_Celgene_2001_107/wwanbing/output")
tumor_2001 <- read_csv("Breast_Celgene_2001_107_tumor.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_SanofiU_2004_135/wwanbing/output")
tumor_2004 <- read_csv("Breast_SanofiU_2004_135_tumor_18MAR2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
Merged_TTFD <- read_csv("Merged_PRO_TTFD_28JUN2025.csv")


#Merge 2001 and 2004
tumor_2001$Study <- "2001"
tumor_2004$Study <- "2004"

all(colnames(tumor_2001) == colnames(tumor_2004))

tumor_2001 <- tumor_2001 %>%
  rename(UID = LUID)

tumor_2001 <- tumor_2001 %>%
  mutate(
    UID = as.character(as.numeric(UID)),
    ID = as.character(as.numeric(ID))
  )

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = as.numeric(SIZE))

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = ifelse(is.na(SIZE), -999, SIZE))


#Merge TTFD
# Select all variables starting with TTFD10, TTFD10P, censored, censored10P
vars_to_merge <- colnames(Merged_TTFD)[grepl("^(TTFD10|TTFD10P|censored|censored10P)_", colnames(Merged_TTFD))]

vars_to_merge <- c("UID", vars_to_merge)

symptoms_TTFD <- Merged_TTFD[, vars_to_merge]

# Merge into tumor_ALL
tumor_2001 <- left_join(tumor_2001, symptoms_TTFD, by = "UID")
tumor_2004 <- left_join(tumor_2004, symptoms_TTFD, by = "UID")

tumor_all <- bind_rows(tumor_2001, tumor_2004)

###Do target or non-target lesions more strongly correlate with symptom deterioration? 
tumor_all$LESTYPE <- ifelse(tumor_all$LESTYPE == 1, 1, 0)
table(tumor_all$LESTYPE)

tumor_all$Organ_category <- ifelse(tumor_all$ORGAN %in% c("Liver"), "Liver",
                                   ifelse(tumor_all$ORGAN %in% c("Lung"), "Lung",
                                          ifelse(tumor_all$ORGAN %in% c("Lymph node", "Lymph nodes"), "Lymph", "Other")))


# Step 2: Extract whether there was ever organ metastasis (by ID)
organ_flags <- tumor_all %>%
  filter(!is.na(ORGAN)) %>%
  mutate(Organ_clean = case_when(
    ORGAN %in% "Liver" ~ "Liver",
    ORGAN %in% "Lung"  ~ "Lung",
    ORGAN %in% c("Lymph node","Lymph nodes") ~ "Lymph",
    TRUE ~ "Other"            # Breast is also merged into Other
  )) %>%
  distinct(ID, Organ_clean) %>%
  mutate(flag = 1) %>%
  pivot_wider(names_from = Organ_clean, values_from = flag, values_fill = 0)

#plot 1: A bar chart showing the number of patients with "several organ metastases", such as 0, 1, 2, 3, and 4 organs.
# Step 1: Count the number of organs involved per patient
organ_flags$Organ_count <- rowSums(organ_flags[, c("Liver", "Lung", "Lymph", "Other")])

# Step 2: Plot
ggplot(organ_flags, aes(x = as.factor(Organ_count))) +
  geom_bar(fill = "#4DBBD5FF") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 4) +  # Add count labels
  labs(
    title = "Distribution of Number of Organs Involved",
    x = "Number of Organs Involved",
    y = "Number of Patients"
  ) +
  theme_minimal(base_size = 13)

#plot 2: By organ combination (Venn-like barplot)
# Step 1: Create organ combination labels
organ_flags$combo <- apply(organ_flags[, c("Liver", "Lung", "Lymph", "Other")], 1, function(row) {
  organs <- c("Liver", "Lung", "Lymph", "Other")[which(row == 1)]
  if (length(organs) == 0) return("None")
  paste(organs, collapse = "+")
})

# Step 2: Plot
ggplot(organ_flags, aes(x = combo)) +
  geom_bar(fill = "#E64B35FF") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 3.8) +  # Add count labels
  labs(
    title = "Organ Involvement Combination",
    x = "Organ Combination",
    y = "Number of Patients"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###each organ -log rank test, KM curve 
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_Celgene_2001_107/wwanbing/output")
tumor_2001 <- read_csv("Breast_Celgene_2001_107_tumor.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_SanofiU_2004_135/wwanbing/output")
tumor_2004 <- read_csv("Breast_SanofiU_2004_135_tumor_18MAR2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
Merged_TTFD <- read_csv("Merged_PRO_TTFD_28JUN2025.csv")


#Merge 2001 and 2004
tumor_2001$Study <- "2001"
tumor_2004$Study <- "2004"

all(colnames(tumor_2001) == colnames(tumor_2004))

tumor_2001 <- tumor_2001 %>%
  rename(UID = LUID)

tumor_2001 <- tumor_2001 %>%
  mutate(
    UID = as.character(as.numeric(UID)),
    ID = as.character(as.numeric(ID))
  )

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = as.numeric(SIZE))

tumor_2004 <- tumor_2004 %>%
  mutate(SIZE = ifelse(is.na(SIZE), -999, SIZE))


#Merge TTFD
# Select all variables starting with TTFD10, TTFD10P, censored, censored10P
vars_to_merge <- colnames(Merged_TTFD)[grepl("^(TTFD10|TTFD10P|censored|censored10P)_", colnames(Merged_TTFD))]

vars_to_merge <- c("UID", vars_to_merge)

symptoms_TTFD <- Merged_TTFD[, vars_to_merge]

# Merge into tumor_ALL
tumor_2001 <- left_join(tumor_2001, symptoms_TTFD, by = "UID")
tumor_2004 <- left_join(tumor_2004, symptoms_TTFD, by = "UID")

tumor_all <- bind_rows(tumor_2001, tumor_2004)

###Do target or non-target lesions more strongly correlate with symptom deterioration? 
tumor_all$LESTYPE <- ifelse(tumor_all$LESTYPE == 1, 1, 0)
table(tumor_all$LESTYPE)

# Organ classification (new variable Organ_clean)
tumor_all$Organ_category <- ifelse(tumor_all$ORGAN %in% c("Liver"), "Liver",
                                   ifelse(tumor_all$ORGAN %in% c("Lung"), "Lung",
                                          ifelse(tumor_all$ORGAN %in% c("Lymph node", "Lymph nodes"), "Lymph", "Other")))


# Extract whether there was ever organ metastasis (by ID)
organ_flags <- tumor_all %>%
  filter(!is.na(ORGAN)) %>%
  mutate(Organ_clean = case_when(
    ORGAN %in% "Liver" ~ "Liver",
    ORGAN %in% "Lung"  ~ "Lung",
    ORGAN %in% c("Lymph node","Lymph nodes") ~ "Lymph",
    TRUE ~ "Other"            # Breast is also merged into Other
  )) %>%
  distinct(ID, Organ_clean) %>%
  mutate(flag = 1) %>%
  pivot_wider(names_from = Organ_clean, values_from = flag, values_fill = 0)


# Merge organ flags back to TTFD data
ttfd_cols <- grep("^TTFD10_", names(tumor_all), value = TRUE)
censor_cols <- grep("^censored_", names(tumor_all), value = TRUE)

# Set symptom name
symptom <- "AP"
ttfd_col <- paste0("TTFD10_", symptom)
censored_col <- paste0("censored_", symptom)

# Prepare data
df <- tumor_all %>%
  select(ID, UID, all_of(ttfd_col), all_of(censored_col)) %>%
  left_join(organ_flags, by = "ID") %>%
  filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
         !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)

# Each organ
organ_list <- c("Liver", "Lung", "Lymph", "Other")

for (org in organ_list) {
  # Fallback
  if (!org %in% names(df)) df[[org]] <- 0
  df[[org]] <- ifelse(is.na(df[[org]]), 0L, as.integer(df[[org]]))
  
  surv_obj <- Surv(df[[ttfd_col]], 1 - df[[censored_col]])
  
  # log-rank test
  surv_diff <- survdiff(surv_obj ~ df[[org]])
  pval <- 1 - pchisq(surv_diff$chisq, df = 1)
  
  cat(paste0("\n===== ", org, " vs ", symptom, " =====\n"))
  cat("Log-rank test p-value: ", round(pval, 4), "\n")
  
  # KM curve plot
  fit <- survfit(surv_obj ~ df[[org]])
  p <- ggsurvplot(
    fit, data = df,
    legend.title = paste0(org, " (0 = No, 1 = Yes)"),
    legend.labs = c("No", "Yes"),
    xlab = paste0("Time to ", symptom, " deterioration"),
    ylab = "Survival Probability",
    title = paste0("KM Curve: ", symptom, " by ", org),
    pval = TRUE, pval.coord = c(1, 0.1),
    ggtheme = theme_minimal()
  )
  print(p)
}

### 15 symptoms
symptoms <- c("AP", "CO", "DI", "DY", "FA", "FI", "GHS", "NV", "PA", "SL", "PF2", "RF2", "EF", "CF", "SF")
organs <- c("Liver", "Lung", "Lymph", "Other")

# Define function: for each symptom, output a 2×2 KM panel plot
plot_km_panel_for_symptom <- function(symptom, tumor_all, organ_flags) {
  ttfd_col <- paste0("TTFD10_", symptom)
  censored_col <- paste0("censored_", symptom)
  
  df <- tumor_all %>%
    select(ID, UID, all_of(ttfd_col), all_of(censored_col)) %>%
    left_join(organ_flags, by = "ID") %>%
    filter(!is.na(.data[[ttfd_col]]), .data[[ttfd_col]] != -999,
           !is.na(.data[[censored_col]]), .data[[censored_col]] != -999)
  
  # Fallback: ensure each organ column exists
  for (org in organs) {
    if (!org %in% names(df)) df[[org]] <- 0L
    df[[org]] <- ifelse(is.na(df[[org]]), 0L, as.integer(df[[org]]))
  }
  
  # Plot for each organ
  plot_list <- list()
  for (org in organs) {
    surv_obj <- Surv(time = df[[ttfd_col]], event = 1 - df[[censored_col]])
    surv_diff <- survdiff(surv_obj ~ df[[org]])
    pval <- 1 - pchisq(surv_diff$chisq, df = 1)
    
    fit <- survfit(surv_obj ~ df[[org]])
    
    g <- ggsurvplot(
      fit, data = df,
      legend.title = paste0(org, " involvement"),
      legend.labs = c("No", "Yes"),
      palette = c("#4DBBD5FF", "#E64B35FF"),
      xlab = paste0("Time to ", symptom, " deterioration"),
      ylab = "Survival Probability",
      title = paste0(org, " (p = ", format.pval(pval, digits = 3), ")"),
      ggtheme = theme_minimal(base_size = 13),
      pval = FALSE
    )
    
    plot_list[[org]] <- g$plot
  }
  
  # Combine into a 2×2 panel
  panel <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2)
  panel <- annotate_figure(panel, top = text_grob(paste0("KM Curves for ", symptom, " by Organ Involvement"), face = "bold", size = 14))
  
  print(panel)
}

for (symptom in symptoms) {
  print(paste("Plotting 2x2 KM panel for:", symptom))
  print(plot_km_panel_for_symptom(symptom, tumor_all, organ_flags))
}