### Updated version for 19_tumor_size_analysis.R

rm(list=ls())
library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(haven)
library(tidyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggsci)
library(reshape2)
library(survival)
library(survminer)
library(RColorBrewer)

###Do tumor size correlate with symptoms deterioration? 
setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/data")

tumor_all<-read_csv("PRO_TTFD_tumor05AUG2025.csv")

Monolix_2001 <- read_csv("Breast_Monolix_2001_tumor.csv")
Monolix_2004 <- read_csv("Breast_Monolix_2004_tumor.csv")

Monolix_2001$Year <- "2001"
Monolix_2004$Year <- "2004"

Monolix_2001$UID <- as.character(Monolix_2001$UID)
Monolix_2004$UID <- as.character(Monolix_2004$UID)

Monolix <- bind_rows(Monolix_2001, Monolix_2004)

PRO_2001 <- read_csv("PRO_Size_2001.csv")
PRO_2004 <- read_csv("PRO_Size_2004.csv")

PRO_2001$Year <- "2001"
PRO_2004$Year <- "2004"

PRO_2001$UID <- as.character(PRO_2001$UID)
PRO_2004$UID <- as.character(PRO_2004$UID)

PRO <- bind_rows(PRO_2001, PRO_2004)

# left_join censored variable

#> length(unique(PRO$UID))
#[1] 356
#> length(unique(tumor_all$UID))
#[1] 436

censor_cols <- grep("^censored", colnames(tumor_all), value = TRUE)
censor_df <- tumor_all %>%
  select(UID, all_of(censor_cols)) %>%
  distinct()  

PRO <- PRO %>%
  left_join(censor_df, by = "UID")

Monolix$UID <- as.character(Monolix$UID)

vars_to_merge <- c("UID", "AGE", "MENOS", "CHEMO", "HORMON", "ERS", "PGRS", "ECOG")

monolix_clinical <- Monolix %>%
  select(all_of(vars_to_merge)) %>%
  distinct()

PRO <- PRO %>%
  left_join(monolix_clinical, by = "UID")

# colnames(PRO)
# [1] "UID"                  "ID"                   "TID"                  "ORGAN2"              
# [5] "TTFD10_AP"            "TTFD10P_AP"           "TTFD10_CO"            "TTFD10P_CO"          
# [9] "TTFD10_DI"            "TTFD10P_DI"           "TTFD10_DY"            "TTFD10P_DY"          
# [13] "TTFD10_FA"            "TTFD10P_FA"           "TTFD10_FI"            "TTFD10P_FI"          
# [17] "TTFD10_GHS"           "TTFD10P_GHS"          "TTFD10_NV"            "TTFD10P_NV"          
# [21] "TTFD10_PA"            "TTFD10P_PA"           "TTFD10_SL"            "TTFD10P_SL"          
# [25] "TTFD10_PF2"           "TTFD10P_PF2"          "TTFD10_RF2"           "TTFD10P_RF2"         
# [29] "TTFD10_EF"            "TTFD10P_EF"           "TTFD10_CF"            "TTFD10P_CF"          
# [33] "TTFD10_SF"            "TTFD10P_SF"           "TS_TTFD10_AP"         "TS_TTFD10P_AP"       
# [37] "TS_TTFD10_CO"         "TS_TTFD10P_CO"        "TS_TTFD10_DI"         "TS_TTFD10P_DI"       
# [41] "TS_TTFD10_DY"         "TS_TTFD10P_DY"        "TS_TTFD10_FA"         "TS_TTFD10P_FA"       
# [45] "TS_TTFD10_FI"         "TS_TTFD10P_FI"        "TS_TTFD10_GHS"        "TS_TTFD10P_GHS"      
# [49] "TS_TTFD10_NV"         "TS_TTFD10P_NV"        "TS_TTFD10_PA"         "TS_TTFD10P_PA"       
# [53] "TS_TTFD10_SL"         "TS_TTFD10P_SL"        "TS_TTFD10_PF2"        "TS_TTFD10P_PF2"      
# [57] "TS_TTFD10_RF2"        "TS_TTFD10P_RF2"       "TS_TTFD10_EF"         "TS_TTFD10P_EF"       
# [61] "TS_TTFD10_CF"         "TS_TTFD10P_CF"        "TS_TTFD10_SF"         "TS_TTFD10P_SF"       
# [65] "base_Tumorsize"       "Ratio_TS_TTFD10_AP"   "Ratio_TS_TTFD10P_AP"  "Ratio_TS_TTFD10_CO"  
# [69] "Ratio_TS_TTFD10P_CO"  "Ratio_TS_TTFD10_DI"   "Ratio_TS_TTFD10P_DI"  "Ratio_TS_TTFD10_DY"  
# [73] "Ratio_TS_TTFD10P_DY"  "Ratio_TS_TTFD10_FA"   "Ratio_TS_TTFD10P_FA"  "Ratio_TS_TTFD10_FI"  
# [77] "Ratio_TS_TTFD10P_FI"  "Ratio_TS_TTFD10_GHS"  "Ratio_TS_TTFD10P_GHS" "Ratio_TS_TTFD10_NV"  
# [81] "Ratio_TS_TTFD10P_NV"  "Ratio_TS_TTFD10_PA"   "Ratio_TS_TTFD10P_PA"  "Ratio_TS_TTFD10_SL"  
# [85] "Ratio_TS_TTFD10P_SL"  "Ratio_TS_TTFD10_PF2"  "Ratio_TS_TTFD10P_PF2" "Ratio_TS_TTFD10_RF2" 
# [89] "Ratio_TS_TTFD10P_RF2" "Ratio_TS_TTFD10_EF"   "Ratio_TS_TTFD10P_EF"  "Ratio_TS_TTFD10_CF"  
# [93] "Ratio_TS_TTFD10P_CF"  "Ratio_TS_TTFD10_SF"   "Ratio_TS_TTFD10P_SF"  "Year"                
# [97] "censored_AP"          "censored10P_AP"       "censored_CO"          "censored10P_CO"      
# [101] "censored_DI"          "censored10P_DI"       "censored_DY"          "censored10P_DY"      
# [105] "censored_FA"          "censored10P_FA"       "censored_FI"          "censored10P_FI"      
# [109] "censored_GHS"         "censored10P_GHS"      "censored_NV"          "censored10P_NV"      
# [113] "censored_PA"          "censored10P_PA"       "censored_SL"          "censored10P_SL"      
# [117] "censored_PF2"         "censored10P_PF2"      "censored_RF2"         "censored10P_RF2"     
# [121] "censored_EF"          "censored10P_EF"       "censored_CF"          "censored10P_CF"      
# [125] "censored_SF"          "censored10P_SF"       "AGE"                  "MENOS"               
# [129] "CHEMO"                "HORMON"               "ERS"                  "PGRS"                
# [133] "ECOG"    


### correlation 
symptoms <- c("AP", "CO", "DI", "DY", "FA", "FI", "GHS", "NV", 
              "PA", "SL", "PF2", "RF2", "EF", "CF", "SF")

symptom_names <- c(
  "FA"  = "Fatigue",
  "NV"  = "Nausea/Vomit",
  "PA"  = "Pain",
  "DY"  = "Dyspnea",
  "SL"  = "Insomnia",
  "AP"  = "Appetite loss",
  "CO"  = "Constipation",
  "DI"  = "Diarrhea",
  "FI"  = "Financial",
  "PF2" = "Physical",
  "RF2" = "Role",
  "EF"  = "Emotional",
  "CF"  = "Cognitive",
  "SF"  = "Social",
  "GHS" = "GHS/QOL"
)

cor_results <- data.frame(
  Symptom = character(),
  Correlation = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (symptom in symptoms) {
  ttfd_col <- paste0("TTFD10_", symptom)
  tumor_col <- paste0("TS_TTFD10_", symptom)
  cens_col <- paste0("censored_", symptom)
  
  df_sub <- PRO %>%
    filter(!is.na(.data[[ttfd_col]]), 
           !is.na(.data[[tumor_col]]),
           .data[[cens_col]] == 0)
  
  if (nrow(df_sub) >= 5) {  # require at least 5 observations to test correlation
    res <- cor.test(df_sub[[ttfd_col]], df_sub[[tumor_col]], method = "spearman")
    
    cor_results <- rbind(cor_results, data.frame(
      Symptom = symptom,
      Correlation = res$estimate,
      P_value = res$p.value
    ))
  }
}

# compute -log10(P)
cor_results$log10_P <- -log10(cor_results$P_value)

# replace abbreviations with full names
cor_results$Symptom <- symptom_names[cor_results$Symptom]

# order by correlation coefficient
cor_results$Symptom <- factor(cor_results$Symptom,
                              levels = cor_results$Symptom[order(cor_results$Correlation)])

# heatmap
ggplot(cor_results, aes(x = "", y = Symptom, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n(p=%.3f)", Correlation, P_value)), size = 3.2) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, name = "Spearman\nCorrelation"
  ) +
  labs(title = "Correlation Between Tumor Size and TTFD by Symptom",
       x = "", y = "") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 11))

### cox time dependent model + KM curve

#AP
PRO_AP <- PRO %>%
  filter(
    !is.na(TTFD10_AP), TTFD10_AP != -999,
    !is.na(Ratio_TS_TTFD10_AP), Ratio_TS_TTFD10_AP != -999,
    !is.na(AGE), AGE != -999,
    
    # ECOG: remove "UNK" or any other NA-like entries
    ECOG %in% c("0", "1", "2-3"),
    
    # MENOS: remove "UNK"
    MENOS != "UNK",
    
    # CHEMO: remove "UNK"
    CHEMO != "UNK",
    
    # HORMON: remove "UNK"
    HORMON != "UNK",
    
    # ERS: remove "UNK"
    ERS != "UNK",
    
    # PGRS: remove "UNK"
    PGRS != "UNK"
  ) %>%
  mutate(event_AP = ifelse(censored_AP == 0, 1, 0))


cox_model <- coxph(Surv(TTFD10_AP, event_AP) ~ Ratio_TS_TTFD10_AP + AGE + MENOS + CHEMO + HORMON + 
                     + ERS + PGRS + ECOG, data = PRO_AP)
summary(cox_model)

# result: the larger the tumor growth ratio, the shorter the TTFD (HR = 2.26, p = 5e-07)

# group by median (high vs low tumor growth)
median_ratio <- median(PRO_AP$Ratio_TS_TTFD10_AP, na.rm = TRUE)
PRO_AP <- PRO_AP %>%
  mutate(GrowthGroup = ifelse(Ratio_TS_TTFD10_AP >= median_ratio, "High Growth", "Low Growth"))

# create Surv object and plot KM curves
surv_obj <- Surv(PRO_AP$TTFD10_AP, PRO_AP$event_AP)

km_fit <- survfit(surv_obj ~ GrowthGroup, data = PRO_AP)

# KM plot
ggsurvplot(
  km_fit,
  data = PRO_AP,
  pval = TRUE,                   # add log-rank p-value
  risk.table = TRUE,            # add risk table
  conf.int = TRUE,              # confidence intervals
  palette = c("#E64B35", "#4DBBD5"),
  legend.title = "Tumor Growth",
  legend.labs = c("High Growth", "Low Growth"),
  xlab = "Time to AP Deterioration (days)",
  ylab = "Survival Probability (Not Deteriorated)",
  title = "Kaplan-Meier Curve: Tumor Growth vs. Symptom AP Deterioration",
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE
)

### 15 symptoms + Cox + KM
# 15 symptoms (abbreviations)
symptoms <- c("AP","CO","DI","DY","FA","FI","GHS","NV","PA","SL","PF2","RF2","EF","CF","SF")

symptom_full <- c(
  "FA"="Fatigue","NV"="Nausea/Vomit","PA"="Pain","DY"="Dyspnea",
  "SL"="Insomnia","AP"="Appetite loss","CO"="Constipation","DI"="Diarrhea",
  "FI"="Financial","PF2"="Physical","RF2"="Role","EF"="Emotional",
  "CF"="Cognitive","SF"="Social","GHS"="GHS/QOL"
)

out_dir <- "~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/Figure_summary/Tumor_Size/KM"
dir.create(path.expand(out_dir), recursive = TRUE, showWarnings = FALSE)

plot_km_one <- function(sym){
  ttfd     <- paste0("TTFD10_", sym)
  cens     <- paste0("censored_", sym)
  ratio    <- paste0("Ratio_TS_TTFD10_", sym)
  full_nm  <- symptom_full[[sym]]
  
  dat <- PRO %>%
    filter(
      !is.na(.data[[ttfd]]), .data[[ttfd]] != -999,
      !is.na(.data[[ratio]]), .data[[ratio]] != -999,
      !is.na(AGE), AGE != -999,
      ECOG %in% c("0", "1", "2-3"),
      MENOS != "UNK", CHEMO != "UNK", HORMON != "UNK",
      ERS   != "UNK", PGRS  != "UNK"
    ) %>%
    mutate(event = ifelse(.data[[cens]] == 0, 1, 0))
  
  # proceed only if each group has ≥5 samples
  if(nrow(dat) < 10) {
    message(sprintf("[%s] insufficient sample size (n=%d), skipped.", full_nm, nrow(dat)))
    return(invisible(NULL))
  }
  
  # median split
  med <- median(dat[[ratio]], na.rm = TRUE)
  dat <- dat %>%
    mutate(GrowthGroup = ifelse(.data[[ratio]] >= med, "High Growth", "Low Growth"),
           GrowthGroup = factor(GrowthGroup, levels = c("Low Growth","High Growth")))
  
  # KM
  surv_obj <- Surv(dat[[ttfd]], dat$event)
  fit <- survfit(surv_obj ~ GrowthGroup, data = dat)
  
  p <- ggsurvplot(
    fit, data = dat,
    pval = TRUE, risk.table = TRUE, conf.int = TRUE,
    palette = c("#4DBBD5", "#E64B35"),  # Low, High
    legend.title = "Tumor Growth",
    legend.labs  = c("Low Growth", "High Growth"),
    xlab = paste0("Time to ", full_nm, " Deterioration (days)"),
    ylab = "Survival Probability (Not Deteriorated)",
    title = paste0("Kaplan–Meier: Tumor Growth vs. ", full_nm, " Deterioration"),
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE
  )
  
  # safe, readable filename
  safe_name <- str_replace_all(full_nm, "[^A-Za-z0-9]+", "_")
  f_pdf <- file.path(path.expand(out_dir), paste0("KM_", sym, "_", safe_name, ".pdf"))
  f_png <- file.path(path.expand(out_dir), paste0("KM_", sym, "_", safe_name, ".png"))
  
  # save (main plot + risk table)
  full_plot <- ggpubr::ggarrange(p$plot, p$table, ncol = 1, heights = c(2, 0.9))
  ggsave(f_pdf, plot = full_plot, width = 9, height = 6.8, device = "pdf")
  ggsave(f_png, plot = full_plot, width = 9, height = 6.8, dpi = 320)
  
  message(sprintf("[%s] saved: %s and %s", full_nm, f_pdf, f_png))
  invisible(p)
}

# batch
invisible(lapply(symptoms, plot_km_one))



### Hazard ratio plot
symptoms <- c("AP", "CO", "DI", "DY", "FA", "FI", "GHS", "NV", 
              "PA", "SL", "PF2", "RF2", "EF", "CF", "SF")

symptom_names <- c(
  "FA"  = "Fatigue",
  "NV"  = "Nausea/Vomit",
  "PA"  = "Pain",
  "DY"  = "Dyspnea",
  "SL"  = "Insomnia",
  "AP"  = "Appetite loss",
  "CO"  = "Constipation",
  "DI"  = "Diarrhea",
  "FI"  = "Financial",
  "PF2" = "Physical",
  "RF2" = "Role",
  "EF"  = "Emotional",
  "CF"  = "Cognitive",
  "SF"  = "Social",
  "GHS" = "GHS/QOL"
)

results <- data.frame()

for (sym in symptoms) {
  ttfd_col <- paste0("TTFD10_", sym)
  cens_col <- paste0("censored_", sym)
  ratio_col <- paste0("Ratio_TS_TTFD10_", sym)
  
  dat <- PRO %>%
    filter(!is.na(.data[[ttfd_col]]),
           !is.na(.data[[ratio_col]]),
           !is.na(AGE),
           !is.na(ECOG),
           .data[[ttfd_col]] != -999,
           .data[[cens_col]] != -999,
           .data[[ratio_col]] != -999,
           AGE != -999,
           !MENOS %in% "UNK",
           !CHEMO %in% "UNK",
           !HORMON %in% "UNK",
           !ERS %in% "UNK",
           !PGRS %in% "UNK",
           !ECOG %in% "UNK") %>%
    mutate(event = ifelse(.data[[cens_col]] == 0, 1, 0))
  
  if (nrow(dat) >= 30) {
    formula_str <- as.formula(paste0(
      "Surv(", ttfd_col, ", event) ~ ", ratio_col,
      " + AGE + MENOS + CHEMO + HORMON + ERS + PGRS + ECOG"
    ))
    
    cox_model <- coxph(formula_str, data = dat)
    sum_model <- summary(cox_model)
    
    hr <- sum_model$coefficients[1, "exp(coef)"]
    lower <- sum_model$conf.int[1, "lower .95"]
    upper <- sum_model$conf.int[1, "upper .95"]
    pval <- sum_model$coefficients[1, "Pr(>|z|)"]
    
    results <- rbind(results, data.frame(
      Symptom = sym,
      HR = hr,
      Lower_CI = lower,
      Upper_CI = upper,
      P_value = pval
    ))
  }
}

results$Symptom <- symptom_names[results$Symptom]

results <- results %>%
  arrange(HR) %>%
  mutate(Symptom = fct_inorder(Symptom),
         Significance = ifelse(P_value < 0.05, "p < 0.05", "ns"),
         Label = sprintf("%.2f (%.2f–%.2f)", HR, Lower_CI, Upper_CI))

ggplot(results, aes(x = Symptom, y = HR, ymin = Lower_CI, ymax = Upper_CI)) +
  geom_pointrange(aes(color = Significance), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  geom_text(aes(label = Label), hjust = -0.1, size = 4) +   # right-side HR (95% CI)
  scale_color_manual(values = c("p < 0.05" = "red", "ns" = "gray40")) +
  labs(
    title = "Hazard Ratio of Tumor Growth for Symptom Deterioration",
    y = "Hazard Ratio (HR)", x = "Symptom",
    color = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(5.5, 80, 5.5, 5.5)) 


###updated forest plot

# abbreviation → full name
symptom_names <- c(
  "FA"="Fatigue","NV"="Nausea/Vomit","PA"="Pain","DY"="Dyspnea",
  "SL"="Insomnia","AP"="Appetite loss","CO"="Constipation","DI"="Diarrhea",
  "FI"="Financial","PF2"="Physical","RF2"="Role","EF"="Emotional",
  "CF"="Cognitive","SF"="Social","GHS"="GHS/QOL"
)

# —— Key tweak: as.character() to avoid factor misalignment —— #
results_viz <- results %>%
  mutate(
    Symptom_chr = as.character(Symptom),
    Symptom_full = dplyr::coalesce(               # if already full name, keep it
      unname(symptom_names[Symptom_chr]),         # map by name
      Symptom_chr
    ),
    SigFlag = ifelse(P_value < 0.05, "p < 0.05", "ns"),
    Stars   = case_when(P_value < 0.001 ~ "***",
                        P_value < 0.01  ~ "**",
                        P_value < 0.05  ~ "*",
                        TRUE ~ ""),
    Label   = sprintf("%.2f (%.2f–%.2f)%s",
                      HR, Lower_CI, Upper_CI,
                      ifelse(Stars=="","", paste0("  ", Stars)))
  ) %>%
  select(-Symptom, -Symptom_chr) %>%
  rename(Symptom = Symptom_full) %>%
  mutate(Symptom = fct_reorder(Symptom, HR, .desc = FALSE))  # larger at the top

# colors
col_sig <- "#D81B60"; col_ns <- "#6B7280"

ggplot(results_viz, aes(y = Symptom, x = HR, color = SigFlag)) +
  geom_vline(xintercept = 1, linetype = "11", linewidth = 0.6, color = "#9aa0a6") +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI),
                 height = 0, linewidth = 1.1, alpha = 0.95) +
  geom_point(size = 3.5) +
  geom_text(aes(label = Label, x = pmax(Upper_CI, HR)),
            hjust = -0.1, size = 4.2, color = "black") +
  scale_color_manual(values = c("p < 0.05" = col_sig, "ns" = col_ns),
                     name = "Significance") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(title = "Hazard Ratio of Tumor Growth for Symptom Deterioration",
       x = "Hazard Ratio (HR)", y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    plot.title.position = "plot",
    axis.line.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    panel.grid.major.x = element_line(color = "#eeeeee"),
    plot.margin = margin(10, 120, 10, 10)
  ) +
  coord_cartesian(clip = "off")


