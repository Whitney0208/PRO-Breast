rm(list=ls())

library(reader)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
Merge_TTFD <- read_csv("Merged_PRO_TTFD_01OCT2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output")

Merge_EQLQ_24APR2025 <- read_csv("Merge_EQLQ_24APR2025.csv")

merge_eqlq <- Merge_EQLQ_24APR2025 %>%
  dplyr::filter(UID %in% Merge_TTFD$UID)

length(unique(merge_eqlq$UID))  # 2738

colnames(merge_eqlq)
# > colnames(merge_eqlq)
# [1] "UID"          "ID"           "STAGE"        "AGE"          "SEX"          "RACE"         "ARM"          "EMPLOY"      
# [9] "MENOS"        "LESION1"      "ECOG"         "CHEMO"        "SURGERY"      "RADIO"        "HORMON"       "ERS"         
# [17] "PGRS"         "HER2"         "MHDEPRESSION" "MDANXIETY"    "HEIGHT"       "WEIGHT"       "DTHDY"        "DTH"         
# [25] "PFSDY"        "PFS"          "QSFLAG"       "FLAG"         "DV"           "TIME"         "B01"          "B02"         
# [33] "B03"          "B04"          "B05"          "B06"          "B07"          "B08"          "B09"          "B10"         
# [41] "B11"          "B12"          "B13"          "B14"          "B15"          "B16"          "B17"          "B18"         
# [49] "B19"          "B20"          "B21"          "B22"          "B23"          "B24"          "B25"          "B26"         
# [57] "B27"          "B28"          "B29"          "B30"          "STDID"        "NCT"  

dat <- merge_eqlq %>%
  select(UID, ID, STAGE, AGE, SEX, RACE, ARM, EMPLOY,
         MENOS, LESION1, ECOG, CHEMO, SURGERY, RADIO, HORMON, ERS,
         PGRS, HER2, MHDEPRESSION, MDANXIETY, HEIGHT, WEIGHT,
         DTHDY, DTH, PFSDY, PFS, STDID) %>%
  distinct()

# 1) Variable lists
cont_vars <- c("AGE", "LESION1", "HEIGHT", "WEIGHT", "DTHDY", "PFSDY")
cat_vars  <- c("SEX", "STAGE", "RACE", "ARM", "EMPLOY", "MENOS", "ECOG",
               "CHEMO", "SURGERY", "RADIO", "HORMON", "ERS", "PGRS", "HER2",
               "MHDEPRESSION", "MDANXIETY", "DTH", "PFS")

# Generate a table for each categorical variable
cat_levels <- lapply(cat_vars, function(v) {
  tab <- table(dat[[v]], useNA = "ifany")
  list(var = v, levels = tab)
})

# Print results
for (x in cat_levels) {
  cat("\n==========", x$var, "==========\n")
  print(x$levels)
}

# Categorical variables: NA -> "-999", keep the original "UNK"
na_to_minus999_chr <- function(x){
  x <- as.character(x)
  x[is.na(x)] <- "-999"
  x
}

# Continuous variables: treat -999 as missing (set to NA) for summaries
minus999_to_na_num <- function(x){
  x <- as.numeric(x)
  x[x == -999] <- NA
  x
}

# Output format for continuous variables
# summarize_cont <- function(x){
#   if (all(is.na(x))) return(NA_character_)
#   sprintf("%.1f (%.1f, %.1f)",
#           stats::median(x, na.rm = TRUE),
#           min(x, na.rm = TRUE),
#           max(x, na.rm = TRUE))
# }

summarize_cont <- function(x){
  if (all(is.na(x))) return(NA_character_)
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.1f (%.1f, %.1f)", q[2], q[1], q[3])
}

# ---- Cleaning ----
dat_clean <- dat %>%
  # Categorical variables: only convert NA -> "-999"
  mutate(across(all_of(cat_vars), na_to_minus999_chr)) %>%
  # Continuous variables: only convert -999 -> NA (exclude from stats)
  mutate(across(all_of(cont_vars), minus999_to_na_num)) %>%
  mutate(STDID = as.character(STDID))

# ---- Continuous variable summaries (-999 already excluded) ----
cont_summary <- dat_clean %>%
  group_by(STDID) %>%
  summarise(across(all_of(cont_vars), summarize_cont), .groups = "drop")

cont_total <- dat_clean %>%
  summarise(across(all_of(cont_vars), summarize_cont)) %>%
  mutate(STDID = "total")

cont_summary <- bind_rows(cont_summary, cont_total)

# ---- Categorical variable summaries (keep "-999" and "UNK" as distinct levels) ----
summarize_cat_one <- function(var, df){
  tab <- table(df[[var]], useNA = "no")  # "-999" is already a character level
  n   <- sum(tab)
  tibble(
    variable = var,
    level    = as.character(names(tab)),
    value    = sprintf("%d (%.1f%%)", as.integer(tab), 100*as.numeric(tab)/n)
  )
}

cat_by_year <- dat_clean %>%
  group_split(STDID) %>%
  map_dfr(function(df){
    yr <- unique(df$STDID)
    bind_rows(lapply(cat_vars, summarize_cat_one, df = df)) %>%
      mutate(STDID = yr)
  })

cat_total <- bind_rows(lapply(cat_vars, summarize_cat_one, df = dat_clean)) %>%
  mutate(STDID = "total")

cat_summary <- bind_rows(cat_by_year, cat_total)

# ---- Combine into a wide table ----
cont_long <- cont_summary %>%
  pivot_longer(-STDID, names_to = "variable", values_to = "value") %>%
  mutate(level = NA_character_)

final_wide <- bind_rows(cont_long, cat_summary) %>%
  select(variable, level, STDID, value) %>%
  pivot_wider(names_from = STDID, values_from = value)

final_wide









