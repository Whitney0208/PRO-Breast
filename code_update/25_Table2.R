rm(list=ls())
library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(haven)
library(tidyr)
library(purrr)
library(readr)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

Merged_TTFD <- read_csv("Merged_PRO_TTFD_01OCT2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output")

Merge_EQLQ_24APR2025 <- read_csv("Merge_EQLQ_24APR2025.csv")

Merged_TTFD <- Merged_TTFD %>%
  left_join(
    Merge_EQLQ_24APR2025 %>% distinct(UID, STDID),
    by = "UID"
  )

df0 <- Merged_TTFD

# Grouping column (if not STDID, change to your grouping column name)
year_var <- "STDID"
stopifnot(year_var %in% names(df0))

# Identify the 15 symptoms: start with TTFD10_ and not the P version (avoid TTFD10P_)
symptoms <- grep("^TTFD10_[A-Z0-9]+$", names(df0), value = TRUE)
symptoms <- sub("^TTFD10_", "", symptoms)

# Metric function (single dataset)
# summ_one_group <- function(d){
#   have_n <- sum(!is.na(d$TTFD))
#   med    <- if (have_n > 0) median(d$TTFD, na.rm = TRUE) else NA_real_
#   mn     <- if (have_n > 0) min(d$TTFD, na.rm = TRUE)    else NA_real_
#   mx     <- if (have_n > 0) max(d$TTFD, na.rm = TRUE)    else NA_real_
#   nc_n   <- sum(d$cens == 0 & !is.na(d$TTFD))
#   pct    <- if (have_n > 0) 100 * nc_n / have_n else NA_real_
#   
#   tibble(
#     row_lab = c("有 TTFD (median, min, max)", "TTFD non-censored N (n%)"),
#     value   = c(
#       ifelse(is.na(med), "NA",
#              sprintf("%d (%d, %d)", round(med), round(mn), round(mx))),
#       ifelse(is.na(pct), "NA",
#              sprintf("%d (%.1f%%)", nc_n, pct))
#     )
#   )
# }

summ_one_group <- function(d){
  have_n <- sum(!is.na(d$TTFD))
  if (have_n > 0) {
    q <- quantile(d$TTFD, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    med <- q[2]; q1 <- q[1]; q3 <- q[3]
  } else {
    med <- q1 <- q3 <- NA_real_
  }
  nc_n <- sum(d$cens == 1 & !is.na(d$TTFD))
  pct  <- if (have_n > 0) 100 * nc_n / have_n else NA_real_
  
  tibble(
    row_lab = c("TTFD (median, Q1, Q3)", "TTFD censored N (n%)"),
    value   = c(
      ifelse(is.na(med), "NA",
             sprintf("%d (%d, %d)", round(med), round(q1), round(q3))),
      ifelse(is.na(pct), "NA",
             sprintf("%d (%.1f%%)", nc_n, pct))
    )
  )
}

# Generate a wide table for a single symptom
make_table_for_symptom <- function(sym){
  time_col   <- paste0("TTFD10_", sym)
  cens_col   <- paste0("censored_", sym)
  if(!(time_col %in% names(df0) && cens_col %in% names(df0))){
    warning(sprintf("Skip %s: %s or %s missing", sym, time_col, cens_col))
    return(NULL)
  }
  
  dat <- df0 %>%
    transmute(
      year  = as.character(.data[[year_var]]),
      TTFD  = na_if(.data[[time_col]], -999),
      cens  = .data[[cens_col]]
    )
  
  # By year
  by_year <- dat %>%
    group_by(year) %>%
    group_modify(~summ_one_group(.x)) %>%
    ungroup()
  
  # total
  total_tab <- summ_one_group(dat) %>% mutate(year = "total")
  
  # Combine to wide table
  wide <- bind_rows(by_year, total_tab) %>%
    mutate(year = factor(year, levels = c("1997","2000","2001","2004","total"))) %>%
    arrange(row_lab, year) %>%
    tidyr::pivot_wider(names_from = year, values_from = value) %>%
    mutate(Symptom = sym, .before = 1)
  
  wide
}

# 1) One table per symptom (list)
tables_per_symptom <- symptoms %>%
  set_names() %>%
  map(make_table_for_symptom)

# 2) Combined master table (two rows × 15 symptoms)
table_all_symptoms <- tables_per_symptom %>%
  bind_rows()

# View a specific symptom (e.g., AP)
tables_per_symptom$AP

# View the combined master table
table_all_symptoms







