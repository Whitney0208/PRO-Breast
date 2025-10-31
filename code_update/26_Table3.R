rm(list = ls())
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")

Merged_TTFD <- read_csv("Merged_PRO_TTFD_01OCT2025.csv")

Merge_EQLQ_24APR2025 <- read_csv(
  "~/Desktop/Zhou Lab/Breast Cancer PRO project/Merge/wwanbing/output/Merge_EQLQ_24APR2025.csv"
)

# Merge STDID
Merged_TTFD <- Merged_TTFD %>%
  left_join(Merge_EQLQ_24APR2025 %>% distinct(UID, STDID), by = "UID")

df0 <- Merged_TTFD

# ==== Set grouping variable and identify symptoms =================================
year_var <- "STDID"
stopifnot(year_var %in% names(df0))

# Identify the 15 symptoms: start with TTFD10_ and exclude the P version (avoid TTFD10P_)
symptoms <- grep("^TTFD10_[A-Z0-9]+$", names(df0), value = TRUE)
symptoms <- sub("^TTFD10_", "", symptoms)

# ==== Metric function (for a single grouped dataset) ===============================
# Outputs three rows:
# 1) Has TTFD (median, min, max)
# 2) TTFD non-censored N (n%)
# 3) TTFD missing (NA or -999)
summ_one_group <- function(d){
  # d: must contain columns year, TTFD (with -999 already converted to NA), cens
  
  have_n <- sum(!is.na(d$TTFD))                         # number with valid TTFD (after removing -999)
  med    <- if (have_n > 0) median(d$TTFD, na.rm = TRUE) else NA_real_
  mn     <- if (have_n > 0) min(d$TTFD, na.rm = TRUE)     else NA_real_
  mx     <- if (have_n > 0) max(d$TTFD, na.rm = TRUE)     else NA_real_
  nc_n   <- sum(d$cens == 0 & !is.na(d$TTFD))           # non-censored and has TTFD
  pct    <- if (have_n > 0) 100 * nc_n / have_n else NA_real_
  
  na_n   <- sum(is.na(d$TTFD))                          # number missing (including -999)
  
  tibble(
    row_lab = c("Has TTFD (median, min, max)",
                "TTFD non-censored N (n%)",
                "TTFD missing (NA or -999)"),
    value   = c(
      ifelse(is.na(med), "NA",
             sprintf("%d (%d, %d)", round(med), round(mn), round(mx))),
      ifelse(is.na(pct), "NA",
             sprintf("%d (%.1f%%)", nc_n, pct)),
      sprintf("%d", na_n)
    )
  )
}

# ---- If you want to separately count “original NA” and “original -999”, use the extended version below to replace the function above ----
# (Disabled by default; only override if you truly need to split them)
# summ_one_group <- function(d_raw){
#   # d_raw: must contain columns year, time_raw (original), cens; -999 not yet converted to NA
#   total_n <- nrow(d_raw)
#   na_only <- sum(is.na(d_raw$time_raw))                     # original NA
#   neg999  <- sum(d_raw$time_raw == -999, na.rm = TRUE)      # original -999
#
#   d <- d_raw %>% mutate(TTFD = na_if(time_raw, -999))
#   have_n <- sum(!is.na(d$TTFD))
#   med    <- if (have_n > 0) median(d$TTFD, na.rm = TRUE) else NA_real_
#   mn     <- if (have_n > 0) min(d$TTFD, na.rm = TRUE)    else NA_real_
#   mx     <- if (have_n > 0) max(d$TTFD, na.rm = TRUE)    else NA_real_
#   nc_n   <- sum(d$cens == 0 & !is.na(d$TTFD))
#   pct    <- if (have_n > 0) 100 * nc_n / have_n else NA_real_
#
#   tibble(
#     row_lab = c("Has TTFD (median, min, max)",
#                 "TTFD non-censored N (n%)",
#                 "TTFD missing (NA or -999)",
#                 "TTFD missing (original NA)",
#                 "TTFD missing (original -999)"),
#     value   = c(
#       ifelse(is.na(med), "NA", sprintf("%d (%d, %d)", round(med), round(mn), round(mx))),
#       ifelse(is.na(pct), "NA", sprintf("%d (%.1f%%)", nc_n, pct)),
#       sprintf("%d", na_only + neg999),
#       sprintf("%d", na_only),
#       sprintf("%d", neg999)
#     )
#   )
# }
# --------------------------------------------------------------------------------

# ==== Build a wide table for a single symptom =====================================
make_table_for_symptom <- function(sym){
  time_col <- paste0("TTFD10_", sym)
  cens_col <- paste0("censored_", sym)
  
  if(!(time_col %in% names(df0) && cens_col %in% names(df0))){
    warning(sprintf("Skip %s: %s or %s missing", sym, time_col, cens_col))
    return(NULL)
  }
  
  dat <- df0 %>%
    transmute(
      year  = as.character(.data[[year_var]]),
      TTFD  = na_if(.data[[time_col]], -999),   # treat -999 as missing
      cens  = .data[[cens_col]]
    )
  
  # Group statistics by year (study)
  by_year <- dat %>%
    group_by(year) %>%
    group_modify(~summ_one_group(.x)) %>%
    ungroup()
  
  # Overall (total)
  total_tab <- summ_one_group(dat) %>% mutate(year = "total")
  
  # Combine into a wide table (rows: three metrics; columns: each year + total)
  # If some study years do not exist, the levels below can be reduced to those present
  all_levels <- c("1997","2000","2001","2004","total")
  present_levels <- intersect(all_levels, unique(c(by_year$year, "total")))
  
  wide <- bind_rows(by_year, total_tab) %>%
    mutate(year = factor(year, levels = present_levels)) %>%
    arrange(row_lab, year) %>%
    tidyr::pivot_wider(names_from = year, values_from = value) %>%
    mutate(Symptom = sym, .before = 1)
  
  wide
}

# ==== Generate tables for all symptoms (list) and the combined big table ==========
tables_per_symptom <- symptoms %>%
  set_names() %>%
  map(make_table_for_symptom)

table_all_symptoms <- tables_per_symptom %>%
  discard(is.null) %>%
  bind_rows()

# Total N per study
study_n <- c("1997" = 735,
             "2000" = 1564,
             "2001" = 226,
             "2004" = 213,
             "total" = 2738)

