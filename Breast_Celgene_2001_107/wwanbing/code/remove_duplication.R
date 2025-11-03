eqlq_long <- eqlq %>%
  #filter(VISIT_ID == 1) %>%  # Keep only rows where VISIT_ID equals 1
  pivot_longer(cols = EQLQ_001:EQLQ_030,  # Convert only EQLQ_001 to EQLQ_030 into long format
               names_to = "QSFLAG", 
               values_to = "DV") %>%
  mutate(
    FLAG = as.numeric(substr(QSFLAG, 6, 8)),  # Extract numerical question number
    TIME = EQLQDAY_000,  # Assign TIME from EQLQDAY_000
    DV = replace_na(DV, -999)  # Replace NA values in DV with -999
  ) %>%
  select(RSUBJ_ID, QSFLAG, FLAG, TIME, DV)%>% # Keep necessary columns
  distinct()  # Remove duplicated observations

eqlq_long_raw <- eqlq %>%
  #filter(VISIT_ID == 1) %>%  # Keep only rows where VISIT_ID equals 1
  pivot_longer(cols = EQLQ_001:EQLQ_030,  # Convert only EQLQ_001 to EQLQ_030 into long format
               names_to = "QSFLAG", 
               values_to = "DV") %>%
  mutate(
    FLAG = as.numeric(substr(QSFLAG, 6, 8)),  # Extract numerical question number
    TIME = EQLQDAY_000,  # Assign TIME from EQLQDAY_000
    DV = replace_na(DV, -999)  # Replace NA values in DV with -999
  ) %>%
  select(RSUBJ_ID, QSFLAG, FLAG, TIME, DV) # Keep necessary columns


df_raw <- eqlq_long_raw %>%
  count(RSUBJ_ID, name = "n_raw")

df_long <- eqlq_long %>%
  count(RSUBJ_ID, name = "n_long")

df_compare <- full_join(df_raw, df_long, by = "RSUBJ_ID") %>%
  mutate(
    n_raw  = replace_na(n_raw, 0),
    n_long = replace_na(n_long, 0),
    diff   = n_raw - n_long
  ) %>%
  arrange(RSUBJ_ID)

df_compare_diff <- df_compare %>%
  filter(diff != 0)


#add one variable for the EQLQ dataset, STAGE, and all assign "IV" in the dataset 
#add the date in the dataset file name, for example, name the dataset as Breast_Celgene_2001_107_16MAR2025.csv

df <- read.csv("Breast_Celgene_2001_107_updated.csv", stringsAsFactors = FALSE)

# add one variable STAGE, all assign "IV"
df$STAGE <- "IV"

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_Celgene_2001_107_", toupper(current_date), ".csv")

write.csv(df, new_filename, quote=FALSE,row.names=FALSE)







