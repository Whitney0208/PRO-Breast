rm(list=ls())
library(plyr)
library(dplyr)
library(haven)
library(tidyr)

ae <- read_sas(data_file='ae.sas7bdat') # AE
aenone <- read_sas(data_file='aenone.sas7bdat') # not useful
antitumo <- read_sas(data_file='antitumo.sas7bdat') #chemo/hormonal therapy
bcrelap <- read_sas(data_file='bcrelap.sas7bdat') # relapse day, regional or distant, which organ
bl_chem <- read_sas(data_file='bl_chem.sas7bdat') # lab test, not included
bloodhdr <- read_sas(data_file='bloodhdr.sas7bdat') # not useful
chf <- read_sas(data_file='chf.sas7bdat') #drug dose
chfhdr <- read_sas(data_file='chfhdr.sas7bdat') #not useful
creat_cl <- read_sas(data_file='creat_cl.sas7bdat') #not useful
creathdr <- read_sas(data_file='creathdr.sas7bdat') #not useful
criteria <- read_sas(data_file='criteria.sas7bdat') # exclusion/inclusion criteria
death <- read_sas(data_file='death.sas7bdat') # death week
delay <- read_sas(data_file='delay.sas7bdat') # not useful
demog <- read_sas(data_file='demog.sas7bdat')
diag1 <- read_sas(data_file='diag1.sas7bdat') #prior surgery
diag2 <- read_sas(data_file='diag2.sas7bdat') #LN resection
diag3 <- read_sas(data_file='diag3.sas7bdat') #baseline primary tumor size, stage
diaghdr <- read_sas(data_file='diaghdr.sas7bdat') #not useful
disctamo <- read_sas(data_file='disctamo.sas7bdat') #not useful
ekg <- read_sas(data_file='ekg.sas7bdat') # not useful
endchemo <- read_sas(data_file='endchemo.sas7bdat') # not useful
euroqol <- read_sas(data_file='euroqol.sas7bdat')
event_dt <- read_sas(data_file='event_dt.sas7bdat') # not useful
evidis <- read_sas(data_file='evidis.sas7bdat')
febneut <- read_sas(data_file='febneut.sas7bdat') #fever
fishlab <- read_sas(data_file='fishlab.sas7bdat') #HER2/NEU
fshlh <- read_sas(data_file='fshlh.sas7bdat')  # not useful
fshlhhdr <- read_sas(data_file='fshlhhdr.sas7bdat') # not useful
hematol <- read_sas(data_file='hematol.sas7bdat') #lab test
hemhdr <- read_sas(data_file='hemhdr.sas7bdat') #not useful
hormrec <- read_sas(data_file='hormrec.sas7bdat') #receptor
hospi <- read_sas(data_file='hospi.sas7bdat') #hospital admission
hospihdr <- read_sas(data_file='hospihdr.sas7bdat') #not useful
labels <- read_sas(data_file='labels.sas7bdat') #not useful
lvef <- read_sas(data_file='lvef.sas7bdat') #not useful
lvefnd <- read_sas(data_file='lvefnd.sas7bdat') #not useful
meno <- read_sas(data_file='meno.sas7bdat') #menopause/hormone treatment
mx <- read_sas(data_file='mx.sas7bdat') # anxiety depression
mxnone <- read_sas(data_file='mxnone.sas7bdat') #not useful
nonsysno <- read_sas(data_file='nonsysno.sas7bdat') #not useful
nsysanti <- read_sas(data_file='nsysanti.sas7bdat') # surgery and radiation
othproc <- read_sas(data_file='othproc.sas7bdat')  #not useful
outcare <- read_sas(data_file='outcare.sas7bdat') #not useful
patstat <- read_sas(data_file='patstat.sas7bdat') #relapse, death, lost follow up, status
pctx <- read_sas(data_file='pctx.sas7bdat') # dose
pctxna <- read_sas(data_file='pctxna.sas7bdat') #not useful
pe <- read_sas(data_file='pe.sas7bdat') #not useful
phneopl <- read_sas(data_file='phneopl.sas7bdat') #histology
pregtest <- read_sas(data_file='pregtest.sas7bdat') # pregnancy test
ptxd <- read_sas(data_file='ptxd.sas7bdat') # drug treatment
ptxdhdr <- read_sas(data_file='ptxdhdr.sas7bdat') # not useful
qol1 <- read_sas(data_file='qol1.sas7bdat') # Q1-28
qol2 <- read_sas(data_file='qol2.sas7bdat') #Q29-30
qol3 <- read_sas(data_file='qol3.sas7bdat')
qolhdr <- read_sas(data_file='qolhdr.sas7bdat') # time
qualife <- read_sas(data_file='qualife.sas7bdat')
radhdr <- read_sas(data_file='radhdr.sas7bdat')
radiscon <- read_sas(data_file='radiscon.sas7bdat') # not useful
radthera <- read_sas(data_file='radthera.sas7bdat')
rando <- read_sas(data_file='rando.sas7bdat')
reptxd <- read_sas(data_file='reptxd.sas7bdat') # not useful
scan <- read_sas(data_file='scan.sas7bdat')
scannd <- read_sas(data_file='scannd.sas7bdat') # not useful
signatur <- read_sas(data_file='signatur.sas7bdat') # not useful
spmalig <- read_sas(data_file='spmalig.sas7bdat') # other cancer
sysnon <- read_sas(data_file='sysnon.sas7bdat') # not useful
tsship <- read_sas(data_file='tsship.sas7bdat') # not useful
vital <- read_sas(data_file='vital.sas7bdat')

#000301-000-901-180 has duplicate AGENO(68, NA)
demog2 <- demog %>%
  transmute(
    UID = RUSUBJID,
    ID  = RSUBJID,
    AGE = AGENO,
    SEX = SEXCD,
    RACE    = -999,
    ETHNIC  = -999,               
    REGION  = -999,
    ARM     = 0,                  
  ) %>%
  distinct()%>%
  mutate(
    SEX = case_when(
      SEX == "Male"   ~ 0,
      SEX == "Female" ~ 1,
      TRUE            ~ -999
    ))%>%
  group_by(ID) %>%        #remove 000301-000-901-180 AGENO(NA)
  filter(!is.na(AGE)) %>%
  slice(1) %>%
  ungroup()

#clean *
diag3_clean <- diag3 %>%
  mutate(
    PATHOT = gsub("\\*", "", PATHOT),
    PATHON = gsub("\\*", "", PATHON),
    PATHOM = gsub("\\*", "", PATHOM)
  )

diag_stage <- diag3_clean %>%
  # DIAGC = "pT2,pN1,M0"
  mutate(
    DIAGC = paste(PATHOT, PATHON, PATHOM, sep = ",")
  ) %>%
  # coding
  mutate(
    STAGE = case_when(
      DIAGC %in% c("pT1,pN1,M1", "pT2,pN1,M1") ~ "IV",
      DIAGC %in% c("pT4,pN1,M0", "pT1,pN2,M0", "pT3,pN2,M0") ~ "III",
      DIAGC %in% c("*pTIS,pN1,M0") ~ "I",
      DIAGC %in% c("pT1,pN1,M0", "pT2,pN1,M0", "pT3,pN1,M0") ~ "II",
      TRUE ~ "-999"
    )
  ) %>%
  transmute(
    UID = RUSUBJID,
    STAGE
  )

euroqol2 <- euroqol %>%
  filter(EVENT_ID == "BASELINE")%>%
  transmute(
    UID = RUSUBJID,
    EMPLOY = case_when(
      ACTLST %in% c("keeping house", "seeking work") ~ 0,
      ACTLST == "employed" ~ 1,
      ACTLST == "retired" ~ 2,
      ACTLST %in% c("student", "other") ~ 3,
      ACTLST %in% c(".", "A", "7","0") | is.na(ACTLST) ~ -999,
      TRUE ~ -999
    )
  )%>%
  distinct()

#euroqol2 %>%
#  group_by(UID) %>%
#  filter(n() > 1) %>%
#  print(n = Inf)  
#
# UID                EMPLOY
# <chr>               <dbl>
# 1 000301-000-901-595      3
# 2 000301-000-901-595   -999
# 3 000301-000-999-083   -999
# 4 000301-000-999-083      1
# 5 000301-000-999-613      1
# 6 000301-000-999-613   -999
# 7 000301-000-999-796   -999
# 8 000301-000-999-796      0

euroqol3 <- euroqol2 %>%
  group_by(UID) %>%
  arrange(EMPLOY == -999) %>%  # Put -999 in the back
  summarise(EMPLOY = first(EMPLOY), .groups = "drop")

pregtest2 <- pregtest %>%
  transmute(
    UID = RUSUBJID,
    MENOS= case_when(
      CBPYN == "Yes" ~ 1,         # Premenopausal
      CBPYN == "No" ~ 3,          # Postmenopausal
      TRUE ~ -999
    )
  )%>%
  distinct()

vital2 <- vital %>%
  filter(EVENT_ID == "BASELINE" & PAG_NAME == "BR1") %>%
  distinct()%>%
  mutate(
    ECOG = case_when(
      PSKAR %in% c(110, 100, 90) ~ 0,
      PSKAR %in% c(80, 70) ~ 1,
      PSKAR %in% c(60, 50, 40, 30, 20, 0) ~ 2,
      PSKAR %in% c(".", "A", "", " ") | is.na(PSKAR) ~ -999,
      TRUE ~ -999
    )
  ) %>%
  select(UID = RUSUBJID, ECOG)

# 000301-000-901-378 NA
nsysanti2 <- nsysanti %>%
  filter(EVENT_ID == "BASELINE") %>%
  group_by(RUSUBJID) %>%
  summarise(
    SURGERY = case_when(
      all(is.na(TYPNON) | TYPNON == "") ~ -999L,
      any(TYPNON == "Surgery", na.rm = TRUE) ~ 1L,
      TRUE ~ 2L
    ),
    RADIO = case_when(
      all(is.na(TYPNON) | TYPNON == "") ~ -999L,
      any(TYPNON == "Radiotherapy", na.rm = TRUE) ~ 1L,
      TRUE ~ 2L
    )
  ) %>%
  rename(UID = RUSUBJID)

#Duplicate ID: 1 000301-000-901-550  -999  -999
#              2 000301-000-901-550     1     1
hormrec2 <- hormrec %>% 
  filter(HRMTEST == "Immunohistochemistry") %>%
  mutate(
    ERS = case_when(
      ERSTA == "Positive" ~ 1,
      ERSTA == "Negative" ~ 2,
      ERSTA %in% c("Not Done", "Not Assessable/Not Done") ~ 3,
      ERSTA %in% c(".", "", NA) ~ -999,
      TRUE ~ -999
    ),
    PGRS = case_when(
      PGRSTA == "Positive" ~ 1,
      PGRSTA == "Negative" ~ 2,
      PGRSTA %in% c("Not Done", "Not Assessable/Not Done") ~ 3,
      PGRSTA %in% c(".", "", NA) ~ -999,
      TRUE ~ -999
    )
  ) %>%
  select(UID = RUSUBJID, ERS, PGRS)

#Remove duplicate ID
hormrec3 <- hormrec2 %>%
  group_by(UID) %>%
  mutate(n_uid = n()) %>%  # 统计每个 UID 的行数
  filter(
    (n_uid == 1) |                                # If there is only one line, keep it
      (n_uid > 1 & !(ERS == -999 & PGRS == -999))   # If there are more than one line, remove all -999 lines
  ) %>%
  slice(1) %>%  # If there are duplicate UIDs, only the first valid record is retained
  ungroup() %>%
  select(-n_uid)  # Remove the auxiliary column

fishlab2 <- fishlab %>%
  mutate(
    HER2 = case_when(
      FISHRES == "Her2/Neu Positive" ~ 1,
      FISHRES == "Her2/Neu Negative" ~ 2,
      FISHRES %in% c(".", "", NA) ~ -999,
      TRUE ~ -999
    )
  ) %>%
  select(UID = RUSUBJID, HER2)

# View duplicate UID lines
#fishlab2 %>%
#  filter(duplicated(UID) | duplicated(UID, fromLast = TRUE))

mx2 <- mx %>%
  group_by(RUSUBJID) %>%
  summarise(
    MHDEPRESSION = case_when(
      any(MXCL == "DEPRESSION", na.rm = TRUE) ~ 1,
      all(MXCL == " " | is.na(MXCL)) ~ -999,
      TRUE ~ 2
    ),
    MDANXIETY = case_when(
      any(MXCL == "ANXIETY", na.rm = TRUE) ~ 1,
      all(MXCL == " " | is.na(MXCL)) ~ -999,
      TRUE ~ 2
    )
  ) %>%
  ungroup()%>%
  rename(UID = RUSUBJID)

# HEIGHT & WEIGHT
vital3 <- vital %>%
  # Mark some value as NA
  mutate(
    HT_chr = ifelse(HT %in% c(".", " ", "A") | is.na(HT), NA, as.character(HT)),
    WT_chr = ifelse(WT %in% c(".", " ", "A") | is.na(WT), NA, as.character(WT))
  ) %>%
  # Convert to numeric
  mutate(
    HT_num = as.numeric(HT_chr),
    WT_num = as.numeric(WT_chr)
  ) %>%
  # Based on UID, take the first valid HEIGHT and WEIGHT of each group
  group_by(RUSUBJID) %>%
  summarise(
    HEIGHT = {
      ht_row <- which(!is.na(HT_num))[1]
      if (is.na(ht_row)) {
        -999
      } else {
        if (HT_U[ht_row] == "cm") HT_num[ht_row]
        else if (HT_U[ht_row] == "in") HT_num[ht_row] * 2.54
        else -999
      }
    },
    WEIGHT = {
      wt_row <- which(!is.na(WT_num))[1]
      if (is.na(wt_row)) {
        -999
      } else {
        if (WT_U[wt_row] == "kg") WT_num[wt_row]
        else if (WT_U[wt_row] == "lb") WT_num[wt_row] * 0.45
        else -999
      }
    },
    .groups = "drop"
  ) %>%
  rename(UID = RUSUBJID)



# DTHDY & DTH
death_vars <- death %>%
  transmute(
    UID = RUSUBJID,
    DTHDY = DTH_WEEK * 7,
    DTH = 1
  )

# PFSDY & PFS
bcrelap_vars <- bcrelap %>%
  transmute(
    UID = RUSUBJID,
    PFSDY = RLPS_DY,
    PFS = 1
  )

patstat2 <- full_join(death_vars, bcrelap_vars, by = "UID") %>%
  mutate(
    DTHDY = ifelse(is.na(DTHDY), -999, DTHDY),
    DTH   = ifelse(is.na(DTH), -999, DTH),
    PFSDY = ifelse(is.na(PFSDY), -999, PFSDY),
    PFS   = ifelse(is.na(PFS), -999, PFS)
  )%>%
  distinct()

#USE first() to remove duplicate ID
patstat3 <- patstat2 %>%
  group_by(UID) %>%
  slice(1) %>%  # Keep the first record of each UID
  ungroup()

#eqlq

# > unique(qol1$QOLRESP)
# [1] "A Little"    "Quite a Bit" "Not at All"  "Very Much"   "."           "G"

qol1_2 <- qol1 %>%
  mutate(QOLRESP = case_when(
    QOLRESP == "Not at All"   ~ 1,
    QOLRESP == "A Little"     ~ 2,
    QOLRESP == "Quite a Bit"  ~ 3,
    QOLRESP == "Very Much"    ~ 4,
    QOLRESP %in% c(".", "G")  ~ -999,
    TRUE                      ~ -999
  )) 

# merge qol1_2 & qol2
combined_qol <- bind_rows(
  qol1_2 %>% rename(QOLRESP_FINAL = QOLRESP),
  qol2 %>% rename(QOLRESP_FINAL = QOLRESP7)
)

# Sort: press RUSUBJID first, then EVENT_ID
combined_qol <- combined_qol %>%
  arrange(RUSUBJID, EVENT_ID, QOL_NO)

# Merge qolhdr
# Extract the unique matching pair from qolhdr
qolhdr_clean <- qolhdr %>%
  select(EVENT_ID, RUSUBJID, QOL_DY) %>%
  distinct(EVENT_ID, RUSUBJID, .keep_all = TRUE)

combined_qol2 <- combined_qol %>%
  left_join(qolhdr_clean, by = c("EVENT_ID", "RUSUBJID"))

combined_qol3 <- combined_qol2 %>%
  arrange(RUSUBJID, EVENT_ID, QOL_NO)

# Step 1: Long format 
qol_long <- combined_qol3 %>%
  mutate(
    QSFLAG = paste0("Q", sprintf("%02d", QOL_NO)),  # Q01 Q02 ...
    FLAG = QOL_NO,                                  # Question number
    DV = QOLRESP_FINAL,                             # answer
    TIME = QOL_DY                                   # time
  ) %>%
  mutate(
    DV = ifelse(is.na(DV), -999, DV),
    TIME
  ) %>%
  select(RUSUBJID, QSFLAG, FLAG, TIME, DV) %>%
  distinct()

qol_long2 <- qol_long %>%
  group_by(RUSUBJID, FLAG, TIME) %>%
  arrange(
    (DV == -999 & any(DV != -999)),  # Put "-999 of groups mixed with -999" at the back
    row_number()                     # The rest of the order does not move
  ) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(RUSUBJID, TIME, FLAG)

# Step 2: Wide format - baseline, total lines=1567(1649)
qol_wide <- qol_long2 %>%
  arrange(RUSUBJID, TIME, FLAG) %>%  
  group_by(RUSUBJID) %>%
  mutate(
    all_na = all(is.na(TIME)),
    TIME_valid = ifelse(is.na(TIME), Inf, TIME),
    baseline_time = ifelse(all_na, NA, min(TIME_valid))
  ) %>%
  # Retain the entire group of 30 items corresponding to the baseline time or retain the first 30 items when all nas are present
  mutate(row_idx = row_number()) %>%
  filter(
    (!all_na & TIME == baseline_time) |
      (all_na & row_idx <= 30)
  ) %>%
  ungroup() %>%
  select(RUSUBJID, FLAG, DV) %>%
  mutate(FLAG = paste0("B", sprintf("%02d", FLAG))) %>%
  tidyr::pivot_wider(names_from = FLAG, values_from = DV) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), -999, .x)))


# Step 3: Merge long and wide formats
qol_final <- qol_long2 %>%
  left_join(qol_wide, by = "RUSUBJID") %>%
  rename(UID = RUSUBJID)

qol_final <- qol_final %>%
  mutate(across(everything(), ~ ifelse(.x == "", NA, .x))) %>%  # Convert empty strings to NA
  mutate(
    across(where(is.numeric), ~ replace_na(.x, -999)),  
    across(where(is.character), ~ replace_na(.x, "-999"))  
  )

#demog2 diag_stage euroqol2 pregtest2 vital2 nsysanti2 hormrec2 fishlab2 mx2 vital3 patstat2
final_data <- demog2 %>%
  left_join(diag_stage,   by = "UID") %>%
  left_join(euroqol3, by = "UID") %>% 
  left_join(pregtest2,   by = "UID") %>%
  left_join(vital2,    by = "UID") %>% 
  left_join(nsysanti2,   by = "UID") %>%
  left_join(hormrec3, by = "UID") %>%
  left_join(fishlab2,   by = "UID") %>%
  left_join(mx2,    by = "UID") %>% 
  left_join(vital3,    by = "UID") %>% 
  left_join(patstat3,   by = "UID") %>%
  distinct()

#qol has no value for 000301-000-901-000
final_data2 <- final_data %>%
  full_join(qol_final, by = "UID")   # Perform a full join to keep all IDs

dataset2 <- final_data2 %>%
  mutate(across(everything(), ~ ifelse(.x == "", NA, .x))) %>%  # Convert empty strings to NA
  mutate(
    across(where(is.numeric), ~ replace_na(.x, -999)),  
    across(where(is.character), ~ replace_na(.x, "-999"))  
  )%>%
  arrange(UID, TIME, FLAG)

# dataset2 SURGERY RADIO -999 -> 2
dataset3 <- dataset2 %>%
  mutate(
    SURGERY = if_else(SURGERY == -999, 2L, SURGERY),
    RADIO = if_else(RADIO == -999, 2L, RADIO)
  )

current_date <- format(Sys.Date(), "%d%b%Y")
new_filename <- paste0("Breast_SanofiU_2000_118_", toupper(current_date), ".csv")

write.csv(dataset3, new_filename, quote=FALSE,row.names=FALSE)


