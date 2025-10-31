# Baseline  survival correlation variable (lasso selection) 

# Covariate: exclude baseline flag... 

library(tidyverse)
library(plyr)
library(dplyr)
library(survival)
library(Matrix)
library(glmnet)
library(GGally)

select = dplyr::select
filter = dplyr::filter
rename = dplyr::rename
mutate = dplyr::mutate
relocate = dplyr::relocate
summarise = dplyr::summarise

#Fatigue=Breast_PRO_Fatigue_updated_19MAY2025

fatigue_baseline <- Fatigue %>%
  group_by(UID) %>%
  filter(TIME == min(TIME, na.rm = TRUE)) %>%
  ungroup()

fatigue_baseline <- fatigue_baseline %>%
  mutate(RACE = as.character(RACE))

OS <- fatigue_baseline %>%
  filter(DTHDY != -999 & DTH != -999) %>%
  select(UID, DTHDY, DTH, STAGE:WEIGHT) %>%
  distinct()

OS[OS == -999] <- NA

# rescale continuous variables
set.seed(1234)
OS_scaled <- OS %>%
  # Convert all character to factor
  mutate(across(where(is.character), as.factor)) %>%
  # Scale numeric columns except DTHDY and DTH
  mutate(across(where(is.numeric) & !any_of(c("DTHDY", "DTH")), scale)) %>%
  # Remove rows with NA
  na.omit()

# corrlation plot(continuous)
con.var <- OS %>% select(where(is.numeric)) %>% select(-DTHDY, -DTH)
pdf("plot/Con_covariate_correlation_fatigue_OS.pdf", width=12, height=12)
print(ggpairs(con.var))
dev.off()

# Step 3: LASSO cross-validation to select covariates
set.seed(1234)

# filter DTHDY <= 0 
OS_scaled_filtered <- OS_scaled %>% filter(DTHDY > 0)

# Delete UID、DTHDY、DTH
x_data <- OS_scaled_filtered %>% 
  select(-c(UID, DTHDY, DTH)) %>% 
  data.matrix()

y_data <- Surv(OS_scaled_filtered$DTHDY, OS_scaled_filtered$DTH)

cvfitOS <- cv.glmnet(
  x = x_data,
  y = y_data,
  family = "cox",
  alpha = 1,
  nfolds = 5
)

png(filename = 'plot/LASSO_OS.png')
plot(cvfitOS)
dev.off()

# Get optimal lambda
best.lambda <- cvfitOS$lambda.min
best.lambda # 0.01435765

Coefficients <- coef(cvfitOS, s = best.lambda)
Coefficients2 <- Coefficients[which(Coefficients != 0),]
Coefficients2
#STAGE         AGE        RACE         ARM      EMPLOY     LESION1        ECOG       CHEMO     SURGERY       RADIO 
#-0.03034113 -0.02391880  0.03033769  0.21959570  0.08957765  0.09147601  0.37994750  1.54833355 -0.15657734  0.38861858 
#HER2      HEIGHT      WEIGHT 
#0.21348912  0.06359369 -0.04421026 

# Cox regression with LASSO selected covariates
res.cox <- coxph(Surv(DTHDY, DTH) ~ STAGE+AGE+RACE+ARM+EMPLOY+LESION1+ECOG+CHEMO+SURGERY+RADIO
                 +HER2+HEIGHT+WEIGHT, data=OS_scaled) #use no filter DTH=0?
summary(res.cox)

####============================================================================####
#### LASSO for PFS

# Generate PFS dataset
PFS<-fatigue_baseline %>%
  filter(PFSDY != -999 & PFS != -999) %>%
  select(UID, PFSDY, PFS, STAGE:WEIGHT) %>%
  distinct()

# 2. rescale continuous variables
set.seed(1234)
PFS_scaled <- PFS %>%
  # Convert all character to factor
  mutate(across(where(is.character), as.factor)) %>%
  # Scale numeric columns except DTHDY and DTH
  mutate(across(where(is.numeric) & !any_of(c("PFSDY", "PFS")), scale)) %>%
  # Remove rows with NA
  na.omit()

# corrlation plot(continuous)
con.var <- PFS %>% select(where(is.numeric)) %>% select(-PFSDY, -PFS)
pdf("plot/Con_covariate_correlation_fatigue_PFS.pdf", width=12, height=12)
print(ggpairs(con.var))
dev.off()

# 3. cross-validation LASSO for covarites on PFS
set.seed(1234)

# filter DTHDY <= 0 
PFS_scaled_filtered <- PFS_scaled %>% filter(PFSDY > 0)

# Delete UID、DTHDY、DTH
x_data <- PFS_scaled_filtered %>% 
  select(-c(UID, PFSDY, PFS)) %>% 
  data.matrix()

y_data <- Surv(PFS_scaled_filtered$PFSDY, PFS_scaled_filtered$PFS)

cvfitPFS <- cv.glmnet(x = x_data,
                      y = y_data,
                      family='cox', alpha=1, nfolds=5)

png(filename = 'plot/LASSO_PFS.png')
plot(cvfitPFS)
dev.off()

# Get optimal lambda
best.lambda <- cvfitPFS$lambda.min
best.lambda # 0.01294616

Coefficients <- coef(cvfitPFS, s = best.lambda)
Coefficients2 <- Coefficients[which(Coefficients != 0),]
Coefficients2
#STAGE         AGE        RACE         ARM       MENOS        ECOG       CHEMO     SURGERY       RADIO         ERS 
#-0.07170791 -0.08764567  0.05148065  0.14171493 -0.02313328  0.23304224  1.82635637 -0.02876714  0.23195497 -0.09336526 
#HER2      WEIGHT 
#0.10286378 -0.01046544 

# 4. Cox regression with LASSO selected covariates
res.cox.PFS <- coxph(Surv(PFSDY, PFS) ~ STAGE + AGE + RACE + ARM + MENOS + ECOG +
                       CHEMO + SURGERY + RADIO + ERS + HER2 + WEIGHT, data=PFS_scaled)
summary(res.cox.PFS)

################################################################################
## Summary of LASSO results
## OS: STAGE、AGE、RACE、ARM、EMPLOY、LESION1、ECOG、CHEMO、SURGERY、RADIO、HER2、HEIGHT、WEIGHT
## PFS: STAGE, AGE, RACE, ARM, MENOS, ECOG, CHEMO, SURGERY, RADIO, ERS, HER2, WEIGHT



















