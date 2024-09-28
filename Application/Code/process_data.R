#=================================================
# Processing PRCS and PROSPECT application data
# Author: Stephanie Wu
# Date created: 2024/08/07
# Date updated: 2024/08/07
#=================================================

rm(list = ls())

library(tidyverse)  # data wrangling
library(readxl)     # read excel files
library(mice)       # missingness pattern

# Directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
# wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/"
data_dir <- "Application/Raw_Data/"  # Raw data directory
res_dir <- "Application/Cleaned_Data/"                  # Cleaned data directory
code_dir <- "Application/Code/"  # Code directory

#============ Clean PRCS probability sample data ===============================

# Read in PRCS person-level data
prcs_p <- read.csv(paste0(wd, data_dir, "prcs_psam_p72.csv"))
dim(prcs_p)
# Read in PRCS household-level data
prcs_h <- read.csv(paste0(wd, data_dir, "prcs_psam_h72.csv"))
dim(prcs_h)
# Read in metropolitan areas
metro <- read_xlsx(paste0(wd, data_dir, "prcs_MSA2013_PUMA2020_crosswalk.xlsx"))
# Restrict to PR and convert metro PUMA code to numeric
metro <- metro %>% 
  filter(`State Name` == "Puerto Rico") %>%
  mutate(PUMA20 = as.numeric(`PUMA Code`))


### Merge in metropolitan area information
# Check no missing 2020 PUMAs other than group quarters (-9)
setdiff(unique(prcs_h$PUMA20), unique(metro$PUMA20))
# Merge in metro area info
prcs_metro <- prcs_h %>%
  left_join(metro, by = join_by(PUMA20), relationship = "many-to-one",
            multiple = "first")  # set this to match PUMA 700 with MSA 25020
# How many PUMAs are MSAs: around 22%
mean(!is.na(prcs_metro$`MSA Code`))
# Create urban variable: 0 = rural, 1 = urban
prcs_metro <- prcs_metro %>%
  mutate(Urban = factor(ifelse(is.na(`MSA Code`), 0, 1), levels = c(1, 0)))

### Obtain variables to be used to model selection
# Age, sex, education, individual income, weights
prcs_vars <- prcs_p %>%
  select(SERIALNO, SPORDER, AGEP, SEX, SCHL, PINCP, PWGTP) %>%
  mutate(
    Sex = factor(SEX - 1, levels = c(0, 1)),  # M, F 
    Educ = factor(case_when(
      SCHL %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) ~ 1, # "<=8th grade", 
      SCHL %in% c(12, 13, 14, 15, 16, 17) ~ 2, # "9th grade to GED", 
      SCHL %in% c(18, 19, 20, 21) ~ 3, #"some college or bachelors",
      SCHL %in% c(22, 23, 24) ~ 4, #"graduate degree", 
      .default = SCHL))) %>%
  rename(
    Age = AGEP,
    Income_indiv = PINCP,
    Weight = PWGTP
  )

### Add in annual household income
# Recode annual household income
prcs_h <- prcs_h %>%
  mutate(Inc_hh = factor(case_when(
    HINCP <= 10000 ~ 1,  # 0-10k
    HINCP > 10000 & HINCP <= 20000 ~ 2, # 10-20k
    HINCP > 20000 ~ 3, # >20k
    .default = HINCP
  )))
prcs_vars <- prcs_vars %>%
  left_join(prcs_h %>% select(SERIALNO, Inc_hh), 
            by = join_by(SERIALNO))
### Add in urban variable
prcs_vars <- prcs_vars %>% 
  left_join(prcs_metro %>% select(SERIALNO, Urban),
            by = join_by(SERIALNO))
# Reorder variables 
prcs_cleaned <- prcs_vars %>% select(SERIALNO, SPORDER, Weight, Age, Sex, Educ, 
                                     Inc_hh, Urban, Income_indiv)
# Drop na
prcs_drop_na <- prcs_cleaned %>% drop_na()

### Save PRCS sample data
#write.csv(prcs_cleaned, file = paste0(wd, res_dir, "prcs_cleaned.csv"), 
#row.names = FALSE)

#============ Clean PROSPECT non-probability sample data =======================
# Read in prospect data
prospect_raw <- read.csv(paste0(wd, data_dir, "prospect_raw.csv"))
summary(prospect_raw)

# Which variables have NAs
count_nas <- sort(apply(prospect_raw, 2, function(x) sum(is.na(x))), 
                  decreasing = TRUE)
count_nas <- as.data.frame(count_nas)

###======== Obtain the dietary behavior variables
# Select variables
dbh_all <- prospect_raw %>%
  select(studyid, dbhi1:dbhv4)
# Recode variables
dbh <- dbh_all %>% mutate(
  # IDs
  studyid = studyid,
  # Purchasing and cooking roles
  purchase_person = case_match(
    dbhi1,
     1 ~ 1, # myself
     c(2:7, 9) ~ 2, # partner, family member, or other people
     8 ~ 3, # usually don't eat at house
     .default = NA), # refused or unknown
  purchase_location = case_match(
    dbhi2,
     1 ~ 1, # supermarket chain
     3 ~ 2, # discount stores like walmart and costco
     c(2, 4) ~ 3, # convenience store in neighborhood, or other
     .default = NA), # refused or unknown
  purchase_freq = case_match(
    dbhi3,
     c(1, 2) ~ 1, # once a week or more
     3 ~ 2, # once every two weeks
     4 ~ 3, # once a month or less
     .default = NA), # refused or unknown
  cook_person = case_match(
    dbhi4, 
     1 ~ 1, # myself
     c(2:7, 9) ~ 2, # partner, family member, or other people
     8 ~ 3, # usually don't eat at house
     .default = NA), # refused or unknown                           
  # Eating frequency and timing
  breakfast_freq = case_match(
    dbhii1,
     c(3, 4) ~ 1, # many times
     2 ~ 2, # sometimes
     1 ~ 3, # rarely or never
     .default = NA), # refused or unknown
  lunch_freq = case_match(
    dbhii3,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  dinner_freq = case_match(
    dbhii5,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  snack_freq = case_match(
    dbhii7,
    1 ~ 1, # rarely or never
    2 ~ 2, # sometimes
    c(3, 4) ~ 3, # many times
    .default = NA), # refused or unknown
  # Eating out
  fast_food_freq = case_match(
    dbhiii2a,
    1 ~ 1, # rarely or never
    2 ~ 2, # sometimes
    c(3, 4) ~ 3, # many times
    .default = NA), # refused or unknown
  restaurant_freq = case_match(
    dbhiii2c,
    1 ~ 1, # rarely or never
    2 ~ 2, # sometimes
    c(3, 4) ~ 3, # many times
    .default = NA), # refused or unknown
  take_out_freq = case_match(
    dbhiii2d,
    1 ~ 1, # rarely or never
    2 ~ 2, # sometimes
    c(3, 4) ~ 3, # many times
    .default = NA), # refused or unknown
  food_truck_freq = case_match(
    dbhiii2f,
    1 ~ 1, # rarely or never
    2 ~ 2, # sometimes
    c(3, 4) ~ 3, # many times
    .default = NA), # refused or unknown
  cafe_bakery_freq = case_match(
    dbhiii2e,
    1 ~ 1, # rarely or never
    2 ~ 2, # sometimes
    c(3, 4) ~ 3, # many times
    .default = NA), # refused or unknown
  party_freq = case_match(
    dbhiii8,
    1 ~ 1, # rarely or never
    2 ~ 2, # sometimes
    c(3, 4) ~ 3, # many times
    .default = NA), # refused or unknown
  # Eating practices
  meal_alone_freq = case_match(
    dbhiv3,
    1 ~ 1, # rarely or never
    2 ~ 2, # sometimes
    c(3, 4) ~ 3, # many times
    .default = NA), # refused or unknown
  meal_tv_freq = case_match(
    dbhiv5,
    1 ~ 1, # rarely or never
    2 ~ 2, # sometimes
    c(3, 4) ~ 3, # many times
    .default = NA), # refused or unknown
  meal_quality = case_match(
    dbhiv9,
    c(1, 2) ~ 1, # excellent or very good
    3 ~ 2, # good
    c(4, 5) ~ 3, # fair or poor
    .default = NA), # refused or unknown
  meal_healthy = case_match(
    dbhiv10,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  meal_accult = case_match(
    dbhiv11,
    c(1, 2) ~ 1, # mostly puerto rican meals
    3 ~ 2, # same amount of puerto rican and american meals
    c(4, 5) ~ 3, # mostly american meals
    .default = NA), # refused or unknown
  # Food restriction practices
  control_salt = case_match(
    dbhiv14a,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  control_fat = case_match(
    dbhiv14b,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  control_carbs_sugars = case_match(
    dbhiv14c,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  control_portions = case_match(
    dbhiv14d,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  count_calories = case_match(
    dbhiv14e,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  # Weight and supplements
  dysfunctional_weight_loss = case_when(  
    (dbhiv16___1 == 1) | (dbhiv16___2 == 1) | (dbhiv16___3 == 1) | 
      (dbhiv16___4 == 1) ~ 2, # yes
    dbhiv16___5 %in% c(0, 1) ~ 1, # no (catches indivs who have 0 for all)
    .default = NA), # refused or unknown
  vitamins = case_match(
    dbhiv17,
    c(1, 2) ~ 1, # once or twice a day
    3 ~ 2, # 2-3 days per week
    0 ~ 3, # no
    .default = NA), # refused or unknown
  probiotics = case_match(
    dbhiv18,
    c(1, 2) ~ 1, # once or twice a day
    c(3, 4) ~ 2, # 3 days per week or less
    0 ~ 3, # no
    .default = NA), # refused or unknown
  check_weight = case_match(
    dbhiv19,
    1 ~ 1, # everyday
    c(2, 3) ~ 2, # 2-4 times per month
    c(4, 5) ~ 3, # once per month or less
    .default = NA), # refused or unknown
  # Nutrition awareness and knowledge
  food_guide = case_match(
    dbhv1,
    1 ~ 1, # yes
    c(0, 2) ~ 2, # no or doesn't know what it is
    .default = NA), # refused or unknown
  nut_panel_freq = case_when(
    dbhv3 == 1 ~ 1, # every day or nearly everyday
    dbhv3 == 2 ~ 2, # a few times per week 
    (dbhv3 %in% c(3, 4)) | (dbhv2 == 0) ~ 3, # a few times per month or less
    .default = NA), # refused or unknown
  nut_info_freq = case_match(
    dbhv4,
    c(0, 1) ~ 1, # every day or nearly everyday
    2 ~ 2, # a few times per week 
    c(3, 4) ~ 3, # a few times per month or less
    .default = NA), # refused or unknown
  # Sustainability
  local_food_freq = case_match(
    dbhiv12,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  buy_organic_freq = case_match(
    dbhiv14f,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  water_type = case_match(
    dbhiv13,
    1 ~ 1, # right from the faucet without filtering
    2 ~ 2, # filtered or boiled water from the faucet
    c(3, 4) ~ 3, # bottled water or other
    .default = NA), # refused or unknown
  eat_vegetarian = case_match(
    dbhiv14g,
    c(3, 4) ~ 1, # many times
    2 ~ 2, # sometimes
    1 ~ 3, # rarely or never
    .default = NA), # refused or unknown
  .keep = "none"
)
# Convert NA to level 4
dbh_4_levels <- as.data.frame(apply(dbh, 2, function(x) ifelse(is.na(x), 4, x)))

# Are all the NAs the same individuals for the dbh variables?
# dbh: n=1700
# dbh_drop_na: n=1550 (150 dropped)
dbh_drop_na <- dbh %>%
  drop_na()
# count_nas
dbh_count_nas <- data.frame(num_nas = sort(apply(dbh, 2, function(x) 
  sum(is.na(x))), decreasing = TRUE))
# All individuals with NAs 
dbh_na <- dbh %>%
  filter(!(studyid %in% dbh_drop_na$studyid))
# Plot missingness pattern
test <- md.pattern(dbh_na %>% select(-studyid), rotate.names = TRUE)



###======== Obtain the outcome variables
# Select variables
outcomes_raw <- prospect_raw %>%
  select(studyid, SR_HTN, SR_TXHTN, SR_DIABETES, SR_TXDIABETES, SR_HCHOLESTEROL,
         SR_TXHCHOLESTEROL, bpsys_avg, bpdias_avg, FBS, HBA1C, CHOLESTEROL, HDL, 
         LDL)
# Recode variables
outcomes <- outcomes_raw %>% 
  mutate(
    hypertension = case_when(
      SR_HTN == 1 | SR_TXHTN == 1 | bpsys_avg >= 130 | bpdias_avg >= 80 ~ 1,
      is.na(SR_HTN) & is.na(SR_TXHTN) & is.na(bpsys_avg) & is.na(bpdias_avg) ~ NA,
      .default = 0),
    htn_aware = case_when(
      SR_HTN == 1 | SR_TXHTN == 1 ~ 1,
      is.na(SR_HTN) ~ NA,
      .default = 0),
    htn_treated = case_when(
      SR_TXHTN == 1 ~ 1,
      is.na(SR_HTN) ~ NA, 
      .default = 0),
    htn_categ = case_when(
      htn_treated == 1 ~ 4,  # treated
      htn_aware == 1 & htn_treated == 0 ~ 3,  # aware but untreated
      hypertension == 1 & htn_aware == 0 ~ 2,  # at risk and unaware
      hypertension == 0 ~ 1,  # not at risk
      .default = NA),
    diabetes = case_when(
      SR_DIABETES == 1 | SR_TXDIABETES == 1 | FBS >= 126 | HBA1C >= 6.5 ~ 1,
      is.na(SR_DIABETES) & is.na(SR_TXDIABETES) & is.na(FBS) & is.na(HBA1C) ~ NA,
      .default = 0),
    t2d_aware = case_when(
      SR_DIABETES == 1 | SR_TXDIABETES == 1 ~ 1,
      is.na(SR_DIABETES) ~ NA,
      .default = 0),
    t2d_treated = case_when(
      SR_TXDIABETES == 1 ~ 1, 
      is.na(SR_DIABETES) ~ NA, 
      .default = 0),
    t2d_categ = case_when(
      t2d_treated == 1 ~ 4,  # treated
      t2d_aware == 1 & t2d_treated == 0 ~ 3,  # aware but untreated
      diabetes == 1 & t2d_aware == 0 ~ 2,  # at risk and unaware
      diabetes == 0 ~ 1,  # not at risk
      .default = NA),
    high_cholesterol = case_when(
      SR_HCHOLESTEROL == 1 | SR_TXHCHOLESTEROL == 1 | CHOLESTEROL > 200 ~ 1,
      is.na(SR_HCHOLESTEROL) & is.na(SR_TXHCHOLESTEROL) & is.na(CHOLESTEROL) ~ NA,
      .default = 0),
    chol_aware = case_when(
      SR_HCHOLESTEROL == 1 | SR_TXHCHOLESTEROL == 1 ~ 1,
      is.na(SR_HCHOLESTEROL) ~ NA,
      .default = 0),
    chol_treated = case_when(
      SR_TXHCHOLESTEROL == 1 ~ 1, 
      is.na(SR_HCHOLESTEROL) ~ NA,
      .default = 0),
    chol_categ = case_when(
      chol_treated == 1 ~ 4,  # treated
      chol_aware == 1 & chol_treated == 0 ~ 3,  # aware but untreated
      high_cholesterol == 1 & chol_aware == 0 ~ 2,  # at risk and unaware
      high_cholesterol == 0 ~ 1,  # not at risk
      .default = NA),
    high_ldl = ifelse(LDL >= 100, 1, 0),
    systolic = bpsys_avg, 
    diastolic = bpdias_avg, 
    fasting_blood_sugar = FBS, 
    HBA1C = HBA1C, 
    cholesterol = CHOLESTEROL, 
    HDL = HDL, 
    LDL = LDL, .keep = "unused")
# Drop na
outcomes_drop_na <- outcomes %>% drop_na
data.frame(num_nas = sort(apply(outcomes_raw, 2, function(x) 
  sum(is.na(x))), decreasing = TRUE))
# Life-threatening hypertensive crisis is >300 for sys and >150 for dias
boxplot(outcomes_raw$bpsys_avg) # okay
boxplot(outcomes_raw$bpdias_avg) # okay
# Diabetic ketoacidosis is >600
boxplot(outcomes_raw$FBS) # okay
# Extreme hyperglycemia highest ever is 25%
boxplot(outcomes_raw$HBA1C) # okay
# Homozygous familial hypercholesterolemia highest ever is 2000
boxplot(outcomes_raw$CHOLESTEROL) # okay
# highest ever is 1000
boxplot(outcomes_raw$LDL) # okay
# highest ever is 150-180
boxplot(outcomes_raw$HDL) # okay

###======== Obtain the selection variables
# Select variables
selection_raw <- prospect_raw %>%
  select(studyid, age, female, rural, educ4cat, income3cat)
# Recode variables
selection <- selection_raw %>% 
  mutate(
    Sex = factor(female, levels = c(0, 1, 2)),  # M, F, O
    Educ = factor(educ4cat, levels = c(1, 2, 3, 4)),
    Age = age,
    Inc_hh = factor(income3cat, levels = c(1, 2, 3)),
    Urban = factor(case_match(rural, 
                                 0 ~ 1,  # urban
                                 c(1, 2) ~ 0,  # rural or neither or in-between
                                 .default = NA), levels = c(1, 0)), 
    .keep = "unused")
# Drop na
selection_drop_na <- selection %>% drop_na()
# Count na's per variable
count_nas <- sort(apply(selection, 2, function(x) sum(is.na(x))), 
                  decreasing = TRUE)


###======== Obtain the additional sociodemographic covariates
# Select variables
demog_raw <- prospect_raw %>%
  select(studyid, ethnicity, smoking_status:GAD_CAT2)
# Recode variables
demog <- demog_raw %>% 
  mutate(
    Ethnicity = factor(ethnicity, levels = c(0, 1)), # almost all Puerto Rican
    Smoking_status = factor(smoking_status, levels = c(0, 1, 2)),  # never, former, current
    Drinking_status = factor(drinking_status, levels = c(0, 1, 2)),
    Physical_activity = factor(PA_CAT, levels = c(0, 1, 2)),  # sedentary, light, moderate/vigorous
    Food_security = factor(fs_status_2cat, levels = c(1, 0)),  # secure, insecure
    WIC_SNAP = factor(WICSNAP, levels = c(0, 1)),  # no, yes (either)
    Social_support = SS_SCORE,
    Perceived_stress = PSS_SCORE,
    # Social_support = factor(case_when(  # CHECK THIS!!!!!
    #   SS_SCORE <= 18 ~ 0, # low support
    #   SS_SCORE >= 19 ~ 1, # high support (reference)
    # ), levels = c(1, 0)),  
    # Perceived_stress = factor(case_when(  # CHECK THIS!!!
    #   PSS_SCORE <= 28 ~ 0, # low stress (reference)
    #   PSS_SCORE >= 29 ~ 1, # high stress
    # ), levels = c(0, 1)),
    Depression = factor(CESD_GE_16, levels = c(0, 1)),  # no, yes (likely)
    Anxiety = factor(GAD_CAT2, levels = c(0, 1)),  # no/mild, yes/moderate/severe
    .keep = "unused")
# Drop na
demog_drop_na <- demog %>% drop_na()


###======== Combine all cleaned PROSPECT variables together
prospect_cleaned <- dbh %>%
  full_join(outcomes, by = join_by(studyid)) %>%
  full_join(selection, by = join_by(studyid)) %>%
  full_join(demog, by = join_by(studyid))
# Save data
#write.csv(prospect_cleaned, file = paste0(wd, res_dir, "prospect_cleaned.csv"),
#          row.names = FALSE)


