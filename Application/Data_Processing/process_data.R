#=================================================
# Processing application data
# Author: Stephanie Wu
# Date created: 2024/08/07
# Date updated: 2024/08/07
#=================================================

rm(list = ls())

library(tidyverse)  # data wrangling
library(readxl)     # read excel files

# Directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
# wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/"
data_dir <- "Application/Data_Processing/"  # Raw data directory
res_dir <- "Application/"                  # Cleaned data directory
code_dir <- "Application/Data_Processing/"  # Code directory

# Read in PRCS person-level data
prcs_p <- read.csv(paste0(wd, data_dir, "psam_p72.csv"))
dim(prcs_p)
# Read in PRCS household-level data
prcs_h <- read.csv(paste0(wd, data_dir, "psam_h72.csv"))
dim(prcs_h)
# Read in metropolitan areas
metro <- read_xlsx(paste0(wd, data_dir, "MSA2013_PUMA2020_crosswalk.xlsx"))
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
# Create rural variable
prcs_metro <- prcs_metro %>%
  mutate(Rural = ifelse(is.na(`MSA Code`), 1, 0))

### Obtain variables to be used to model selection
# Age, sex, education, individual income, weights
prcs_vars <- prcs_p %>%
  select(SERIALNO, SPORDER, AGEP, SEX, SCHL, PINCP, PWGTP) %>%
  mutate(
    Sex = factor(SEX, levels = c(1, 2), labels = c("M", "F")),
    Educ = case_when(
      SCHL %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) ~ 0, # "<=8th grade", 
      SCHL %in% c(12, 13, 14, 15, 16, 17) ~ 1, # "9th grade to GED", 
      SCHL %in% c(18, 19, 20, 21) ~ 2, #"some college or bachelors",
      SCHL %in% c(22, 23, 24) ~ 3, #"graduate degree", 
      .default = SCHL)) %>%
  rename(
    Age = AGEP,
    Income_indiv = PINCP,
    Weight = PWGTP
  )

### Add in annual household income
# Recode annual household income
prcs_h <- prcs_h %>%
  mutate(Inc_hh = case_when(
    HINCP <= 10000 ~ 1,  # 0-10k
    HINCP > 10000 & HINCP <= 20000 ~ 2, # 10-20k
    HINCP > 20000 ~ 3, # >20k
    .default = HINCP
  ))
prcs_vars <- prcs_vars %>%
  left_join(prcs_h %>% select(SERIALNO, Inc_hh), 
            by = join_by(SERIALNO))
### Add in rural variable
prcs_vars <- prcs_vars %>% 
  left_join(prcs_metro %>% select(SERIALNO, Rural),
            by = join_by(SERIALNO))

# Number of individuals with complete data on selected variables
prcs_comp <- prcs_vars %>% drop_na()
# Reorder variables 
prcs_comp <- prcs_comp %>% select(SERIALNO, SPORDER, Weight, Age, Sex, Educ, 
                                  Inc_hh, Urban, Income_indiv)

### Save PRCS sample data
#write.csv(prcs_comp, file = paste0(wd, res_dir, "prcs.csv"))
