#=================================================
# Plot and display results
# Author: Stephanie Wu
# Date created: 2024/09/30
# Date updated: 2024/09/30
#=================================================

library(baysc)
library(furniture)
library(survey)
library(knitr)
library(kableExtra)
library(tidyverse)

# Directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
# wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/"
data_dir <- "Application/Cleaned_Data/"  # Data directory
res_dir <- "Application/Results/"        # Results directory
code_dir <- "Application/Code/"  # Code directory

# Load functions
source(paste0(wd, code_dir, "app_functions.R"))

# Load DBH results
load(paste0(wd, res_dir, "cc_dbh_wolcan_results.RData"))
load(paste0(wd, res_dir, "cc_dbh_wolcan_weights.RData"))

#==================== Read in and prepare data =================================
# Read in PRCS data
prcs_cleaned <- read.csv(paste0(wd, data_dir, "prcs_cleaned.csv"))
# Restrict to ages 30 to 75
prcs_cleaned <- prcs_cleaned %>% filter(Age >= 30 & Age < 76)

# Read in PROSPECT data
prospect_cleaned <- read.csv(paste0(wd, data_dir, "prospect_cleaned.csv"))
# Only the dietary behavior variables
prospect_dbh <- prospect_cleaned %>% 
  select(studyid, purchase_person:eat_vegetarian) 
# Drop NA: from 1690 to 1544 (146 removed)
prospect_dbh_drop_na <- prospect_dbh %>% drop_na()
# Restrict to DBH complete cases
restrict_data <- prospect_cleaned %>% 
  filter(studyid %in% prospect_dbh_drop_na$studyid)

# # Plot missingness pattern
# sort(apply(prospect_cleaned, 2, function(x) sum(is.na(x))), decreasing = TRUE)
# test <- md.pattern(prospect_cleaned, rotate.names = TRUE)

item_labels <- colnames(prospect_dbh)[-1]
class_title <- "Dietary Behavior Pattern"
categ_title <- "Risk Level"
categ_labels <- c("Low", "Med", "High")
cov_df <- restrict_data %>%
  select(Sex, Educ, Age, Inc_hh, Urban, Ethnicity, Physical_activity, 
         Smoking_status, Drinking_status, Food_security, WIC_SNAP, Social_support, 
         Perceived_stress, Depression, Anxiety) %>%
  mutate_at(c("Educ", "Inc_hh", "Physical_activity", "Smoking_status", 
              "Drinking_status"), as.factor)


#================ Results plots for the estimated weights ======================
# Median weight
summary(est_weights$wts)

# Create dataframe with weights and adding label names for factors
restrict_weights <- restrict_data %>%
  mutate(weight = est_weights$wts) %>%
  mutate(Sex = factor(Sex, levels = c(0,1), labels = c("Male", "Female")),
         Educ = factor(Educ, levels = c(1,2,3,4), 
                       labels = c("<HS", "HS", "College", "Graduate")),
         Inc_hh = factor(Inc_hh, levels = c(1,2,3),
                         labels = c("0-10k", "10-20k", ">20k")),
         Smoking_status = factor(Smoking_status, levels = c(0,1,2),
                                 labels = c("Never", "Former", "Current")),
         Drinking_status = factor(Drinking_status, levels = c(0,1,2),
                                  labels = c("Never", "Former", "Current")),
         Physical_activity = factor(Physical_activity, levels = c(0,1,2),
                                    labels = c("Sedentary", "Light", "Moderate")),
         Urban = factor(Urban, levels = c(0,1), labels = c("Rural", "Urban")),
         Ethnicity = factor(Ethnicity, levels = c(0,1), 
                            labels = c("Puerto Rican", "Other")),
         Food_security = factor(Food_security, levels = c(0,1),
                                labels = c("Insecure", "Secure")),
         WIC_SNAP = factor(WIC_SNAP, levels=c(0,1), labels = c("No", "Yes")),
         Depression = factor(Depression, levels = c(0,1), 
                             labels = c("None/Mild", "Moderate/Severe")),
         Anxiety = factor(Anxiety, levels = c(0,1),
                          labels = c("None/Mild", "Moderate/Severe")))
# Sample version with all weights set to 1
samp_data <- restrict_weights %>%
  mutate(weight = 1)


var_vec <- c("Age", "Sex", "Educ", "Inc_hh", "Ethnicity", "Urban",  
             "Smoking_status", "Drinking_status", "Physical_activity", 
             "Food_security", "WIC_SNAP", "Depression", "Anxiety", 
             "Social_support", "Perceived_stress")
# Proportions in sample
wts_tb_samp <- get_mean_props_wolcan(var_vec = var_vec, wts = samp_data$weight, 
                              data = samp_data, digits = 1)
# Estimated proportions in population, using weights
wts_tb_bart <- get_mean_props_wolcan(var_vec = var_vec, wts = restrict_weights$weight, 
                              data = restrict_weights, digits = 1)
# Restrict to those with higher-than-median weight
ind_high50 <- which(restrict_weights$weight >= median(restrict_weights$weight))
wts_tb_bart_high50 <- get_mean_props_wolcan(var_vec = var_vec, 
                                     wts = rep(1, length(ind_high50)),
                                     # wts = restrict_weights$weight[ind_high50], 
                                     data = restrict_weights[ind_high50, ], 
                                     digits = 1)
# Restrict to those with higher-than-75-percentile weight
ind_high75 <- which(restrict_weights$weight >= 
                      quantile(restrict_weights$weight, 0.75))
wts_tb_bart_high75 <- get_mean_props_wolcan(var_vec = var_vec, 
                                     wts = rep(1, length(ind_high75)),
                                     # wts = restrict_weights$weight[ind_high75], 
                                     data = restrict_weights[ind_high75, ], 
                                     digits = 1)
# Restrict to those with higher-than-90-percentile weight
ind_high90 <- which(restrict_weights$weight >= 
                      quantile(restrict_weights$weight, 0.9))
wts_tb_bart_high90 <- get_mean_props_wolcan(var_vec = var_vec, 
                                     wts = rep(1, length(ind_high90)),
                                     # wts = restrict_weights$weight[ind_high90], 
                                     data = restrict_weights[ind_high90, ], 
                                     digits = 1)
# Get "population" values using PRCS weights
# Strata and clustering are not needed here b/c no variance calculation
# Relabel factors
prcs_new <- prcs_cleaned %>%
  mutate(Sex = factor(Sex, levels = c(0,1), labels = c("Male", "Female")),
         Educ = factor(Educ, levels = c(1,2,3,4), 
                       labels = c("<HS", "HS", "College", "Graduate")),
         Inc_hh = factor(Inc_hh, levels = c(1,2,3),
                         labels = c("0-10k", "10-20k", ">20k")),
         Ethnicity = factor(Ethnicity, levels = c(0,1), 
                            labels = c("Puerto Rican", "Other")))
wts_tb_prcs <- get_mean_props_wolcan(var_vec = var_vec[1:5], wts = prcs_new$Weight, 
                              data = prcs_new, digits = 1)

# Combine all together into one table
wts_all <- wts_tb_samp %>% 
  dplyr::left_join(wts_tb_bart_high50, by = join_by(Variable, Level),
                   suffix = c("_samp", "_wts_q50")) %>%
  dplyr::left_join(wts_tb_bart_high90, by = join_by(Variable, Level)) %>%
  dplyr::left_join(wts_tb_bart, by = join_by(Variable, Level),
                   suffix = c("_wts_q90", "_pop"))
# wts_all <- wts_tb_samp %>% 
#   dplyr::left_join(wts_tb_bart, by = join_by(Variable, Level),
#                    suffix = c("_samp", "_bart_all")) %>%
#   dplyr::left_join(wts_tb_bart_high50, by = join_by(Variable, Level)) %>%
#   dplyr::left_join(wts_tb_bart_high75, by = join_by(Variable, Level),
#                    suffix = c("_bart_high50", "_bart_high75")) %>%
#   dplyr::left_join(wts_tb_bart_high90, by = join_by(Variable, Level)) %>%
#   dplyr::left_join(wts_tb_prcs, by = join_by(Variable, Level),
#                    suffix = c("_bart_high90", "_pop"))

knitr::kable(wts_all) %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

# colnames(wts_all) <- c("Variable", "Level", "Mean in\nSample", 
#                        "Mean in\nWeights >P50", "Mean in\nWeights >P90")
knitr::kable(wts_all, format = "latex", booktabs = TRUE)

### Variable importance plot
bart_mod <- est_weights$fit_pi_B
# Compute row percentages
percent_varcount = bart_mod$varcount / apply(bart_mod$varcount, 1, sum)
# Mean of row percentages
mean_percent_varcount = apply(percent_varcount, 2, mean)
# Quantiles of row percentages
quant_percent_varcount = apply(percent_varcount, 2, quantile, probs=c(.025,.975))
# Create plot
rgy = range(quant_percent_varcount)
p <- ncol(percent_varcount)
plot(c(1,p),rgy,type="n",xlab="variable",ylab="post mean, percent var use",axes=FALSE)
axis(1,at=1:p,labels=names(mean_percent_varcount),cex.lab=0.7,cex.axis=0.7)
axis(2,cex.lab=1.2,cex.axis=1.2)
lines(1:p,mean_percent_varcount,col="black",lty=4,pch=4,type="b",lwd=1.5)
for(i in 1:p) {
  lines(c(i,i),quant_percent_varcount[,i],col="blue",lty=3,lwd=1.0)
}

### Partial dependence plot: more useful for continuous variables
# pred_covs_B <- selection_covs <- c("Sex", "Educ", "Age", "Inc_hh", "Urban")
# dat_B <- restrict_data %>%  # Covariates for NPS
#   select(all_of(selection_covs)) %>%
#   mutate_at(c("Sex", "Educ", "Inc_hh", "Urban"), as.factor)
# dat_R <- prcs_cleaned %>%  # Covariates for reference
#   select(all_of(selection_covs)) %>%
#   mutate_at(c("Sex", "Educ", "Inc_hh", "Urban"), as.factor)
# n1 <- nrow(dat_B)
# n0 <- nrow(dat_R)
# samp_comb <- rbind(dat_B, dat_R)
# z <- rep(1:0, c(n1, n0))
# pdb1 <- dbarts::pdbart(x.train = samp_comb[, pred_covs_B, drop = FALSE],
#                      y.train = z, ntree = 50L, nskip = 100L, ndpost = num_post)


# # Table 1 draft version
# tb1 <- furniture::table1(restrict_weights, Sex, Age, Educ, Inc_hh, Urban, Ethnicity, 
#                          Smoking_status, Drinking_status, Physical_activity, 
#                          Food_security, WIC_SNAP, Social_support, Perceived_stress,
#                          Depression, Anxiety, splitby = ~high_weight, row_wise = FALSE, 
#                          na.rm = FALSE, total = TRUE, type = "condensed")
# tb1 <- as.data.frame(tb1$Table1)
# knitr::kable(tb1) %>% 
#   kableExtra::kable_classic(full_width = F, html_font = "Cambria")
  


#================ Figure for dietary behavior patterns =========================
baysc::plot_class_dist(res)
# Traceplot for pi
plot(res$estimates_adjust$pi_red[, 1], type = "l", ylim = c(0, 1))
hist(res$data_vars$sampling_wt, breaks = 30)

plot_pattern_probs_wolcan(res = res, item_labels = item_labels, 
                          class_title = class_title, 
                          categ_title = categ_title,
                          categ_labels = categ_labels)
res$estimates_adjust$pi_med

baysc::plot_pattern_probs(res = res, item_labels = item_labels, 
                             class_title = class_title, 
                             categ_title = categ_title,
                             categ_labels = categ_labels) + 
  ggplot2::facet_wrap(factor(Item, labels = item_labels) ~ ., nrow = 4)

### Unweighted results
load(paste0(wd, res_dir, "Sept_20/agePRCS_unwt_cc_dbh_wolca_results.RData"))
res_unwt <- res
res_unwt <- baysc::reorder_classes(res = res_unwt, new_order = c(3, 4, 2, 1))
plot_pattern_probs_wolcan(res = res_unwt, item_labels = item_labels, 
                          class_title = class_title, 
                          categ_title = categ_title,
                          categ_labels = categ_labels)
res_unwt$estimates$pi_med

baysc::plot_pattern_probs(res = res_unwt, item_labels = item_labels, 
                          class_title = class_title, 
                          categ_title = categ_title,
                          categ_labels = categ_labels)
baysc::plot_pattern_profiles(res = res_unwt, item_labels = item_labels, 
                          class_title = class_title, 
                          categ_title = categ_title,
                          categ_labels = categ_labels)

#================ Table of dietary behavior patterns by covariates =============

vars_cov_df <- restrict_data %>%
  select(Sex, Educ, Age, Inc_hh, Ethnicity, Urban, Physical_activity,
         Smoking_status, Drinking_status, Food_security, WIC_SNAP, Social_support, 
         Perceived_stress, Depression, Anxiety) %>%
  mutate(Educ = factor(Educ, levels = c(1,2,3,4), 
                       labels = c("<HS", "HS", "College", "Graduate")),
         Inc_hh = factor(Inc_hh, levels = c(1,2,3),
                         labels = c("0-10k", "10-20k", ">20k")),
         Smoking_status = factor(Smoking_status, levels = c(0,1,2),
                                 labels = c("Never", "Former", "Current")),
         Drinking_status = factor(Drinking_status, levels = c(0,1,2),
                                 labels = c("Never", "Former", "Current")),
         Physical_activity = factor(Physical_activity, levels = c(0,1,2),
                                    labels = c("Sedentary", "Light", "Moderate")))
dbh_vars_dist <- vars_across_class_wolcan(c_all = as.factor(res$estimates_adjust$c_all),
                                   cov_df = vars_cov_df, 
                                   sampling_wt = est_weights$wts, res = res,
                                   digits = 2)
colnames(dbh_vars_dist) <- c("Variable", "Level", "DBP1", "DBP2", "DBP3", "DBP4",
                             "Overall")
knitr::kable(dbh_vars_dist) %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")




  