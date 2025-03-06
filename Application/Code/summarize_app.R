#=================================================
# Plot and display results
# Author: Stephanie Wu
# Date created: 2024/09/30
# Date updated: 2024/11/15
#=================================================

library(baysc)
library(furniture)
library(survey)
library(knitr)
library(kableExtra)
library(tidyverse)
library(ggpubr)

# Directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
# wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/"
data_dir <- "Application/Cleaned_Data/"  # Data directory
res_dir <- "Application/Results/"        # Results directory
code_dir <- "Application/Code/"  # Code directory

# Load functions
source(paste0(wd, code_dir, "app_functions.R"))

# Load DBH results
# load(paste0(wd, res_dir, "cc_dbh_wolcan_results.RData"))
# load(paste0(wd, res_dir, "cc_dbh_wolcan_weights.RData"))
load(paste0(wd, res_dir, "no_varadj_cc_dbh_wolcan_results.RData"))
load(paste0(wd, res_dir, "no_varadj_cc_dbh_wolcan_weights.RData"))

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
item_labels[25] <- "dysfunct_weight_loss"
item_title <- "Risk Level Probability"
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
  mutate(weight = est_weights$wts,
         Physical_activity = as.factor(ifelse(Physical_activity == 2, 1, 
                                              Physical_activity))) %>%
  mutate(Sex = factor(Sex, levels = c(0,1), labels = c("Male", "Female")),
         Educ = factor(Educ, levels = c(1,2,3,4), 
                       labels = c("<HS", "HS", "College", "Graduate")),
         Inc_hh = factor(Inc_hh, levels = c(1,2,3),
                         labels = c("0-10k", "10-20k", ">20k")),
         Smoking_status = factor(Smoking_status, levels = c(0,1,2),
                                 labels = c("Never", "Former", "Current")),
         Drinking_status = factor(Drinking_status, levels = c(0,1,2),
                                  labels = c("Never", "Former", "Current")),
         Physical_activity = factor(Physical_activity, levels = c(0,1),
                                    labels = c("Sedentary", "Active")),
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
# Display table
knitr::kable(wts_all) %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

# Latex version
knitr::kable(wts_all, format = "latex", booktabs = TRUE)


#================ Figure for dietary behavior patterns =========================
    # # Set variance adjustment to NULL if it isn't already
    # res$estimates_adjust <- NULL

# Plot distribution of class prevalences
baysc::plot_class_dist(res)
# ggsave(filename = paste0(wd, "Tables_Figures/", "class_dist_plot.png"),
#        width = 3000, height = 3500, units = "px", dpi = 700)

# Traceplot for pi
plot(res$estimates_adjust$pi_red[, 1], type = "l", ylim = c(0, 1))
plot(res$estimates$pi_red[, 1], type = "l", ylim = c(0, 1))

# Traceplot for theta
plot(res$estimates_adjust$theta_red[, 1, 1, 1], type = "l", ylim = c(0, 1))
plot(res$estimates$theta_red[, 1, 1, 1], type = "l", ylim = c(0, 1))

# # Plot pattern probabilities horizontally
# plot_pattern_probs_wolcan(res = res, item_labels = item_labels, 
#                           class_title = class_title, 
#                           categ_title = categ_title,
#                           categ_labels = categ_labels)

# Plot pattern probabilities vertically
baysc::plot_pattern_probs(res = res, item_labels = item_labels, 
                          class_title = class_title, y_title = item_title,
                          categ_title = categ_title,
                          categ_labels = categ_labels) + 
  ggplot2::facet_wrap(factor(Item, labels = item_labels) ~ ., nrow = 5) + 
  theme(strip.text = element_text(size = 11),
        strip.background = element_rect(fill = "aliceblue"))
# # Save pattern probabilities
# ggsave(filename = paste0(wd, "Tables_Figures/", "pattern_probs_plot.png"),
#        width = 8200, height = 6700, units = "px", dpi = 700)

# Plot pattern profiles (USE NEW VERSION)
plot_pattern_profiles(res = res, item_labels = item_labels, 
                      class_title = class_title, 
                      class_labels = c("Nutrition\nSensitive",
                                                 "Social\nEating",
                                                 "Out-of-Home\nEating",
                                                 "Nutrition\nInsensitive"),
                      y_title = item_title,
                      categ_title = categ_title,
                      categ_labels = categ_labels)
# Save pattern profiles
# ggsave(filename = paste0(wd, "Tables_Figures/", "pattern_profiles_plot.png"),
#        width = 3800, height = 5600, units = "px", dpi = 700)


#================ Table of dietary behavior patterns by covariates =============
# Create table showing distribution of dietary behavior patterns by 
# sociobehavioral covariates
vars_cov_df <- restrict_data %>%
  mutate(# Combine moderate/vigorous activity into light activity
    Physical_activity = as.factor(ifelse(Physical_activity == 2, 1, 
                                         Physical_activity))) %>%
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
         Physical_activity = factor(Physical_activity, levels = c(0,1),
                                    labels = c("Sedentary", "Active")))
dbh_vars_dist <- vars_across_class_wolcan(c_all = as.factor(res$estimates_adjust$c_all),
                                          cov_df = vars_cov_df, 
                                          sampling_wt = est_weights$wts, res = res,
                                          digits = 2)
colnames(dbh_vars_dist) <- c("Variable", "Level", "DBP1", "DBP2", "DBP3", "DBP4",
                             "Overall")
# Display table
knitr::kable(dbh_vars_dist) %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")
# Latex version
knitr::kable(dbh_vars_dist, format = "latex", booktabs = TRUE)


#=============== Create outcome regression table and figure ====================

### Load in full models with all covariates, no interactions, no subsetting
# List of model paths to load
mod_names_list <- list(paste0(wd, res_dir, "t2d_fulllogreg.RData"),
                       paste0(wd, res_dir, "htn_fulllogreg.RData"),
                       paste0(wd, res_dir, "chol_fulllogreg.RData"))
# List of loaded models
mod_res_list <- list()
for (i in 1:length(mod_names_list)) {
  load(mod_names_list[[i]])
  mod_res_list[[i]] <- wtd_logreg_res
}
full_no_subset_mod_res_list <- mod_res_list
# List of summarized models, exponentiated to OR scale
summ_mod_list <- lapply(mod_res_list, function(x) 
  summarize_parms(x, quant_lb = 0.025, quant_ub = 0.975, exponentiate = TRUE))
# Plot models
create_mod_plot(summ_mod_list) 


### Load in full models with all covariates, no interactions, subsetting to 
### exclude self-reported diagnosis 
# List of model paths to load
mod_names_list <- list(paste0(wd, res_dir, "t2d_unaware_fulllogreg.RData"),
                       paste0(wd, res_dir, "htn_unaware_fulllogreg.RData"),
                       paste0(wd, res_dir, "chol_unaware_fulllogreg.RData"))
# List of loaded models
mod_res_list <- list()
for (i in 1:length(mod_names_list)) {
  load(mod_names_list[[i]])
  mod_res_list[[i]] <- wtd_logreg_res
}
full_mod_res_list <- mod_res_list
# List of summarized models
summ_mod_list <- lapply(mod_res_list, function(x) 
  summarize_parms(x, quant_lb = 0.025, quant_ub = 0.975, exponentiate = TRUE))
# Plot models
create_mod_plot(summ_mod_list)



# ### Core, no interactions, no subset
# # List of model paths to load
# mod_names_list <- list(paste0(wd, res_dir, "t2d_corelogreg.RData"),
#                        paste0(wd, res_dir, "htn_corelogreg.RData"),
#                        paste0(wd, res_dir, "chol_corelogreg.RData"))
# # List of loaded models
# mod_res_list <- list()
# for (i in 1:length(mod_names_list)) {
#   load(mod_names_list[[i]])
#   mod_res_list[[i]] <- wtd_logreg_res
# }
# core_no_subset_mod_res_list <- mod_res_list
# # List of summarized models
# summ_mod_list <- lapply(mod_res_list, function(x) 
#   summarize_parms(x, quant_lb = 0.025, quant_ub = 0.975, exponentiate = TRUE))
# # Plot models
# create_mod_plot(summ_mod_list)
# 
# 
# ### Core, no interactions
# # List of model paths to load
# mod_names_list <- list(paste0(wd, res_dir, "t2d_unaware_corelogreg.RData"),
#                        paste0(wd, res_dir, "htn_unaware_corelogreg.RData"),
#                        paste0(wd, res_dir, "chol_unaware_corelogreg.RData"))
# # List of loaded models
# mod_res_list <- list()
# for (i in 1:length(mod_names_list)) {
#   load(mod_names_list[[i]])
#   mod_res_list[[i]] <- wtd_logreg_res
# }
# core_mod_res_list <- mod_res_list
# # List of summarized models
# summ_mod_list <- lapply(mod_res_list, function(x) 
#   summarize_parms(x, quant_lb = 0.025, quant_ub = 0.975, exponentiate = TRUE))
# # Plot models
# create_mod_plot(summ_mod_list)
# 
# 
# ### Marginal
# # List of model paths to load
# mod_names_list <- list(paste0(wd, res_dir, "t2d_unaware_marglogreg.RData"),
#                        paste0(wd, res_dir, "htn_unaware_marglogreg.RData"),
#                        paste0(wd, res_dir, "chol_unaware_marglogreg.RData"))
# # List of loaded models
# mod_res_list <- list()
# for (i in 1:length(mod_names_list)) {
#   load(mod_names_list[[i]])
#   mod_res_list[[i]] <- wtd_logreg_res
# }
# # List of summarized models
# summ_mod_list <- lapply(mod_res_list, function(x) 
#   summarize_parms(x, quant_lb = 0.025, quant_ub = 0.975, exponentiate = TRUE))
# # Plot models
# create_mod_plot(summ_mod_list)
# marg_mod_res_list <- mod_res_list



### Create facetted figure displaying the non-subsetted and subsetted models for 
### the full model with all covariates
summ_full_no_subset <- lapply(full_no_subset_mod_res_list, function(x) 
  summarize_parms(x, quant_lb = 0.025, quant_ub = 0.975, exponentiate = TRUE))
summ_full <- lapply(full_mod_res_list, function(x) 
  summarize_parms(x, quant_lb = 0.025, quant_ub = 0.975, exponentiate = TRUE))
create_mod_plot_multiple(
  summ_mod_list_multiple = list(summ_full_no_subset, summ_full),
  mod_labels = c("Not Subsetted", "Subsetted to Exclude \nSelf-Reported Diagnosis"),
  legend_labels = c("Social Eating", "Out-of-Home Eating",
                    "Nutrition-Insensitive")) + 
  theme(strip.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))  # add space between facets
# ggsave(filename = paste0(wd, "Tables_Figures/", "full_outcome_plot.png"), 
#        width = 5000, height = 4800, units = "px", dpi = 700)


### Create table displaying the non-subsetted and subsetted models for the full 
### model with all covariates
# Clean output for models
digits <- 2
df_full_no_subset <- lapply(summ_full_no_subset, function(x) {
  x %>% mutate(Est_CI = paste0(format(round(Mean, digits), digits = digits), " (", 
                               format(round(`2.5%`, digits), digits = digits), ", ",
                               format(round(`97.5%`, digits), digits = digits), ")"),
               Post_Prob = format(round(ifelse(`P(xi>0)` > 0.5, `P(xi>0)`,
                                               `P(xi<0)`), digits), digits = digits)) 
})
df_full <- lapply(summ_full, function(x) {
  x %>% mutate(Est_CI = paste0(format(round(Mean, digits), digits = digits), " (", 
                               format(round(`2.5%`, digits), digits = digits), ", ",
                               format(round(`97.5%`, digits), digits = digits), ")"),
               Post_Prob = format(round(ifelse(`P(xi>0)` > 0.5, `P(xi>0)`,
                                               `P(xi<0)`), digits), digits = digits)) 
})
df_both <- append(df_full_no_subset, df_full)
names(df_both) <- c("t2d_no_subset", "htn_no_subset", "chol_no_subset",
                    "t2d", "htn", "chol")
# Create output table
pattern_labels <- c("Nutrition-Sensitive (Ref)", "Social Eating",
                    "Out-of-Home Eating", "Nutrition-Insensitive")
tab_dbh_disease <- as.data.frame(matrix(NA, nrow = 10, ncol = 7))
tab_dbh_disease[, 1] <- c("Full Sample Model",
                          pattern_labels,
                          "Subsetted to Exclude Self-Reported Diagnosis",
                          pattern_labels)
colnames(tab_dbh_disease) <- c("Covariate", "Est_CI_t2d", "Post_prob_t2d", 
                               "Est_CI_htn", "Post_prob_htn", 
                               "Est_CI_chol", "Post_prob_chol")
# No subset results
tab_dbh_disease[2:5, c(2, 4, 6)] <- sapply(df_full_no_subset, function(x) 
  x[c(28, 1:3), "Est_CI"]) 
tab_dbh_disease[2:5, c(3, 5, 7)] <- sapply(df_full_no_subset, function(x) 
  x[c(28, 1:3), "Post_Prob"]) 
# Subset results
tab_dbh_disease[7:10, c(2, 4, 6)] <- sapply(df_full, function(x) 
  x[c(28, 1:3), "Est_CI"]) 
tab_dbh_disease[7:10, c(3, 5, 7)] <- sapply(df_full, function(x) 
  x[c(28, 1:3), "Post_Prob"]) 

# Display table
knitr::kable(tab_dbh_disease, format = "latex", booktabs = TRUE) %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")


# # With Intercept instead of b_Intercept
# tab_dbh_disease[2:5, c(2, 4, 6)] <- sapply(df_full_no_subset, function(x) 
#   x[c(26, 1:3), "Est_CI"]) 
# tab_dbh_disease[2:5, c(3, 5, 7)] <- sapply(df_full_no_subset, function(x) 
#   x[c(26, 1:3), "Post_Prob"]) 
# # Subset results
# tab_dbh_disease[7:10, c(2, 4, 6)] <- sapply(df_full, function(x) 
#   x[c(26, 1:3), "Est_CI"]) 
# tab_dbh_disease[7:10, c(3, 5, 7)] <- sapply(df_full, function(x) 
#   x[c(26, 1:3), "Post_Prob"]) 
# knitr::kable(tab_dbh_disease, format = "latex", booktabs = TRUE) %>% 
#   kableExtra::kable_classic(full_width = F, html_font = "Cambria")


### Create table displaying all covariate values for the non-subsetted and 
### subsetted models for the full model with all covariates
# Clean output for models
digits <- 2
df_full_no_subset <- lapply(summ_full_no_subset, function(x) {
  x %>% mutate(Est_CI = paste0(format(round(Mean, digits), digits = digits), " (", 
                               format(round(`2.5%`, digits), digits = digits), ", ",
                               format(round(`97.5%`, digits), digits = digits), ")"),
               Post_Prob = format(round(ifelse(`P(xi>0)` > 0.5, `P(xi>0)`,
                                               `P(xi<0)`), digits), digits = digits)) 
})
df_full <- lapply(summ_full, function(x) {
  x %>% mutate(Est_CI = paste0(format(round(Mean, digits), digits = digits), " (", 
                               format(round(`2.5%`, digits), digits = digits), ", ",
                               format(round(`97.5%`, digits), digits = digits), ")"),
               Post_Prob = format(round(ifelse(`P(xi>0)` > 0.5, `P(xi>0)`,
                                               `P(xi<0)`), digits), digits = digits)) 
})
df_both <- append(df_full_no_subset, df_full)
names(df_both) <- c("t2d_no_subset", "htn_no_subset", "chol_no_subset",
                    "t2d", "htn", "chol")
# Create output table
pattern_labels <- c("Nutrition-Sensitive (Ref)", "Social Eating",
                    "Out-of-Home Eating", "Nutrition-Insensitive")
tab_dbh_disease <- as.data.frame(matrix(NA, nrow = 10, ncol = 7))
tab_dbh_disease[, 1] <- c("Full Sample Model",
                          pattern_labels,
                          "Subsetted to Exclude Self-Reported Diagnosis",
                          pattern_labels)
colnames(tab_dbh_disease) <- c("Covariate", "Est_CI_t2d", "Post_prob_t2d", 
                               "Est_CI_htn", "Post_prob_htn", 
                               "Est_CI_chol", "Post_prob_chol")
# No subset results
tab_dbh_disease[2:5, c(2, 4, 6)] <- sapply(df_full_no_subset, function(x) 
  x[c(28, 1:3), "Est_CI"]) 
tab_dbh_disease[2:5, c(3, 5, 7)] <- sapply(df_full_no_subset, function(x) 
  x[c(28, 1:3), "Post_Prob"]) 
# Subset results
tab_dbh_disease[7:10, c(2, 4, 6)] <- sapply(df_full, function(x) 
  x[c(28, 1:3), "Est_CI"]) 
tab_dbh_disease[7:10, c(3, 5, 7)] <- sapply(df_full, function(x) 
  x[c(28, 1:3), "Post_Prob"]) 

tab_cov_full_t2d <- df_full_no_subset[[1]]
tab_cov_full_t2d <- tab_cov_full_t2d[c(28, 1:25), 
                                     c("Variable", "Est_CI", "Post_Prob")]
tab_cov_full_t2d$Variable[1] <- "Intercept"
rownames(tab_cov_full_t2d) <- NULL
# Display table
knitr::kable(tab_cov_full_t2d, format = "html", booktabs = TRUE) %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")

# Display table
knitr::kable(tab_dbh_disease, format = "latex", booktabs = TRUE) %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")


#================ UNWEIGHTED dietary behavior patterns =========================

reorder_classes <- function (res, new_order) {
  if (!inherits(res, c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting \n         from a call to one of these functions")
  } else if ((inherits(res, "wolca")) & !is.null(res$estimates_svyglm)) {
    warning(paste0("For WOLCA, reordering of classes should be done before ", 
                   "calling wolca_svyglm(). res$estimates_svyglm should be NULL ", 
                   "prior to running this function."))
  }
  ### CHANGED
  if (any(sort(unique(new_order)) != sort(unique(res$estimates$c_all)))) {
    stop("New ordering must have same number of classes and unique values as old ordering.")
  }
  res_new <- res
  K <- length(new_order)
  ### END CHANGED
  if (!is.null(res$estimates_adjust)) {
    ### CHANGED
    res_new$estimates_adjust$pi_red <- res$estimates_adjust$pi_red[, 
                                                                   new_order, drop = FALSE]
    res_new$estimates_adjust$theta_red <- res$estimates_adjust$theta_red[, 
                                                                         , new_order, , drop = FALSE]
    res_new$estimates_adjust$pi_med <- res$estimates_adjust$pi_med[new_order, drop = FALSE]
    res_new$estimates_adjust$theta_med <- res$estimates_adjust$theta_med[, 
                                                                         new_order, , drop = FALSE]
    
    for (i in 1:K) {
      res_new$estimates_adjust$c_all[res$estimates_adjust$c_all == 
                                       new_order[i]] <- i
    }
    ### END CHANGED
    if (is(res, "swolca")) {
      ### CHANGED
      res_new$estimates_adjust$xi_red <- res$estimates_adjust$xi_red[, 
                                                                     new_order, , drop = FALSE]
      res_new$estimates_adjust$xi_med <- res$estimates_adjust$xi_med[new_order, , drop = FALSE]
    }
  }
  else {
    res_new$estimates$pi_red <- res$estimates$pi_red[, new_order, drop = FALSE]
    res_new$estimates$theta_red <- res$estimates$theta_red[, 
                                                           , new_order, , drop = FALSE]
    res_new$estimates$pi_med <- res$estimates$pi_med[new_order, drop = FALSE]
    res_new$estimates$theta_med <- res$estimates$theta_med[, 
                                                           new_order, , drop = FALSE]
    for (i in 1:K) {
      res_new$estimates$c_all[res$estimates$c_all == new_order[i]] <- i
    }
    if (is(res, "swolca")) {
      res_new$estimates$xi_red <- res$estimates$xi_red[, 
                                                       new_order, , drop = FALSE]
      res_new$estimates$xi_med <- res$estimates$xi_med[new_order, , drop = FALSE]
      ### END CHANGED
    }
  }
  return(res_new)
}

plot_pattern_profiles <- function(res, item_labels = NULL, item_title = "Item",
                                  categ_labels = NULL, 
                                  categ_title = "Consumption Level",
                                  class_labels = NULL, 
                                  class_title = "Dietary Pattern", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain theta estimates
  if (!is.null(res$estimates_adjust)) {
    # Adjusted estimates
    est_item_probs <- res$estimates_adjust$theta_med
  } else {
    # Unadjusted estimates 
    est_item_probs <- res$estimates$theta_med
  }
  
  # Obtain theta modes as the category with highest probability for each item
  mode_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), which.max))
  
  # Define item, latent class, and category names if not specified
  if (is.null(item_labels)) {
    item_labels <- 1:res$data_vars$J
  } else if (length(item_labels) != res$data_vars$J) {
    stop(paste0("length of item_labels must equal the number of exposure items, J = ",
                res$data_vars$J))
  }
  K <- dim(mode_item_probs)[2]
  if (is.null(class_labels)) {
    class_labels <- 1:K
  } else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", K))
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  } else if (length(categ_labels) != res$data_vars$R) {
    stop(paste0("length of categ_labels must equal the number of exposure categories, R = ", 
                res$data_vars$R))
  }
  rownames(mode_item_probs) <- item_labels
  colnames(mode_item_probs) <- class_labels
  mode_item_probs$Item <- rownames(mode_item_probs)
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Class <- Item <- Level <- NULL
  
  # Create plot
  mode_plot <- mode_item_probs %>% 
    tidyr::gather("Class", "Level", -Item) %>%
    mutate(Class = factor(Class, levels = class_labels))  ## ADDED!
  mode_plot %>% ggplot2::ggplot(ggplot2::aes(x=Class, 
                                             y=factor(Item, 
                                                      levels = rev(item_labels)), 
                                             fill = factor(Level))) + 
    ggplot2::geom_tile(color="black", linewidth = 0.3) + 
    ggplot2::scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
                               name = categ_title, labels = categ_labels) +
    ggplot2::labs(x = class_title, y = item_title) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete() + 
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "top")
}


# Load unweighted DBH results
load(paste0(wd, res_dir, "agePRCS_unwt_cc_dbh_wolca_results.RData"))
res_unwt <- res
# Reorder classes to match weighted model ordering
res_unwt <- reorder_classes(res = res_unwt, new_order = c(3, 4, 2, 1))

# Plot distribution of class prevalences
baysc::plot_class_dist(res = res_unwt)
# ggsave(filename = paste0(wd, "Tables_Figures/", "unwt_class_dist_plot.png"), 
#        width = 3000, height = 3500, units = "px", dpi = 700)

# Traceplot for pi
plot(res_unwt$estimates$pi_red[, 1], type = "l", ylim = c(0, 1))

# Plot pattern probabilities vertically
baysc::plot_pattern_probs(res = res_unwt, item_labels = item_labels, 
                          class_title = class_title, y_title = item_title,
                          categ_title = categ_title,
                          categ_labels = categ_labels) + 
  ggplot2::facet_wrap(factor(Item, labels = item_labels) ~ ., nrow = 5) + 
  theme(strip.text = element_text(size = 11),
        strip.background = element_rect(fill = "aliceblue"))
# Save pattern probabilities
# ggsave(filename = paste0(wd, "Tables_Figures/", "unwt_pattern_probs_plot.png"),
#        width = 8200, height = 6700, units = "px", dpi = 700)

# res_unwt$estimates$c_all <- factor(res_unwt$estimates$c_all, 
#                                    levels = c(1, 2, 3, 4))

# Plot pattern profiles
plot_pattern_profiles(res = res_unwt, item_labels = item_labels, 
                      class_title = class_title, 
                      class_labels = c("Nutrition\nSensitive",
                                       "Social\nEating",
                                       "Out-of-Home\nEating",
                                       "Nutrition\nInsensitive"),
                      y_title = item_title,
                      categ_title = categ_title,
                      categ_labels = categ_labels)
# Save pattern profiles
# ggsave(filename = paste0(wd, "Tables_Figures/", "unwt_pattern_profiles_plot.png"),
#        width = 3800, height = 5600, units = "px", dpi = 700)


#=============== UNWEIGHTED outcome regression table and figure ================

### Load in full models with all covariates, no interactions, no subsetting
load(paste0(wd, res_dir, "unwt_t2d_htn_chol_fulllogreg.RData"))

full_no_subset_mod_res_list <- mod_res_list

# List of summarized models, exponentiated to OR scale
summ_full_no_subset <- lapply(mod_res_list, function(x) 
  summarize_parms_unwt(x, quant_lb = 0.025, quant_ub = 0.975, 
                       exponentiate = TRUE))
# Plot models
create_mod_plot(summ_full_no_subset) 


### Load in full models with all covariates, no interactions, subsetting to 
### exclude self-reported diagnosis 
load(paste0(wd, res_dir, "unwt_t2d_htn_chol_unaware_fulllogreg.RData"))

full_mod_res_list <- mod_res_list_unaware

# List of summarized models, exponentiated to OR scale
summ_full <- lapply(mod_res_list_unaware, function(x) 
  summarize_parms_unwt(x, quant_lb = 0.025, quant_ub = 0.975, 
                       exponentiate = TRUE))
# Plot models
create_mod_plot(summ_full) 



### Create facetted figure displaying the non-subsetted and subsetted models for 
### the full model with all covariates
create_mod_plot_multiple(
  summ_mod_list_multiple = list(summ_full_no_subset, summ_full),
  mod_labels = c("Not Subsetted", "Subsetted to Exclude \nSelf-Reported Diagnosis"),
  legend_labels = c("Social Eating", "Out-of-Home Eating",
                    "Nutrition-Insensitive")) + 
  theme(strip.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines"))  # add space between facets
# ggsave(filename = paste0(wd, "Tables_Figures/", "unwt_full_outcome_plot.png"), 
#        width = 5000, height = 4800, units = "px", dpi = 700)


### Create table displaying the non-subsetted and subsetted models for the full 
### model with all covariates
# Clean output for models
digits <- 2
df_full_no_subset <- lapply(summ_full_no_subset, function(x) {
  x %>% mutate(Est_CI = paste0(format(round(Mean, digits), digits = digits), " (", 
                               format(round(`2.5%`, digits), digits = digits), ", ",
                               format(round(`97.5%`, digits), digits = digits), ")"),
               `P-value` = format(round(`P-value`, digits), digits = digits)) 
})
df_full <- lapply(summ_full, function(x) {
  x %>% mutate(Est_CI = paste0(format(round(Mean, digits), digits = digits), " (", 
                               format(round(`2.5%`, digits), digits = digits), ", ",
                               format(round(`97.5%`, digits), digits = digits), ")"),
               `P-value` = format(round(`P-value`, digits), digits = digits)) 
})
df_both <- append(df_full_no_subset, df_full)
names(df_both) <- c("t2d_no_subset", "htn_no_subset", "chol_no_subset",
                    "t2d", "htn", "chol")
# Create output table
pattern_labels <- c("Nutrition-Sensitive (Ref)", "Social Eating",
                    "Out-of-Home Eating", "Nutrition-Insensitive")
tab_dbh_disease <- as.data.frame(matrix(NA, nrow = 10, ncol = 7))
tab_dbh_disease[, 1] <- c("Full Sample Model",
                          pattern_labels,
                          "Subsetted to Exclude Self-Reported Diagnosis",
                          pattern_labels)
colnames(tab_dbh_disease) <- c("Covariate", "Est_CI_t2d", "P_val_t2d", 
                               "Est_CI_htn", "P_val_htn", 
                               "Est_CI_chol", "P_val_chol")
# No subset results
tab_dbh_disease[2:5, c(2, 4, 6)] <- sapply(df_full_no_subset, function(x) 
  x[c(26, 1:3), "Est_CI"]) 
tab_dbh_disease[2:5, c(3, 5, 7)] <- sapply(df_full_no_subset, function(x) 
  x[c(26, 1:3), "P-value"]) 
# Subset results
tab_dbh_disease[7:10, c(2, 4, 6)] <- sapply(df_full, function(x) 
  x[c(26, 1:3), "Est_CI"]) 
tab_dbh_disease[7:10, c(3, 5, 7)] <- sapply(df_full, function(x) 
  x[c(26, 1:3), "P-value"]) 

# Display table
knitr::kable(tab_dbh_disease, format = "latex", booktabs = TRUE) %>% 
  kableExtra::kable_classic(full_width = F, html_font = "Cambria")



#======================== Correlation between weight draws =====================

# `get_cor_compare` compares the average pairwise correlation between different 
# ways of sampling from the posterior predictive distribution of the weights
# Inputs:
#   wts_post: nxM matrix of posterior weights, where n is the sample size of the 
# NPS and M is the number of MCMC samples that are preserved
#   D: Number of draws to keep from the posterior weight distribution when 
# propagating the uncertainty of the weights throughout the model
# Outputs: A list containing the following elements:
#   avg_cor_all: Average pairwise correlation using all samples from the 
# posterior predictive weight distribution
#   cor_avgs: Vector of the average pairwise correlations across 100 SRSs 
# of D draws
#   avg_cor_SRS: Average across `cor_avgs`
#   dist_cor_SRS: Summary of `cor_avgs`
#   avg_cor_quants: Average pairwise correlation for a quantile-based method to 
# obtain the D draws 
get_cor_compare <- function(wts_post, D) {
  set.seed(1)
  
  # Obtain the correlation distribution for all posterior samples
  cor_mat_all <- cor(wts_post)
  dim(cor_mat_all)
  cor_mat_all[1:5, 1:5]
  cor_unique_all <- cor_mat_all[upper.tri(cor_mat_all)]
  avg_cor_all <- round(mean(cor_unique_all), 3)
  print(paste0("Average correlation across all samples: ", avg_cor_all))
  
  # Obtain the mean correlation using SRS to select D draws
  # Try this 100 times
  cor_avgs <- numeric(100)
  for (i in 1:100) {
    ind_SRS <- sample(1:ncol(wts_post), size = D)
    wts_SRS <- wts_post[, ind_SRS]
    cor_mat_SRS <- cor(wts_SRS)
    cor_unique_SRS <- cor_mat_SRS[upper.tri(cor_mat_SRS)]
    cor_avgs[i] <- mean(cor_unique_SRS)
  }
  avg_cor_SRS <- round(mean(cor_avgs), 3)
  dist_cor_SRS <- summary(cor_avgs)
  print(paste0("Average correlation for SRS for 100 simulations: ", avg_cor_SRS))
  print("Distribution of correlations across 100 simulations:")
  print(dist_cor_SRS)
  
  
  # Obtain the mean correlation using stratified quantiles to select D draws
  # Get quantiles of medians of weights posterior
  col_meds <- c(apply(wts_post, 2, median))
  cutoffs <- (seq(1, D, length.out = D) - 0.5) / D
  med_quants <- stats::quantile(x = col_meds, probs = cutoffs)
  wts_quants <- matrix(NA, nrow = nrow(wts_post), ncol = D)
  for (d in 1:D) {
    # Draw from weights posterior closest to quantile
    draw <- which.min(abs(col_meds - med_quants[d]))
    wts_quants[, d] <- c(wts_post[, draw])
  }
  cor_mat_quants <- cor(wts_quants)
  cor_unique_quants <- cor_mat_quants[upper.tri(cor_mat_quants)]
  avg_cor_quants <- round(mean(cor_unique_quants), 3)
  print(paste0("Average correlation for quantiles: ", avg_cor_quants))
  
  return(list(avg_cor_all = avg_cor_all, cor_avgs = cor_avgs, 
              avg_cor_SRS = avg_cor_SRS, dist_cor_SRS = dist_cor_SRS, 
              avg_cor_quants = avg_cor_quants))
}

# Obtain correlations for D = 5, 10, 20, 100
d5 <- get_cor_compare(wts_post = est_weights$wts_post, D = 5)
d10 <- get_cor_compare(wts_post = est_weights$wts_post, D = 10)
d20 <- get_cor_compare(wts_post = est_weights$wts_post, D = 20)
d100 <- get_cor_compare(wts_post = est_weights$wts_post, D = 100)

# Create dataframe including all vectors for the SRS scenario
d_all <- data.frame(d5$cor_avgs, d10$cor_avgs, d20$cor_avgs, d100$cor_avgs)
colnames(d_all) <- c("D=5", "D=10", "D=20", "D=100")

# Convert to long format
d_all_long <- d_all %>%
  pivot_longer(cols = everything(), names_to = "Num_Draws", 
               values_to = "Correlation")
# Add in correlation when using all posterior samples
d_all_long$cor_all <- d5$avg_cor_all
# Add in correlation for the average SRS scenario, and for the quantile scenario
d_all_long <- d_all_long %>%
  mutate(
    Num_Draws = factor(Num_Draws, levels = c("D=5", "D=10", "D=20", "D=100")),
    # SRS average across 100 iterations
    mean_SRS = case_when(
      Num_Draws == "D=5" ~ d5$avg_cor_SRS,
      Num_Draws == "D=10" ~ d10$avg_cor_SRS,
      Num_Draws == "D=20" ~ d20$avg_cor_SRS,
      Num_Draws == "D=100" ~ d100$avg_cor_SRS,
      .default = NA
    ),
    # Quantile method
    mean_quants = case_when(
      Num_Draws == "D=5" ~ d5$avg_cor_quants,
      Num_Draws == "D=10" ~ d10$avg_cor_quants,
      Num_Draws == "D=20" ~ d20$avg_cor_quants,
      Num_Draws == "D=100" ~ d100$avg_cor_quants,
      .default = NA
    )
  )

# Create plot displaying results
d_all_long %>%
  ggplot(aes(x = Correlation, fill = Num_Draws)) + 
  geom_histogram(bins = 20) +
  facet_wrap(. ~ Num_Draws, nrow = 2) + 
  geom_vline(aes(xintercept = mean_SRS, color = "Mean across SRS Draws"), 
             linetype = "dashed") + 
  geom_vline(aes(xintercept = mean_quants, color = "Quantile Draws"), 
             linetype = "dashed") + 
  geom_vline(aes(xintercept = cor_all, color = "All Posterior Samples"), 
             linetype = "dashed") + 
  facet_wrap(~ Num_Draws, nrow = 2) + 
  scale_color_manual(values = c("Mean across SRS Draws" = "red", 
                                "Quantile Draws" = "orange",
                                "All Posterior Samples" = "black")) +
  theme_bw() +
  labs(color = "Type of Sampling for Draws") + 
  xlab("Correlation for 100 SRS Draws") + ylab("Frequency") 


# plot(wts_post[, 1], wts_post[, 1000])
# cor(wts_post[, 1], wts_post[, 1000])
# plot(wts_post[, 900], wts_post[, 1000])
# cor(wts_post[, 900], wts_post[, 1000])

#======================== Old code =============================================
# ### Variable importance plot
# bart_mod <- est_weights$fit_pi_B
# # Compute row percentages
# percent_varcount = bart_mod$varcount / apply(bart_mod$varcount, 1, sum)
# # Mean of row percentages
# mean_percent_varcount = apply(percent_varcount, 2, mean)
# # Quantiles of row percentages
# quant_percent_varcount = apply(percent_varcount, 2, quantile, probs=c(.025,.975))
# # Create plot
# rgy = range(quant_percent_varcount)
# p <- ncol(percent_varcount)
# plot(c(1,p),rgy,type="n",xlab="variable",ylab="post mean, percent var use",axes=FALSE)
# axis(1,at=1:p,labels=names(mean_percent_varcount),cex.lab=0.7,cex.axis=0.7)
# axis(2,cex.lab=1.2,cex.axis=1.2)
# lines(1:p,mean_percent_varcount,col="black",lty=4,pch=4,type="b",lwd=1.5)
# for(i in 1:p) {
#   lines(c(i,i),quant_percent_varcount[,i],col="blue",lty=3,lwd=1.0)
# }

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
#
#
# #============ Interaction models
# ### Core, interactions
# # List of model paths to load
# mod_names_list <- list(paste0(wd, res_dir, "t2d_unaware_core_intlogreg.RData"),
#                        paste0(wd, res_dir, "htn_unaware_core_intlogreg.RData"),
#                        paste0(wd, res_dir, "chol_unaware_core_intlogreg.RData"))
# # List of loaded models
# mod_res_list <- list()
# for (i in 1:length(mod_names_list)) {
#   load(mod_names_list[[i]])
#   mod_res_list[[i]] <- wtd_logreg_res
# }
# core_int_mod_res_list <- mod_res_list
# 
# 
# #int_categs <- c("Female_Active", "Female_Sedentary", "Male_Active", "Male_Sedentary")
# int_categs <- c("Female", "Male")
# class_levels <- 2:4
# dbh_classes <- paste0("dbh_class", class_levels)
# 
# quant_lb <- 0.025
# quant_ub <- 0.975
# 
# int_df_list <- list()
# for (i in 1:length(mod_res_list)) {
#   fit_obj <- mod_res_list[[i]]
#   fit_samples <- fit_obj$all_adj_parms
#   # Parameter names
#   mod_mat <- brms::make_standata(formula = fit_obj$brms_mod$formula, 
#                                  data = fit_obj$data, family = "bernoulli")
#   parm_names <- colnames(mod_mat$X)
#   parm_names_cs <- parm_names[-1]  # Remove intercept
#   colnames(fit_samples)[1:length(parm_names_cs)] <- parm_names_cs
#   
#   outcome_df_list <- list()
#   for (k in 1:length(dbh_classes)) {
#     dbh_class <- dbh_classes[k]
#     coef_names <- list(dbh_class, 
#                        c(paste0(dbh_class, c("", ":SexMale"))))
#     # coef_names <- list(dbh_class, 
#     #                    c(paste0(dbh_class, c("", ":Physical_activitySedentary"))),
#     #                    c(paste0(dbh_class, c("", ":SexMale"))),
#     #                    c(paste0(dbh_class, c("", ":SexMale", ":Physical_activitySedentary"))))
#     coef_inds <- lapply(coef_names, function(x) which(parm_names_cs %in% x))
#     int_samples <- lapply(coef_inds, function(x) 
#       apply(as.matrix(fit_samples[, x]), 1, prod))
#     int_summ <- lapply(int_samples, function(x) 
#       c(mean(x, na.rm = TRUE), 
#         quantile(x, c(quant_lb, quant_ub), na.rm = TRUE),
#         mean(x > 0, na.rm = TRUE),
#         mean(x < 0, na.rm = TRUE)))
#     # Exponentiate
#     int_summ_df <- as.data.frame(do.call("rbind", int_summ))
#     int_summ_exp <- exp(int_summ_df)
#     colnames(int_summ_exp) <- c("Mean", paste0(quant_lb*100, "%"), 
#                                paste0(quant_ub*100, "%"), "P(xi>0)", "P(xi<0)")
#     int_summ_exp$int_categs <- int_categs
#     outcome_df_list[[k]] <- int_summ_exp
#   }
#   
#   outcome_df <- do.call("rbind", outcome_df_list)
#   outcome_df$dbh_class <- rep(class_levels, each = nrow(outcome_df_list[[1]]))
#   
#   int_df_list[[i]] <- outcome_df
# }
# int_df <- do.call("rbind", int_df_list)
# int_df$Outcome <- rep(c("Type 2 Diabetes", "Hypertension", "High Cholesterol"),
#                       each = nrow(int_df_list[[1]]))
# 
# # Create Plot
# plot_df <- int_df
# colnames(plot_df) <- c("OR", "Lower", "Upper", "P(xi>0)", "P(xi<0)",
#                        "Variable", "Pattern", "Outcome")
# plot_df %>%
#   mutate(Variable = factor(Variable, levels = unique(Variable)),
#          Pattern = factor(Pattern)) %>%
#   ggplot(aes(x=fct_rev(Variable), y=OR, ymin=Lower, ymax=Upper,col=fct_rev(Pattern),
#              fill=fct_rev(Pattern))) + 
#   scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
#   facet_grid(. ~ Outcome) + 
#   theme_bw() + 
#   #specify position here
#   geom_linerange(size = 1, position=position_dodge(width = 0.5)) +
#   geom_hline(yintercept=1, lty=2) +
#   #specify position here too
#   geom_point(size=2, shape=21, colour="white", stroke = 0.5,
#              position=position_dodge(width = 0.5)) +
#   xlab("Variable") + ylab("Odds Ratio Scale") + 
#   geom_hline(yintercept = 1, linetype = "dashed", col = "gray") + 
#   guides(fill = guide_legend(reverse = TRUE),
#          color = guide_legend(reverse = TRUE)) + 
#   labs(fill = "Dietary\nBehavior\nClass", col = "Dietary\nBehavior\nClass") +
#   ylim(0, 30)+ 
#   coord_flip() 
# 
# 
# 
# 
# 
# #============ DIAGNOSIS => DBH
# ### Full, no interactions
# # List of model paths to load
# mod_names_list <- list(paste0(wd, res_dir, "t2d_diag_fulllogreg.RData"),
#                        paste0(wd, res_dir, "htn_diag_fulllogreg.RData"),
#                        paste0(wd, res_dir, "chol_diag_fulllogreg.RData"))
# # List of loaded models
# mod_res_list <- list()
# for (i in 1:length(mod_names_list)) {
#   load(mod_names_list[[i]])
#   mod_res_list[[i]] <- wtd_logreg_res
# }
# # List of summarized models
# summ_mod_list <- lapply(mod_res_list, function(x) 
#   summarize_parms(x, quant_lb = 0.025, quant_ub = 0.975, exponentiate = TRUE,
#                   diag = TRUE))
# # Plot models
# # create_mod_plot(summ_mod_list)
# diag_mod_res_list <- mod_res_list
# 
# 
# n_mods <- length(summ_mod_list)
# n_vars <- nrow(summ_mod_list[[1]])
# summ_reshape <- do.call("rbind", summ_mod_list)
# summ_reshape$Outcome <- rep(c("Type 2 Diabetes", "Hypertension", "High Cholesterol"),
#                             each = n_vars)
# summ_reshape <- summ_reshape %>%
#   drop_na(dbh_class)  # remove the lprior and b_mu terms
# 
# plot_df <- summ_reshape %>%
#   filter(Variable != "Intercept")
# colnames(plot_df) <- c("OR", "Lower", "Upper", "P(xi>0)", "P(xi<0)", "Pattern",
#                        "Variable", "Outcome")
# p_t2d <- plot_df %>% 
#   filter(Outcome == "Type 2 Diabetes") %>%
#   mutate(Variable = factor(Variable, levels = unique(Variable))) %>%
#   ggplot(aes(x=fct_rev(Variable), y=OR, ymin=Lower, ymax=Upper,col=fct_rev(Pattern),
#              fill=fct_rev(Pattern))) + 
#   scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
#   facet_grid(. ~ Pattern) + 
#   theme_bw() + 
#   #specify position here
#   geom_linerange(size = 1, position=position_dodge(width = 0.5)) +
#   geom_hline(yintercept=1, lty=2) +
#   #specify position here too
#   geom_point(size=2, shape=21, colour="white", stroke = 0.5,
#              position=position_dodge(width = 0.5)) +
#   xlab("Variable") + ylab("Odds Ratio Scale") + 
#   geom_hline(yintercept = 1, linetype = "dashed", col = "gray") + 
#   guides(fill = guide_legend(reverse = TRUE),
#          color = guide_legend(reverse = TRUE)) + 
#   labs(fill = "Dietary\nBehavior\nClass", col = "Dietary\nBehavior\nClass") +
#   coord_flip() 
# p_htn <- plot_df %>% 
#   filter(Outcome == "Hypertension") %>%
#   mutate(Variable = factor(Variable, levels = unique(Variable))) %>%
#   ggplot(aes(x=fct_rev(Variable), y=OR, ymin=Lower, ymax=Upper,col=fct_rev(Pattern),
#              fill=fct_rev(Pattern))) + 
#   scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
#   facet_grid(. ~ Pattern) + 
#   theme_bw() + 
#   #specify position here
#   geom_linerange(size = 1, position=position_dodge(width = 0.5)) +
#   geom_hline(yintercept=1, lty=2) +
#   #specify position here too
#   geom_point(size=2, shape=21, colour="white", stroke = 0.5,
#              position=position_dodge(width = 0.5)) +
#   xlab("Variable") + ylab("Odds Ratio Scale") + 
#   geom_hline(yintercept = 1, linetype = "dashed", col = "gray") + 
#   guides(fill = guide_legend(reverse = TRUE),
#          color = guide_legend(reverse = TRUE)) + 
#   labs(fill = "Dietary\nBehavior\nClass", col = "Dietary\nBehavior\nClass") +
#   coord_flip() 
# p_chol <- plot_df %>%
#   filter(Outcome == "High Cholesterol") %>%
#   mutate(Variable = factor(Variable, levels = unique(Variable))) %>%
#   ggplot(aes(x=fct_rev(Variable), y=OR, ymin=Lower, ymax=Upper,col=fct_rev(Pattern),
#              fill=fct_rev(Pattern))) + 
#   scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
#   facet_grid(. ~ Pattern) + 
#   theme_bw() + 
#   #specify position here
#   geom_linerange(size = 1, position=position_dodge(width = 0.5)) +
#   geom_hline(yintercept=1, lty=2) +
#   #specify position here too
#   geom_point(size=2, shape=21, colour="white", stroke = 0.5,
#              position=position_dodge(width = 0.5)) +
#   xlab("Variable") + ylab("Odds Ratio Scale") + 
#   geom_hline(yintercept = 1, linetype = "dashed", col = "gray") + 
#   guides(fill = guide_legend(reverse = TRUE),
#          color = guide_legend(reverse = TRUE)) + 
#   labs(fill = "Dietary\nBehavior\nClass", col = "Dietary\nBehavior\nClass") +
#   coord_flip() 
# ggpubr::ggarrange(p_t2d, p_htn, p_chol, nrow = 3, ncol = 1, 
#                   common.legend = TRUE, legend = "right")
# 
# 
# 
# 
# 
# 
# # load(paste0(wd, res_dir, "t2d_unaware_marglogreg.RData"))
# # t2d_unaware_marg <- summarize_parms(wtd_logreg_res, quant_lb = 0.05, 
# #                                     quant_ub = 0.95, exponentiate = TRUE)
# # load(paste0(wd, res_dir, "htn_unaware_marglogreg.RData"))
# # htn_unaware_marg <- summarize_parms(wtd_logreg_res, quant_lb = 0.05, 
# #                                     quant_ub = 0.95, exponentiate = TRUE)
# # load(paste0(wd, res_dir, "chol_unaware_marglogreg.RData"))
# # chol_unaware_marg <- summarize_parms(wtd_logreg_res, quant_lb = 0.05, 
# #                                      quant_ub = 0.95, exponentiate = TRUE)
# # load(paste0(wd, res_dir, "t2d_unaware_corelogreg.RData"))
# # t2d_unaware_core <- summarize_parms(wtd_logreg_res, quant_lb = 0.05, 
# #                                     quant_ub = 0.95, exponentiate = TRUE)
# # load(paste0(wd, res_dir, "htn_unaware_corelogreg.RData"))
# # htn_unaware_core <- summarize_parms(wtd_logreg_res, quant_lb = 0.05, 
# #                                     quant_ub = 0.95, exponentiate = TRUE)
# # load(paste0(wd, res_dir, "chol_unaware_corelogreg.RData"))
# # chol_unaware_core <- summarize_parms(wtd_logreg_res, quant_lb = 0.05, 
# #                                      quant_ub = 0.95, exponentiate = TRUE)
# # load(paste0(wd, res_dir, "t2d_unaware_fulllogreg.RData"))
# # t2d_unaware_full <- summarize_parms(wtd_logreg_res, quant_lb = 0.05, 
# #                                     quant_ub = 0.95, exponentiate = TRUE)
# # load(paste0(wd, res_dir, "htn_unaware_fulllogreg.RData"))
# # htn_unaware_full <- summarize_parms(wtd_logreg_res, quant_lb = 0.05, 
# #                                     quant_ub = 0.95, exponentiate = TRUE)
# # load(paste0(wd, res_dir, "chol_unaware_fulllogreg.RData"))
# # chol_unaware_full <- summarize_parms(wtd_logreg_res, quant_lb = 0.05, 
# #                                      quant_ub = 0.95, exponentiate = TRUE)
#   
#   tab_dbh_disease <- as.data.frame(matrix(NA, nrow = 15, ncol = 4))
# tab_dbh_disease[, 1] <- c("Unadjusted Model", "Reference", "DBP2", 
#                           "DBP3", "DBP4", "Model 1", "Reference", "DBP2", 
#                           "DBP3", "DBP4", "Model 2", "Reference", "DBP2", 
#                           "DBP3", "DBP4")
# colnames(tab_dbh_disease) <- c("Covariate", "Type 2 Diabetes", "Hypertension", 
#                                "High Cholesterol")
# tab_dbh_disease[3:5, 2] <- apply(t2d_unaware_marg[1:3, ], 1, function(x)
#   paste0(x[1], " (", x[2], ", ", x[3], ") [", x[4], "]"))
# tab_dbh_disease[3:5, 3] <- apply(htn_unaware_marg[1:3, ], 1, function(x)
#   paste0(x[1], " (", x[2], ", ", x[3], ") [", x[4], "]"))
# tab_dbh_disease[3:5, 4] <- apply(chol_unaware_marg[1:3, ], 1, function(x)
#   paste0(x[1], " (", x[2], ", ", x[3], ") [", x[4], "]"))
# tab_dbh_disease[8:10, 2] <- apply(t2d_unaware_core[1:3, ], 1, function(x)
#   paste0(x[1], " (", x[2], ", ", x[3], ") [", x[4], "]"))
# tab_dbh_disease[8:10, 3] <- apply(htn_unaware_core[1:3, ], 1, function(x)
#   paste0(x[1], " (", x[2], ", ", x[3], ") [", x[4], "]"))
# tab_dbh_disease[8:10, 4] <- apply(chol_unaware_core[1:3, ], 1, function(x)
#   paste0(x[1], " (", x[2], ", ", x[3], ") [", x[4], "]"))
# tab_dbh_disease[13:15, 2] <- apply(t2d_unaware_full[1:3, ], 1, function(x)
#   paste0(x[1], " (", x[2], ", ", x[3], ") [", x[4], "]"))
# tab_dbh_disease[13:15, 3] <- apply(htn_unaware_full[1:3, ], 1, function(x)
#   paste0(x[1], " (", x[2], ", ", x[3], ") [", x[4], "]"))
# tab_dbh_disease[13:15, 4] <- apply(chol_unaware_full[1:3, ], 1, function(x)
#   paste0(x[1], " (", x[2], ", ", x[3], ") [", x[4], "]"))
# 
# #write.csv(tab_dbh_disease,
# #          file = paste0(wd, res_dir, "tab_dbh_disease_90CI.csv"), row.names = FALSE)
# 
# 
# 
# #=================== DIAGNOSIS -> DBH ==========================================
# load("~/Documents/GitHub/WOLCAN/Application/Results/t2d_diag_fulllogreg.RData")
# sum_temp <- as.data.frame(round(t(apply(wtd_logreg_res$all_adj_parms, 2, function(x) 
#   c(mean(x, na.rm = TRUE), 
#     quantile(x, c(0.05, 0.95), na.rm = TRUE),
#     mean(x > 0, na.rm = TRUE),
#     mean(x < 0, na.rm = TRUE)))), 3))
# mod_mat <- brms::make_standata(formula = wtd_logreg_res$brms_mod$formula, 
#                                data = wtd_logreg_res$data, family = "categorical")
# parm_names <- colnames(mod_mat$X_mu2)
# sum_temp$names <- c(rep(c(parm_names[-1], parm_names[1]), 3), 
#                     "lprior", "bmu2", "bmu3", "bmu4")