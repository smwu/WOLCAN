#=================================================
# Running PROSPECT application data
# Author: Stephanie Wu
# Date created: 2024/09/13
# Date updated: 2024/09/13
#=================================================


rm(list = ls())

library(tidyverse)  # data wrangling
library(readxl)     # read excel files
library(mice)       # missingness pattern
library(BART)  # BART
library(survey)  # survey functions
library(parallel)  # parallel processing
library(abind)
library(baysc)
library(furniture)  # table1
library(corrplot)
library(brms)  # outcome regression
library(bayesplot)
library(csSampling)  # variance adjustment
library(rstan)
# rstan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
# wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/"
data_dir <- "Application/Cleaned_Data/"  # Data directory
res_dir <- "Application/Results/"        # Results directory
code_dir <- "Application/Code/"  # Code directory

# Source R model functions
source(paste0(wd, "Model_Code/", "model_functions.R"))
source(paste0(wd, "Summary_Code/", "summary_functions.R"))

#==================== Read in and prepare data =================================
# Read in PRCS data
prcs_cleaned <- read.csv(paste0(wd, data_dir, "prcs_cleaned.csv"))
# Restrict to ages 30 to 75
prcs_cleaned <- prcs_cleaned %>% filter(Age >= 30 & Age < 76)
# Drop NA: from 79858 to 77907 (1951 removed)
prcs_drop_na <- prcs_cleaned %>% drop_na()

# Read in PROSPECT data
prospect_cleaned <- read.csv(paste0(wd, data_dir, "prospect_cleaned.csv"))
# Drop NA: from 1690 to 1018 (672 removed)
prospect_drop_na <- prospect_cleaned %>% drop_na()


# Plot missingness pattern
sort(apply(prospect_cleaned, 2, function(x) sum(is.na(x))), decreasing = TRUE)
test <- md.pattern(prospect_cleaned, rotate.names = TRUE)

# REMEMBER TO RELEVEL DEPENDING ON NECESSARY REFERENCE
# Social_support = factor(case_when(  # CHECK THIS!!!!!
#   SS_SCORE <= 18 ~ 0, # low support
#   SS_SCORE >= 19 ~ 1, # high support (reference)
# ), levels = c(1, 0)),  
# Perceived_stress = factor(case_when(  # CHECK THIS!!!
#   PSS_SCORE <= 28 ~ 0, # low stress (reference)
#   PSS_SCORE >= 29 ~ 1, # high stress
# ), levels = c(0, 1)),

# Only the dietary behavior variables
prospect_dbh <- prospect_cleaned %>% 
  select(studyid, purchase_person:eat_vegetarian) 
  # %>%mutate_all(as.factor)
# Drop NA: from 1690 to 1544 (146 removed)
prospect_dbh_drop_na <- prospect_dbh %>% drop_na()

# Restrict to DBH complete cases
restrict_data <- prospect_cleaned %>% 
  filter(studyid %in% prospect_dbh_drop_na$studyid)


item_labels <- colnames(prospect_dbh)[-1]
class_title <- "Dietary Behavior Pattern"
categ_title <- "Risk Level"
categ_labels <- c("Low", "Med", "High")
# Display plots of results
res_plots <- function(res) {
  print(baysc::plot_pattern_profiles(res, item_labels = item_labels, 
                                     class_title = class_title, 
                                     categ_title = categ_title,
                                     categ_labels = categ_labels))
  print(baysc::plot_pattern_probs(res, item_labels = item_labels, 
                                  class_title = class_title, 
                                  categ_title = categ_title,
                                  categ_labels = categ_labels))
  print(baysc::plot_class_dist(res))
  # Traceplot for pi
  plot(res$estimates_adjust$pi_red[, 1], type = "l", ylim = c(0, 1))
  hist(res$data_vars$sampling_wt, breaks = 30)
}

cov_df <- restrict_data %>%
  select(Sex, Educ, Age, Inc_hh, Urban, Ethnicity, Smoking_status, 
         Drinking_status, Food_security, WIC_SNAP, Social_support, 
         Perceived_stress, Depression, Anxiety) %>%
  mutate_at(c("Educ", "Inc_hh", "Smoking_status", "Drinking_status"), as.factor)

#==================== Correlation analysis =====================================

# Use spearman for ordinal data
cor <- round(cor(prospect_dbh_drop_na[, -1], method = "spearman"), 1)
corrplot(cor, order = "AOE")
corrplot(cor)
corrplot(cor(prospect_dbh_drop_na[, -1], method = "kendall"))
# png(paste0(wd, "Tables_Figures/", "dbh_corrplot.png"), width = 900, height = 900)
corrplot.mixed(cor, number.cex = 0.8, tl.cex = 1.1, 
               tl.pos = "lt", upper = "ellipse", number.digits = 1)
# dev.off()

#==================== Run WOLCAN removing missing DBH ==========================

### Define data and parameters
# Estimating pseudo-weights
x_mat <- as.matrix(prospect_dbh_drop_na[, -1])  # Multivariate categorical variables
selection_covs <- c("Sex", "Educ", "Age", "Inc_hh", "Urban")
dat_B <- restrict_data %>%  # Covariates for NPS
  select(all_of(selection_covs)) %>%
  mutate_at(c("Sex", "Educ", "Inc_hh", "Urban"), as.factor)
nrow(dat_B %>% drop_na()) # 414 dropped
dat_R <- prcs_cleaned %>%  # Covariates for reference
  select(all_of(selection_covs)) %>%
  # mutate(Age = Age - mean(Age, na.rm = TRUE)) %>% # Center Age
  mutate_at(c("Sex", "Educ", "Inc_hh", "Urban"), as.factor)
nrow(dat_R %>% drop_na()) # 6185 dropped
pred_model <- "bart"  # Prediction model for selection
pred_covs_B <- selection_covs  # Covariates to predict NPS selection
pred_covs_R <- selection_covs  # Covariates to predict RS selection
pi_R <- 1 / prcs_cleaned$Weight  # RS selection probabilites for those in RS
# Set pi_R = 1 to 0.999999 to prevent Inf/Nan in BART
pi_R <- ifelse(pi_R == 1, 0.999999, pi_R)
hist(pi_R, breaks = 30)
summary(pi_R)
hat_pi_R <- NULL  # RS selection probabilities for those in NPS
num_post <- 1000  # Number of posterior BART draws for estimating weights
frame_B <- 1  # Coverage probability of NPS frame
frame_R <- 1  # Coverage probability of RS frame
trim_method = "t2"  # Trimming method using IQR
trim_c = 20         # Trimming constant

# Model estimation
D <- 20            # Number of sets of MI pseudo-weights
parallel <- FALSE   # Whether to parallelize for MI
n_cores <- 4       # Number of cores to use for parallelization
wts_adj <- "MI"    # Adjustment method to account for weights uncertainty
adjust <- TRUE     # Whether to adjust variance for pseudo-likelihood
tol <- 1e-8        # Underflow tolerance
num_reps <- 100    # Number of bootstrap replicates for WS adjustment
run_adapt <- TRUE  # Whether to run adaptive sampler to get K
K_max <- 30        # Max number of latent classes for adaptive sampler
adapt_seed <- 1    # Seed for adaptive sampler
fixed_seed <- 1    # Seed for fixed sampler
n_runs <- 20000    # Number of MCMC iterations
burn <- 10000      # Burn-in
thin <- 5          # Thinning
update <- 5000     # Display update frequency
class_cutoff <- 0.05


# Save as complete case DBH (note: this was used for the restricted age run)
save_res <- TRUE
save_path <- paste0(wd, res_dir, "cc_dbh")

# Run weighted model
res <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R, 
              pred_model = pred_model, pred_covs_B = pred_covs_B, 
              pred_covs_R = pred_covs_R, 
              pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post, 
              frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
              trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
              wts_adj = wts_adj, adjust = adjust, tol = tol, 
              num_reps = num_reps, run_adapt = run_adapt, 
              K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
              fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
              n_runs = n_runs, burn = burn, thin = thin, update = update, 
              save_res = save_res, save_path = save_path, 
              overwrite = TRUE)
class(res) <- "wolca"

baysc::plot_pattern_profiles(res)
baysc::plot_pattern_probs(res)
baysc::plot_class_dist(res)
sum(res$estimates_adjust$messages_draws == "no error")
# Traceplot for pi
plot(res$estimates_adjust$pi_red[, 1], type = "l", ylim = c(0, 1))
temp <- res$estimates_adjust$pi_red[, 1]
# Troubleshooting
# Check draw 19
res$estimates_adjust$K_red_draws_all
pi_draw19 <- res$estimates_adjust$pi_red[36001:38000, ]

hist(est_weights$wts)
summary(est_weights$wts)


### Re-run without variance adjustment
adjust <- FALSE
save_path <- paste0(wd, res_dir, "no_varadj_cc_dbh")

res <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R, 
              pred_model = pred_model, pred_covs_B = pred_covs_B, 
              pred_covs_R = pred_covs_R, 
              pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post, 
              frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
              trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
              wts_adj = wts_adj, adjust = adjust, tol = tol, 
              num_reps = num_reps, run_adapt = run_adapt, 
              K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
              fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
              n_runs = n_runs, burn = burn, thin = thin, update = update, 
              save_res = save_res, save_path = save_path, 
              overwrite = FALSE)
load(paste0(wd, res_dir, "no_varadj_cc_dbh", "_wolcan_results.RData"))
res_no_varadj_cc_dbh <- res

res_plots(res = res_no_varadj_cc_dbh)
vars_across_class_temp(c_all = as.factor(res_no_varadj_cc_dbh$estimates_adjust$c_all), 
                         cov_df = cov_df, 
                         sampling_wt = res_no_varadj_cc_dbh$data_vars$sampling_wt,
                         digits = 1, cluster_id = 1:nrow(restrict_data),
                         stratum_id = rep(1, nrow(restrict_data)), 
                         col_props = TRUE, res = res_no_varadj_cc_dbh)
# baysc::vars_across_class(c_all = as.factor(res_no_varadj_cc_dbh$estimates_adjust$c_all), 
#                          cov_df = cov_df, 
#                          sampling_wt = res_no_varadj_cc_dbh$data_vars$sampling_wt,
#                          digits = 1, cluster_id = 1:nrow(restrict_data),
#                          stratum_id = rep(1, nrow(restrict_data)), 
#                          col_props = TRUE, res = res_no_varadj_cc_dbh)


### Re-run without variance adjustment and with age-subsetted PRCS data
adapt_seed <- 1
adjust <- FALSE
parallel <- FALSE
save_path <- paste0(wd, res_dir, "agePRCS_no_varadj_cc_dbh")
save_res <- TRUE

res_agePRCS <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R, 
              pred_model = pred_model, pred_covs_B = pred_covs_B, 
              pred_covs_R = pred_covs_R, 
              pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post, 
              frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
              trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
              wts_adj = wts_adj, adjust = adjust, tol = tol, 
              num_reps = num_reps, run_adapt = run_adapt, 
              K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
              fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
              n_runs = n_runs, burn = burn, thin = thin, update = update, 
              save_res = save_res, save_path = save_path, 
              overwrite = FALSE)
load(paste0(save_path, "_wolcan_results.RData"))
res_agePRCS <- res
res_plots(res = res_agePRCS)


### Run separately for males and females, w/ subsetted age, w/out variance adjustment
# Should be okay because PRCS post-stratifies by sex, so the subtotals should 
# be matched to the census totals by sex
# Run weighted model for females
adjust <- FALSE
females_prosp <- which(restrict_data$Sex == 1)  # 1155/1550 (75%)
females_prcs <- which(prcs_cleaned$Sex == 1)  # 43134/79858 (54%)
save_path <- paste0(wd, res_dir, "f_no_varadj_cc_dbh")
res_f <- wolcan(x_mat = x_mat[females_prosp, ], dat_B = dat_B[females_prosp, ], 
                dat_R = dat_R[females_prcs, ], 
                pred_model = pred_model, pred_covs_B = pred_covs_B, 
                pred_covs_R = pred_covs_R, 
                pi_R = pi_R [females_prcs], hat_pi_R = hat_pi_R, num_post = num_post, 
                frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
                trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
                wts_adj = wts_adj, adjust = adjust, tol = tol, 
                num_reps = num_reps, run_adapt = run_adapt, 
                K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
                fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
                n_runs = n_runs, burn = burn, thin = thin, update = update, 
                save_res = save_res, save_path = save_path, 
                overwrite = FALSE)

# Run weighted model for males
adjust <- FALSE
males_prosp <- which(restrict_data$Sex == 0)  # 395/1550 (25%)
males_prcs <- which(prcs_cleaned$Sex == 0)  # 36724/79858 (46%)
save_path <- paste0(wd, res_dir, "m_no_varadj_cc_dbh")
res_m <- wolcan(x_mat = x_mat[males_prosp, ], dat_B = dat_B[males_prosp, ], 
                dat_R = dat_R[males_prcs, ], 
                pred_model = pred_model, pred_covs_B = pred_covs_B, 
                pred_covs_R = pred_covs_R, 
                pi_R = pi_R [males_prcs], hat_pi_R = hat_pi_R, num_post = num_post, 
                frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
                trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
                wts_adj = wts_adj, adjust = adjust, tol = tol, 
                num_reps = num_reps, run_adapt = run_adapt, 
                K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
                fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
                n_runs = n_runs, burn = burn, thin = thin, update = update, 
                save_res = save_res, save_path = save_path, 
                overwrite = FALSE)
res_plots(res = res_f)
res_plots(res = res_m)


### Re-run with age-subsetted and class cutoff = 0.01
class_cutoff <- 0.01
adjust <- FALSE
save_path <- paste0(wd, res_dir, "cutoff_agePRCS_no_varadj_cc_dbh")
res_cutoff_agePRCS <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R, 
              pred_model = pred_model, pred_covs_B = pred_covs_B, 
              pred_covs_R = pred_covs_R, 
              pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post, 
              frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
              trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
              wts_adj = wts_adj, adjust = adjust, tol = tol, 
              num_reps = num_reps, run_adapt = run_adapt, 
              K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
              fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
              n_runs = n_runs, burn = burn, thin = thin, update = update, 
              save_res = save_res, save_path = save_path, 
              overwrite = FALSE)
load(paste0(save_path, "_wolcan_results.RData"))
res_cutoff_agePRCS <- res
res_plots(res = res_cutoff_agePRCS)
# Only 3 of the iterations had K=5. The rest had more


### Re-run with c = 10, without variance adjustment and with age-subsetted PRCS 
trim_c = 10        
adapt_seed <- 1
adjust <- FALSE
parallel <- FALSE
save_path <- paste0(wd, res_dir, "trim10_agePRCS_no_varadj_cc_dbh")
save_res <- TRUE

res_trim10 <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R, 
                      pred_model = pred_model, pred_covs_B = pred_covs_B, 
                      pred_covs_R = pred_covs_R, 
                      pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post, 
                      frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
                      trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
                      wts_adj = wts_adj, adjust = adjust, tol = tol, 
                      num_reps = num_reps, run_adapt = run_adapt, 
                      K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
                      fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
                      n_runs = n_runs, burn = burn, thin = thin, update = update, 
                      save_res = save_res, save_path = save_path, 
                      overwrite = FALSE)
res_plots(res = res_trim10)

# Troubleshooting 
res_all <- list()
for (d in 1:D) {
  load(paste0(save_path, "draw_", d, "_results.RData"))
  res_all[[d]] <- est_d
}


### Re-run with c = 25, without variance adjustment and with age-subsetted PRCS 
trim_c = 25        
adapt_seed <- 1
adjust <- FALSE
parallel <- FALSE
save_path <- paste0(wd, res_dir, "trim25_agePRCS_no_varadj_cc_dbh")
save_res <- TRUE

res_trim25 <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R, 
                     pred_model = pred_model, pred_covs_B = pred_covs_B, 
                     pred_covs_R = pred_covs_R, 
                     pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post, 
                     frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
                     trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
                     wts_adj = wts_adj, adjust = adjust, tol = tol, 
                     num_reps = num_reps, run_adapt = run_adapt, 
                     K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
                     fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
                     n_runs = n_runs, burn = burn, thin = thin, update = update, 
                     save_res = save_res, save_path = save_path, 
                     overwrite = FALSE)
res_plots(res = res_trim25)

### Run weighted model with subsetted age, no variance adjustment, and removing
### correlation redundancy
x_mat <- as.matrix(prospect_dbh_drop_na %>% # Multivariate categorical variables
                     select(-c(studyid, control_fat, control_portions, 
                               take_out_freq, restaurant_freq))) 
adapt_seed <- 1
adjust <- FALSE
parallel <- FALSE
save_res <- TRUE
save_path <- paste0(wd, res_dir, "redun_agePRCS_no_varadj_cc_dbh")
res_redun <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R, 
              pred_model = pred_model, pred_covs_B = pred_covs_B, 
              pred_covs_R = pred_covs_R, 
              pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post, 
              frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
              trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
              wts_adj = wts_adj, adjust = adjust, tol = tol, 
              num_reps = num_reps, run_adapt = run_adapt, 
              K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
              fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
              n_runs = n_runs, burn = burn, thin = thin, update = update, 
              save_res = save_res, save_path = save_path, 
              overwrite = FALSE)
baysc::plot_pattern_profiles(res = res_redun, item_labels = colnames(x_mat),
                             categ_labels = categ_labels, 
                             categ_title = categ_title, 
                             class_title = class_title)
res_plots(res = res_redun)

#===================== Unweighted model ========================================
### Run unweighted model 
class_cutoff <- 0.05
save_path <- paste0(wd, res_dir, "agePRCS_unwt_cc_dbh")
res_agePRCS_unwt <- wolca(x_mat = x_mat, sampling_wt = rep(1, nrow(x_mat)), 
                  run_sampler = "both", 
                   K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
                   fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
                   n_runs = n_runs, burn = burn, thin = thin, update = update, 
                   save_res = save_res, save_path = save_path)
load(paste0(save_path, "_wolca_results.RData"))
res_agePRCS_unwt <- res
res_plots(res = res_agePRCS_unwt)

### Run unweighted model with full PRCS dataset WILL BE SAME B/C DOESN'T USE PRCS
class_cutoff <- 0.05
save_path <- paste0(wd, res_dir, "unwt_cc_dbh")
res_unwt <- wolca(x_mat = x_mat, sampling_wt = rep(1, nrow(x_mat)), 
                  run_sampler = "both", 
                  K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
                  fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
                  n_runs = n_runs, burn = burn, thin = thin, update = update, 
                  save_res = save_res, save_path = save_path)
load(paste0(save_path, "_wolca_results.RData"))
res_unwt <- res
res_plots(res = res_unwt)

# CHECK UNWT plot_class_probs ERROR!!!!!

#==================== Examine weights ==========================================
load("~/Documents/GitHub/WOLCAN/Application/Results/agePRCS_no_varadj_cc_dbh_wolcan_weights.RData")
hist(est_weights$wts, breaks = 30)
summary(est_weights$wts)
boxplot(est_weights$wts)
# Sequence of cutoffs depending on c value
sapply(c(1, 5, 10, 15, 20), 
       function(c) median(est_weights$wts) + c*IQR(est_weights$wts))
# All posterior weight draws
hist(c(est_weights$wts_post), breaks = 30)
summary(c(est_weights$wts_post))
sapply(c(1, 5, 10, 15, 20), 
       function(c) median(c(est_weights$wts_post)) + c*IQR(c(est_weights$wts_post)))

# Upweighted individuals
restrict_weights <- restrict_data %>%
  mutate(weight = est_weights$wts, 
         high_weight = ifelse(weight >= 4000, 1, 0)) %>%
  mutate_at(c("Sex", "Educ", "Inc_hh", "Urban", "Ethnicity", "Smoking_status", 
              "Drinking_status", "Physical_activity", "Food_security", "WIC_SNAP", 
              "Depression", "Anxiety"), as.factor)

tb1 <- furniture::table1(restrict_weights, Sex, Age, Educ, Inc_hh, Urban, Ethnicity, 
                  Smoking_status, Drinking_status, Physical_activity, 
                  Food_security, WIC_SNAP, Social_support, Perceived_stress,
                  Depression, Anxiety, splitby = ~high_weight, row_wise = FALSE, 
                  na.rm = FALSE, total = TRUE, type = "condensed")
tb1 <- as.data.frame(tb1$Table1)
tb1
#==================== Comparison with sociodemog vars ==========================



#==================== Examination of outcome variables =========================
df <- restrict_data %>%
  select(hypertension:high_ldl)
# Table of prevalences
outcome_prev <- as.data.frame(matrix(NA, nrow = 3, ncol = 8))
rownames(outcome_prev) <- c("Hypertension", "Type 2 Diabetes", "High Cholesterol")
colnames(outcome_prev) <- c("Num At Risk", "% At Risk", "Num Aware", 
                            "% At Risk Aware", "Num Treated", "% Aware Treated", 
                            "% At Risk Treated", "Num NA")
outcome_prev[1, ] <- c(sum(df$hypertension, na.rm = TRUE), 
                       round(mean(df$hypertension, na.rm = TRUE) * 100, 1),
                       sum(df[df$hypertension == 1, "htn_aware"], na.rm = TRUE),
                       round(mean(df[df$hypertension == 1, "htn_aware"], 
                                  na.rm = TRUE) * 100, 1),
                       sum(df[df$htn_aware == 1, "htn_treated"], na.rm = TRUE),
                       round(mean(df[df$htn_aware == 1, "htn_treated"], 
                                  na.rm = TRUE) * 100, 1),
                       round(mean(df[df$hypertension == 1, "htn_treated"], 
                                  na.rm = TRUE) * 100, 1),
                       sum(is.na(df$hypertension)))
outcome_prev[2, ] <- c(sum(df$diabetes, na.rm = TRUE), 
                       round(mean(df$diabetes, na.rm = TRUE) * 100, 1),
                       sum(df[df$diabetes == 1, "t2d_aware"], na.rm = TRUE),
                       round(mean(df[df$diabetes == 1, "t2d_aware"], 
                                  na.rm = TRUE) * 100, 1),
                       sum(df[df$t2d_aware == 1, "t2d_treated"], na.rm = TRUE),
                       round(mean(df[df$t2d_aware == 1, "t2d_treated"], 
                                  na.rm = TRUE) * 100, 1),
                       round(mean(df[df$diabetes == 1, "t2d_treated"], 
                                  na.rm = TRUE) * 100, 1),
                       sum(is.na(df$diabetes)))
outcome_prev[3, ] <- c(sum(df$high_cholesterol, na.rm = TRUE), 
                       round(mean(df$high_cholesterol, na.rm = TRUE) * 100, 1),
                       sum(df[df$high_cholesterol == 1, "chol_aware"], na.rm = TRUE),
                       round(mean(df[df$high_cholesterol == 1, "chol_aware"], 
                                  na.rm = TRUE) * 100, 1),
                       sum(df[df$chol_aware == 1, "chol_treated"], na.rm = TRUE),
                       round(mean(df[df$chol_aware == 1, "chol_treated"], 
                                  na.rm = TRUE) * 100, 1),
                       round(mean(df[df$high_cholesterol == 1, "chol_treated"], 
                                  na.rm = TRUE) * 100, 1),
                       sum(is.na(df$high_cholesterol)))


outcome_long <- outcome_prev %>%
  mutate(Total = nrow(restrict_data),
         Condition = rownames(outcome_prev)) %>%
  pivot_longer(cols = c(Total, `Num At Risk`, `Num Aware`, `Num Treated`),
               names_to = "Group", 
               values_to = "Number") %>%
  mutate(Group = factor(Group, levels = c("Total", "Num At Risk", "Num Aware", 
                                          "Num Treated")),
         Condition = factor(Condition, levels = c("Hypertension", "Type 2 Diabetes",
                                                  "High Cholesterol")))
  
outcome_long %>% ggplot(aes(x = Condition, y = Number, fill = Group)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
  scale_fill_brewer(palette = "Set2")

num_disease <- apply(restrict_data %>% select(hypertension, diabetes, high_cholesterol),
                     1, sum)
table(num_disease, useNA = "always")

#==================== Outcome regression =======================================

### Unweighted logistic regression
regr_dat <- as.data.frame(cbind(
  restrict_data, 
  dbh_class = factor(res_unwt$estimates$c_all)))
summary(glm(hypertension ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban, 
            data = regr_dat, family = binomial()))

summary(glm(diabetes ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban, 
            data = regr_dat, family = binomial()))

summary(glm(high_cholesterol ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban, 
            data = regr_dat, family = binomial()))


### Survey-weighted logistic regression
wtd_logreg <- function(res, data, formula_y) {
  svy_data <- data.frame(dbh_class = as.factor(res$estimates_adjust$c_all),
                         restrict_data, 
                         wts = res$data_vars$sampling_wt)
  svy_des <- survey::svydesign(ids = ~1, weights = ~wts, data = svy_data)
  log_reg <- survey::svyglm(formula = as.formula(formula_y), 
                            design = svy_des, 
                            family = stats::quasibinomial(link = "logit"))
  return(summary(log_reg))
}
wtd_logreg(res = res_no_varadj_cc_dbh, data = restrict_data, 
           formula_y = "hypertension ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban")
wtd_logreg(res = res_cc_dbh, data = restrict_data, 
           formula_y = "hypertension ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban")
wtd_logreg(res = res_agePRCS, data = restrict_data, 
           formula_y = "hypertension ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban")
wtd_logreg(res = res_cutoff_agePRCS, data = restrict_data, 
           formula_y = "hypertension ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban")


wtd_logreg(res = res_no_varadj_cc_dbh, data = restrict_data, 
           formula_y = "diabetes ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban")
wtd_logreg(res = res_cc_dbh, data = restrict_data, 
           formula_y = "diabetes ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban")

wtd_logreg(res = res_no_varadj_cc_dbh, data = restrict_data, 
           formula_y = "high_cholesterol ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban")
wtd_logreg(res = res_cc_dbh, data = restrict_data, 
           formula_y = "high_cholesterol ~ dbh_class + Age + Sex + Educ + Inc_hh + Urban")

# Main covariate set: food security, physical activity, smoking status, age, sex, 
# education, income, urban
# Other covariates: Ethnicity, drinking status, physical activity, food security, 
# social support, perceived stress, depression, anxiety, WIC/SNAP

# Create survey data using estimated dietary pattern assignments
svy_data <- data.frame(dbh_class = as.factor(res$estimates_adjust$c_all),
                       restrict_data, 
                       wts = res$data_vars$sampling_wt)
svy_data <- svy_data %>%
  mutate_at(c("Educ", "Inc_hh", "Smoking_status", "Drinking_status", 
              "Physical_activity"), as.factor) %>%
  mutate(htn_categ = as.factor(ifelse(htn_categ == 4, 3, htn_categ)),
         t2d_categ = as.factor(ifelse(t2d_categ == 4, 3, t2d_categ)),
         chol_categ = as.factor(ifelse(chol_categ == 4, 3, chol_categ)))
svy_des <- survey::svydesign(ids = ~1, weights = ~wts, data = svy_data)

survey::svyglm(formula = as.formula(formula_y), 
               design = svy_des, 
               family = stats::quasibinomial(link = "logit"))

library(brms)
# mod1 <- brms::brm(dbh_class ~ chol_categ)

summary(survey::svyglm(formula = 
                 as.formula(paste0("hypertension ~ dbh_class + Age + Sex + Educ",
                 "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
                 "+ Food_security")), 
               design = svy_des, 
               family = stats::quasibinomial(link = "logit")))

summary(survey::svyglm(formula = 
   as.formula(paste0("hypertension ~ dbh_class * htn_aware + Age + Sex + Educ",
                     "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
                     "+ Food_security")), 
   design = svy_des, 
   family = stats::quasibinomial(link = "logit")))

# ### Measuring disease category as an exposure for DBH class
# # Hypertension
# fit1 <- survey::svyglm(formula = 
#                          as.formula(paste0("dbh_class ~ htn_categ + Age + Sex + Educ",
#                                            "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
#                                            "+ Food_security")), 
#                        design = svy_des, 
#                        family = stats::quasibinomial(link = "logit"))
# summary(fit1)
# dim(fit1$model)
# 
# # Type 2 Diabetes
# fit2 <- survey::svyglm(formula = 
#                          as.formula(paste0("dbh_class ~ t2d_categ + Age + Sex + Educ",
#                                            "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
#                                            "+ Food_security")), 
#                        design = svy_des, 
#                        family = stats::quasibinomial(link = "logit"))
# summary(fit2)
# dim(fit2$model)
# 
# # High Cholesterol
# fit3 <- survey::svyglm(formula = 
#                          as.formula(paste0("dbh_class ~ chol_categ + Age + Sex + Educ",
#                                            "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
#                                            "+ Food_security")), 
#                        design = svy_des, 
#                        family = stats::quasibinomial(link = "logit"))
# summary(fit3)
# dim(fit3$model)
# 
# survey::svyglm(formula = 
#                  as.formula(paste0("hypertension ~ dbh_class + Age + Sex + Educ",
#                                    "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
#                                    "+ Food_security + Ethnicity + Drinking_status",
#                                    "+ WIC_SNAP + Social_support + Perceived_stress",
#                                    "+ Depression + Anxiety")), 
#                design = svy_des, 
#                family = stats::quasibinomial(link = "logit"))


### Restrict to those not at risk or at risk and unaware (dbh -> disease)

# Hypertension: 845 out of 1550 (55%)
svy_data_unaware_htn <- svy_data %>%
  filter(htn_categ != 3)
svy_des_unaware_htn <- survey::svydesign(ids = ~1, weights = ~wts, 
                                     data = svy_data_unaware_htn)
fit11 <- survey::svyglm(formula = 
                          as.formula(paste0("hypertension ~ dbh_class + Age + Sex + Educ",
                                            "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
                                            "+ Food_security")), 
                        design = svy_des_unaware_htn, 
                        family = stats::quasibinomial(link = "logit"))
summary(fit11)
dim(fit11$model)


# T2D: 1280 out of 1550 (83%)
svy_data_unaware_t2d <- svy_data %>%
  filter(t2d_categ != 3)
svy_des_unaware_t2d <- survey::svydesign(ids = ~1, weights = ~wts, 
                                     data = svy_data_unaware_t2d)
fit12 <- survey::svyglm(formula = 
                          as.formula(paste0("diabetes ~ dbh_class + Age + Sex + Educ",
                                            "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
                                            "+ Food_security")), 
                        design = svy_des_unaware_t2d, 
                        family = stats::quasibinomial(link = "logit"))
summary(fit12)
dim(fit12$model)

# Chol: 952 out of 1550 (61%)
svy_data_unaware_chol <- svy_data %>%
  filter(chol_categ != 3)
svy_des_unaware_chol <- survey::svydesign(ids = ~1, weights = ~wts, 
                                     data = svy_data_unaware_chol)
fit13 <- survey::svyglm(formula = 
                          as.formula(paste0("high_cholesterol ~ dbh_class + Age + Sex + Educ",
                                            "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
                                            "+ Food_security")), 
                        design = svy_des_unaware_chol, 
                        family = stats::quasibinomial(link = "logit"))
summary(fit13)
dim(fit13$model)


### Restrict to those with the disease (diagnosis -> dbh)

# Hypertension: 845 out of 1550 (55%)
svy_data_disease_htn <- svy_data %>%
  filter(htn_categ != 1) %>% 
  mutate(htn_categ = factor(htn_categ, levels = c("2", "3"), 
                            labels = c("Unaware", "Aware")))
svy_data_disease_htn <- survey::svydesign(ids = ~1, weights = ~wts, 
                                     data = svy_data_disease_htn)
fit21 <- survey::svyglm(formula = 
                          as.formula(paste0("dbh_class ~ htn_categ * (Age + Sex + Educ + t2d_aware + chol_aware)",
                                            "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
                                            "+ Food_security")), 
                        design = svy_data_disease_htn, 
                        family = stats::quasibinomial(link = "logit"))
summary(fit21)
dim(fit21$model)

# T2D: 395 out of 1550 (25%)
svy_data_disease_t2d <- svy_data %>%
  filter(t2d_categ != 1) %>% 
  mutate(t2d_categ = factor(t2d_categ, levels = c("2", "3"), 
                            labels = c("Unaware", "Aware")))
svy_data_disease_t2d <- survey::svydesign(ids = ~1, weights = ~wts, 
                                      data = svy_data_disease_t2d)
fit22 <- survey::svyglm(formula = 
                          as.formula(paste0("dbh_class ~ t2d_categ * (Age + Sex + Educ)",
                                            "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
                                            "+ Food_security")), 
                        design = svy_data_disease_t2d, 
                        family = stats::quasibinomial(link = "logit"))
summary(fit22)
dim(fit22$model)

# Chol:  out of 1550 (61%)
svy_data_disease_chol <- svy_data %>%
  filter(chol_categ != 1) %>% 
  mutate(chol_categ = factor(chol_categ, levels = c("2", "3"), 
                            labels = c("Unaware", "Aware")))
svy_data_disease_chol <- survey::svydesign(ids = ~1, weights = ~wts, 
                                      data = svy_data_disease_chol)
fit23 <- survey::svyglm(formula = 
                          as.formula(paste0("dbh_class ~ chol_categ * (Age + Sex + Educ)",
                                            "+ Inc_hh + Urban + Physical_activity + Smoking_status ",
                                            "+ Food_security")), 
                        design = svy_data_disease_chol, 
                        family = stats::quasibinomial(link = "logit"))
summary(fit23)
dim(fit23$model)



#========== Run regression models using csSampling package =====================
# Create survey data using estimated dietary pattern assignments
svy_data <- data.frame(dbh_class = as.factor(res$estimates_adjust$c_all),
                       restrict_data, 
                       wts = res$data_vars$w_all)
svy_data <- svy_data %>%
  mutate_at(c("Educ", "Inc_hh", "Smoking_status", "Drinking_status", 
              "Physical_activity"), as.factor) %>%
  mutate(htn_categ = as.factor(ifelse(htn_categ == 4, 3, htn_categ)),
         t2d_categ = as.factor(ifelse(t2d_categ == 4, 3, t2d_categ)),
         chol_categ = as.factor(ifelse(chol_categ == 4, 3, chol_categ)))

# Restrict to those not at risk or at risk and unaware (dbh -> disease)
# T2D: n=1446
svy_data_unaware_t2d <- svy_data %>%
  filter(t2d_categ != 3) %>%
  drop_na(diabetes, Age, Sex, Educ, Inc_hh, Urban, Physical_activity, 
          Smoking_status, Food_security)
svy_des_unaware_t2d <- survey::svydesign(ids = ~1, weights = ~wts, 
                                         data = svy_data_unaware_t2d)
model_formula <- as.formula(paste0("diabetes | weights(wts) ~ dbh_class ",
                                   "+ Age + Sex + Educ+ Inc_hh + Urban",
                                   " + Physical_activity + Smoking_status ",
                                   "+ Food_security")) 
brms_mod <- brms::brmsformula(model_formula, center = TRUE)
set.seed(1)
fit_htn_risk <- csSampling::cs_sampling_brms(svydes = svy_des_unaware_t2d, 
                                             brmsmod = brms_mod, 
                                              data = svy_data_unaware_t2d, 
                                              family = bernoulli(link = "logit"),
                                             ctrl_stan = 
                                               list(chains = 3, iter = 2000, 
                                                    warmup = 1000, thin = 5))
plot(fit_htn_risk)
temp <- summary(fit_htn_risk$stan_fit, probs = c(0.025, 0.975))
round(temp$summary, 3)



wts_draws <- res$estimates_adjust$wts_draws
D <- ncol(wts_draws)
adj_parms_draws <- list()
sampled_parms_draws <- list()
stan_fit_draws <- list()
for (d in 1:D) {
  set.seed(d)
  svy_data_d <- svy_data %>%
    mutate(wts_d = wts_draws[, d])
  svy_data_unaware_t2d_d <- svy_data_d %>%
    filter(t2d_categ != 3) %>%
    drop_na(diabetes, Age, Sex, Educ, Inc_hh, Urban, Physical_activity, 
            Smoking_status, Food_security)
  svy_des_unaware_t2d_d <- survey::svydesign(ids = ~1, weights = ~wts_d, 
                                           data = svy_data_unaware_t2d_d)
  fit_d <- csSampling::cs_sampling_brms(svydes = svy_des_unaware_t2d_d, 
                                        brmsmod = brms_mod, 
                                        data = svy_data_unaware_t2d_d, 
                                        family = bernoulli(link = "logit"),
                                        ctrl_stan = list(chains = 3, iter = 2000, 
                                                         warmup = 1000, thin = 5))
  sampled_parms_draws[[d]] <- fit_d$sampled_parms
  adj_parms_draws[[d]] <- fit_d$adjusted_parms
  stan_fit_draws[[d]] <- fit_d$stan_fit
}

# Function to run weighted logistic regression, incorporating variability from 
# estimated weights
# Inputs:
#   wts_draws: nxD matrix of the D draws from the weights posterior distribution
#   subset: String specifying whether to subset to those unaware ("unaware") to 
# measure impact of diet behavior pattern on disease, or to subset to those with 
# the disease ("disease") to measure impact of diagnosis on diet behavior. Must 
# be one of "unaware" or "disease".
#   condition: String specifying the condition to focus on. Must be one of "htn"
# (hypertension), "t2d" (type 2 diabetes), or "chol" (high cholesterol)
# Add variability of the weights
wtd_logreg_wolcan <- function(wts_draws, subset, condition, save_res = TRUE,
                              save_path) {
  # Check input arguments
  if (subset == "unaware") {
    filter_categ <- 3
  } else if (subset == "disease") {
    filter_categ <- 1
  } else {
    stop("Input argument 'subset' must be one of 'unaware' or 'disease'")
  }
  
  if (condition == "htn") {
    cond_categ <- "htn_categ"
    cond <- "hypertension"
  } else if (condition == "t2d") {
    cond_categ <- "t2d_categ"
    cond <- "diabetes"
  } else if (condition == "chol") {
    cond_categ <- "chol_categ"
    cond <- "cholesterol"
  } else {
    stop("Input argument 'condition' must be one of 'htn', 't2d', or 'chol'")
  }
  
  # Number of draws
  D <- ncol(wts_draws)  
  # Initialize lists
  adj_parms_draws <- list()
  sampled_parms_draws <- list()
  stan_fit_draws <- list()
  
  # Define brms formula
  model_formula <- as.formula(paste0(cond, " | weights(wts_d) ~ dbh_class ",
                                     "+ Age + Sex + Educ+ Inc_hh + Urban",
                                     " + Physical_activity + Smoking_status ",
                                     "+ Food_security")) 
  brms_mod <- brms::brmsformula(model_formula, center = TRUE)
  
  # For each draw, run the variance-adjusted survey-weighted regression and 
  # store results
  for (d in 1:D) {
    print("d: ", d)
    set.seed(d)
    svy_data_d <- svy_data %>%
      mutate(wts_d = wts_draws[, d])
    svy_data_subset_d <- svy_data_d %>%
      filter(cond_categ != filter_categ) %>% 
      drop_na(cond, Age, Sex, Educ, Inc_hh, Urban, Physical_activity, 
              Smoking_status, Food_security)
    svy_des_subset_d <- survey::svydesign(ids = ~1, weights = ~wts_d,
                                          data = svy_data_subset_d)
    fit_d <- csSampling::cs_sampling_brms(svydes = svy_des_subset_d, 
                                          brmsmod = brms_mod, 
                                          data = svy_data_subset_d, 
                                          family = bernoulli(link = "logit"),
                                          ctrl_stan = list(chains = 3, iter = 2000, 
                                                           warmup = 1000, thin = 5))
    sampled_parms_draws[[d]] <- fit_d$sampled_parms
    adj_parms_draws[[d]] <- fit_d$adjusted_parms
    stan_fit_draws[[d]] <- fit_d$stan_fit
  }
  
  all_adj_parms <- as.data.frame(do.call(rbind, adj_parms_draws))
  
  if (save_res) {
    wtd_logreg_res <- list(all_adj_parms = all_adj_parms, 
                           sampled_parms_draws = sampled_parms_draws,
                           stan_fit_draws = stan_fit_draws)
    save(wtd_logreg_res, file = paste0(save_path, "logreg.RData"))
  }
}

fit_t2d_unaware <- wtd_logreg_wolcan(wts_draws = wts_draws, subset = "unaware", 
                                     condition = "t2d", 
                                     save_path = paste0(wd, res_dir, "t2d_unaware"))
fit_htn_unaware <- wtd_logreg_wolcan(wts_draws = wts_draws, subset = "unaware", 
                                     condition = "htn", 
                                     save_path = paste0(wd, res_dir, "htn_unaware"))
fit_chol_unaware <- wtd_logreg_wolcan(wts_draws = wts_draws, subset = "unaware", 
                                     condition = "chol", 
                                     save_path = paste0(wd, res_dir, "chol_unaware"))

fit_t2d_disease <- wtd_logreg_wolcan(wts_draws = wts_draws, subset = "disease", 
                                     condition = "t2d", 
                                     save_path = paste0(wd, res_dir, "t2d_disease"))
fit_htn_disease <- wtd_logreg_wolcan(wts_draws = wts_draws, subset = "disease", 
                                     condition = "htn", 
                                     save_path = paste0(wd, res_dir, "htn_disease"))
fit_chol_disease <- wtd_logreg_wolcan(wts_draws = wts_draws, subset = "disease", 
                                      condition = "chol", 
                                      save_path = paste0(wd, res_dir, "chol_disease"))


all_adj_parms <- as.data.frame(do.call(rbind, adj_parms_draws))
summary_adj_parms <- t(apply(all_adj_parms, 2, function(x) c(mean(x, na.rm = TRUE), 
                                      quantile(x, c(0.025, 0.975), na.rm = TRUE),
                                      mean(x > 0, na.rm = TRUE),
                                      mean(x < 0, na.rm = TRUE))))
colnames(summary_adj_parms) <- c("Mean", "2.5%", "97.5%", "P(xi>0)", "P(xi<0)")

family <- bernoulli(link = "logit")
stancode <- do.call(make_stancode, c(list(brms_mod, data = svy_data_unaware_t2d,
                                          family = family, prior = NULL,
                                          stanvars = NULL,
                                          knots = NULL), stancode_args = list()))
mod_brms <- stan_model(model_code = stancode)
data_brms <- do.call(make_standata, c(list(brms_mod, data = svy_data_unaware_t2d,
                                           family = family, prior = NULL,
                                           stanvars = NULL,
                                           knots = NULL), standata_args = list()))
fit_htn_risk <- csSampling::cs_sampling(svydes = svy_des, mod_stan = mod_brms,
                                        data_stan = data_brms, 
                                        ctrl_stan = list(chains = 3, 
                                        iter = 2000, warmup = 1000, thin = 5))
out_stan <- do.call(sampling, c(list(object = mod_brms, data = data_brms, 
                                     chains = 3, iter = 4000, 
                                     warmup = 2000, thin = 5), sampling_args = list()))
# Just weighted brms version
brms_out <- brms::brm(formula = brms_mod, data = svy_data_unaware_t2d, 
                      family = family)
prior_summary(brms_out)
summary(brms_out)


# Try unweighted brms
formula_unwt <- as.formula(paste0("diabetes ~ dbh_class ",
                                  "+ Age + Sex + Educ+ Inc_hh + Urban",
                                  " + Physical_activity + Smoking_status ",
                                  "+ Food_security")) 
brms_unwt <- brms::brm(formula = brmsformula(formula_unwt, center = TRUE), 
                       data = svy_data_unaware_t2d, family = family)
prior_summary(brms_unwt)
summary(brms_unwt)
#==================== Run WOLCAN using NAs in DBH ==============================

# Convert NA to level 4
dbh_4_levels <- as.data.frame(apply(prospect_dbh, 2, 
                                    function(x) ifelse(is.na(x), 4, x)))
# Restrict to DBH cases
restrict_data <- prospect_cleaned %>% 
  filter(studyid %in% dbh_4_levels$studyid)

### Define data and parameters
# Estimating pseudo-weights
x_mat <- as.matrix(dbh_4_levels[, -1])  # Multivariate categorical variables
selection_covs <- c("Sex", "Educ", "Age", "Inc_hh", "Urban")
dat_B <- restrict_data %>%  # Covariates for NPS
  select(all_of(selection_covs)) %>%
  mutate_at(c("Sex", "Educ", "Inc_hh", "Urban"), as.factor)
nrow(dat_B %>% drop_na()) # 414 dropped
dat_R <- prcs_cleaned %>%  # Covariates for reference
  select(all_of(selection_covs)) %>%
  mutate_at(c("Sex", "Educ", "Inc_hh", "Urban"), as.factor)
nrow(dat_R %>% drop_na()) # 6185 dropped
pred_model <- "bart"  # Prediction model for selection
pred_covs_B <- selection_covs  # Covariates to predict NPS selection
pred_covs_R <- selection_covs  # Covariates to predict RS selection
pi_R <- 1 / prcs_cleaned$Weight  # RS selection probabilites for those in RS
# Set pi_R = 1 to 0.999999 to prevent Inf/Nan in BART
pi_R <- ifelse(pi_R == 1, 0.999999, pi_R)
hist(pi_R, breaks = 30)
summary(pi_R)
hat_pi_R <- NULL  # RS selection probabilities for those in NPS
num_post <- 1000  # Number of posterior BART draws for estimating weights
frame_B <- 1  # Coverage probability of NPS frame
frame_R <- 1  # Coverage probability of RS frame
trim_method = "t2"  # Trimming method using IQR
trim_c = 20         # Trimming constant

# Model estimation
D <- 20            # Number of sets of MI pseudo-weights
parallel <- FALSE   # Whether to parallelize for MI
n_cores <- 4       # Number of cores to use for parallelization
wts_adj <- "MI"    # Adjustment method to account for weights uncertainty
adjust <- TRUE     # Whether to adjust variance for pseudo-likelihood
tol <- 1e-8        # Underflow tolerance
num_reps <- 100    # Number of bootstrap replicates for WS adjustment
run_adapt <- TRUE  # Whether to run adaptive sampler to get K
K_max <- 30        # Max number of latent classes for adaptive sampler
adapt_seed <- 1    # Seed for adaptive sampler
fixed_seed <- 1    # Seed for fixed sampler
n_runs <- 20000    # Number of MCMC iterations
burn <- 10000      # Burn-in
thin <- 5          # Thinning
update <- 5000     # Display update frequency
class_cutoff <- 0.05


# Save as complete case DBH
save_res <- TRUE
save_path <- paste0(wd, res_dir, "cat4_dbh")
# Source R model functions
source(paste0(wd, "Model_Code/", "model_functions.R"))

### Run weighted model
res <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R, 
              pred_model = pred_model, pred_covs_B = pred_covs_B, 
              pred_covs_R = pred_covs_R, 
              pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post, 
              frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
              trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
              wts_adj = wts_adj, adjust = adjust, tol = tol, 
              num_reps = num_reps, run_adapt = run_adapt, 
              K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
              fixed_seed = fixed_seed, class_cutoff = class_cutoff, 
              n_runs = n_runs, burn = burn, thin = thin, update = update, 
              save_res = save_res, save_path = save_path, 
              overwrite = TRUE)

#================= Troubleshooting =============================================

M <- floor(n_runs / thin) - floor(burn / thin) 
n <- nrow(x_mat)
J <- ncol(x_mat)
# Keep only the draws with K_red = K_fixed
K_red_draws_all <- sapply(1:D, function(d) res_all[[d]]$K_red)
keep_d <- which(K_red_draws_all == K)
D_red <- length(keep_d)
# Initialize parameter and variables across draws
K_red_draws <- numeric(D_red)  # restricted to draws with K_red = K_fixed
pi_red_draws <- vector("list", length = D_red)
theta_red_draws <- vector("list", length = D_red)
c_all_draws <- matrix(NA, nrow = D_red, ncol = n)
wts_draws <- matrix(NA, nrow = n, ncol = D_red)
messages_draws <- vector("list", length = D_red)  # var_adjust error messages  
# Extract results for draws with okay K
for (d in 1:D_red) {
  # for (d in 1:D) {
  # Extract results
  est_d <- res_all[[keep_d[d]]]
  
  K_red_draws[d] <- est_d$K_red
  c_all_draws[d, ] <- est_d$c_all
  pi_red_draws[[d]] <- est_d$pi_red
  theta_red_draws[[d]] <- est_d$theta_red
  wts_draws[, d] <- est_d$wts
  messages_draws[[d]] <- ifelse(is.null(est_d$message), "no error", est_d$message)
}


#=============== Try weights with glm ==========================================
# Using full PRCS data
wts_glm <- get_weights_bart(dat_B = dat_B, dat_R = dat_R %>% 
                              mutate(Educ = relevel(Educ, ref = "3")), 
                                pred_covs_B = pred_covs_B, 
                                pred_covs_R = pred_covs_R, 
                                num_post = num_post, pi_R = pi_R, 
                                hat_pi_R = hat_pi_R,
                                frame_B = frame_B, frame_R = frame_R,
                                trim_method = trim_method, trim_c = trim_c,
                                pred_model = "glm")
summary(wts_glm$fit_pi_B)
hist(wts_glm$wts, breaks = 30)
wts_bart <- get_weights_bart(dat_B = dat_B, dat_R = dat_R %>% 
                              mutate(Educ = relevel(Educ, ref = "3")), 
                            pred_covs_B = pred_covs_B, 
                            pred_covs_R = pred_covs_R, 
                            num_post = num_post, pi_R = pi_R, 
                            hat_pi_R = hat_pi_R,
                            frame_B = frame_B, frame_R = frame_R,
                            trim_method = trim_method, trim_c = trim_c,
                            pred_model = "bart")
summary(wts_bart$wts)
hist(wts_bart$wts, breaks = 30)

# blue
hist(full_wts$wts, breaks = 50, xlim = c(0, 11000), ylim = c(0, 1100), 
     col = rgb(0,0,1,1/5))
# red
hist(age_wts$wts, breaks = 50, xlim = c(0, 11000), col=rgb(1,0,0,1/5), add = TRUE)
