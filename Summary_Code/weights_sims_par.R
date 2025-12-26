#===============================================================================
# Run and summarize scenarios for predicting pseudo-weights: parallel version
# Author: Stephanie Wu
# Date created: 2025/11/19
# Date updated: 2025/11/19
#===============================================================================

# Read in two command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_scen <- args[[1]]   # Data scenario
model_scen_num <- args[[2]]  # Model type

# Data scenarios:
# High overlap: 0, 8, 13, 14, 24
# Low overlap: 10, 19, 20, 21, 25
scenario <- data_scen

# Specify model type
if (model_scen_num == 1) {
  model_scen <- list("none")           # No weights (weights all set to 1)
} else if (model_scen_num == 2) {
  model_scen <- list("logistic_true")  # Log reg w/ true model covariate forms
} else if (model_scen_num == 3) { 
  model_scen <- list("logistic")       # Log reg w/ main effects and basic interaction
} else if (model_scen_num == 4) {
  model_scen <- list("bart_500")       # BART 500 samples w/ main effects and basic interaction
} else if (model_scen_num == 5) {
  model_scen <- list("bart_1000")      # BART 1000 samples w/ main effects and basic interaction
} else if (model_scen_num == 6) {
  model_scen <- list("bart_2000")      # BART 2000 samples w/ main effects and basic interaction
} else if (model_scen_num == 7) {
  model_scen <- list("logistic_cov")   # Log reg w/ main effects and missing covs
} else if (model_scen_num == 8) {
  model_scen <- list("bart_1000_cov")  # BART 1000 samples w/ main effects and missing covs
} else if (model_scen_num == 9) {      # All models
  model_scen <- list("none", "logistic_true", "logistic", 
                     "bart_500", "bart_1000", "bart_2000", 
                     "logistic_cov", "bart_1000_cov")
} 

# Load libraries
library(dplyr)
library(gridExtra)
library(knitr)
library(kableExtra)
library(BART)  # BART

# Set directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
wd <- "/n/netscratch/stephenson_lab/Lab/stephwu18/WOLCAN/"
data_dir <- "Data/"    # Data directory
res_dir <- "Results/"  # Results directory
code_dir <- "Summary_Code/"  # Model code directory
sum_dir <- "Summary_Results/"  # Summary results directory

# Source R model functions
source(paste0(wd, "Model_Code/model_functions.R"))

#============== Run and store weights ==========================================

### Default settings
hat_pi_R <- NULL  # RS selection probabilities for those in NPS
frame_B <- 1  # Coverage probability of NPS frame
frame_R <- 1  # Coverage probability of RS frame
trim_method = "t2"  # Trimming method using IQR
trim_c = 20         # Trimming constant
ci_level <- 0.95

### Models and scenarios

# Total number of scenarios (specified scenario, so just number of models)
num_scen <- length(model_scen)
num_samps <- 100 # Number of samples per scenario

# Initialize results list for all models
wts_res_all_models <- vector(mode = "list", length = num_scen)
names(wts_res_all_models) <- unlist(model_scen)

    # Initialize results data frame
    # weights_res <- as.data.frame(matrix(NA, nrow = num_scen, ncol = num_samps))
    # pi_res <- as.data.frame(matrix(NA, nrow = num_scen, ncol = num_samps))
    # se_w_B_mean_res <- as.data.frame(matrix(NA, nrow = num_scen, ncol = num_samps))
    # ci_width_w_B_mean_res <- as.data.frame(matrix(NA, nrow = num_scen, ncol = num_samps))
    # se_pi_mean_res <- as.data.frame(matrix(NA, nrow = num_scen, ncol = num_samps))
    # ci_width_pi_mean_res <- as.data.frame(matrix(NA, nrow = num_scen, ncol = num_samps))
    # rownames(weights_res) <- rownames(pi_res) <- rownames(se_w_B_mean_res) <-
    #   rownames(ci_width_w_B_mean_res) <- rownames(se_pi_mean_res) <- 
    #   rownames(ci_width_pi_mean_res) <- model_names


# Iterate over the model types
for (j in 1:length(model_scen)) {  # Prediction model for weights
  model <- model_scen[[j]]
  wts_res <- as.data.frame(matrix(NA, nrow = 6, ncol = num_samps))
  rownames(wts_res) <- c("weights", "pi_B", "se_w_B_mean", "ci_width_w_B_mean", 
                         "se_pi_B_mean", "ci_width_pi_B_mean")
  colnames(wts_res) <- 1:num_samps
  
  if (scenario %in% c(24, 25)) {  # Application-based weights (supplement)
    # Covariates to include in prediction model
    if (model == "logistic_true") {
      # True covariate forms
      pred_covs_B <- c("age_cent", "age_cent2", "sex_f_num", "educ_num", 
                       "educ_num2", "hhinc_num_low", "hhinc_num_high",
                        "ethn_num", "urban", "disability", 
                        "sex_f_num_urban", "age_cent_disability")
      pred_covs_R <- c("age_cent", "age_cent2", "sex_f_num", 
                       "educ_num", "hhinc_num_low", "hhinc_num_high",
                       "ethn_num", "urban", "disability", 
                       "sex_f_num_urban", "age_cent_disability")
    } else if (model %in% c("bart_1000_cov", "logistic_cov")) {
      # Only main effects, and missing urban and disability
      pred_covs_B <- pred_covs_R <- c("age", "sex_f_num", "educ_num", 
                                      "hhinc_num_low", "hhinc_num_high", 
                                      "ethn_num") 
    } else {
      # Only main effects
      pred_covs_B <- pred_covs_R <- c("age", "sex_f_num", "educ_num", 
                                      "hhinc_num_low", "hhinc_num_high",
                                      "ethn_num", "urban", "disability") 
    }
    
  } else {  # Default simulation settings
    # Covariates to include in prediction model
    if (model == "logistic_true") {
      # True covariate forms
      pred_covs_B <- pred_covs_R <- c("A1", "A12", "A2", "logA2", "A1A2", 
                                      "sinA1A2", "A3")
    } else if (model %in% c("bart_1000_cov", "logistic_cov")) {
      # Missing A3 and only main effects of A1 and A2
      pred_covs_B <- pred_covs_R <- c("A1", "A2")
    } else {
      # Main effects for A1, A2, A3, and interaction A1:A2
      pred_covs_B <- pred_covs_R <- c("A1", "A2", "A1A2", "A3")
    }
  }
  
  # Number of BART samples to use
  if (model == "bart_500") {
    num_post <- 500
  } else if (model == "bart_2000") {
    num_post <- 2000
  } else {
    num_post <- 1000
  }
    
  # Type of prediction model
  if (model %in% c("logistic_true", "logistic", "logistic_cov")) {  # Logistic regression
    pred_model <- "glm"
  } else {  # BART model
    pred_model <- "bart"
  }
  
  # Get weights for each sample
  for (k in 1:num_samps) {
    samp_i <- k
    # Load sample data
    load(paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", samp_i, "_B_wolcan.RData"))
    load(paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", samp_i, "_R_wolcan.RData"))
    
    # Other variables
    x_mat <- sim_samp_B$X_data  # Multivariate categorical variables
    dat_B <- data.frame(sim_samp_B$covs)  # Covariates for NPS
    dat_R <- data.frame(sim_samp_R$covs)  # Covariates for reference
    pi_R <- sim_samp_R$pi_R  # RS selection probabilites for those in RS
    
    # Prepare sample data for more complex variables in true formula for 
    # application-based weights
    if (scenario %in% c(24, 25)) {
      dat_B$age_cent <- dat_B$age - 53
      dat_B$age_cent2 <- dat_B$age_cent^2
      dat_B$educ_num2 <- dat_B$educ_num^2  # Simulated data treats as numeric
      dat_B$hhinc_num_low <- ifelse(dat_B$hhinc_num == 0, 1, 0)
      dat_B$hhinc_num_high <- ifelse(dat_B$hhinc_num == 2, 1, 0)
      dat_B$sex_f_num_urban <- dat_B$sex_f_num * dat_B$urban
      dat_B$age_cent_disability <- dat_B$age_cent * dat_B$disability
      
      dat_R$age_cent <- dat_R$age - 53
      dat_R$age_cent2 <- dat_R$age_cent^2
      dat_R$educ_num2 <- dat_R$educ_num^2  # Simulated data treats as numeric
      dat_R$hhinc_num_low <- ifelse(dat_R$hhinc_num == 0, 1, 0)
      dat_R$hhinc_num_high <- ifelse(dat_R$hhinc_num == 2, 1, 0)      
      dat_R$sex_f_num_urban <- dat_R$sex_f_num * dat_R$urban
      dat_R$age_cent_disability <- dat_R$age_cent * dat_R$disability
    }
    
    # set seed
    set.seed(samp_i)  
    
    if (model == "none") {  # no weights
      wts <- rep(1, nrow(x_mat))
      res <- mean(abs(wts - 1/(sim_samp_B$true_pi_B)))
      pi <- mean(abs(1/wts - sim_samp_B$true_pi_B))
      se_w_B_mean <- ci_width_w_B_mean <- se_pi_mean <- ci_width_pi_mean <- 1
      
    } else {
      # Get weights
      invisible(capture.output(
        est_weights <- get_weights_bart(dat_B = dat_B, dat_R = dat_R, 
                                        pred_covs_B = pred_covs_B, 
                                        pred_covs_R = pred_covs_R, 
                                        num_post = num_post, pi_R = pi_R, 
                                        hat_pi_R = hat_pi_R,
                                        frame_B = frame_B, frame_R = frame_R,
                                        trim_method = trim_method, trim_c = trim_c,
                                        pred_model = pred_model, 
                                        ci_level = ci_level)))
      
      # Mean absolute difference of weights compared with true
      res <- mean(abs(est_weights$wts - (1/sim_samp_B$true_pi_B)))
      # Mean absolute different of inclusion probabilities compared with true
      pi <- mean(abs(est_weights$hat_pi_B - sim_samp_B$true_pi_B))
      
      # SE and CI-width for weights
      se_w_B_mean <- mean(est_weights$se_w_B)
      ci_width_w_B_mean <- mean(est_weights$w_B_high - est_weights$w_B_low)
      
      # SE and CI-width for inclusion probabilities
      se_pi_mean <- mean(est_weights$se_pi_B)
      ci_width_pi_mean <- mean(est_weights$pi_B_high - est_weights$pi_B_low)
    }
    
    wts_res[, k] <- c(res, pi, se_w_B_mean, ci_width_w_B_mean, se_pi_mean,
                      ci_width_pi_mean)
        # weights_res[counter, k] <- res
        # pi_res[counter, k] <- pi
        # se_w_B_mean_res[counter, k] <- se_w_B_mean
        # ci_width_w_B_mean_res[counter, k] <- ci_width_w_B_mean
        # se_pi_mean_res[counter, k] <- se_pi_mean
        # ci_width_pi_mean_res[counter, k] <- ci_width_pi_mean
  }
  
  # Add scenario results to list of all model scenarios
  wts_res_all_models[[j]] <- wts_res
  
  # Save specific model results
  save(wts_res, file = paste0(wd, sum_dir, "scen_", scenario, "/model_", model, 
                                  "_wts_res.RData"))
  
      # temp2 <- pi_res[counter, ]
      # save(temp2, file = paste0(wd, sum_dir, "scen_", scenario, "/model_", model, 
      #                           "_pi.RData"))
      # temp3 <- se_w_B_mean_res[counter, ]
      # save(temp3, file = paste0(wd, sum_dir, "scen_", scenario, "/model_", model, 
      #                           "_se_w_B_mean.RData"))
      # temp4 <- ci_width_w_B_mean_res[counter, ]
      # save(temp4, file = paste0(wd, sum_dir, "scen_", scenario, "/model_", model, 
      #                           "_ci_width_w_B_mean.RData"))
      # temp5 <- se_pi_mean_res[counter, ]
      # save(temp5, file = paste0(wd, sum_dir, "scen_", scenario, "/model_", model, 
      #                           "_se_pi_mean.RData"))
      # temp6 <- ci_width_pi_mean_res[counter, ]
      # save(temp6, file = paste0(wd, sum_dir, "scen_", scenario, "/model_", model, 
      #                           "_ci_width_pi_mean.RData"))
  
  # Print progress
  print(paste0("Scenario ", scenario, " model ", model, " done!"))
}

# If running multiple models, save all model results for the scenario together
if (model_scen_num == 9) {  # All models
  
  # Save all model results together
  save(wts_res_all_models, 
       file = paste0(wd, sum_dir, "scen_", scenario, "/wts_res_all_models.RData"))
  
      # rownames(weights_res) <- unlist(model_scen)
      # rownames(pi_res) <- unlist(model_scen)
      # rownames(se_w_B_mean_res) <- unlist(model_scen)
      # rownames(ci_width_w_B_mean_res) <- unlist(model_scen)
      # rownames(se_pi_mean_res) <- unlist(model_scen)
      # rownames(ci_width_pi_mean_res) <- unlist(model_scen)
      # 
      # save(weights_res, file = paste0(wd, sum_dir, "scen_", scenario, "/weights_res_all.RData"))
      # save(pi_res, file = paste0(wd, sum_dir, "scen_", scenario, "/pi_res_all.RData"))
      # save(se_w_B_mean_res, file = paste0(wd, sum_dir, "scen_", scenario, "/se_w_B_mean_res_all.RData"))
      # save(ci_width_w_B_mean_res, file = paste0(wd, sum_dir, "scen_", scenario, "/ci_width_w_B_mean_res_all.RData"))
      # save(se_pi_mean_res, file = paste0(wd, sum_dir, "scen_", scenario, "/se_pi_mean_res_all.RData"))
      # save(ci_width_pi_mean_res, file = paste0(wd, sum_dir, "scen_", scenario, "/ci_width_pi_mean_res_all.RData"))
}

