#=================================================
# Create pseudo-weights and run WOLCAN model
# Author: Stephanie Wu
# Date created: 2024/04/16
# Date updated: 2024/07/10
#=================================================

# Read in two command line arguments
args <- commandArgs(trailingOnly = TRUE)
scenario <- args[[1]]   # Simulation scenario
samp_i <- args[[2]]     # Sample number
samp_i <- as.numeric(samp_i)

# Whether or not baysc package has been installed. If FALSE, local functions 
# will be used instead
baysc_package = TRUE

# Load libraries
library(BART)  # BART
library(abind)  # binding arrays
library(parallel)  # parallel processing
library(plyr)
library(dplyr)
library(LaplacesDemon)  # distributions
library(truncnorm)
library(fastDummies)
library(matrixStats)
library(Matrix)
library(gtools)
library(e1071)
library(rstan)  # stan
library(survey)  # survey functions
library(Rcpp)  # Rcpp interface
library(RcppArmadillo)
library(RcppTN)

# Set directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/"
data_dir <- "Data/"    # Data directory
res_dir <- "Results/"  # Results directory
code_dir <- "Model_Code/"  # Model code directory

# Create scenario results folder if it doesn't exist
dir_path <- paste0(wd, res_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

# Define path to save results
save_path <- paste0(wd, res_dir, "scen_", scenario, "/samp_", samp_i)

# Check if results already exist
already_done <- file.exists(paste0(save_path, "_wolcan_results.RData"))
if (already_done) {
  print(paste0('Scenario ', scenario, ' samp ', samp_i, ' already exists.'))
} else {
  
  print(paste0('Running scenario ', scenario, ' samp ', samp_i, '...'))
  
  # Source R model functions
  source(paste0(wd, code_dir, "model_functions.R"))
  
  # Source baysc functions
  if (baysc_package) {
    library(baysc)
    mod_stan <- NULL
  } else {
    # Source Rcpp functions from baysc package
    Rcpp::sourceCpp(paste0(wd, code_dir, "baysc_functions/", "main_Rcpp_functions.cpp"))
    
    # Source functions from baysc package
    file_list <- c("mcmc_functions_wolca.R", "mcmc_functions.R", "utilities.R", 
                   "wolca_var_adjust.R", "wolca.R", "wolca_svyglm.R")
    invisible(lapply(file_list, function(x) 
      source(paste0(wd, code_dir, "baysc_functions/", x))))
    
    # Source stan models
    mod_stan <- stan_model(paste0(wd, code_dir, "baysc_functions/", "WOLCA_main.stan"))
  }
  
  ### Read in population and sample data
  load(paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData"))
  load(paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", samp_i, "_B_wolcan.RData"))
  load(paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", samp_i, "_R_wolcan.RData"))
    
  # # Number of overlapping individuals
  # length(intersect(sim_samp_R$ind_R, sim_samp_B$ind_B))
  
  ### Define data and parameters
  # Estimating pseudo-weights
  x_mat <- sim_samp_B$X_data  # Multivariate categorical variables
  dat_B <- data.frame(sim_samp_B$covs)  # Covariates for NPS
  dat_R <- data.frame(sim_samp_R$covs)  # Covariates for reference
  pred_covs_B <- c("A1", "A2", "A1A2", "A3")  # Covariates to predict NPS selection
  pred_covs_R <- c("A1", "A2", "A1A2", "A3")  # Covariates to predict RS selection
  pi_R <- sim_samp_R$pi_R  # RS selection probabilites for those in RS
  hat_pi_R <- NULL  # RS selection probabilities for those in NPS
  num_post <- 1000  # Number of posterior BART draws for estimating weights
  frame_B <- 1  # Coverage probability of NPS frame
  frame_R <- 1  # Coverage probability of RS frame
  trim_method = "t2"  # Trimming method using IQR
  trim_c = 20         # Trimming constant
  
  # Model estimation
  D <- 20            # Number of sets of MI pseudo-weights
  parallel <- TRUE   # Whether to parallelize for MI
  n_cores <- 8       # Number of cores to use for parallelization
  wts_adj <- "MI"    # Adjustment method to account for weights uncertainty
  adjust <- TRUE     # Whether to adjust variance for pseudo-likelihood
  tol <- 1e-8        # Underflow tolerance
  num_reps <- 100    # Number of bootstrap replicates for WS adjustment
  run_adapt <- TRUE  # Whether to run adaptive sampler to get K
  K_max <- 30          # Max number of latent classes for adaptive sampler
  adapt_seed <- samp_i # Seed for adaptive sampler
  fixed_seed <- samp_i # Seed for fixed sampler
  n_runs <- 20000    # Number of MCMC iterations
  burn <- 10000      # Burn-in
  thin <- 5          # Thinning
  update <- 5000     # Display update frequency
  
  # # For tests
  # n_runs <- 1000
  # burn <- 500
  # thin <- 5
  # update <- 500
  
  ### Modifications based on scenario
  if (scenario == 2) {  # WS one-step adjustment with draws
    wts_adj = "WS all"
  } else if (scenario == 3) {  # WS one-step adjustment with draws
    wts_adj = "WS mean"
  } else if (scenario == 4) {  # No pseudo-likelihood adjustment
    adjust = FALSE
  } else if (scenario == 5) {  # No propagation of weights uncertainty
    wts_adj = "none"
  } else if (scenario == 14) {  # RS weights known for NPS
    hat_pi_R = sim_samp_B$true_pi_R
  } 
  
  ### Run weighted model
  res <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R, 
                pred_covs_B = pred_covs_B, pred_covs_R = pred_covs_R, 
                pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post, 
                frame_B = frame_B, frame_R = frame_R, trim_method = trim_method, 
                trim_c = trim_c, D = D, parallel = parallel, n_cores = n_cores,
                wts_adj = wts_adj, adjust = adjust, tol = tol, 
                num_reps = num_reps, run_adapt = run_adapt, 
                K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL, 
                fixed_seed = fixed_seed, class_cutoff = 0.05, 
                n_runs = n_runs, burn = burn, thin = thin, update = update, 
                save_res = TRUE, save_path = save_path)
  
  ### Run unweighted model
  if (scenario %in% c(0, 1, 7:11)) {
    res_unwt <- wolca(x_mat = x_mat, sampling_wt = rep(1, nrow(x_mat)), 
                             run_sampler = "both", K_max = K_max, 
                             adapt_seed = adapt_seed, class_cutoff = 0.05, 
                             n_runs = n_runs, burn = burn, thin = thin, 
                             update = update, save_res = TRUE, 
                             save_path = save_path)
  }
}


# #================== Local version ==============================================
# # Don't need to run if using cluster!
# # Run model locally for all samples
# num_samps <- 100
# for (samp_i in 1:num_samps) {
#   
#   # Define path to save results
#   save_path <- paste0(wd, res_dir, "scen_", scenario, "/samp_", samp_i)
#   
#   print(paste0('Running scenario ', scenario, ' samp ', samp_i, '...'))
# 
#   # Source functions
#   source(paste0(wd, code_dir, "model_functions.R"))
# 
#   ### Read in population and sample data
#   load(paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData"))
#   load(paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", samp_i, "_B_wolcan.RData"))
#   load(paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", samp_i, "_R_wolcan.RData"))
# 
#   # # Number of overlapping individuals
#   # length(intersect(sim_samp_R$ind_R, sim_samp_B$ind_B))
# 
#   ### Define data and parameters
#   # Estimating pseudo-weights
#   x_mat <- sim_samp_B$X_data  # Multivariate categorical variables
#   dat_B <- data.frame(sim_samp_B$covs)  # Covariates for NPS
#   dat_R <- data.frame(sim_samp_R$covs)  # Covariates for reference
#   pred_covs_B <- c("A1", "A2", "A1A2")  # Covariates to predict NPS selection
#   pred_covs_R <- c("A1", "A2", "A1A2")  # Covariates to predict RS selection
#   pi_R <- sim_samp_R$pi_R  # RS selection probabilites for those in RS
#   hat_pi_R <- NULL  # RS selection probabilities for those in NPS
#   num_post <- 1000  # Number of posterior BART draws for estimating weights
#   frame_B <- 1  # Coverage probability of NPS frame
#   frame_R <- 1  # Coverage probability of RS frame
# 
#   # Model estimation
#   D <- 10            # Number of sets of MI pseudo-weights
#   parallel <- TRUE   # Whether to parallelize for MI
#   MI <- TRUE         # Whether to run MI procedure for variance estimation
#   adjust <- TRUE     # Whether to adjust variance for pseudo-likelihood
#   tol <- 1e-8        # Underflow tolerance
#   run_adapt <- TRUE  # Whether to run adaptive sampler to get K
#   K_max <- 30          # Max number of latent classes for adaptive sampler
#   adapt_seed <- samp_i # Seed for adaptive sampler
#   fixed_seed <- samp_i # Seed for fixed sampler
#   n_runs <- 10000    # Number of MCMC iterations
#   burn <- 5000      # Burn-in
#   thin <- 5          # Thinning
#   update <- 1000     # Display update frequency
#   
#   # # For tests
#   # n_runs <- 1000
#   # burn <- 500
#   # thin <- 5
#   # update <- 500
# 
#   ### Modifications based on scenario
#   if (scenario == 2) {  # No propagation of weights uncertainty
#     MI <- FALSE
#   }
#   if (scenario == 3) {  # No pseudo-likelihood adjustment
#     adjust = FALSE
#   }
#   if (scenario == 4) {  # RS weights known for NPS
#     hat_pi_R = sim_samp_B$true_pi_R
#   }
# 
#   ### Run weighted model
#   res <- wolcan(x_mat = x_mat, dat_B = dat_B, dat_R = dat_R,
#                 pred_covs_B = pred_covs_B, pred_covs_R = pred_covs_R,
#                 pi_R = pi_R, hat_pi_R = hat_pi_R, num_post = num_post,
#                 frame_B = frame_B, frame_R = frame_R, D = D, parallel = parallel,
#                 MI = MI, adjust = adjust, tol = tol, run_adapt = run_adapt,
#                 K_max = K_max, adapt_seed = adapt_seed, K_fixed = NULL,
#                 fixed_seed = fixed_seed, class_cutoff = 0.05,
#                 n_runs = n_runs, burn = burn, thin = thin, update = update,
#                 save_res = TRUE, save_path = save_path)
# 
#   ### Run unweighted model
#   if (scenario %in% c(1, 7:11)) {
#     res_unwt <- wolca(x_mat = x_mat, sampling_wt = rep(1, nrow(x_mat)),
#                              run_sampler = "both", K_max = K_max,
#                              adapt_seed = adapt_seed, class_cutoff = 0.05,
#                              n_runs = n_runs, burn = burn, thin = thin,
#                              update = update, save_res = TRUE,
#                              save_path = save_path)
#   }
# }


