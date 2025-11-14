#=================================================
# Run outcome model simulations
# Author: Stephanie Wu
# Date created: 2024/09/02
# Date updated: 2024/09/02
#=================================================

rm(list = ls())

library(baysc)
library(abind)  # binding arrays
library(flextable)
library(dplyr)
library(bayesplot)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(knitr)
library(kableExtra)
library(gt)
library(parallel)
library(stringr)
library(Hmisc)  # plotting
library(ggbeeswarm)  # beeswarm plot
library(survey)
library(rstan)

# Set directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
wd <- "/n/netscratch/stephenson_lab/Lab/stephwu18/WOLCAN/"
data_dir <- "Data/"    # Data directory
res_dir <- "Results/"  # Results directory
code_dir <- "Summary_Code/"  # Model code directory
sum_dir <- "Summary_Results/"  # Summary results directory

# Source summary functions
source(paste0(wd, code_dir, "summary_functions.R"))

# rstan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#rstan_options(javascript = FALSE)

#================== Run simulations and summary for outcome model ==============
# library(survey)
# pop$Y <- Y_data
# 
# test <- glm(as.formula(Y ~ c_all + A1 + A2 + c_all:A1), data = pop, 
#             family = binomial(link = "logit"))
# round(summary(test)$coefficients[, 1], 3)
# xi_vec_y
# 
# samp_data <- data.frame(c_all = sim_samp_B$c_all, Y = sim_samp_B$Y_data, 
#                         sim_samp_B$covs, wts = 1/sim_samp_B$true_pi_B, 
#                         clus = 1:nrow(sim_samp_B$covs))
# svy_des <- survey::svydesign(ids = ~clus, weights = ~wts, data = samp_data)
# svy_test <- survey::svyglm(formula = as.formula(Y ~ c_all + A1 + A2 + c_all:A1),
#                            design = svy_des, 
#                            family = stats::binomial(link = "logit"))
# round(summary(svy_test)$coefficients[, 1], 3)



#formula_y <- "y_all ~ c_all + A1 + A2 + c_all:A1"
#formula_y <- "y_all ~ c_all"
formula_y <- "y_all ~ c_all + A1 + A2"
dist_type = "mean_abs"
samp_i_seq <- 1:100
#scenarios <- c(0, 6:10, 15, 18:19)
scenarios <- 0  # where model results are
scenarios_data <- 16  # where data is
model <- "wolcan"
save <- FALSE
subset <- TRUE
informative <- FALSE

xi_bias_all <- xi_var_all <- xi_cover_all <- matrix(NA, nrow = length(scenarios), 
                                                    ncol = length(samp_i_seq))
for (i in 1:length(scenarios)) {
  scenario <- scenarios[i]
  scenario_data <- scenarios_data[i]
  print(paste0("Starting scenario ", scenario, "..."))
  # Load simulated population data
  pop_data_path <- paste0(wd, data_dir, "scen_", scenario_data, "/sim_pop_wolcan.RData")
  load(pop_data_path)
  true_params <- get_true_params_wolcan(sim_pop = sim_pop)
  true_K <- length(true_params$true_pi)
  
  # True parameters
  # pop <- data.frame(c_all = sim_pop$c_all, sim_pop$pop, y_all = sim_pop$Y_data)
  # pop_logreg <- glm(as.formula(formula_y), data = pop, 
  #                               family = binomial(link = "logit"))
  # xi_vec_y <- pop_logreg$coefficients
  xi_vec_y <- sim_pop$true_xi
  
  # Get bias, variance, and coverage for each sample
  for (j in 1:length(samp_i_seq)) {
    samp_i <- samp_i_seq[j]
    if (samp_i %% 10 == 0) {
      print(samp_i)
    }
    # Load data
    # Check that sample data file exists
    sim_data_path <- paste0(wd, data_dir, "scen_", scenario_data, "/sim_samp_", 
                            samp_i, "_B_wolcan.RData")
    load(sim_data_path)
    
    # Load model output
    sim_res_path <- paste0(wd, res_dir, "scen_", scenario, "/samp_", samp_i, "_", 
                           model, "_results.RData")
    
    if (!file.exists(sim_res_path)) {
      print(paste0("File does not exist: ", sim_res_path))
    } else {
      load(sim_res_path)
      
      if (informative) {  # New weights for informative sampling
        
        if (model == "wolcan") {
          source(paste0(wd, "Model_Code/model_functions.R"))
          load(paste0(wd, data_dir, "scen_", scenario_data, "/sim_samp_",
                      samp_i, "_R_wolcan.RData"))
          # Re-estimate weights w/ y
          est_weights <- get_weights_bart(dat_B = data.frame(sim_samp_B$covs, 
                                                             y_all = sim_samp_B$Y_data), 
                                          dat_R = data.frame(sim_samp_R$covs, 
                                                             y_all = sim_samp_R$Y_data), 
                                          pred_model = "bart",
                                          pred_covs_B = c("A1", "A2", "A1A2", "A3", "y_all"), 
                                          pred_covs_R = c("A1", "A2", "A1A2", "A3", "y_all"), 
                                          num_post = 1000, pi_R = sim_samp_R$pi_R , 
                                          frame_B = 1, frame_R = 1,
                                          trim_method = "t2", trim_c = 20)
          # Mean estimated weights
          wts <- est_weights$wts
          # Posterior distribution of weights
          wts_post <- est_weights$wts_post
          res$data_vars$sampling_wt <- wts
          
           
          # res$estimates_adjust$c_all <- sim_samp_B$c_all
          
        } else if (model == "wolca") {
          # res$estimates$c_all <- sim_samp_B$c_all
        }
      }
      
      # Use true c_all
      true_c_all <- sim_samp_B$c_all
    
      y_metrics_i <- get_y_metrics_i(res = res, true_xi = sim_pop$true_xi, 
                                     true_params = true_params, 
                                     true_K = true_K, model = model,
                                     dist_type = dist_type, sim_samp_B = sim_samp_B, 
                                     subset = subset, method = "svyglm", 
                                     true_c_all = true_c_all)
    }
    
    xi_bias_all[i, j] <- y_metrics_i$xi_dist
    xi_var_all[i, j] <- y_metrics_i$xi_var_mean
    xi_cover_all[i, j] <- y_metrics_i$xi_cover_mean
  }
  
  if (save) {
    # Save results for each scenario
    scen_res <- list(xi_bias = xi_bias_all[i, ],
                     xi_var = xi_var_all[i, ],
                     xi_cover = xi_cover_all[i, ])
    dir_path <- paste0(wd, sum_dir, "scen_", scenario_data, "/")
    if (!dir.exists(dir_path)) {
      dir.create(file.path(dir_path))
    }
    save(scen_res, file = paste0(dir_path, "summary_y_", model, ".RData"))
  }
}

rowMeans(xi_bias_all, na.rm = TRUE)
rowMeans(xi_var_all, na.rm = TRUE)
rowMeans(xi_cover_all, na.rm = TRUE)

summary_all <- list(xi_bias_all = xi_bias_all, 
                    xi_var_all = xi_var_all,
                    xi_cover_all = xi_cover_all)
# save(summary_all, file = paste0(wd, sum_dir, "summary_y_", model, 
#                                 "_all_scens.RData"))



## Test code
formula_y <- "y_all ~ c_all + A1 + A2"
load(paste0(wd, data_dir, "scen_0/sim_samp_1_B_wolcan.RData"))
load(paste0(wd, data_dir, "scen_0/sim_samp_1_R_wolcan.RData"))
load(paste0(wd, res_dir, "scen_0/samp_1_wolcan_results.RData"))
n <- res$data_vars$n
y_all <- sim_samp_B$Y_data
w_all <- res$data_vars$w_all
svy_data <- data.frame(c_all = as.factor(res$estimates_adjust$c_all), 
                       y_all = y_all,
                       A1 = sim_samp_B$covs$A1,
                       A2 = sim_samp_B$covs$A2,
                       w_all = w_all, 
                       cluster_id = 1:n)
V <- model.matrix(as.formula(formula_y), data = svy_data)
q <- ncol(V)
data_stan <- list(n = n, q = q, y = y_all, V = V, weights = w_all)
par_stan <- c('xi')  # subset of parameters interested in
# Compile stan model
mod_stan <- rstan::stan_model(file = paste0(wd, "Model_Code/wtd_logistic.stan"))
out_stan <- rstan::sampling(object = mod_stan, data = data_stan, 
                            pars = par_stan)
out_stan <- rstan::sampling(object = mod_stan, data = data_stan, 
                            pars = par_stan, chains = 4, iter = 2000, 
                            warmup = 1000, thin = 5)



bayes_glm <- outcome_var_adj(D = 20, wts_post = res$estimates_adjust$wts_draws, 
                             num_reps = 100, save_res = FALSE, adjust_seed = 1,
                             res = res)

# Compile stan model
mod_stan <- rstan::stan_model(file = paste0(wd, "Model_Code/wtd_logistic.stan"))

# wts_post: (n1)xM posterior distribution of weights for individuals in the NPS, 
# with each column corresponding to a posterior draw
outcome_var_adj <- function(res, sim_samp_B, D = 20, wts_post, num_reps = 100, 
                        save_res = TRUE, save_path = NULL, 
                        adjust_seed = NULL, mod_stan) {
  
  # run stan model
  print("stan fitting")
  out_stan <- sampling(object = mod_stan)
  
  # Begin runtime tracker
  start_time <- Sys.time()
  
  # Set seed
  if (!is.null(adjust_seed)) {
    set.seed(adjust_seed)
  }
  
  #================= Extract dimensions and catch errors =======================
  # Check num_reps
  if ((num_reps %% 1 != 0) | num_reps < 1) {
    stop("num_reps must be a whole number greater than 0, recommended to be at least 50.
    More replicates will lead to more accurate results but will take longer to run.")
  }
  
  # Extract data elements into the global environment
  n <- res$data_vars$n
  y_all <- sim_samp_B$Y_data
  w_all <- res$data_vars$w_all
  
  # Survey data frame for specifying survey design
  # Clusters are simply individuals
  svy_data <- data.frame(c_all = as.factor(res$estimates_adjust$c_all), 
                         y_all = y_all,
                         A1 = sim_samp_B$covs$A1,
                         A2 = sim_samp_B$covs$A2,
                         w_all = w_all, 
                         cluster_id = 1:n)
  V <- model.matrix(as.formula(formula_y), data = svy_data)
  q <- ncol(V)
  
  #=============== Run variance adjustment for each set of weights =============
  print("Running variance adjustment")
  
  # Normalize posterior weight distribution to sum to the sample size
  w_all_post <- apply(wts_post, 2, function(wts_l) wts_l / (sum(wts_l) / n))
  
  # Get quantiles of medians of weights posterior
  col_meds <- c(apply(w_all_post, 2, median))
  cutoffs <- (seq(1, D, length.out = D) - 0.5) / D
  med_quants <- stats::quantile(x = col_meds, probs = cutoffs)
  
  # Initialize parameter and variables across draws
  pi_red_draws <- vector("list", length = D)
  theta_red_draws <- vector("list", length = D)
  wts_draws <- matrix(NA, nrow = n, ncol = D)
  
  for (d in 1:D) {  # For each posterior weight set
    print(paste0("Draw ", d))
    # Draw from weights posterior closest to quantile
    draw <- which.min(abs(col_meds - med_quants[d]))
    w_all <- c(w_all_post[, draw])
    wts_draws[, d] <- wts_post[, draw]
    
    #=============== Run Stan model ==============================================
    # Define data for Stan model
    data_stan <- list(n = n, q = q, y = y_all, V = V, weights = w_all)
    
    # Stan parameters of interest
    par_stan <- c('xi')  # subset of parameters interested in
    
    # Run Stan model
    # Stan will pass warnings from calling 0 chains, but will still create an 
    # out_stan object for the 'grad_log_prob()' method
    out_stan <- rstan::sampling(object = mod_stan, data = data_stan, 
                                pars = par_stan, chains = 4, iter = 2000, 
                                warmup = 1000, thin = 5)
    
    # Get posterior mean across all chains 
    par_samps <- extract(out_stan, pars = par_stan, permuted = FALSE)
    par_hat <- colMeans(par_samps, dim = 2)  # dim = 2 across chains
    
    #=============== Post-processing adjustment in unconstrained space ===========
    # Estimate Hessian
    H_hat <- -1*stats::optimHess(par_hat, 
                                 gr = function(x){rstan::grad_log_prob(out_stan, x)})
    
    ### Create survey replicates 
    
    # Specify survey design
    svydes <- survey::svydesign(ids = ~cluster_id, weights = ~w_all,
                                data = svy_data)
    # Create svrepdesign
    svyrep <- survey::as.svrepdesign(design = svydes, type = "mrbbootstrap",
                                     replicates = num_reps)
    # Get survey replicates and gradient for each replicate
    rep_temp <- survey::withReplicates(design = svyrep, theta = grad_par,
                                       stan_mod = mod_stan, stan_data = data_stan,
                                       par_stan = par_stan, u_pars = par_hat)
    # Empirical variance of gradient across replicates
    J_hat <- stats::vcov(rep_temp)
    
    # Compute adjustment
    H_inv <- solve(H_hat)
    V1 <- H_inv %*% J_hat %*% H_inv
    
    # Check for issues with negative diagonals
    if (min(diag(V1)) < 0) {
      print("V1 has negative variances")
    }
    if (min(diag(H_inv)) < 0) {
      print("H_inv has negative variances")
    }
    # If matrices are not p.d. due to rounding issues, convert to nearest p.d. 
    # matrix using method proposed in Higham (2002)
    if (min(Re(eigen(V1)$values)) < 0) { 
      V1_pd <- Matrix::nearPD(V1)
      R1 <- chol(V1_pd$mat)
      V1_pd_diff <- sum(abs(eigen(V1)$values - eigen(V1_pd$mat)$values))
      print(paste0("V1: absolute eigenvalue difference to nearest p.d. matrix: ", 
                   V1_pd_diff))
    } else {
      R1 <- chol(V1)
    }
    if (min(Re(eigen(H_inv)$values)) < 0) {
      H_inv_pd <- Matrix::nearPD(H_inv)
      R2_inv <- chol(H_inv_pd$mat)
      H_inv_pd_diff <- sum(abs(eigen(H_inv)$values - eigen(H_inv_pd$mat)$values))
      print(paste0("H_inv: absolute eigenvalue difference to nearest p.d. matrix: ", 
                   H_inv_pd_diff))
      if (H_inv_pd_diff > 5) {
        stop("NaNs created during variance adjustment, likely due to lack of 
      smoothness in the posterior. Please run the sampler for more iterations or 
      do not run variance adjustment.")
      }
    } else {
      R2_inv <- chol(H_inv)
    }
    # Obtain the variance adjustment matrix
    R2 <- solve(R2_inv)
    R2R1 <- R2 %*% R1
    
    # Apply variance adjustment to parameters
    par_adj <- apply(par_samps, 1, DEadj, par_hat = par_hat, R2R1 = R2R1, 
                     simplify = FALSE)
    par_adj <- matrix(unlist(par_adj), byrow=TRUE, nrow=M)
    
    #=============== Output adjusted parameters ==================================
    # Store adjusted parameters
    xi_red_draws[[d]] <- par_adj
  } 
  
  # Stack together the draws
  xi_red <- do.call("abind", c(xi_red_draws, along = 1))

  # Get adjustment parameters by averaging across draws
  # Get posterior median estimates
  xi_med <- apply(xi_red, 2, stats::median, na.rm = TRUE)
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  
  xi_adjust <- list(xi_red = xi_red, xi_med = xi_med, runtime = runtime)
  
  return(xi_adjust)
}