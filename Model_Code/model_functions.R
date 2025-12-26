#=================================================
# Functions for running WOLCAN model
# Author: Stephanie Wu
# Date created: 2024/05/29
# Date updated: 2024/05/29
#=================================================

# Runs the WOLCAN model
# Run Weighted Overfitted Latent Class Analysis for Non-probability samples (WOLCAN)
# `wolcan` runs a weighted overfitted latent class analysis for non-probability 
# samples (WOLCAN), then saves and returns the results.
# Inputs:
# Inherit from wolca
#   x_mat: Matrix of multivariate categorical exposures for those in the 
# non-probability sample (NPS). (n_B)xJ
#   dat_B: Dataframe of covariates for those in the NPS. nxq
#   dat_R: Dataframe of covariates for those in the reference sample (RS). nxq
#   pred_covs_B: String vector of covariates to be used to predict selection 
# into the non-probability sample.
#   pred_covs_R: String vector of covariates to be used to predict selection 
# into the reference sample
#   pi_R: Vector of known selection probabilities for those in the reference sample. 
#   hat_pi_R: Vector of reference sample selection probabilities for those in 
# the non-probability sample. If `NULL` (default), assumed to be unknown and 
# subsequently estimated using BART. Otherwise, assumed to be known if provided.
#   num_post: Number of posterior BART draws for estimating weights. Default is 1000. 
#   frame_B: Coverage probability of the non-probability sample sampling frame. 
# Default is 1, indicating complete overlap with the population of interest.
#   frame_R: Coverage probability of the reference sample sampling frame. 
# Default is 1, indicating complete overlap with the population of interest.
#   D: Number of posterior draws for the multiple imputation variance calculation.
# Default is 10.
#   parallel: Boolean specifying if parallelization is to be used for the 
# multiple imputation variance calculation (default is `TRUE`).
#   n_cores: Number of cores to run in parallel. Default is 4.
#   wts_adj: String specifying how variability in the pseudo-weights should be 
# accounted for. Default ("MI") uses multiple imputation and re-runs the WOLCA 
# model with different pseudo-weights. Other options are: "WS all", where WOLCA 
# is run once but the Williams & Savitsky variance adjustment uses different 
# pseudo-weights; "WS mean", where WOLCA is run once and the variance adjustment 
# uses the mean pseudo-weights; and "none", where only the mean pseudo-weights 
# are used and no uncertainty is accounted for. 
#   wts_quantile: Boolean specifying whether quantiles of the pseudo-weight 
# distribution should be used to improve efficiency, or if these should be 
# random draws. Default is `TRUE` and uses quantiles.
#   adjust: Boolean specifying whether to incorporate the post-processing variance 
# adjustment to ensure proper uncertainty intervals. Default is `TRUE`.
#   adapt_seed: Default is 1.
#   fixed_seed: Default is 1.
#   run_adapt: Boolean specifying whether the adaptive sampler should be run to 
# determine the number of classes, K. Default is `TRUE`. If `FALSE`, `K_fixed`
# must be specified.
#   alpha_fixed, eta_fixed: only specify if K_fixed is not `NULL`. Otherwise, 
# defaults to alpha = K, eta = 1.
#   num_reps: Number of bootstrap replicates for the WS variance adjustment
#   pred_model: String specifying type of prediction model to use to create the 
# pseudo-weights. Must be one of `"bart"` (default) or `"glm"`.
# overwrite: Whether to overwrite existing results. Default is TRUE.
# weights_only: Whether only the weights should be obtained. Default is FALSE.
#   save_res_d: Boolean indicating if the fixed sampler results for each draw d, 
# used when MI = TRUE, should be saved. Default is `FALSE`.
# Outputs:
# 
wolcan <- function(x_mat, dat_B, dat_R, pred_model = c("bart", "glm"), 
                   pred_covs_B, pred_covs_R, pi_R, hat_pi_R = NULL, 
                   num_post = 1000, frame_B = 1, frame_R = 1, 
                   trim_method = "t2", trim_c = 20, 
                   D = 10, parallel = TRUE, n_cores = 4, 
                   wts_adj = c("MI", "WS all", "WS mean", "none"), 
                   wts_quantile = TRUE, 
                   adjust = TRUE, tol = 1e-8, num_reps = 100,
                   run_adapt = TRUE, K_max = 30, adapt_seed = 1, 
                   K_fixed = NULL, fixed_seed = 1, class_cutoff = 0.05,
                   n_runs = 20000, burn = 10000, thin = 5, update = 1000,
                   save_res = TRUE, save_res_d = FALSE, save_path = NULL,
                   alpha_adapt = NULL, eta_adapt = NULL,
                   alpha_fixed = NULL, eta_fixed = NULL, 
                   overwrite = TRUE, weights_only = FALSE) {
  
  # Check errors
  if (!run_adapt) {
    if (is.null(K_fixed)) {
      stop("K_fixed must be specified if run_adapt is set to FALSE")
    }
  }
  if (n_cores > detectCores()) {
    stop("n_cores must not be larger than the number of cores available, obtained by detectCores()")
  }
  
  # Begin runtime tracker
  start_time <- Sys.time()
  
  #================= Read in data ==============================================
  print("Read in data")
  
  # Obtain dimensions
  n <- dim(x_mat)[1]        # Number of individuals
  J <- dim(x_mat)[2]        # Number of exposure items
  R_j <- apply(x_mat, 2,    # Number of exposure categories for each item
               function(x) length(unique(x)))  
  R <- max(R_j)             # Maximum number of exposure categories across items
  
  if (!is.null(adapt_seed)) {
    set.seed(adapt_seed)
  }
  
  #================ Get pseudo-weights =========================================
  
  if (file.exists(paste0(save_path, "_wolcan_weights.RData")) & !overwrite) {
    # Load weights if they already exist
    load(paste0(save_path, "_wolcan_weights.RData"))
    
  } else {
    print("Getting pseudo-weights...")
    est_weights <- get_weights_bart(dat_B = dat_B, dat_R = dat_R, 
                                    pred_model = pred_model,
                                    pred_covs_B = pred_covs_B, 
                                    pred_covs_R = pred_covs_R, 
                                    num_post = num_post, pi_R = pi_R, 
                                    hat_pi_R = hat_pi_R,
                                    frame_B = frame_B, frame_R = frame_R,
                                    trim_method = trim_method, trim_c = trim_c)
  } 
  
  # Mean estimated weights
  wts <- est_weights$wts
  # Posterior distribution of weights
  wts_post <- est_weights$wts_post
  
  #================ Run WOLCA model with mean weights ==========================
  if (weights_only) {

    # Stop runtime tracker
    runtime <- Sys.time() - start_time
    est_weights$runtime <- runtime
    
    # Save weights
    if (save_res & overwrite) {
      save(est_weights, file = paste0(save_path, "_wolcan_weights.RData"))
    }
    
    return(est_weights)
    
  } else {
    if (run_adapt) {
      # Load adaptive results if they already exist
      if (file.exists(paste0(save_path, "_wolcan_adapt.RData")) & !overwrite) {
        load(paste0(save_path, "_wolcan_adapt.RData"))
      } else {
        ### Using mean estimated weights, run WOLCA model adaptive sampler to obtain 
        ### number of latent classes K
        print("Obtaining number of latent classes...")
        res <- wolca(x_mat = x_mat, sampling_wt = wts, run_sampler = "adapt", 
                     K_max = K_max, adapt_seed = adapt_seed, 
                     class_cutoff = class_cutoff, 
                     n_runs = n_runs, burn = burn, thin = thin, 
                     update = update, save_res = FALSE, 
                     alpha_adapt = alpha_adapt, eta_adapt = eta_adapt)
        
        # Post-processing to recalibrate labels and remove extraneous empty classes
        res$post_MCMC_out <- post_process_wolca(MCMC_out = res$MCMC_out, J = J, 
                                                R = R, class_cutoff = class_cutoff)
        
        # Obtain posterior estimates, reduce number of classes, analyze results
        res$estimates <- get_estimates_wolca(MCMC_out = res$MCMC_out, 
                                             post_MCMC_out = res$post_MCMC_out, 
                                             n = n, J = J, x_mat = x_mat)
        # Save adaptive sampler results
        if (save_res) {
          save(res, file = paste0(save_path, "_wolcan_adapt.RData"))
        }
      }
      # Save space
      res$MCMC_out <- NULL
      
      # Define K, depending on if adaptive sampler results exist
      K <- res$estimates$K_red
      
    } else {
      # Define K
      K <- K_fixed
      # Initialize results list
      res <- list()
    }
    print(paste0("K: ", K))
    # if (run_sampler == "adapt") {
    #   K <- res$K_fixed
    # } else {
    #   K <- res$estimates$K_red
    # }
    
    #================ Calculate variance for WOLCA using MI ======================
    # Hyperparameter for prior for pi
    if (is.null(alpha_fixed)) {
      alpha_fixed <- rep(K, K)
    }
    # Hyperparameter for prior for theta
    # Unviable categories have value 0.01 to prevent rank deficiency issues
    if (is.null(eta_fixed)) {
      eta_fixed <- matrix(0.01, nrow = J, ncol = R) 
      for (j in 1:J) {
        eta_fixed[j, 1:R_j[j]] <- rep(1, R_j[j]) 
      }
    }
    
    if (wts_adj == "MI") {
      # Get estimates by running the WOLCA fixed sampler for multiple draws from the 
      # weights posterior, using parallelization if desired. 
      # Apply post-processing variance adjustment for proper variance estimation.
      print("Running estimation and variance...")
      comb_estimates <- get_var_dir(D = D, wts_post = wts_post, 
                                    wts_quantile = wts_quantile, K = K, 
                                    tol = tol, parallel = parallel, 
                                    n_cores = n_cores, adjust = adjust, 
                                    x_mat = x_mat, fixed_seed = fixed_seed,
                                    class_cutoff = class_cutoff, n_runs = n_runs, 
                                    burn = burn, thin = thin, update = update, 
                                    alpha_fixed = alpha_fixed, eta_fixed = eta_fixed,
                                    save_res_d = save_res_d, save_path = save_path)
      
      # Add variance estimates to output list
      res$estimates_adjust <- comb_estimates
      
    } else {    # Run fixed sampler
      
      print("Running fixed sampler...")
      
      # Set seed
      if (!is.null(fixed_seed)) {
        set.seed(fixed_seed)
      }
      
      # Initialize OLCA model using fixed number of classes. Obtain pi, theta, c_all
      OLCA_params <- init_OLCA(alpha = alpha_fixed, eta = eta_fixed, n = n,
                               K = K, J = J, R = R)
      # Run MCMC algorithm using fixed number of classes
      # Obtain pi_MCMC, theta_MCMC, c_all_MCMC
      w_all <- c(wts / (sum(wts) / n)) # Mean weights normalized to sum to n
      res$MCMC_out <- run_MCMC_Rcpp_wolca(OLCA_params = OLCA_params, n_runs = n_runs, 
                                          burn = burn, thin = thin, K = K, J = J, 
                                          R = R, n = n, w_all = w_all, x_mat = x_mat, 
                                          alpha = alpha_fixed, eta = eta_fixed,
                                          update = update)
      # Post-processing to recalibrate labels and remove extraneous empty classes
      # Obtain K_med, pi, theta
      res$post_MCMC_out <- post_process_wolca(MCMC_out = res$MCMC_out, J = J, R = R,
                                              class_cutoff = class_cutoff)
      # Obtain posterior estimates, reduce number of classes, analyze results
      # Obtain K_red, pi_red, theta_red, pi_med, theta_med, c_all, pred_class_probs
      res$estimates <- get_estimates_wolca(MCMC_out = res$MCMC_out, 
                                           post_MCMC_out = res$post_MCMC_out, n = n, J = J,
                                           x_mat = x_mat)
      K_red <- res$estimates$K_red
      print(paste0("K_red: ", K_red))
      
      if (wts_adj == "WS all") {
        print("Accounting for weights uncertainty...")
        res$estimates_adjust <- suppressMessages(
          get_var_adj(D = D, K = K_red, res = res, wts_post = wts_post, tol = tol,
                      num_reps = num_reps, save_res = save_res, 
                      save_path = save_path, adjust_seed = fixed_seed))
      } else if (wts_adj == "WS mean") {
        print("Using mean weights...")
        res$estimates_adjust <- suppressMessages(
          get_var_adj_mean(K = K_red, res = res, wts = wts, tol = tol, 
                           num_reps = num_reps, save_res = save_res, 
                           save_path = save_path, adjust_seed = fixed_seed))
      } else {
        print("Not accounting for weights uncertainty...")
      }
    }
    #================= Save and return output ==================================
    # Stop runtime tracker
    runtime <- Sys.time() - start_time
    res$runtime <- runtime
    
    # # Store estimated pseudo-weights
    # res$est_weights <- est_weights
    
    # Store additional data variables used for estimating pseudo-weights
    res$data_vars <- c(res$data_vars, list(dat_B = dat_B, dat_R = dat_R, 
                                           pred_covs_B = pred_covs_B, 
                                           pred_covs_R = pred_covs_R, pi_R = pi_R))
    
    class(res) <- "wolca"
    
    # Save output
    if (save_res) {
      save(res, file = paste0(save_path, "_wolcan_results.RData"))
    }
    
    return(res)
  }
} 

# Only use the mean posterior weights
# wts: mean weights
get_var_adj_mean <- function(K, res, wts, tol, num_reps = 100, 
                             alpha = NULL, eta = NULL,  
                        save_res = TRUE, save_path = NULL, 
                        adjust_seed = NULL) {
  
  # Begin runtime tracker
  start_time <- Sys.time()
  
  # Set seed
  if (!is.null(adjust_seed)) {
    set.seed(adjust_seed)
  }
  
  #================= Extract dimensions and catch errors =======================
  # Check num_reps
  if ((num_reps %% 100 != 0) | num_reps < 1) {
    stop("num_reps must be a whole number greater than 0, recommended to be at least 50.
    More replicates will lead to more accurate results but will take longer to run.")
  }
  
  # Extract data elements into the global environment
  J <- res$data_vars$J
  R_j <- res$data_vars$R_j
  R <- res$data_vars$R
  n <- res$data_vars$n
  x_mat <- res$data_vars$x_mat
  w_all <- wts / (sum(wts) / n)

  #================= Initialize hyperparameters ================================
  # Default hyperparameters for pi and theta
  if (is.null(alpha)) {
    alpha <- rep(1, K) / K   # Hyperparameter for prior for pi
  }
  if (is.null(eta)) {
    # Hyperparameter for prior for theta
    # Unviable categories have value 0.01 to prevent rank deficiency issues
    eta <- matrix(0.01, nrow = J, ncol = R) 
    for (j in 1:J) {
      eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
    }
  }
  
  #=============== Run Stan model ==============================================
  print("Running variance adjustment")
  
  # Define data for Stan model
  data_stan <- list(K = K, J = J, R = R, n = n, X = x_mat, weights = w_all, 
                    alpha = alpha, eta = eta)
  
  # Stan parameters of interest
  par_stan <- c('pi', 'theta')  # subset of parameters interested in
  
  # Stan model
  mod_stan <- stanmodels$WOLCA_main
  
  # Run Stan model
  # Stan will pass warnings from calling 0 chains, but will still create an 
  # out_stan object for the 'grad_log_prob()' method
  out_stan <- rstan::sampling(object = mod_stan, data = data_stan, 
                              pars = par_stan, chains = 0, iter = 0, refresh = 0)
  
  #=============== Convert to unconstrained parameters =========================
  # Convert params from constrained space to unconstrained space
  unc_par_hat <- rstan::unconstrain_pars(out_stan, 
                                         list("pi" = res$estimates$pi_med,
                                              "theta" = res$estimates$theta_med))
  # Get posterior MCMC samples in unconstrained space for all parameters
  M <- dim(res$estimates$pi_red)[1]
  unc_par_samps <- lapply(1:M, unconstrain_wolca, stan_model = out_stan, K = K, 
                          pi = res$estimates$pi_red, 
                          theta = res$estimates$theta_red)
  unc_par_samps <- matrix(unlist(unc_par_samps), byrow = TRUE, nrow = M)
  
  #=============== Post-processing adjustment in unconstrained space ===========
  # Estimate Hessian
  H_hat <- -1*stats::optimHess(unc_par_hat, 
                               gr = function(x){rstan::grad_log_prob(out_stan, x)})


  ### Create survey replicates 
  # Survey data frame for specifying survey design
  # Clusters are simply individuals
  svy_data <- data.frame(cluster_id = 1:nrow(x_mat), x_mat = x_mat, 
                         w_all = w_all)
  # Specify survey design
  svydes <- survey::svydesign(ids = ~cluster_id, weights = ~w_all,
                              data = svy_data)
  # Create svrepdesign
  svyrep <- survey::as.svrepdesign(design = svydes, type = "mrbbootstrap",
                                   replicates = num_reps)
  # Get survey replicates and gradient for each replicate
  rep_temp <- survey::withReplicates(design = svyrep, theta = grad_par,
                                     stan_mod = mod_stan, stan_data = data_stan,
                                     par_stan = par_stan, u_pars = unc_par_hat)
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
      stop(paste0("NaNs created during variance adjustment, likely due to lack of ",
      "smoothness in the posterior. Please run the sampler for more iterations or ", 
      "do not run variance adjustment."))
    }
  } else {
    R2_inv <- chol(H_inv)
  }
  # Obtain the variance adjustment matrix
  R2 <- solve(R2_inv)
  R2R1 <- R2 %*% R1
  
  # Apply variance adjustment to parameters
  par_adj <- apply(unc_par_samps, 1, DEadj, par_hat = unc_par_hat, R2R1 = R2R1, 
                   simplify = FALSE)
  par_adj <- matrix(unlist(par_adj), byrow=TRUE, nrow=M)
  
  #=============== Convert adjusted to constrained space =======================
  # Constrained adjusted parameters for all MCMC samples
  pi_red_adj <- matrix(NA, nrow=M, ncol=K)
  theta_red_adj <- array(NA, dim=c(M, J, K, R))
  for (i in 1:M) {
    ##### FIX WITH CUSTOMIZED ERROR
    constr_pars <- rstan::constrain_pars(out_stan, par_adj[i,])
    pi_red_adj[i, ] <- constr_pars$pi
    theta_red_adj[i,,,] <- constr_pars$theta
  }
  
  #=============== Output adjusted parameters ==================================
  # Re-normalize pi and theta for each iteration
  pi_red_adj = pi_red_adj / rowSums(pi_red_adj)  
  theta_red_adj <- plyr::aaply(theta_red_adj, c(1, 2, 3), function(x) x / sum(x),
                               .drop = FALSE) 
  
  # Get posterior median estimates
  pi_med_adj <- apply(pi_red_adj, 2, stats::median)
  theta_med_adj <- apply(theta_red_adj, c(2,3,4), stats::median)
  
  # Renormalize posterior median estimates for pi and theta to sum to 1
  pi_med_adj <- pi_med_adj / sum(pi_med_adj)  
  theta_med_adj <- plyr::aaply(theta_med_adj, c(1, 2), function(x) x / sum(x),
                               .drop = FALSE)  # Re-normalize
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  # Add variance adjustment runtime to overall runtime
  sum_runtime <- runtime + res$runtime
  res$runtime <- sum_runtime
  
  estimates_adjust <- list(pi_red = pi_red_adj, theta_red = theta_red_adj, 
                           pi_med = pi_med_adj, theta_med = theta_med_adj, 
                           c_all = res$estimates$c_all,
                           pred_class_probs = res$estimates$pred_class_probs)
  
  return(estimates_adjust)
}


# wts_post: (n1)xM posterior distribution of weights for individuals in the NPS, 
# with each column corresponding to a posterior draw
get_var_adj <- function(D, K, res, wts_post, wts_quantile, tol, num_reps = 100, 
                             alpha = NULL, eta = NULL,  
                             save_res = TRUE, save_path = NULL, 
                             adjust_seed = NULL) {
  
  # Begin runtime tracker
  start_time <- Sys.time()
  
  # Set seed
  if (!is.null(adjust_seed)) {
    set.seed(adjust_seed)
  }
  
  #================= Extract dimensions and catch errors =======================
  # Check num_reps
  if ((num_reps %% 100 != 0) | num_reps < 1) {
    stop("num_reps must be a whole number greater than 0, recommended to be at least 50.
    More replicates will lead to more accurate results but will take longer to run.")
  }
  
  # Extract data elements into the global environment
  J <- res$data_vars$J
  R_j <- res$data_vars$R_j
  R <- res$data_vars$R
  n <- res$data_vars$n
  x_mat <- res$data_vars$x_mat
  
  #================= Initialize hyperparameters ================================
  # Default hyperparameters for pi and theta
  if (is.null(alpha)) {
    alpha <- rep(1, K) / K   # Hyperparameter for prior for pi
  }
  if (is.null(eta)) {
    # Hyperparameter for prior for theta
    # Unviable categories have value 0.01 to prevent rank deficiency issues
    eta <- matrix(0.01, nrow = J, ncol = R) 
    for (j in 1:J) {
      eta[j, 1:R_j[j]] <- rep(1, R_j[j]) 
    }
  }
  
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
    
    # Get pseudo-weight draw 
    if (wts_quantile) {  # Use quantiles for efficiency
      
      # Draw from weights posterior closest to quantile
      draw <- which.min(abs(col_meds - med_quants[d]))
      w_all <- c(w_all_post[, draw])
      wts_draws[, d] <- wts_post[, draw]
      
    } else {  # Use random sample
      # Draw from weights posterior
      set.seed(d)
      draw <- sample(1:ncol(wts_post), size = 1)
      wts_d <- c(wts_post[, draw])
    }
    
    #=============== Run Stan model ==============================================
    
    # Define data for Stan model
    data_stan <- list(K = K, J = J, R = R, n = n, X = x_mat, weights = w_all, 
                      alpha = alpha, eta = eta)
    
    # Stan parameters of interest
    par_stan <- c('pi', 'theta')  # subset of parameters interested in
    
    # Stan model
    mod_stan <- stanmodels$WOLCA_main
    
    # Run Stan model
    # Stan will pass warnings from calling 0 chains, but will still create an 
    # out_stan object for the 'grad_log_prob()' method
    out_stan <- rstan::sampling(object = mod_stan, data = data_stan, 
                                pars = par_stan, chains = 0, iter = 0, refresh = 0)
    
    #=============== Convert to unconstrained parameters =========================
    # Convert params from constrained space to unconstrained space
    unc_par_hat <- rstan::unconstrain_pars(out_stan, 
                                           list("pi" = res$estimates$pi_med,
                                                "theta" = res$estimates$theta_med))
    # Get posterior MCMC samples in unconstrained space for all parameters
    M <- dim(res$estimates$pi_red)[1]
    unc_par_samps <- lapply(1:M, unconstrain_wolca, stan_model = out_stan, K = K, 
                            pi = res$estimates$pi_red, 
                            theta = res$estimates$theta_red)
    unc_par_samps <- matrix(unlist(unc_par_samps), byrow = TRUE, nrow = M)
    
    #=============== Post-processing adjustment in unconstrained space ===========
    # Estimate Hessian
    H_hat <- -1*stats::optimHess(unc_par_hat, 
                                 gr = function(x){rstan::grad_log_prob(out_stan, x)})
    
    ### Create survey replicates 
    # Survey data frame for specifying survey design
    # Clusters are simply individuals
    svy_data <- data.frame(cluster_id = 1:nrow(x_mat), x_mat = x_mat, 
                           w_all = w_all)
    # Specify survey design
    svydes <- survey::svydesign(ids = ~cluster_id, weights = ~w_all,
                                data = svy_data)
    # Create svrepdesign
    svyrep <- survey::as.svrepdesign(design = svydes, type = "mrbbootstrap",
                                     replicates = num_reps)
    # Get survey replicates and gradient for each replicate
    rep_temp <- survey::withReplicates(design = svyrep, theta = grad_par,
                                       stan_mod = mod_stan, stan_data = data_stan,
                                       par_stan = par_stan, u_pars = unc_par_hat)
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
    par_adj <- apply(unc_par_samps, 1, DEadj, par_hat = unc_par_hat, R2R1 = R2R1, 
                     simplify = FALSE)
    par_adj <- matrix(unlist(par_adj), byrow=TRUE, nrow=M)
    
    #=============== Convert adjusted to constrained space =======================
    # Constrained adjusted parameters for all MCMC samples
    pi_red_adj <- matrix(NA, nrow=M, ncol=K)
    theta_red_adj <- array(NA, dim=c(M, J, K, R))
    for (i in 1:M) {
      ##### FIX WITH CUSTOMIZED ERROR
      constr_pars <- rstan::constrain_pars(out_stan, par_adj[i,])
      pi_red_adj[i, ] <- constr_pars$pi
      theta_red_adj[i,,,] <- constr_pars$theta
    }
    
    #=============== Output adjusted parameters ==================================
    # Re-normalize pi and theta for each iteration
    pi_red_adj = pi_red_adj / rowSums(pi_red_adj)  
    theta_red_adj <- plyr::aaply(theta_red_adj, c(1, 2, 3), function(x) x / sum(x),
                                 .drop = FALSE) 
    
    # Store adjusted parameters
    pi_red_draws[[d]] <- pi_red_adj
    theta_red_draws[[d]] <- theta_red_adj
  } 
  
  # Stack together the draws
  pi_red <- do.call("abind", c(pi_red_draws, along = 1))
  theta_red <- do.call("abind", c(theta_red_draws, along = 1))
  
  # Get adjustment parameters by averaging across draws
  # Get posterior median estimates
  pi_med <- apply(pi_red, 2, stats::median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), stats::median, na.rm = TRUE)
  theta_med <- ifelse(theta_med < tol, tol, theta_med) # prevent underflow
  theta_med <- plyr::aaply(theta_med, c(1, 2), function(x) x / sum(x),
                           .drop = FALSE)  # Re-normalize
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  # Add variance adjustment runtime to overall runtime
  sum_runtime <- runtime + res$runtime
  res$runtime <- sum_runtime
  
  estimates_adjust <- list(pi_red = pi_red, theta_red = theta_red, 
                           pi_med = pi_med, theta_med = theta_med, 
                           c_all = res$estimates$c_all,
                           pred_class_probs = res$estimates$pred_class_probs)
  
  return(estimates_adjust)
}



# Obtain posterior parameter estimates and class assignments after propagating 
# uncertainty from the posterior weight distribution
# Inputs:
#   D: number of draws from the posterior weight distribution. WOLCA fixed sampler 
# will be run for each of these draws. Default is 10.
#   wts_post: Posterior distribution of estimated pseudo-weights
#   K: Number of classes, determined by the adaptive sampler
#   tol: Underflow tolerance
#   adjust: Boolean indicating if post-processing variance adjustment should be 
# applied. Default is TRUE.
#   alpha_fixed, eta_fixed. Only non-NULL if K came from K_fixed
get_var_dir <- function(D = 10, wts_post, K, wts_quantile, x_mat, fixed_seed = 1, 
                        tol = 1e-8, parallel = TRUE, n_cores = 4,
                        adjust = TRUE, class_cutoff = 0.05, n_runs = 20000, 
                        burn = 10000, thin = 5, update = 1000,
                        alpha_fixed = NULL, eta_fixed = NULL, 
                        save_res_d = FALSE, save_path = NULL) {
  
  # Get dimensions
  M <- floor(n_runs / thin) - floor(burn / thin) 
  n <- nrow(x_mat)
  J <- ncol(x_mat)
  
  if (parallel) {  # Parallel version
    # # Number of cores available
    # n_cores <- detectCores()  
    
    # Create a cluster using sockets
    cluster <- parallel::makeCluster(n_cores)
        # # Export global functions to be available to all workers
        # clusterExport(cluster, varlist = c("wolca", "wolca_var_adjust", "catch_errors", 
        #                                    "init_OLCA", "run_MCMC_Rcpp_wolca", 
        #                                    "post_process_wolca", "get_estimates_wolca"))
    
    # Run wolca_d in parallel
    res_all <- parLapply(cluster, 1:D, wolca_d, wts_post = wts_post, 
                         wts_quantile = wts_quantile, 
                         x_mat = x_mat, K = K, adjust = adjust, 
                         class_cutoff = class_cutoff, n_runs = n_runs, 
                         burn = burn, thin = thin, update = update,
                         fixed_seed = fixed_seed, alpha_fixed = alpha_fixed, 
                         eta_fixed = eta_fixed, D = D,
                         save_res_d = save_res_d, save_path = save_path)
    # Shutdown cluster
    stopCluster(cluster)
    
  } else {  # Serial version
    # Run wolca_d serially
    res_all <- lapply(1:D, wolca_d, wts_post = wts_post, 
                      wts_quantile = wts_quantile, x_mat = x_mat, K = K, 
                      adjust = adjust, class_cutoff = class_cutoff, 
                      n_runs = n_runs, burn = burn, thin = thin, update = update,
                      fixed_seed = fixed_seed, alpha_fixed = alpha_fixed, 
                      eta_fixed = eta_fixed, D = D, 
                      save_res_d = save_res_d, save_path = save_path)
  }
  
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
    # Convert theta to numeric array
    theta_red_numer <- array(as.numeric(est_d$theta_red), dim = dim(est_d$theta_red))
    theta_red_draws[[d]] <- theta_red_numer
    wts_draws[, d] <- est_d$wts
    messages_draws[[d]] <- ifelse(is.null(est_d$message), "no error", est_d$message)
    
    # # Add dimensions if necessary
    # if (K_d < K) {
    #   # Add NA's to full MCMC outputs for the missing classes
    #   miss_dims <- K - K_d
    #   pi_red_draws[[d]] <- abind(pi_red_draws[[d]], 
    #                              array(NA, dim=c(M, miss_dims)), along=2)
    #   theta_red_draws[[d]] <- abind(theta_red_draws[[d]], 
    #                                 array(NA, dim=c(M, J, miss_dims, 
    #                                           dim(theta_red_draws[[d]])[4])), 
    #                                 along = 3)
    # }
  }
  
  
  # Reorder the classes across draws using minimum theta distance
  # Initialize lists across draws
  pi_red_new <- vector("list", length = D_red)
  theta_red_new <- vector("list", length = D_red)
  # Set base order using draw 1
  theta_1 <- apply(theta_red_draws[[1]], c(2, 3, 4), stats::median)
  pi_red_new[[1]] <- pi_red_draws[[1]]
  theta_red_new[[1]] <- theta_red_draws[[1]]
  # Harmonize each draw to match draw 1 order
  if (D_red > 1) {
    for (d in 2:D_red) {
      theta_d <- apply(theta_red_draws[[d]], c(2, 3, 4), stats::median)
      # Minimum theta distance using mean abs distance
      theta_perm <- get_theta_dist_wolcan(est_theta = theta_d,
                                          true_theta = theta_1,
                                          est_K = K, true_K = K,
                                          subset = FALSE,
                                          dist_type = "mean_abs")
      theta_dist <- theta_perm$theta_dist
      # New order
      order <- theta_perm$order
      # Relabel pi and theta classes for draw d using new order
      pi_red_new[[d]] <- pi_red_draws[[d]][, order]
      theta_red_new[[d]] <- theta_red_draws[[d]][, , order, ]
    }
  }
  
  
      # Clustering to harmonize classes
      # Cluster individuals into reduced number of classes using agglomerative clustering
      # Calculate pairwise distance matrix using Hamming distance: proportion of
      # iterations where two individuals have differing class assignments
      distMat <- e1071::hamming.distance(t(c_all_draws))
      # Hierarchical clustering dendrogram
      dendrogram <- stats::hclust(stats::as.dist(distMat), method = "complete")
      # # Group individuals into K classes
      # red_c_all <- stats::cutree(dendrogram, k = K)
      # # Modify classes if any classes are less than the cutoff percentage
      # class_prop <- prop.table(table(red_c_all))
      # if (any(class_prop < class_cutoff)) {
      #   # Get classes that are too small
      #   small <- which(class_prop < class_cutoff)
      #   # Group individuals into a larger number of classes
      #   red_c_all_temp <- stats::cutree(dendrogram, k = K + length(small))
      #   red_c_all <- red_c_all_temp
      #   class_prop_temp <- prop.table(table(red_c_all_temp))
      #   # Get updated classes that are too small
      #   small_temp <- sort(which(class_prop_temp < class_cutoff))
      #   for (small_c in 1:length(small_temp)) {
      #     c_ind <- small_temp[small_c]
      #     class_small <- which(red_c_all_temp == c_ind)
      #     # Get nearest class
      #     inds <- 1:length(class_prop_temp)
      #     class_dist <- sapply(inds, function(x)
      #       mean(distMat[class_small, which(red_c_all_temp == x)]))
      #     # Set small class distance to Inf
      #     class_dist[small_temp] <- Inf
      #     nearest <- which.min(class_dist[-c_ind])
      #     red_c_all[red_c_all_temp == c_ind] <- nearest
      #   }
      #   class_prop <- prop.table(table(red_c_all))
      # }
      # 
      # # Get unique reduced classes to aid relabeling
      # unique_red_classes <- unique(red_c_all)
      # 
      # # Relabel so that the class numbering doesn't exceed the number of classes
      # red_c_all_old <- red_c_all
      # for (k in 1:length(unique_red_classes)) {
      #   red_c_all[red_c_all == unique_red_classes[k]] <- k
      # }
      # 
      # # For each draw, relabel new classes using the most common old class assignment
      # pi_red_new <- vector("list", length = D_red)
      # theta_red_new <- vector("list", length = D_red)
      # for (d in 1:D_red) {
      #   # Initialize pi and theta
      #   pi <- matrix(NA, nrow = M, ncol = K)
      #   theta <- array(NA, dim = c(M, J, K, dim(theta_red_draws[[d]])[4]))
      # 
      #   # Get new order
      #   new_order <- numeric(K)
      #   for (k in 1:K) {
      #     # Individuals who were assigned to old class k
      #     old_k_inds <- which(c_all_draws[d, ] == k)
      #     # Most common new class assignment among those individuals
      #     mode_new <- get_mode(red_c_all[old_k_inds])
      #     new_order[k] <- mode_new
      #     # Reorder classes for all pi and theta estimates
      #     if (!is.na(mode_new)) {
      #       pi[, new_order[k]] <- pi_red_draws[[d]][, k]  # reorder
      #       theta[, , new_order[k], ] <- theta_red_draws[[d]][, , k, ]
      #     }
      #   }
      # 
      #   # pi <- pi_red_draws[[d]][, new_order, drop = FALSE]  # reorder
      #   pi <- pi / rowSums(pi, na.rm = TRUE)
      #   pi_red_new[[d]] <- pi
      #   # theta <- theta_red_draws[[d]][, , new_order, , drop = FALSE]  # reorder
      #   theta <- ifelse(theta < tol, tol, theta)  # prevent underflow
      #   for (m in 1:dim(theta)[1]) {  # normalize
      #     for (j in 1:dim(theta)[2]) {
      #       for (k in 1:dim(theta)[3]) {
      #         theta[m, j, k, ] <- theta[m, j, k, ] / sum(theta[m, j, k, ])
      #       }
      #     }
      #   }
      #   theta_red_new[[d]] <- theta
      # }
  
  # Stack together the draws
  pi_red_stack <- do.call("abind", c(pi_red_new, along = 1))
  theta_red_stack <- do.call("abind", c(theta_red_new, along = 1))
  # theta_red_stack <- apply(theta_red_new, c(2, 3, 4), c)
  
  #============= Check for duplicate classes using modal exposure categories ===
  # Posterior median estimate for theta across iterations
  theta_med_temp <- apply(theta_red_stack, c(2, 3, 4), stats::median)
  # Posterior modal exposure categories for each exposure item and reduced class
  theta_modes <- apply(theta_med_temp, c(1, 2), which.max)
  # Find NA columns
  theta_modes <- apply(theta_modes, 2, as.numeric)
  na_cols <- is.na(colSums(theta_modes))
  # Identify unique classes, removing duplicates and NA columns
  unique_classes <- which(!duplicated(theta_modes, MARGIN = 2))
  unique_classes <- unique_classes[!na_cols]
  # Number of unique classes
  K_red <- length(unique_classes) 
  
  # Combine duplicated classes and re-normalize pi to sum to 1
  pi_red <- pi_red_stack[, unique_classes, drop = FALSE] # Initialize pi for unique classes
      # if (K_red < dim(pi_red_stack)[2]) {  # Check if there are duplicated classes
      #   for (k in 1:K_red) {
      #     # Find duplicated modal theta patterns
      #     dupes_k <- apply(theta_modes, 2, function(x)  
      #       identical(x,theta_modes[, unique_classes[k]]))
      #     # Combine class proportions for all duplicated patterns together
      #     pi_red[, k] <- apply(as.matrix(pi_red_stack[, dupes_k]), 1, sum)  
      #   }
      # }
  # Re-normalize to ensure pi sums to 1 for each iteration
  pi_red = pi_red / rowSums(pi_red)  
  
  # Get posterior parameter samples for unique classes for theta 
  theta_red <- theta_red_stack[, , unique_classes, ]
  theta_red <- ifelse(theta_red < 1e-8, 1e-8, theta_red) # prevent underflow
  theta_red <- plyr::aaply(theta_red, c(1, 2, 3), function(x) x / sum(x),
                           .drop = FALSE) # Re-normalize
  
  # Posterior median estimates
  pi_med <- apply(pi_red, 2, stats::median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red, c(2, 3, 4), stats::median, na.rm = TRUE)
  theta_med <- ifelse(theta_med < tol, tol, theta_med) # prevent underflow
  theta_med <- plyr::aaply(theta_med, c(1, 2), function(x) x / sum(x),
                           .drop = FALSE)  # Re-normalize
  
  #============== Update c using unique classes and posterior estimates ========
  c_all <- numeric(n)  # Placeholder class assignments
  pred_class_probs <- matrix(NA, nrow = n, ncol = K_red) # Posterior class membership probabilities
  log_cond_c <- matrix(NA, nrow = n, ncol = K_red)       # Individual log-likelihood for each class
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    for (k in 1:K_red) {
      # Calculate theta component of individual log-likelihood assuming class k
      log_theta_comp_k <- 0
      for (j in 1:J) {
        log_theta_comp_k <- log_theta_comp_k + log(theta_med[j, k, x_mat[i, j]])
      }
      # Individual log-likelihood for class k
      log_cond_c[i, k] <- log(pi_med[k]) + log_theta_comp_k
    }
    # Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs[i, ] <- exp(log_cond_c[i, ] - matrixStats::logSumExp(log_cond_c[i, ]))
    # Update class assignment using the posterior probabilities
    c_all[i] <- LaplacesDemon::rcat(n = 1, p = pred_class_probs[i, ])
  }
  
  comb_estimates <- list(K_red_draws_all = K_red_draws_all, 
                         K_red_draws = K_red_draws, c_all_draws = c_all_draws, 
                         pi_red = pi_red, theta_red = theta_red,
                         pi_med = pi_med, theta_med = theta_med, 
                         c_all = c_all, pred_class_probs = pred_class_probs,
                         dendrogram = dendrogram, wts_draws = wts_draws, 
                         K_red = K_red, K = K, messages_draws = messages_draws)
  return(comb_estimates)
}


# Function to run WOLCA using weights from one posterior draw
# Inputs:
#   d: Numeric index for posterior draw
#   wts_post: Posterior distribution of  estimated pseudo-weights
#   adjust: Boolean indicating if post-processing variance adjustment should be applied
# inherit parameters from get_var_dir
wolca_d <- function(d, wts_post, wts_quantile, x_mat, K, adjust, class_cutoff, 
                    n_runs, burn, thin, update, fixed_seed, alpha_fixed, 
                    eta_fixed, D, save_res_d, save_path) {
  
  # # Source Rcpp functions from baysc package
  # Rcpp::sourceCpp(rcpp_path)
  
  print(paste0("Draw ", d))
  
  # Get pseudo-weight draw
  if (wts_quantile) {  # Use quantiles for efficiency
    print("Using quantiles for efficiency...")
    # Get quantiles of medians of weights posterior
    col_meds <- c(apply(wts_post, 2, median))
    cutoffs <- (seq(1, D, length.out = D) - 0.5) / D
    med_quants <- stats::quantile(x = col_meds, probs = cutoffs)
    
    # Draw from weights posterior closest to quantile
    draw <- which.min(abs(col_meds - med_quants[d]))
    wts_d <- c(wts_post[, draw])
  
  } else { # Random sample
    print("Randomly drawing weights...")
    # Draw from weights posterior
    set.seed(d)
    draw <- sample(1:ncol(wts_post), size = 1)
    wts_d <- c(wts_post[, draw])
  }
  
  # Initialize fixed sampler hyperparameters
  if (is.null(alpha_fixed)) {
    alpha_fixed = rep(K, K)  # not sparsity-inducing
  }
  # Run fixed sampler
  res_d <- baysc::wolca(x_mat = x_mat, sampling_wt = wts_d, run_sampler = "fixed", 
                        K_fixed = K, fixed_seed = fixed_seed, 
                        alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
                        class_cutoff = class_cutoff, n_runs = n_runs, burn = burn, 
                        thin = thin, update = update, save_res = FALSE)
  # Apply post-processing variance adjustment for pseudo-likelihood if desired
  if (adjust) {
    try_d <- tryCatch(baysc::wolca_var_adjust(res = res_d, num_reps = 100, 
                                              save_res = FALSE, 
                                              adjust_seed = fixed_seed), 
                      error = function(e) e)
    if (is.null(try_d$message)) {  # variance adjustment worked
      res_d <- try_d
      # Get estimates
      est_d <- res_d$estimates_adjust
    } else {  # variance adjustment failed
      # Get estimates
      est_d <- res_d$estimates
      est_d$message <- try_d$message
    }
     
    # Get number of classes
    est_d$K_red <- res_d$estimates$K_red
  } else {
    # Get estimates
    est_d <- res_d$estimates
  }
  # Return the sampled weights
  est_d$wts <- wts_d
  
  # Save output
  if (save_res_d) {
    save(est_d, file = paste0(save_path, "draw_", d, "_results.RData"))
  }
  
  return(est_d)
}


# Obtain estimated weights using BART
# Inherit inputs from wolcan
# Outputs: List `est_weights` containing the following objects:
#   wts: Vector of mean posterior weights
#   wts_post: Matrix of weights from the BART posterior draws. Each column is 
# the set of weights from a single draw. Nx(num_post).
#   hat_pi_B: Vector of mean estimated NPS selection probabilities. Nx1.
#   hat_pi_B_dist: Matrix of estimated NPS selection probabilities from the 
# BART posterior draws. Each column corresonds to a single draw. Nx(num_post).
#   hat_pi_R: Vector of mean estimated reference sample selection probabilities. Nx1.
#   hat_pi_R_dist: Matrix of estimated reference sample selection probabilities 
# from the BART posterior draws. Each column corresonds to a single draw. Nx(num_post).
#   pred_model: String specifying type of prediction model to use to create the 
# pseudo-weights. Must be one of `"bart"` (default) or `"glm"`.
#   ci_level: Level for confidence interval or posterior interval. Default is 0.95.
get_weights_bart <- function(dat_B, dat_R, pred_covs_B, pred_covs_R, 
                             num_post = 1000, pi_R, hat_pi_R = NULL, 
                             frame_B = 1, frame_R = 1, trim_method = "t2", 
                             trim_c = 20, pred_model = c("bart", "glm"),
                             ci_level = 0.95) {
  # Initialize output
  est_weights <- list()
  
  # Create stacked sample
  n1 <- nrow(dat_B)
  n0 <- nrow(dat_R)
  samp_comb <- rbind(dat_B, dat_R)
  # Get known sampling probabilities for i in S_R and convert to logit
  logit_pi_R <- log(pi_R / (1 - pi_R))
  # Indicator for NPS given inclusion in stacked sample
  z <- rep(1:0, c(n1, n0))
  
  # If pi_R known for all (not typically used)
  if (!is.null(hat_pi_R)) {
    # Predict pi_B for i in S_B
    # First predict pi_z
    fit_pi_B <- pbart(x.train = samp_comb[, pred_covs_B, drop = FALSE],
                      y.train = z, ntree = 50L, nskip = 100L, ndpost = num_post)
    hat_pi_z <- fit_pi_B$prob.train.mean
    # Then convert to pi_B using CRISP
    hat_pi_B <- hat_pi_z[z == 1] * hat_pi_R * frame_R / 
      (frame_B * (1 - hat_pi_z[z == 1]))
    
    # Get distribution of hat_pi_B
    hat_pi_z_dist <- fit_pi_B$prob.train
    hat_pi_B_dist <- matrix(NA, nrow = num_post, ncol = n1)
    # hat_pi_B_dist <- vector(mode = "list", length = num_post)
    for (m in 1:num_post) {
      pi_z_m <- hat_pi_z_dist[m, ]
      pi_B_m <- pi_z_m[z == 1] * hat_pi_R * frame_R / 
        (frame_B * (1 - pi_z_m[z == 1]))
      # Restrict to [0, 1]
      hat_pi_B_dist[m, ] <- sapply(1:n1, function(x) max(min(pi_B_m[x], 1), 0))
    }
    # Add BART model to output
    est_weights$fit_pi_B <- fit_pi_B
    
    # Get pseudo-weights
    w_B <- 1 / hat_pi_B
    w_B_post <- 1 / hat_pi_B_dist
    
  # If pi_R not known for all, predict pi_R for those in S_B 
  } else {
    
    if (pred_model == "bart") {
      
      ### Predict logit(pi_R) for i in S_B using wbart
      
      fit_pi_R <- wbart(x.train = samp_comb[z == 0, pred_covs_R, drop = FALSE],
                        y.train = logit_pi_R,
                        x.test = samp_comb[z == 1, pred_covs_R, drop = FALSE],
                        ntree = 50L, nskip = 100L, ndpost = num_post)
      hat_logit_pi_R <- fit_pi_R$yhat.test.mean  # length n1
      # Convert to pi_R scale using expit
      hat_pi_R <- exp(hat_logit_pi_R) / (1 + exp(hat_logit_pi_R))  # length n1
      # Get distribution of hat_pi_R using expit: num_post x n1
      hat_logit_pi_R_dist <- fit_pi_R$yhat.test
      hat_pi_R_dist <- exp(hat_logit_pi_R_dist) / (1 + exp(hat_logit_pi_R_dist))
      se_pi_R   <- apply(hat_pi_R_dist, 2, sd)
        
            # # Binary bart model for predicting pi_R for i in S_B
            # fit_pi_R <- pbart(x.train = samp_comb[z == 0, pred_covs_R, drop = FALSE],
            #                   y.train = pi_R, 
            #                   x.test = samp_comb[z == 1, pred_covs_R, drop = FALSE],
            #                   ntree = 50L, nskip = 100L, ndpost = num_post)
            # hat_pi_R <- colMeans(fit_pi_R$yhat.test)
            # # Get distribution of hat_pi_R
            # hat_pi_R_dist <-  fit_pi_R$yhat.test
      
      
      ### Predict pi_B for i in S_B using pbart
      
      fit_pi_B <- pbart(x.train = samp_comb[, pred_covs_B, drop = FALSE],
                        y.train = z, 
                        ntree = 50L, nskip = 100L, ndpost = num_post)
      hat_pi_z <- fit_pi_B$prob.train.mean  # length n1 + n0
      
      
      ### Get hat_pi_B for i in S_B using CRISP formula 
      
      # Get distribution of hat_pi_B using joint draws of pi_z and pi_R
      hat_pi_z_dist <- fit_pi_B$prob.train  # num_post x (n1 + n0)
      hat_pi_B_dist <- matrix(NA, nrow = num_post, ncol = n1)
      for (m in 1:num_post) {
        pi_z_m <- hat_pi_z_dist[m, ]  # length n1 + n0
        pi_B_m <- pi_z_m[z == 1] * hat_pi_R_dist[m, ] * frame_R / 
          (frame_B * (1 - pi_z_m[z == 1]))
        # Restrict to [0, 1]
        hat_pi_B_dist[m, ] <- sapply(1:n1, function(x) max(min(pi_B_m[x], 1), 0))
      }
      
      # Get mean inclusion probability estimate for i in S_B
      hat_pi_B <- colMeans(hat_pi_B_dist)  # length n1
      # Get SE and CI for inclusion probabilities for each indiv across BART samples
      alpha <- 1 - ci_level
      se_pi_B   <- apply(hat_pi_B_dist, 2, sd)
      pi_B_low  <- apply(hat_pi_B_dist, 2, quantile, probs = alpha / 2)
      pi_B_high <- apply(hat_pi_B_dist, 2, quantile, probs = 1 - alpha / 2)
      
      # Posterior mean pseudo-weights
      w_B <- 1 / hat_pi_B 
      
      # Pseudo-weight distribution: note this is not centered at w_B due to inverse
      w_B_post <- 1 / hat_pi_B_dist          # num_post x n1
      se_w_B <- apply(w_B_post, 2, sd)
      w_B_low <- apply(w_B_post, 2, quantile, probs = alpha / 2)
      w_B_high <- apply(w_B_post, 2, quantile, probs = 1 - alpha / 2)

      
    } else if (pred_model == "glm") {
      
      ### Predict logit(pi_R) for i in S_B using regression
      
      # Add logit(pi_R) to data for those in S_R
      glm_data <- cbind(samp_comb[z==0, , drop = FALSE], logit_pi_R)
      form_R <- as.formula(paste0("logit_pi_R ~ ", 
                                   paste0(pred_covs_R, collapse = " + ")))
      fit_pi_R <- lm(form_R, data = glm_data)
      
      # Get predictions on logit scale, including SE
      pred_R <- predict(fit_pi_R, 
                        newdata = samp_comb[z == 1, pred_covs_R, drop = FALSE],
                        se.fit = TRUE)
      hat_logit_pi_R <- as.numeric(pred_R$fit)    
      se_logit_pi_R <- as.numeric(pred_R$se.fit)  
      
      # Convert to pi_R scale using expit
      hat_pi_R <- exp(hat_logit_pi_R) / (1 + exp(hat_logit_pi_R))
      # Use delta method to get SE of pi_R
      se_pi_R <- hat_pi_R * (1 - hat_pi_R) * se_logit_pi_R
      # Distribution of hat_pi_R set to mean
      hat_pi_R_dist <- hat_pi_R
      
      
      ### Predict pi_B for i in S_B using logistic regression instead of BART
      
      # Add z to data
      glm_data2 <- cbind(samp_comb, z)
      form_B <- as.formula(paste0("z ~ ", 
                                  paste0(pred_covs_B, collapse = " + ")))
      # Get hat_logit_pi_z w/ SE
      fit_pi_B <- glm(form_B, data = glm_data2, family = binomial())
      # Convert to pi_z scale, and get SE of pi_z (uses Delta method internally)
      pred_B <- predict(fit_pi_B, type = "response", se.fit = TRUE)
      hat_pi_z <- pred_B$fit
      se_pi_z <- pred_B$se.fit
      # Restrict to S_B
      hat_pi_z_B <- hat_pi_z[z == 1]
      se_pi_z_B  <- se_pi_z[z == 1]
      
      # Get pi_B using CRISP: mean inclusion probabilities for i in S_B
      hat_pi_B <- hat_pi_z[z == 1] * hat_pi_R * frame_R / 
        (frame_B * (1 - hat_pi_z[z == 1]))
      # Dummy distribution of hat_pi_B set to mean: for BART comparison
      hat_pi_B_dist <- hat_pi_B
      
      # Get SE of pi_B using delta method
      k  <- frame_R / frame_B
      pz <- hat_pi_z_B
      pR <- hat_pi_R
      # Gradients wrt pz and pR
      # d g / d pz = k * pR / (1 - pz)^2
      # d g / d pR = k * pz / (1 - pz)
      dg_dpz <- k * pR / (1 - pz)^2
      dg_dpR <- k * pz / (1 - pz)
      # Get SE of pi_B
      se_pi_B <- sqrt((dg_dpz^2) * (se_pi_z_B^2) + (dg_dpR^2) * (se_pi_R^2))
      # 95% CI for pi_B, truncated to [0, 1]
      alpha <- 1 - ci_level
      pi_B_low  <- pmax(0, hat_pi_B - abs(qnorm(alpha/2)) * se_pi_B)
      pi_B_high <- pmin(1, hat_pi_B + abs(qnorm(alpha/2)) * se_pi_B)
      
      # Get pseudo-weights and their SE and CI (using Delta method)
      w_B <- 1 / hat_pi_B   # Mean pseudo-weights
      w_B_post <- 1 / hat_pi_B_dist  # dummy, uses mean: for BART comparison
      # d w / d pi_B = -1 / pi_B^2
      se_w_B <- sqrt((1 / hat_pi_B^2)^2 * se_pi_B^2)
      # Because 1/x is decreasing, CI for w_B uses reversed pi_B bounds
      w_B_low  <- 1 / pi_B_high
      w_B_high <- 1 / pi_B_low
      
    } else {
      stop("pred_model must be one of `'bart'` or `'glm'`")
    }
    
    
    ### Add model to output
    est_weights$fit_pi_R <- fit_pi_R
    est_weights$fit_pi_B <- fit_pi_B
    
    # pi_R
    est_weights$hat_pi_R   <- hat_pi_R
    est_weights$se_pi_R    <- se_pi_R
    est_weights$hat_pi_R_dist <- hat_pi_R_dist
    
    # pi_B
    est_weights$hat_pi_B   <- hat_pi_B
    est_weights$se_pi_B    <- se_pi_B
    est_weights$pi_B_low   <- pi_B_low
    est_weights$pi_B_high  <- pi_B_high
    est_weights$hat_pi_B_dist <- hat_pi_B_dist
    
    # weights 
    est_weights$w_B      <- w_B
    est_weights$se_w_B   <- se_w_B
    est_weights$w_B_low  <- w_B_low
    est_weights$w_B_high <- w_B_high
  }
  
  #===== Form full weights + weight trimming
  
  # Weights for stacked sample
  weights <- numeric(length = nrow(samp_comb))
  weights[z == 1] <- w_B  # Pseudo-weights for i in S_B (NPS)
  weights[z == 0] <- 1 / pi_R  # (Known) weights for i in S_R (PS)
  # Transpose posterior samples for weights
  wts_post <- t(w_B_post)
  
  ### Weight trimming
  wts <- trim_w(w_B, m = trim_method, c = trim_c, max = 10)
  wts_post <- apply(wts_post, 2, function(x) trim_w(x, m = trim_method, 
                                                    c = trim_c, max = 10)) 
  
  # # Normalize weights to sum to N
  # pi_B_m <- pi_B_m * N / sum(1 / pi_B_m)
  # # Normalize to sum to sample size??????????
  # wts_post <- apply(wts_post, 2, function(x) x / sum(x / n1))
  
  est_weights <- c(list(wts = wts, wts_post = wts_post, 
                        hat_pi_B = hat_pi_B, hat_pi_B_dist = hat_pi_B_dist,
                        hat_pi_R = hat_pi_R), est_weights)
  
  
  return(est_weights)
}


### Weight trimming function adapted from Rafei et al. (2020)
# Inputs:
#   w: Vector of untrimmed weights
#   m: Character specifying trimming method. Either 't1' for the entropy 
# procedure or 't2' for the median plus a multiple of the interquartile range
#   c: Integer constant indicating multiple for the cutoff, with larger values 
# allowing for more extreme weights.
#   max: Maximum number of recursive iterations allowed for checking that all 
# weights fall below the cutoff
trim_w <- function(w, m='t1', c=10, max=Inf){
  if(m=='t1'){
    kn <- sqrt(c*mean(w^2, na.rm=T))  # cutoff
    i <- 0
    while(sum(w>kn & !is.na(w))>=1 & i<max){
      s <- sum(w, na.rm=T)
      w[w>kn & !is.na(w)] <- kn
      w[!is.na(w)] <- w[!is.na(w)]*s/sum(w, na.rm=T)
      kn <- sqrt(c*sum(w^2, na.rm=T)/sum(!is.na(w), na.rm=T))
      i <- i+1
    }
  } else if(m=='t2'){
    kn <- median(w, na.rm=T)+c*IQR(w, na.rm=T)  # cutoff
    i <- 0
    while(sum(w>kn, na.rm=T)>=1 & i<max){
      s <- sum(w, na.rm=T)
      w[w>kn & !is.na(w)] <- kn
      w[!is.na(w)] <- w[!is.na(w)]*s/sum(w, na.rm=T)
      kn <- median(w, na.rm=T)+c*IQR(w, na.rm=T)
      i <- i+1
    }
  } 
  # else if(m=='t3'){
  #   a <- 2+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T))
  #   b <- (sum(w, na.rm=T)-1)*(1+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T)))
  #   kn <- 1/(length(w[!is.na(w)])*qbeta(0.99, a, b))
  #   i <- 0
  #   while(sum(w>kn, na.rm=T)>=1 & i<max){
  #     s <- sum(w, na.rm=T)
  #     w[w>kn & !is.na(w)] <- kn
  #     w[!is.na(w)] <- w[!is.na(w)]*s/sum(w, na.rm=T)
  #     a <- 2+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T))
  #     b <- (sum(w, na.rm=T)-1)*(1+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T)))
  #     kn <- 1/(length(w[!is.na(w)])*qbeta(0.99, a, b))
  #     i <- i+1
  #   }
  # }
  w	
}


#' Get mode
#' 
#' `get_mode` is a helper function that obtains the modal value given an input 
#' vector.
#' @param v Input vector
#' @return Outputs most common value found in input vector `v`
#' @keywords internal
#' @export
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


get_theta_dist_wolcan <- function(est_theta, true_theta, est_K, true_K, 
                                  subset, dist_type) {
  
  ### First get minimum distance using full vectors (with additional 0's)
  ### to get optimal ordering
  # Find all permutations of est_theta and true_theta with filler 0's
  all_perms <- gtools::permutations(n = dim(est_theta)[2], r = dim(true_theta)[2])
  # Obtain vector of mean absolute distance between est and true theta, 
  # calculated for each permutation
  dist_all_perms <- numeric(nrow(all_perms))
  for (i in 1:nrow(all_perms)) {
    est_theta_perm <- est_theta[,all_perms[i, ],]
    dist_all_perms[i] <- get_dist_wolcan(est_theta_perm, true_theta, 
                                         dist_type = dist_type) 
  }
  # Obtain optimal ordering of classes
  order <- all_perms[which.min(dist_all_perms), ]
  est_theta_perm <- est_theta[ , order, ]
  
  # Initialize subset ordering of classes
  order_sub_est <- order
  order_sub_true <- 1:true_K
  # Lowest dist out of all permutations
  theta_dist <- min(dist_all_perms)
  
  ### Option to use this minimum distance, or calculate minimum distance after
  ### subsetting (subset == TRUE)
  if (subset) {   # Calculate distance after subsetting
    if (est_K < true_K) {  # If missing a true class
      theta_sub <- get_subset_dist_wolcan(large_par = true_theta[, 1:true_K, ], 
                                          small_par = est_theta[, 1:est_K, ], 
                                          param_name = "theta", dist_type = dist_type)
      theta_dist <- theta_sub$par_dist         # Distance
      order_sub_true <- theta_sub$order_large  # Subset order for true_theta
      order_sub_est <- 1:est_K                 
    } else if (true_K < est_K) {  # If extra class
      theta_sub <- get_subset_dist_wolcan(large_par = est_theta[, 1:est_K, ], 
                                          small_par = true_theta[, 1:true_K, ],
                                          param_name = "theta", dist_type = dist_type)
      theta_dist <- theta_sub$par_dist
      order_sub_est <- theta_sub$order_large  # Subset order for est_theta
      order_sub_true <- 1:true_K
    } else {  # true_K == est_K
      # Lowest dist out of all permutations
      theta_dist <- min(dist_all_perms)
      order_sub_est <- order
      order_sub_true <- 1:true_K
    }
  }
  
  ### Return dist, ordering, reordered estimate, and subsetted orderings
  return(list("theta_dist" = theta_dist, "order" = order, 
              "est_theta_perm" = est_theta_perm, 
              "order_sub_est" = order_sub_est,
              "order_sub_true" = order_sub_true))
}

get_dist_wolcan <- function(par1, par2, dist_type = "mean_abs") {
  if (dist_type == "mean_abs") {  # Mean absolute error
    dist <- mean(abs(par1 - par2))
  } else if (dist_type == "sum_sq") {  # Frobenius norm / squared Euclidean norm
    dist <- sum((par1 - par2)^2)
  } else if (dist_type == "mean_sq") {  # MSE
    dist <- mean((par1 - par2)^2)
  } else {
    stop("Error: dist_type must be 'mean_abs', 'sum_sq', or 'mean_sq' ")
  }
  return(dist)
}

get_subset_dist_wolcan <- function(large_par, small_par, param_name, dist_type) {
  if (param_name == "theta") {
    large_K <- dim(large_par)[2]   ## Change to handle 0's
    sum(large_par[1, , 1] != 0)
    small_K <- dim(small_par)[2]
  } else if (param_name == "pi") {
    large_K <- length(large_par)
    small_K <- length(small_par)
  } else if (param_name == "xi") {
    large_K <- dim(large_par)[1] 
    small_K <- dim(small_par)[1]
  } else {
    stop("Error: 'param_name' must be either 'theta', 'pi', or 'xi'")
  }
  # Find all subsets of large_K with size equal to small_K
  sub_perms <- gtools::permutations(n = large_K, r = small_K)
  # Obtain dist (Frobenius norm) between large_par and small_par per permutation
  dist_sub_perms <- numeric(nrow(sub_perms))
  for (i in 1:nrow(sub_perms)) {
    if (param_name == "theta") {
      large_par_sub <- large_par[ , sub_perms[i, ], ]
    } else if (param_name == "pi") {
      large_par_sub <- large_par[sub_perms[i, ]]
    } else if (param_name == "xi") {
      large_par_sub <- large_par[sub_perms[i, ], ]
    }
    dist_sub_perms[i] <- get_dist_wolcan(small_par, large_par_sub, "mean_abs")
  }
  # Lowest dist out of all permutations
  par_dist <- min(dist_sub_perms)
  # Ordering corresponding to lowest dist
  order_large <- sub_perms[which.min(dist_sub_perms), ]
  
  return(list(par_dist = par_dist, order_large = order_large))
}