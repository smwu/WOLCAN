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
#   adjust: Boolean specifying whether to incorporate the post-processing variance 
# adjustment to ensure proper uncertainty intervals. Default is `TRUE`.
#   adapt_seed: Default is 1.
#   fixed_seed: Default is 1.
#   run_adapt: Boolean specifying whether the adaptive sampler should be run to 
# determine the number of classes, K. Default is `TRUE`. If `FALSE`, `K_fixed`
# must be specified.
#   alpha_fixed, eta_fixed: only specify if K_fixed is not `NULL`. Otherwise, 
# defaults to alpha = K, eta = 1.
#   MI: Boolean indicating if multiple imputation procedure should be performed
# for variance estimation that incorporates variability in the pseudo-weights. 
# Default is `TRUE`.
#   num_reps: Number of bootstrap replicates for the WS variance adjustment
# Outputs:
# 
wolcan <- function(x_mat, dat_B, dat_R, pred_covs_B, pred_covs_R, pi_R, 
                   hat_pi_R = NULL, num_post = 1000, frame_B = 1, frame_R = 1,
                   trim_method = "t2", trim_c = 20, 
                   D = 10, parallel = TRUE, n_cores = 4, 
                   wts_adj = c("MI", "WS all", "WS mean", "none"), 
                   adjust = TRUE, tol = 1e-8, num_reps = 100,
                   run_adapt = TRUE, K_max = 30, adapt_seed = 1, 
                   K_fixed = NULL, fixed_seed = 1, class_cutoff = 0.05,
                   n_runs = 20000, burn = 10000, thin = 5, update = 1000,
                   save_res = TRUE, save_path = NULL,
                   alpha_adapt = NULL, eta_adapt = NULL,
                   alpha_fixed = NULL, eta_fixed = NULL) {
  
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
  print("Getting pseudo-weights...")
  est_weights <- get_weights_bart(dat_B = dat_B, dat_R = dat_R, 
                                  pred_covs_B = pred_covs_B, 
                                  pred_covs_R = pred_covs_R, 
                                  num_post = num_post, pi_R = pi_R, 
                                  hat_pi_R = hat_pi_R,
                                  frame_B = frame_B, frame_R = frame_R,
                                  trim_method = trim_method, trim_c = trim_c)
  # Mean estimated weights
  wts <- est_weights$wts
  # Posterior distribution of weights
  wts_post <- est_weights$wts_post
  
  # Save weights
  if (save_res) {
    save(est_weights, file = paste0(save_path, "_wolcan_weights.RData"))
  }
  
  #================ Run WOLCA model with mean weights ==========================
  if (run_adapt) {
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
    # Save space
    res$MCMC_out <- NULL
  }
  
  # Define K, depending on if adaptive sampler results exist
  if (run_adapt) {
    K <- res$estimates$K_red
  } else {
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
    comb_estimates <- get_var_dir(D = D, wts_post = wts_post, K = K, 
                                  tol = tol, parallel = parallel, 
                                  n_cores = n_cores, adjust = adjust, 
                                  x_mat = x_mat, fixed_seed = fixed_seed,
                                  class_cutoff = class_cutoff, n_runs = n_runs, 
                                  burn = burn, thin = thin, update = update, 
                                  alpha_fixed = alpha_fixed, eta_fixed = eta_fixed)
    
    # Add variance estimates to output list
    res$estimates_adjust <- comb_estimates
  
  } else {    # Run fixed sampler
    
    print("Running fixed sampler...")
    
    # Set seed
    if (!is.null(fixed_seed)) {
      set.seed(fixed_seed)
    }
    set.seed(18)
    
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
  
  #================= Save and return output ====================================
  # Stop runtime tracker
  runtime <- Sys.time() - start_time
  res$runtime <- runtime
  
  # # Store estimated pseudo-weights
  # res$est_weights <- est_weights
  
  # Store additional data variables used for estimating pseudo-weights
  res$data_vars <- c(res$data_vars, list(dat_B = dat_B, dat_R = dat_R, 
                                         pred_covs_B = pred_covs_B, 
                                         pred_covs_R = pred_covs_R, pi_R = pi_R))
  
  class(res) <- "wolcan"
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_wolcan_results.RData"))
  }
  
  return(res)
  
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
get_var_adj <- function(D, K, res, wts_post, tol, num_reps = 100, 
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
    # Draw from weights posterior closest to quantile
    draw <- which.min(abs(col_meds - med_quants[d]))
    w_all <- c(w_all_post[, draw])
    wts_draws[, d] <- wts_post[, draw]
    
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
get_var_dir <- function(D = 10, wts_post, K, x_mat, fixed_seed = 1, 
                        tol = 1e-8, parallel = TRUE, n_cores = 4,
                        adjust = TRUE, class_cutoff = 0.05, n_runs = 20000, 
                        burn = 10000, thin = 5, update = 1000,
                        alpha_fixed = NULL, eta_fixed = NULL) {
  
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
                         x_mat = x_mat, K = K, adjust = adjust, 
                         class_cutoff = class_cutoff, n_runs = n_runs, 
                         burn = burn, thin = thin, update = update,
                         fixed_seed = fixed_seed, alpha_fixed = alpha_fixed, 
                         eta_fixed = eta_fixed, D = D)
    # Shutdown cluster
    stopCluster(cluster)
    
  } else {  # Serial version
    # Run wolca_d serially
    res_all <- lapply(1:D, wolca_d, wts_post = wts_post, x_mat = x_mat, K = K, 
                      adjust = adjust, class_cutoff = class_cutoff, 
                      n_runs = n_runs, burn = burn, thin = thin, update = update,
                      fixed_seed = fixed_seed, alpha_fixed = alpha_fixed, 
                      eta_fixed = eta_fixed, D = D)
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
  
  # Clustering to harmonize classes
  # Cluster individuals into reduced number of classes using agglomerative clustering
  # Calculate pairwise distance matrix using Hamming distance: proportion of
  # iterations where two individuals have differing class assignments
  distMat <- e1071::hamming.distance(t(c_all_draws))
  # Hierarchical clustering dendrogram
  dendrogram <- stats::hclust(stats::as.dist(distMat), method = "complete") 
  # Group individuals into K classes
  red_c_all <- stats::cutree(dendrogram, k = K)     
  # Modify classes if any classes are less than the cutoff percentage
  class_prop <- prop.table(table(red_c_all))
  if (any(class_prop < class_cutoff)) {
    # Get classes that are too small
    small <- which(class_prop < class_cutoff)
    # Group individuals into a larger number of classes 
    red_c_all_temp <- stats::cutree(dendrogram, k = K + length(small))
    red_c_all <- red_c_all_temp
    class_prop_temp <- prop.table(table(red_c_all_temp))
    # Get updated classes that are too small
    small_temp <- sort(which(class_prop_temp < class_cutoff))
    for (small_c in 1:length(small_temp)) {
      c_ind <- small_temp[small_c]
      class_small <- which(red_c_all_temp == c_ind)
      # Get nearest class
      inds <- 1:length(class_prop_temp)
      class_dist <- sapply(inds, function(x) 
        mean(distMat[class_small, which(red_c_all_temp == x)]))
      # Set small class distance to Inf
      class_dist[small_temp] <- Inf
      nearest <- which.min(class_dist[-c_ind])
      red_c_all[red_c_all_temp == c_ind] <- nearest
    }
    class_prop <- prop.table(table(red_c_all))
  }
  
  # Get unique reduced classes to aid relabeling
  unique_red_classes <- unique(red_c_all)
  
  # For each draw, relabel new classes using the most common old class assignment
  pi_red_new <- vector("list", length = D_red)
  theta_red_new <- vector("list", length = D_red)
  for (d in 1:D_red) {
    # Initialize pi and theta
    pi <- matrix(NA, nrow = M, ncol = K)
    theta <- array(NA, dim = c(M, J, K, dim(theta_red_draws[[d]])[4]))
    
    # Get new order
    new_order <- numeric(K)
    for (k in 1:K) {
      # Individuals who were assigned to old class k
      old_k_inds <- which(c_all_draws[d, ] == k)
      # Most common new class assignment among those individuals
      mode_new <- get_mode(red_c_all[old_k_inds])
      new_order[k] <- mode_new
      # Reorder classes for all pi and theta estimates
      if (!is.na(mode_new)) {
        pi[, new_order[k]] <- pi_red_draws[[d]][, k]  # reorder
        theta[, , new_order[k], ] <- theta_red_draws[[d]][, , k, ]
      }
    }
    
    # pi <- pi_red_draws[[d]][, new_order, drop = FALSE]  # reorder
    pi <- pi / rowSums(pi, na.rm = TRUE)
    pi_red_new[[d]] <- pi
    # theta <- theta_red_draws[[d]][, , new_order, , drop = FALSE]  # reorder
    theta <- ifelse(theta < tol, tol, theta)  # prevent underflow
    for (m in 1:dim(theta)[1]) {  # normalize
      for (j in 1:dim(theta)[2]) {
        for (k in 1:dim(theta)[3]) {
          theta[m, j, k, ] <- theta[m, j, k, ] / sum(theta[m, j, k, ])
        }
      }
    }
    theta_red_new[[d]] <- theta
  }
  
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
                         K_red = K_red, K = K)
  return(comb_estimates)
}


# Function to run WOLCA using weights from one posterior draw
# Inputs:
#   d: Numeric index for posterior draw
#   wts_post: Posterior distribution of  estimated pseudo-weights
#   adjust: Boolean indicating if post-processing variance adjustment should be applied
# inherit parameters from get_var_dir
wolca_d <- function(d, wts_post, x_mat, K, adjust, class_cutoff, n_runs, burn, 
                    thin, update, fixed_seed, alpha_fixed, eta_fixed, D) {
  
  # # Source Rcpp functions from baysc package
  # Rcpp::sourceCpp(rcpp_path)
  
  print(paste0("Draw ", d))
  
  # Get quantiles of medians of weights posterior
  col_meds <- c(apply(wts_post, 2, median))
  cutoffs <- (seq(1, D, length.out = D) - 0.5) / D
  med_quants <- stats::quantile(x = col_meds, probs = cutoffs)
  
  # Draw from weights posterior closest to quantile
  draw <- which.min(abs(col_meds - med_quants[d]))
  wts_d <- c(wts_post[, draw])
      # draw <- sample(1:ncol(wts_post), size = 1)
      # wts_d <- c(wts_post[, draw])
  
  # Initialize fixed sampler hyperparameters
  if (is.null(alpha_fixed)) {
    alpha_fixed = rep(1, K)  # not sparsity-inducing
  }
  # Run fixed sampler
  res_d <- baysc::wolca(x_mat = x_mat, sampling_wt = wts_d, run_sampler = "fixed", 
                        K_fixed = K, fixed_seed = fixed_seed, 
                        alpha_fixed = alpha_fixed, eta_fixed = eta_fixed, 
                        class_cutoff = class_cutoff, n_runs = n_runs, burn = burn, 
                        thin = thin, update = update, save_res = FALSE)
  # Apply post-processing variance adjustment for pseudo-likelihood if desired
  if (adjust) {
    res_d <- baysc::wolca_var_adjust(res = res_d, num_reps = 100, save_res = FALSE, 
                              adjust_seed = fixed_seed)
    # Get estimates
    est_d <- res_d$estimates_adjust
    est_d$K_red <- res_d$estimates$K_red
  } else {
    # Get estimates
    est_d <- res_d$estimates
  }
  # Return the sampled weights
  est_d$wts <- wts_d
  
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
get_weights_bart <- function(dat_B, dat_R, pred_covs_B, pred_covs_R, 
                             num_post = 1000, pi_R, hat_pi_R = NULL, 
                             frame_B = 1, frame_R = 1, trim_method = "t2", 
                             trim_c = 20) {
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
  
  # If pi_R known for all
  if (!is.null(hat_pi_R)) {
    # Predict pi_B for i in S_B
    # First predict pi_z
    fit_pi_B <- pbart(x.train = samp_comb[, pred_covs_B, drop = FALSE],
                      y.train = z, ntree = 50L, nskip = 100L, ndpost = num_post)
    hat_pi_z <- fit_pi_B$prob.train.mean
    # Then convert to pi_B using CRISP
    hat_pi_B <- hat_pi_z[z == 1] * hat_pi_R * frame_R / 
      (frame_B * (1 - hat_pi_z[z == 1]))
  
    # # Logistic regression instead of BART
    # temp <- glm(z ~ samp_comb$A1 + samp_comb$A2, family = binomial())
    # temp_pi_z <- temp$fitted.values
    # temp_pi_B <- temp_pi_z[z == 1] * hat_pi_R * frame_R / 
    #   (frame_B * (1 - temp_pi_z[z == 1]))
    # mean(abs(temp_pi_B - sim_samp_B$true_pi_B))
    
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
    
  # If pi_R not known for all, predict pi_R for those in S_B 
  } else {
    # Predict logit(pi_R) for i in S_B
    # Change so that it only predicts for those in S_B rather than all
    fit_pi_R <- wbart(x.train = samp_comb[z == 0, pred_covs_R, drop = FALSE],
                      y.train = logit_pi_R,
                      x.test = samp_comb[z == 1, pred_covs_R, drop = FALSE],
                      ntree = 50L, nskip = 100L, ndpost = num_post)
    hat_logit_pi_R <- fit_pi_R$yhat.test.mean
    # Convert to pi_R scale using expit
    hat_pi_R <- exp(hat_logit_pi_R) / (1 + exp(hat_logit_pi_R))
    # Get distribution of hat_pi_R
    hat_logit_pi_R_dist <- fit_pi_R$yhat.test
    hat_pi_R_dist <- exp(hat_logit_pi_R_dist) / (1 + exp(hat_logit_pi_R_dist))
    
    # Predict pi_B for i in S_B
    fit_pi_B <- pbart(x.train = samp_comb[, pred_covs_B, drop = FALSE],
                      y.train = z, ntree = 50L, nskip = 100L, ndpost = num_post)
    hat_pi_z <- fit_pi_B$prob.train.mean
    hat_pi_B <- hat_pi_z[z == 1] * hat_pi_R * frame_R / 
      (frame_B * (1 - hat_pi_z[z == 1]))
    
    # Get distribution of hat_pi_B, incorporating distribution of hat_pi_R
    hat_pi_z_dist <- fit_pi_B$prob.train
    hat_pi_B_dist <- matrix(NA, nrow = num_post, ncol = n1)
    # hat_pi_B_dist <- vector(mode = "list", length = num_post)
    for (m in 1:num_post) {
      pi_z_m <- hat_pi_z_dist[m, ]
      pi_B_m <- pi_z_m[z == 1] * hat_pi_R_dist[m, ] * frame_R / 
        (frame_B * (1 - pi_z_m[z == 1]))
      # Restrict to [0, 1]
      hat_pi_B_dist[m, ] <- sapply(1:n1, function(x) max(min(pi_B_m[x], 1), 0))
    }
    
    # Add estimated hat_pi_R distribution to output list
    est_weights$hat_pi_R_dist <- hat_pi_R_dist
  }
  
  # Form pseudo-weights
  weights <- numeric(length = nrow(samp_comb))
  weights[z == 1] <- 1 / hat_pi_B
  weights[z == 0] <- 1 / pi_R
  wts <- weights[z == 1]
  
  ### Get distribution of weights for variance estimation
  wts_post <- t(1 / hat_pi_B_dist)
  
  ### Weight trimming
  wts <- trim_w(wts, m = trim_method, c = trim_c, max = 10)
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

