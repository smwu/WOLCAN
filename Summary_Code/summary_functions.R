#=================================================
# Functions for summarizing model output
# Author: Stephanie Wu
# Date created: 2024/05/29
# Date updated: 2024/05/29
#=================================================



save_scen_metrics <- function(scenario, samp_i_seq, WOLCAN = TRUE, WOLCA = TRUE, 
                              save_path, wd, data_dir, res_dir, subset = FALSE, 
                              dist_type = "mean_abs", parallel = FALSE) {
  
  # Get metrics for models
  metrics_all <- list()
  if (WOLCAN) {
    print("Getting WOLCAN results...")
    model <- "wolcan"
    metrics_wolcan <- get_metrics_wolcan(wd = wd, data_dir = data_dir, 
                                         res_dir = res_dir, scenario = scenario, 
                                         model = model, samp_i_seq = samp_i_seq,
                                         subset = subset, dist_type = dist_type,
                                         parallel = parallel, 
                                         save_path = save_path)
    metrics_all$metrics_wolcan <- metrics_wolcan
  } 
  if (WOLCA) {
    print("Getting WOLCA results...")
    model <- "wolca"
    metrics_wolca <- get_metrics_wolcan(wd = wd, data_dir = data_dir, 
                                         res_dir = res_dir, scenario = scenario, 
                                         model = model, samp_i_seq = samp_i_seq,
                                         subset = subset, dist_type = dist_type,
                                        parallel = parallel, 
                                        save_path = save_path)
    metrics_all$metrics_wolca <- metrics_wolca
  }
  
  # Save summary metrics
  save(metrics_all, 
       file = paste0(save_path, "summary.RData"))
  
}

get_metrics_wolcan <- function(wd, data_dir, res_dir, sum_dir, 
                               scenario, model, samp_i_seq, 
                               subset = FALSE, dist_type = "mean_abs",
                               parallel = FALSE, save_path) {
  
  #============== Load data and initialize variables ===========================
  # Load simulated population data
  pop_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData")
  load(pop_data_path)

  # Obtain true observed population parameters
  true_params <- get_true_params_wolcan(sim_pop = sim_pop)
  true_K <- as.vector(sim_pop$K)
  
  # Initialize variables
  L <- length(samp_i_seq)  # Number of samples
  J <- dim(sim_pop$true_global_thetas)[1]        # Number of exposure items
  runtime_all <- numeric(L)
  # Bias squared using posterior median
  K_all <- K_dist <- pi_dist <- theta_dist <- wts_dist <- rep(NA, L) 
  # Posterior variance
  pi_var_all <- theta_var_all <- rep(NA, L) 
  # Coverage variables
  pi_cover <- matrix(0, nrow=L, ncol=length(true_params$true_pi))
  theta_cover <- array(0, c(L, dim(true_params$true_theta)[c(1,2)]))
  # MSE for all iterations
  pi_mse_all <- theta_mse_all <- rep(NA, L)
  
  # Initialize plotting structures
  pi_all <- matrix(NA, nrow=L, ncol=true_K)
  theta_mode_all <- array(NA, dim=c(L, dim(sim_pop$true_global_patterns)))
  mode_mis_all <- rep(NA, L)
  
  
  #============== Get performance metrics for each iteration ===================
  if (parallel) {
    # Number of cores available
    n_cores <- detectCores()  
    # Create a cluster using sockets
    cluster <- parallel::makeCluster(n_cores)
    # Export global functions to be available to all workers
    clusterExport(cluster, varlist = c("get_true_params_wolcan", "get_dist_wolcan", 
                                       "get_theta_dist_wolcan", "get_pi_dist_wolcan",
                                       "get_subset_dist_wolcan", "abind"))
    # Get metrics in parallel
    summ_samp_all <- parLapply(cluster, samp_i_seq, get_metrics_wolcan_i, 
                               sim_pop = sim_pop, wd = wd, data_dir = data_dir, 
                               res_dir = res_dir, scenario = scenario, 
                               model = model, dist_type = dist_type, 
                               subset = subset, save_path = save_path)
    
    # Shutdown cluster
    stopCluster(cluster)
    
    # Get performance metrics for each sample iteration
    for (l in 1:L) { 
      
      summ_l <- summ_samp_all[[l]]
      # print(l)
      
      if (!is.null(summ_l)) {
        runtime_all[l] <- summ_l$runtime
        wts_dist[l] <- summ_l$wts_dist
        K_all[l] <- summ_l$K
        K_dist[l] <- summ_l$K_dist
        theta_dist[l] <- summ_l$theta_dist
        pi_dist[l] <- summ_l$pi_dist
        pi_cover[l, ] <- summ_l$pi_cover
        pi_var_all[l] <- summ_l$pi_var
        pi_mse_all[l] <- summ_l$pi_mse
        theta_cover[l, , 1:(dim(summ_l$theta_cover)[2])] <- summ_l$theta_cover
        theta_var_all[l] <- summ_l$theta_var
        theta_mse_all[l] <- summ_l$theta_mse
        mode_mis_all[l] <- summ_l$mode_mis
        
        # Handle extra estimated classes if necessary
        K_l <- length(summ_l$pi)
        if (true_K < K_l) { 
          # Expand estimated matrix and array sizes
          extra <- K_l - true_K
          # Expand pi_all
          filler_pi <- array(NA, dim=c(L, extra))
          pi_all <- abind::abind(pi_all, filler_pi, along = 2)
          # Expand theta_mode_all
          filler_theta <- array(NA, dim=c(L, J, extra))
          theta_mode_all <- abind::abind(theta_mode_all, filler_theta, along = 3)
        }
        
        theta_mode_all[l, , 1:(dim(summ_l$theta_mode)[2])] <- summ_l$theta_mode
        pi_all[l, 1:length(summ_l$pi)] <- summ_l$pi
      }
    }
  } else {
    
    # Get performance metrics for each sample iteration
    for (l in 1:L) { 
      samp_i <- samp_i_seq[l]
      summ_l <- get_metrics_wolcan_i(samp_i = samp_i,  sim_pop = sim_pop, 
                                     wd = wd, data_dir = data_dir, 
                                     res_dir = res_dir, scenario = scenario, 
                                     model = model, dist_type = dist_type, 
                                     subset = subset, save_path = save_path, 
                                     true_params = true_params)
      
      if (!is.null(summ_l)) {
        runtime_all[l] <- summ_l$runtime
        wts_dist[l] <- summ_l$wts_dist
        K_all[l] <- summ_l$K
        K_dist[l] <- summ_l$K_dist
        theta_dist[l] <- summ_l$theta_dist
        pi_dist[l] <- summ_l$pi_dist
        pi_cover[l, ] <- summ_l$pi_cover
        pi_var_all[l] <- summ_l$pi_var
        pi_mse_all[l] <- summ_l$pi_mse
        theta_cover[l, , 1:(dim(summ_l$theta_cover)[2])] <- summ_l$theta_cover
        theta_var_all[l] <- summ_l$theta_var
        theta_mse_all[l] <- summ_l$theta_mse
        mode_mis_all[l] <- summ_l$mode_mis
        
        # Handle extra estimated classes if necessary
        K_l <- length(summ_l$pi)
        if (true_K < K_l) { 
          # Expand estimated matrix and array sizes
          extra <- K_l - true_K
          # Expand pi_all
          filler_pi <- array(NA, dim=c(L, extra))
          pi_all <- abind::abind(pi_all, filler_pi, along = 2)
          # Expand theta_mode_all
          filler_theta <- array(NA, dim=c(L, J, extra))
          theta_mode_all <- abind::abind(theta_mode_all, filler_theta, along = 3)
        }
        
        theta_mode_all[l, , 1:(dim(summ_l$theta_mode)[2])] <- summ_l$theta_mode
        pi_all[l, 1:length(summ_l$pi)] <- summ_l$pi  
      }
    }
  }
 
  
  #============== Calculate bias^2 averaged over sample iterations =============
  wts_bias <- mean(wts_dist, na.rm = TRUE)
  K_bias <- mean(K_dist, na.rm = TRUE)
  pi_bias <- mean(pi_dist, na.rm = TRUE)
  theta_bias <- mean(theta_dist, na.rm = TRUE)
  
  # Calculated CI width, averaged across iterations
  pi_var <- mean(pi_var_all, na.rm = TRUE)
  theta_var <- mean(theta_var_all, na.rm = TRUE)
  
  # Calculate class-specific coverage, averaged across iterations
  # Coverage for pi
  pi_cover_avg <- colMeans(pi_cover, na.rm = TRUE)
  # Coverage for theta: average over food items
  theta_cover_avg <- colMeans(colMeans(theta_cover, na.rm = TRUE), na.rm = TRUE)
  
  runtime_avg <- mean(runtime_all, na.rm = TRUE)
  
  #============== Return results ===============================================
  ret_list <- list(wts_bias = wts_bias, K_bias = K_bias, pi_bias = pi_bias, 
                   pi_var = pi_var, theta_bias = theta_bias, theta_var = theta_var, 
                   pi_cover_avg = pi_cover_avg, theta_cover_avg = theta_cover_avg, 
                   runtime_avg = runtime_avg, wts_dist = wts_dist, 
                   K_dist = K_dist, pi_dist = pi_dist, theta_dist = theta_dist, 
                   pi_mse_all = pi_mse_all, theta_mse_all = theta_mse_all)
  
  ret_list[["pi_all"]] <- pi_all
  theta_mode <- apply(theta_mode_all, c(2,3), function(x) mean(x, na.rm = TRUE))
  mode_mis <- mean(mode_mis_all, na.rm = TRUE)
  ret_list[["theta_mode"]] <- theta_mode
  ret_list[["mode_mis"]] <- mode_mis
  ret_list[["mode_mis_all"]] <- mode_mis_all
  ret_list[["K_all"]] <- K_all
  
  return(ret_list)
}


# Function to get performance metrics for one iteration
# Inputs:
#   samp_i: sample iteration index
get_metrics_wolcan_i <- function(samp_i, sim_pop, wd, data_dir, res_dir, 
                                 scenario, model, dist_type, subset, save_path,
                                 true_params) {
  
  # Initialize return
  summ_i <- NULL
  
  # Check if summary already exists
  if (file.exists(paste0(save_path, "samp_", samp_i, "_", model, ".RData"))) {
    # Read in summary if it already exists
    print(paste0("Summary for samp ", samp_i, " already exists."))
    load(paste0(save_path, "samp_", samp_i, "_", model, ".RData"))
    
  } else {
    # Obtain true observed population parameters
          # # Need to reload so that filler dimensions do not keep adding over iterations
          # true_params <- get_true_params_wolcan(sim_pop = sim_pop)
    true_K <- as.vector(sim_pop$K)
    
    # Check that sample data file exists
    sim_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", 
                            samp_i, "_B_wolcan.RData")
    # Check that results file exists
    sim_res_path <- paste0(wd, res_dir, "scen_", scenario, "/samp_", samp_i, "_", 
                           model, "_results.RData")
    if (!file.exists(sim_data_path)) {
      print(paste0("File does not exist: ", sim_data_path))
    } else if (!file.exists(sim_res_path)) {
      print(paste0("File does not exist: ", sim_res_path))
    } else {
      print(samp_i)
      
      # Read in sample data
      load(sim_data_path)
      
      # Read in results data
      load(sim_res_path)
      runtime <- res$runtime
      
      # Check model
      if (model == "wolcan") {
        if (scenario == 5) {  # No weights adjustment procedure (2)
          estimates <- res$estimates  # From adaptive sampler
        } else {  
          estimates <- res$estimates_adjust
        }
      } else if (model == "wolca") {
        estimates <- res$estimates
      } else {
        stop("Error: model must be one of 'wolcan' or 'wolca'")
      }
      
      M <- dim(estimates$theta_red)[1]        # Number of MCMC iterations (stacked if wolcan)
      J <- dim(estimates$theta_red)[2]        # Number of exposure items
      R <- dim(estimates$theta_red)[4]        # Number of exposure levels
      K <- length(estimates$pi_med)           # Number of classes
      
      # If number of classes is incorrect, fill remaining components with 0's
      if (K > true_K) {
        # If there are extra estimated classes, add 0s to true parameters
        extra <- K - true_K
        true_params$true_pi <- c(true_params$true_pi, rep(0, extra))
        filler <- array(0, dim=c(dim(estimates$theta_med)[1], extra, 
                                 dim(estimates$theta_med)[3]))
        true_params$true_theta <- abind::abind(true_params$true_theta, filler, along = 2)
      } else if (K < true_K) {
        # If there are missing estimated classes, add 0s to estimated parameters
        missing <- true_K - K
        estimates$pi_med <- c(estimates$pi_med, rep(0, missing))
        filler <- array(0, dim=c(dim(estimates$theta_med)[1], missing, 
                                 dim(estimates$theta_med)[3]))
        estimates$theta_med <- abind::abind(estimates$theta_med, filler, along = 2)  
        
        # Add 0's to full MCMC outputs for the missing classes
        estimates$pi_red <- abind::abind(estimates$pi_red, array(0, dim=c(M, missing)), 
                                         along=2)
        estimates$theta_red <- abind::abind(estimates$theta_red, 
                                            array(0, dim=c(M, J, missing, R)), along = 3)
      }
      
      #============== Calculated mean absolute distance (abs bias) ===============
      ##### Posterior mean weights
      wts_dist <- get_dist_wolcan(par1 = res$data_vars$sampling_wt, 
                                  par2 = 1 / sim_samp_B$true_pi_B, 
                                  dist_type = dist_type)
      
      ##### Number of classes, K
      K_dist <- get_dist_wolcan(K, true_K, dist_type = dist_type)
      
      ##### theta: get dist (Eucl norm) and optimal ordering
      theta_perm <- get_theta_dist_wolcan(est_theta = estimates$theta_med, 
                                          true_theta = true_params$true_theta, 
                                          est_K = K, true_K = true_K, subset = subset,
                                          dist_type = dist_type)
      theta_dist <- theta_perm$theta_dist
      order <- theta_perm$order
      if (subset) {
        order_sub_est <- theta_perm$order_sub_est
        order_sub_true <- theta_perm$order_sub_true
        K_min <- length(order_sub_est)
      } else {
        order_sub_est <- order
        order_sub_true <- 1:length(order)
        K_min <- true_K
      }
      
      ##### pi 
      pi_perm <- get_pi_dist_wolcan(est_pi = estimates$pi_med, 
                                    true_pi = true_params$true_pi, order = order, 
                                    est_K = K, true_K = true_K, subset = subset,
                                    order_sub_est = order_sub_est, 
                                    order_sub_true = order_sub_true,
                                    dist_type = dist_type)
      pi_dist <- pi_perm$pi_dist
      
      #============== Calculate coverage and CI widths ===========================
      ##### pi
      # Obtain credible intervals for each of the K true clusters
      pi_CI <- apply(estimates$pi_red[, order_sub_est], 2, 
                     function(x) quantile(x, c(0.025, 0.975)))
      # Assign 1 if interval covers true value, 0 if not
      # If a class is missing, defaults to 0 (not covered)
      pi_cover <- numeric(length(order_sub_true))
      pi_cover[order_sub_true] <- ifelse(
        (true_params$true_pi[order_sub_true] >= pi_CI[1,]) & 
          (true_params$true_pi[order_sub_true] <= pi_CI[2,]), 1, 0)
      # Subset to the true number of classes for pi_cover
      pi_cover <- pi_cover[1:true_K]
      # CI width averaged over the components
      pi_var <- mean(apply(pi_CI, 2, diff))
      # MSE
      pi_mse <- mean(apply(estimates$pi_red, 1, function(x) 
        get_dist_wolcan(x[order_sub_est], true_params$true_pi[order_sub_true], 
                        "mean_sq")))
      
      ##### theta
      # Theta mode consumption levels for each item and class (pxK)
      est_modes <- apply(estimates$theta_med[, order_sub_est, ], c(1,2), which.max)
      true_modes <- apply(true_params$true_theta[, order_sub_true, ], c(1,2), 
                          which.max)
      # True modal probabilities for each item and class (pxK)
      true_theta_modal <- apply(true_params$true_theta[ , order_sub_true, ], 
                                c(1,2), max) 
      theta_var_temp <- numeric(K_min)
      # Initialize theta_cover
      theta_cover <- array(NA, dim = dim(true_params$true_theta)[c(1,2)])
      for (k in 1:K_min) {
        # Subset theta for cluster k
        est_theta_k <- estimates$theta_red[, , order_sub_est[k], ]
        # Each row provides the indices for one row of modal probabilities
        modal_idx <- cbind(rep(1:M, each = J), rep(1:J, times = M), 
                           rep(est_modes[, k], times = M))
        # estimated probabilities for the mode for cluster k (Mxp)
        est_theta_k_modal <- matrix(est_theta_k[modal_idx], ncol = J, byrow = TRUE)
        # Obtain credible intervals for each item 
        # Margins of apply are the dimensions that should be preserved
        theta_CI <- apply(est_theta_k_modal, 2, 
                          function(x) quantile(x, c(0.025, 0.975)))
        theta_cover[, order_sub_true[k]] <- ifelse(
          (true_theta_modal[, k] >= theta_CI[1, ]) &
            (true_theta_modal[, k] <= theta_CI[2, ]), 1, 0)
        # CI width measures variation in estimating the modes for each k,
        # averaged over the items
        theta_var_temp[k] <- mean(apply(theta_CI, 2, diff))
      }
      # Subset to the true number of classes for theta_cover
      theta_cover <- theta_cover[, 1:true_K]
      # CI width averaged over the classes
      theta_var <- mean(theta_var_temp)
      # MSE
      theta_mse <- mean(apply(estimates$theta_red, 1, function(x) 
        get_dist_wolcan(x[, order_sub_est, ], 
                        true_params$true_theta[, order_sub_true, ], "mean_sq")))
      
      
      #============== Parameter estimate plot outputs ============================
      ##### theta
      # Theta mode consumption levels for each item and class (pxK)
      est_modes <- apply(estimates$theta_med[, order, ], c(1,2), which.max)
      true_modes <- apply(true_params$true_theta[, order, ], c(1,2), 
                          which.max)
      # Get theta mode
      theta_mode <- array(NA, dim = c(J, max(K, true_K)))
      theta_mode[, 1:length(order)] <- est_modes
      # Mode mismatches
      mode_mis <- sum(abs(est_modes[, 1:true_K] - sim_pop$true_global_patterns))
      
      ##### pi
      pi <- numeric(length(order))
      pi[1:length(order)] <- estimates$pi_med[order]
      
      # Return performance metrics for the iteration
      summ_i <- list(runtime = runtime, wts_dist = wts_dist, K_dist = K_dist, 
                     theta_dist = theta_dist, pi_dist = pi_dist, 
                     pi_cover = pi_cover, pi_var = pi_var,  pi_mse = pi_mse, 
                     theta_cover = theta_cover, theta_var = theta_var, 
                     theta_mse = theta_mse, theta_mode = theta_mode, 
                     mode_mis = mode_mis, pi = pi, K = K)
      
      # Save summary metrics
      save(summ_i, file = paste0(save_path, "samp_", samp_i, "_", model, ".RData"))
    }
  }
  
  # Return results
  return(summ_i)
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

get_pi_dist_wolcan <- function(est_pi, true_pi, order, est_K, true_K, subset = TRUE,
                               order_sub_est, order_sub_true, dist_type) {
  ### Use input optimal ordering
  ### Option to subset when calculating minimum distance (subset == TRUE)
  if (!subset) {  # No subsetting
    pi_dist <- get_dist_wolcan(est_pi[order], true_pi, dist_type = dist_type)
  } else {  # Calculate distance after subsetting
    pi_dist <- get_dist_wolcan(est_pi[order_sub_est], true_pi[order_sub_true], 
                        dist_type = dist_type)
  }
  
  ### Return dist, ordering, and reordered estimate
  return(list("pi_dist" = pi_dist, "order" = order, "est_pi_perm" = est_pi[order]))
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

get_true_params_wolcan <- function(sim_pop) {
  # Get true pi using population data
  true_pi <- tabulate(sim_pop$c_all) / length(sim_pop$c_all)
  # Get true theta
  theta_dim <- dim(sim_pop$true_global_thetas)
  true_theta <- array(NA, dim=theta_dim)
  for (j in 1:theta_dim[1]) {
    for (k in 1:theta_dim[2]) {
      for (r in 1:theta_dim[3]) {
        true_theta[j,k,r] <- sum((sim_pop$X_data[,j]==r) & 
                                   (sim_pop$c_all==k)) / sum(sim_pop$c_all==k) 
      }
    }
  }
  return(list(true_pi = true_pi, true_theta = true_theta))
}

#==================== Tables ===================================================
create_app_tables_wolcan <- function(save_paths, scenarios, scen_names, 
                                     overall_name, format = "latex", 
                                     digits = 3, WOLCA = TRUE, WOLCAN = TRUE) {
  num_scen <- length(scenarios)
  # models depending on WOLCA and WOLCAN true/false
  model <- list()
  if (WOLCA) {
    model <- c(model, "Unweighted")
  }
  if (WOLCAN) {
    model <- c(model, "WOLCAN")
  }
  model <- unlist(model)
  # multiplier for number of rows depending on models
  mult <- length(model)
  
  metrics_wolcan_df <- as.data.frame(matrix(NA, nrow = mult*length(scenarios), 
                                       ncol = 10))
  colnames(metrics_wolcan_df) <- c(overall_name, "Model", "Weights Abs Bias",
                              "$K$ Abs Bias", "$\\pi$ Abs Bias",  
                              "$\\theta$ Abs Bias", 
                              "$\\pi$ CI Width", "$\\theta$ CI Width", 
                              "$\\pi$ Coverage","$\\theta$ Coverage")
  metrics_wolcan_df[, 1] <- rep(scen_names, each = mult)
  metrics_wolcan_df[, 2] <- rep(model, num_scen)  
  # output_inds <- 1:7
  output_inds <- c(1, 2, 3, 5, 4, 6)
  row_ind <- 1
  for (i in 1:num_scen) {
    save_path <- save_paths[i]
    load(paste0(save_path, "summary.RData"))
    if (WOLCA) {
      metrics_wolcan_df[row_ind, -c(1,2)] <- 
        c(metrics_all$metrics_wolca[output_inds], 
          mean(metrics_all$metrics_wolca$pi_cover_avg), 
          mean(metrics_all$metrics_wolca$theta_cover_avg))
      row_ind <- row_ind + 1
    }
    if (WOLCAN) {
      metrics_wolcan_df[row_ind, -c(1,2)] <- 
        c(metrics_all$metrics_wolcan[output_inds], 
          mean(metrics_all$metrics_wolcan$pi_cover_avg), 
          mean(metrics_all$metrics_wolcan$theta_cover_avg))
      row_ind <- row_ind + 1
    }
  }
  
  metrics_wolcan_df %>% 
    kbl(digits = digits, align = "rrrrrrrrrrrr", booktabs = TRUE, format = format,
        caption = "Summary of mean absolute bias, 95% credible interval width, and coverage for simulations based on posterior samples.") %>%
    kable_classic() %>%
    kable_styling(full_width = FALSE)
}


#========================= Plot theta ==========================================

plot_theta_patterns_wolcan <- function(wd, data_dir, scenario, save_path) {
  # Load summary
  load(paste0(save_path, "summary.RData"))
  
  # Load true population data
  pop_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData")
  load(pop_data_path)
  # Obtain true observed population parameters
  true_params <- get_true_params_wolcan(sim_pop = sim_pop)
  
  p_true <- theta_mode_plot(sim_pop$true_global_patterns, "True Classes") + 
    guides(fill = guide_legend(reverse = FALSE)) +
    labs(fill = "Modal θ Level")
  p_unwt <- theta_mode_plot(metrics_all$metrics_wolca$theta_mode, "Unweighted Classes")
  p_WOLCAN <- theta_mode_plot(metrics_all$metrics_wolcan$theta_mode, "WOLCAN Classes")
  p_comb <- ggarrange(p_true, 
                      p_unwt + theme(axis.title.y = element_blank()), 
                      p_WOLCAN + theme(axis.title.y = element_blank()), 
                      nrow = 1, common.legend = TRUE, legend = "top")
  return(p_comb)
}

theta_mode_plot <- function(theta_plot_data, x_label) {
  p <- dim(theta_plot_data)[1]
  K <- dim(theta_plot_data)[2]
  Item <- factor(as.character(1:p), levels = as.character(p:1))
  theta_plot <- data.frame(theta_plot_data, Item)
  colnames(theta_plot) <- c(1:K, "Item")
  theta_plot <- theta_plot %>% gather("Class", "Level", 1:K) 
  patterns <- ggplot(theta_plot, aes(x=Class, y=Item, fill=Level)) + 
    theme_classic() +
    xlab(x_label) +
    geom_tile(color="gray") + 
    geom_text(aes(label = round(Level,2)), col="white", cex=2.5) +
    scale_fill_gradient(trans = "reverse")
  return(patterns)
}

#===================== Plot pi =================================================

plot_pi_patterns_wolcan <- function(wd, data_dir, scenario, samp_i_seq, 
                                    save_path, y_lim = c(0,1)) {
  # Load true population data
  pop_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData")
  load(pop_data_path)
  # Obtain true observed population parameters
  true_params <- get_true_params_wolcan(sim_pop = sim_pop)
  true_pi <- true_params$true_pi
  true_K <- as.vector(sim_pop$K)
  samp_pi <- get_avg_over_samps_wolcan(wd = wd, data_dir = data_dir, 
                                       scenario = scenario, 
                                       samp_i_seq = samp_i_seq)$avg_samp_pi
  
  # Load summary
  load(paste0(save_path, "summary.RData"))
  
  # Load simulated sample data
  L <- length(metrics_all$metrics_wolcan$K_dist)
  
  pi_plot_data <- as.data.frame(rbind(metrics_all$metrics_wolca$pi_all[, 1:true_K],
                                      metrics_all$metrics_wolcan$pi_all[, 1:true_K]))
  colnames(pi_plot_data) <- 1:ncol(pi_plot_data)
  pi_plot_data$Model <- c(rep("Unweighted", times=L),
                          rep("WOLCAN", times=L))
  pi_plot_data <- pi_plot_data %>% gather("pi_component", "value", -Model)
  ggplot(pi_plot_data, aes(x=pi_component, y=value, fill=Model)) +
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() +
    geom_segment(mapping=aes(x=0.5, xend=1.5, y=true_params$true_pi[1],
                             yend=true_params$true_pi[1]),color="red") +
    geom_segment(mapping=aes(x=1.5, xend=2.5, y=true_params$true_pi[2],
                             yend=true_params$true_pi[2]),color="red") +
    geom_segment(mapping=aes(x=2.5, xend=3.5, y=true_params$true_pi[3],
                             yend=true_params$true_pi[3]),color="red") +
    geom_segment(mapping=aes(x=0.5, xend=1.5, y=samp_pi[1], yend=samp_pi[1]),
                 color="red", linetype = "dashed") +
    geom_segment(mapping=aes(x=1.5, xend=2.5, y=samp_pi[2], yend=samp_pi[2]),
                 color="red", linetype = "dashed") +
    geom_segment(mapping=aes(x=2.5, xend=3.5, y=samp_pi[3], yend=samp_pi[3]),
                 color="red", linetype = "dashed") +
    # ggtitle(paste0("Parameter estimation for π across samples")) + 
    xlab("Latent Class") + ylab("π Value") + 
    theme(legend.position="top") 
}


get_avg_over_samps_wolcan <- function(wd, data_dir, scenario, samp_i_seq) {
  # Load true population data
  pop_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData")
  load(pop_data_path)
  # Obtain true observed population parameters
  true_params <- get_true_params_wolcan(sim_pop = sim_pop)
  
  L <- length(samp_i_seq)
  samp_pis <- array(NA, dim = c(length(samp_i_seq), length(sim_data$true_pi)))
  for (l in 1:L) { # For each sample iteration
    samp_i = samp_i_seq[l]
    sim_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", 
                            samp_i, "_B_wolcan.RData")
    load(sim_data_path)
    samp_pis[l, ] <- tabulate(sim_samp_B$c_all) / length(sim_samp_B$c_all)
  }
  avg_samp_pi <- apply(samp_pis, 2, function(x) mean(x, na.rm = TRUE))

  return(list(avg_samp_pi = avg_samp_pi))
}


# # Compare distributions
# hist(sim_samp_B$true_pi_B, breaks = 30, col = rgb(0,0,1,1/4), ylim = c(0, 100),
#      xlab = "true_pi_B (blue) vs. hat_pi_B (red)")
# hist(hat_pi_B, breaks = 30, col = rgb(1,0,0,1/4), add = TRUE)
# hist(hat_pi_z, breaks = 20)
# mean(abs(hat_pi_B - sim_samp_B$true_pi_B))
# 
# hist(hat_pi_B_dist, breaks = 30)
# hist(apply(hat_pi_B_dist, 1, max), breaks = 30)
# hist(apply(hat_pi_B_dist, 1, function(x) sum(x > 1)), breaks = 30)
# 
# # Compare distributions
# hist(sim_samp_B$true_pi_B, breaks = 30, col = rgb(0,0,1,1/4), ylim = c(0, 100),
#      xlab = "true_pi_B (blue) vs. hat_pi_B (red)")
# hist(hat_pi_B, breaks = 30, col = rgb(1,0,0,1/4), add = TRUE)
# 
# 
# 
# get_summ_stats <- function(res, res_true = NULL) {
#   if (!is.null(res$estimates_adjust)) {
#     est <- res$estimates_adjust
#   } else {
#     est <- res$estimates
#   }
#   # This is a rough fix right now. Should be updated to proper version using theta
#   est_pi <- sort(est$pi_med)
#   #true_pi <- sort(res_true$estimates$pi_med)
#   true_pi <- sort(prop.table(table(sim_pop$c_all)))
#   ci_est <- apply(est$pi_red, 2, function(x) quantile(x, c(0.025, 0.975)))
#   ci_est <- t(apply(ci_est, 1, sort))
#   
#   abs_bias <- mean(abs(est_pi - true_pi))
#   ci_width <- mean(ci_est[2, ] - ci_est[1, ])
#   cover <- mean(sapply(1:ncol(ci_est), function(x) 
#     ci_est[1, x] <= true_pi[x] & true_pi[x] <= ci_est[2, x]))
#   
#   return(list(abs_bias = abs_bias, ci_width = ci_width, cover = cover))
# }
