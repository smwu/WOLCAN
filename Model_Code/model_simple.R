#======================================================
# Create pseudo-weights and run WOLCA model
# Simplified logistic regression and known pi_R setting
# Author: Stephanie Wu
# Date created: 2024/04/16
# Date updated: 2024/04/16
#======================================================

library(BART)
library(baysc)
library(abind)
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
data_dir <- "Data/"                 # Data directory
res_dir <- "Results/"               # Results directory
code_dir <- "Model_Code/"           # Model code directory
set.seed(1)

load(paste0(wd, data_dir, "sim_pop_simple_wolca.RData"))

# Scenario settings
num_post <- 300
pred_covs_R <- c("A1", "A2")
pred_covs_B <- c("A1", "A2")
frame_R <- 1
frame_B <- 1

i <- 1
for (i in 1:10) {
  print(i)
  run_samp(i, N = nrow(sim_pop$X_data), num_post = num_post, 
           pred_covs_R = pred_covs_R, pred_covs_B = pred_covs_B, 
           frame_R = frame_R, frame_B = frame_B, pi_R_known = TRUE)
}

for (i in 1:10) {
  print(i)
  run_samp2(i, N = nrow(sim_pop$X_data), num_post = num_post, 
           pred_covs_R = pred_covs_R, pred_covs_B = pred_covs_B, 
           frame_R = frame_R, frame_B = frame_B, pi_R_known = TRUE,
           n_runs = 5000, burn = 2500, thin = 5, update = 2500, D = 10)
}


run_samp <- function(i, N, num_post, pred_covs_R, pred_covs_B, frame_R, frame_B,
                     pi_R_known = FALSE) {
  load(paste0(wd, data_dir, "sim_samp_simple", i, "_B_wolca.RData"))
  load(paste0(wd, data_dir, "sim_samp_simple", i, "_R_wolca.RData"))
  
  dat_B <- data.frame(sim_samp_B$covs)
  dat_R <- data.frame(sim_samp_R$covs)

  # Number of overlapping individuals
  length(intersect(sim_samp_R$ind_R, sim_samp_B$ind_B))
  
  ### Get pseudo-weights for PS for both those in NPS and those in PS
  pi_R <- sim_samp_R$true_pi_R
  # Assume pi_R KNOWN for everyone
  if (pi_R_known) {
    hat_pi_R <- sim_samp_B$true_pi_R
  # Need to estimate pi_R for those in S_B
  } else {
    hat_pi_R <- NULL
  }
  weights_temp <- get_weights(dat_B = dat_B, dat_R = dat_R, ind_B = ind_B, 
                              num_post = num_post, pi_R = pi_R, hat_pi_R = hat_pi_R,
                              pred_covs_B = pred_covs_B, pred_covs_R = pred_covs_R, 
                              frame_B = frame_B, frame_R = frame_R, N = N)
  wts <- weights_temp$wts
  wts_post <- weights_temp$wts_post
  
  ### Run WOLCA model
  X_data <- sim_samp_B$X_data
  res <- wolca(x_mat = X_data, sampling_wt = wts, cluster_id = NULL, 
               stratum_id = NULL, run_sampler = "both", K_max = 30, adapt_seed = 1,
               class_cutoff = 0.05, n_runs = 10000, burn = 5000, thin = 5, 
               update = 5000, save_res = TRUE, 
               save_path = paste0(wd, res_dir, "test", i))
  res$data_vars$cluster_id <- 1:nrow(X_data)
  res_adj <- var_adjust(res = res, wts_post = wts_post, 
                        save_res = TRUE, 
                        save_path = paste0(wd, res_dir, "test_adj", i))
  # res_adj <- wolca_var_adjust(res = res, num_reps = 100, save_res = FALSE)
  
  res_unwt <- wolca(x_mat = X_data, sampling_wt = rep(1, nrow(X_data)), cluster_id = NULL, 
                    stratum_id = NULL, run_sampler = "both", K_max = 30, adapt_seed = 1,
                    class_cutoff = 0.05, n_runs = 10000, burn = 5000, thin = 5, 
                    update = 5000, save_res = TRUE, 
                    save_path = paste0(wd, res_dir, "test_unwt", i))
  
  # res_true <- wolca(x_mat = X_data, sampling_wt = (1/sim_samp_B$true_pi_B), 
  #                   cluster_id = 1:nrow(X_data), 
  #                   stratum_id = NULL, run_sampler = "both", K_max = 30, adapt_seed = 1,
  #                   class_cutoff = 0.05, n_runs = 20000, burn = 10000, thin = 5, 
  #                   update = 5000, save_res = TRUE, save_path = paste0(wd, res_dir, "test_true", i))
  # res_true <- wolca_var_adjust(res = res_true, save_res = TRUE, 
  #                              save_path = paste0(wd, res_dir, "test_true_adj", i))
  
  # summ_unadj <- get_summ_stats(res = res, res_true = res_true)
  # summ_adj <- get_summ_stats(res = res_adj, res_true = res_true)
  # summ_unwt <- get_summ_stats(res = res_unwt, res_true = res_true)
  # summ_df <- as.data.frame(rbind(summ_unadj, summ_adj, summ_unwt))
  
  summ_unadj <- get_summ_stats(res = res)
  summ_adj <- get_summ_stats(res = res_adj)
  summ_unwt <- get_summ_stats(res = res_unwt)
  summ_df <- as.data.frame(rbind(summ_unadj, summ_adj, summ_unwt))
  
  print(summ_df)
}


get_weights <- function(dat_B, dat_R, ind_B, ind_R, num_post, pi_R, hat_pi_R = NULL,
                        pred_covs_B, pred_covs_R, frame_B, frame_R, N) {
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
    # # Compare distributions
    # hist(sim_samp_B$true_pi_B, breaks = 30, col = rgb(0,0,1,1/4), ylim = c(0, 100), 
    #      xlab = "true_pi_B (blue) vs. hat_pi_B (red)")
    # hist(hat_pi_B, breaks = 30, col = rgb(1,0,0,1/4), add = TRUE)
    # hist(hat_pi_z, breaks = 20)
    # mean(abs(hat_pi_B - sim_samp_B$true_pi_B))
    
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
    
    # hist(hat_pi_B_dist, breaks = 30)
    # hist(apply(hat_pi_B_dist, 1, max), breaks = 30)
    # hist(apply(hat_pi_B_dist, 1, function(x) sum(x > 1)), breaks = 30)
    
    
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
    # # Compare distributions
    # hist(sim_samp_B$true_pi_B, breaks = 30, col = rgb(0,0,1,1/4), ylim = c(0, 100), 
    #      xlab = "true_pi_B (blue) vs. hat_pi_B (red)")
    # hist(hat_pi_B, breaks = 30, col = rgb(1,0,0,1/4), add = TRUE)
    
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
  }
  
  # Form pseudo-weights
  weights <- numeric(length = nrow(samp_comb))
  weights[z == 1] <- 1 / hat_pi_B
  weights[z == 0] <- 1 / pi_R
  wts <- weights[z == 1]
  
  ### Get distribution of weights for variance estimation
  wts_post <- t(1 / hat_pi_B_dist)
  
  ### Weight trimming
  wts <- trimw(wts, m = "t2", c = 10, max = 10)
  wts_post <- apply(wts_post, 2, function(x) trimw(x, m = "t2", c = 10, max = 10)) 
  
  # # Normalize weights to sum to N
  # pi_B_m <- pi_B_m * N / sum(1 / pi_B_m)
  # 
  # 
  # # Normalize to sum to sample size??????????
  # wts_post <- apply(wts_post, 2, function(x) x / sum(x / n1))
  
  return(list(wts = wts, wts_post = wts_post))
}

pi_R_known <- TRUE
run_samp2 <- function(i, N, num_post, pred_covs_R, pred_covs_B, frame_R, frame_B,
                     pi_R_known = FALSE, n_runs = 10000, burn = 5000, thin = 5, 
                     update = 5000, D = 10) {
  load(paste0(wd, data_dir, "sim_samp_simple", i, "_B_wolca.RData"))
  load(paste0(wd, data_dir, "sim_samp_simple", i, "_R_wolca.RData"))
  
  dat_B <- data.frame(sim_samp_B$covs)
  dat_R <- data.frame(sim_samp_R$covs)
  
  # Number of overlapping individuals
  length(intersect(sim_samp_R$ind_R, sim_samp_B$ind_B))
  
  ### Get pseudo-weights for PS for both those in NPS and those in PS
  pi_R <- sim_samp_R$true_pi_R
  # Assume pi_R KNOWN for everyone
  if (pi_R_known) {
    hat_pi_R <- sim_samp_B$true_pi_R
    # Need to estimate pi_R for those in S_B
  } else {
    hat_pi_R <- NULL
  }
  weights_temp <- get_weights(dat_B = dat_B, dat_R = dat_R, ind_B = ind_B, 
                              num_post = num_post, pi_R = pi_R, hat_pi_R = hat_pi_R,
                              pred_covs_B = pred_covs_B, pred_covs_R = pred_covs_R, 
                              frame_B = frame_B, frame_R = frame_R, N = N)
  wts <- weights_temp$wts
  wts_post <- weights_temp$wts_post
  
  ### Run WOLCA model
  X_data <- sim_samp_B$X_data
  res <- baysc::wolca(x_mat = X_data, sampling_wt = wts, cluster_id = NULL, 
               stratum_id = NULL, run_sampler = "both", K_max = 30, adapt_seed = 1,
               class_cutoff = 0.05, n_runs = n_runs, burn = burn, thin = thin, 
               update = update, save_res = FALSE)
  comb_estimates <- get_var_dir(D = D, wts_post = wts_post, 
                                K = res$estimates$K_red, 
                                x_mat = X_data, res = res, n_runs = n_runs, 
                                burn = burn, thin = thin, update = update)
  res$estimates_adjust <- list(pi_red = comb_estimates$pi_red_stack, 
                               theta_red = comb_estimates$theta_red_stack, 
                               pi_med = comb_estimates$pi_med, 
                               theta_med = comb_estimates$theta_med, 
                               c_all = comb_estimates$c_all,
                               pred_class_probs = comb_estimates$pred_class_probs)
  
  summ_df <- as.data.frame(get_summ_stats(res = res))
  
  print(summ_df)
}

# Obtain posterior parameter estimates and class assignments after propagating 
# uncertainty from the posterior weight distribution
# D: number of draws from the posterior weight distribution. WOLCA fixed sampler 
# will be run for each of these draws
# wts_post: Posterior distribution of estimated pseudo-weights
# K: Number of classes, determined by the adaptive sampler
# tol: Underflow tolerance
get_var_dir <- function(D = 5, wts_post, K, x_mat, res, cluster_id = NULL, 
                        stratum_id = NULL, fixed_seed = 1, alpha_fixed = NULL, 
                        class_cutoff = 0.05, n_runs = 10000, burn = 5000, 
                        thin = 5, update = 1000, save_res = FALSE, 
                        tol = 1e-8) {
  M <- dim(res$estimates$pi_red)[1]
  n <- nrow(x_mat)
  J <- ncol(x_mat)
  
  K_red_draws <- numeric(D)
  pi_red_draws <- vector("list", length = D)
  theta_red_draws <- vector("list", length = D)
  # pi_red_draws <- array(NA, dim = c(D, dim(res$estimates$pi_red)))
  # theta_red_draws <- array(NA, dim = c(D, dim(res$estimates$theta_red)))
  c_all_draws <- matrix(NA, nrow = D, ncol = n)
  for (d in 1:D) {
    print(paste0("Draw ", d))
    # Draw from weights posterior
    draw <- sample(1:ncol(wts_post), size = 1)
    wts_d <- c(wts_post[, draw])
    
    # Run fixed sampler
    if (is.null(alpha_fixed)) {
      alpha_fixed = rep(1, K)  # not sparsity-inducing
    }
    res_d <- baysc::wolca(x_mat = x_mat, sampling_wt = wts_d, cluster_id = cluster_id, 
                   stratum_id = stratum_id, run_sampler = "fixed", 
                   K_fixed = K, fixed_seed = d, alpha_fixed = alpha_fixed, 
                   class_cutoff = class_cutoff, n_runs = n_runs, burn = burn, 
                   thin = thin, update = update, save_res = save_res)
    # Get estimates
    K_red_draws[d] <- res_d$estimates$K_red
    c_all_draws[d, ] <- res_d$estimates$c_all
    pi_red_draws[[d]] <- res_d$estimates$pi_red
    theta_red_draws[[d]] <- res_d$estimates$theta_red
    # pi_red_draws[d, , ] <- res_d$estimates$pi_red
    # theta_red_draws[d, , , , ] <- res_d$estimates$theta_red
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
  pi_red_new <- vector("list", length = D)
  theta_red_new <- vector("list", length = D)
  for (d in 1:D) {
    new_order <- numeric(K)
    for (k in 1:K) {
      # Individuals who were assigned to old class k
      old_k_inds <- which(c_all_draws[d, ] == k)
      # Most common new class assignment among those individuals
      new_order[k] <- baysc::get_mode(red_c_all[old_k_inds])
    }
    # Reorder classes for all pi and theta estimates
    pi <- pi_red_draws[[d]][, new_order, drop = FALSE]  # reorder
    pi <- pi / rowSums(pi, na.rm = TRUE)
    pi_red_new[[d]] <- pi
    theta <- theta_red_draws[[d]][, , new_order, , drop = FALSE]  # reorder
    theta <- ifelse(theta < tol, tol, theta)  # prevent underflow
    theta <- plyr::aaply(theta, c(1, 2, 3),  # normalize
                         function(x) x / sum(x, na.rm = TRUE), .drop = FALSE)
    theta_red_new[[d]] <- theta
  }
  
  # Stack together the draws
  pi_red_stack <- do.call("abind", c(pi_red_new, along = 1))
  theta_red_stack <- do.call("abind", c(theta_red_new, along = 1))
  # theta_red_stack <- apply(theta_red_new, c(2, 3, 4), c)

  # Posterior median estimates
  pi_med <- apply(pi_red_stack, 2, stats::median, na.rm = TRUE)
  pi_med <- pi_med / sum(pi_med)  # Re-normalize
  theta_med <- apply(theta_red_stack, c(2, 3, 4), stats::median, na.rm = TRUE)
  theta_med <- ifelse(theta_med < tol, tol, theta_med) # prevent underflow
  theta_med <- plyr::aaply(theta_med, c(1, 2), function(x) x / sum(x),
                           .drop = FALSE)  # Re-normalize
  
  #============== Update c using unique classes and posterior estimates ========
  c_all <- res$estimates$c_all                     # Placeholder class assignments
  pred_class_probs <- matrix(NA, nrow = n, ncol = K) # Posterior class membership probabilities
  log_cond_c <- matrix(NA, nrow = n, ncol = K)       # Individual log-likelihood for each class
  # Calculate posterior class membership, p(c_i=k|-), for each class k
  for (i in 1:n) {
    for (k in 1:K) {
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
  
  comb_estimates <- list(K_red_draws = K_red_draws, c_all_draws = c_all_draws, 
                         pi_red_stack = pi_red_stack,
                         theta_red_stack = theta_red_stack,
                         pi_med = pi_med, theta_med = theta_med, 
                         c_all = c_all, pred_class_probs = pred_class_probs)
  return(comb_estimates)
}



# wts_post: (n1)xM posterior distribution of weights for individuals in the NPS, 
# with each column corresponding to a posterior draw
var_adjust <- function(res, wts_post, alpha = NULL, eta = NULL,  
                       save_res = TRUE, save_path = NULL, 
                       adjust_seed = NULL) {
  
  # Begin runtime tracker
  start_time <- Sys.time()
  
  # Set seed
  if (!is.null(adjust_seed)) {
    set.seed(adjust_seed)
  }
  
  #================= Extract dimensions and catch errors =======================
  # Check object class and estimates
  if (!inherits(res, "wolca")) {
    stop("res must be an object of class `wolca`, resulting from a call to the 
         `wolca()` function that includes results from the fixed sampler")
  } else if (is.null(res$estimates)) {
    stop("res must include results from the fixed sampler in the `wolca()` function")
  }
  # Check variance adjustment has not already been performed
  if ("estimates_adjust" %in% names(res)) {
    stop("variance adjustment has already been performed, since res$estimates is not NULL")
  }
  # # Check num_reps
  # if ((num_reps %% 100 != 0) | num_reps < 1) {
  #   stop("num_reps must be a whole number greater than 0, recommended to be at least 50. 
  #   More replicates will lead to more accurate results but will take longer to run.")
  # }
  
  # Extract data elements into the global environment
  K <- res$estimates$K_red
  J <- res$data_vars$J
  R_j <- res$data_vars$R_j
  R <- res$data_vars$R
  n <- res$data_vars$n
  x_mat <- res$data_vars$x_mat
  w_all <- res$data_vars$w_all
  stratum_id <- res$data_vars$stratum_id
  cluster_id <- res$data_vars$cluster_id
  
  # Get final number of classes (usually fewer classes than K_fixed)
  K <- res$estimates$K_red
  
  # Check hyperparameter dimensions match K
  if (!is.null(alpha)) {
    if (length(alpha) != K) {
      stop("length of alpha must be the same as K")
    }
  }
  if (!is.null(eta)) {
    if ((nrow(eta) != J) | (ncol(eta) != R)) {
      stop("eta must be a matrix with J rows and R columns, where J is
             the number of exposure items and R is the maximum number of 
             exposure categories")
    }
    if (any(eta == 0)) {
      warning("eta has 0 values and may result in rank-difficiency issues
                during the Hessian calculation in the var_adjust() function")
    }
  }
  
  # Check saving parameters
  if (!is.null(save_res)) {
    if (!is.logical(save_res)) {
      stop("save_res must be a boolean specifying if results should be saved")
    }
    if (save_res) {
      if (is.null(save_path) | !is.character(save_path)) {
        stop("save_path must be a string specifying a path and file name, such as '~/Documents/run'")
      } else {
        last_slash_ind <- regexpr("\\/[^\\/]*$", save_path)
        if (!dir.exists(substr(save_path, start = 1, stop = last_slash_ind))) {
          stop("directory specified in save_path does not exist")
        }
        if (last_slash_ind == length(save_path)) {
          stop("please append the start of a file name to the end of save_path. 
            For example, '~/Documents/run' can produce a saved file named 
            'run_wolca_results.RData'")
        }
      }
    }
  }
  
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
  
  # Normalize posterior weight distribution to sum to the sample size
  w_all_post <- apply(wts_post, 2, function(wts_l) wts_l / (sum(wts_l) / n))
  # Number of posterior weight sets
  L <- ncol(wts_post)
  
    # # Create survey design
    # if (!is.null(stratum_id)) {  # Include stratifying variable
    #   # Survey data frame for specifying survey design
    #   svy_data <- data.frame(stratum_id = stratum_id, cluster_id = cluster_id,
    #                          x_mat = x_mat, w_all = w_all)
    #   # Specify survey design
    #   svydes <- survey::svydesign(ids = ~cluster_id, strata = ~factor(stratum_id), 
    #                               weights = ~w_all, data = svy_data)
    # } else { # No stratifying variable
    #   # Survey data frame for specifying survey design
    #   svy_data <- data.frame(cluster_id = cluster_id, x_mat = x_mat, w_all = w_all)
    #   # Specify survey design
    #   svydes <- survey::svydesign(ids = ~cluster_id, weights = ~w_all, 
    #                               data = svy_data)    
    # }
  
  ### CHANGED TO USE THE POSTERIOR WEIGHT DISTRIBUTION
  ### Compute gradient of log-posterior for each posterior weight set
  j_l <- matrix(NA, nrow = L, ncol = length(unc_par_hat))
  for (l in 1:L) {
    # Survey data frame for specifying survey design
    svy_data <- data.frame(x_mat = x_mat, w_all = w_all_post[, l])
    # Get gradient of log-posterior for each posterior weight set
    j_l[l, ] <- grad_par(pwts = w_all_post[, l], svydata = svy_data, 
                         stan_mod = mod_stan, stan_data = data_stan, 
                         par_stan = par_stan, u_pars = unc_par_hat)
  }
  
  # ### TEST: case where no variability in weights
  # j_l <- matrix(NA, nrow = L, ncol = length(unc_par_hat))
  # for (l in 1:L) {
  #   # Survey data frame for specifying survey design
  #   svy_data <- data.frame(x_mat = x_mat, w_all = w_all)
  #   # Get gradient of log-posterior for each posterior weight set
  #   j_l[l, ] <- grad_par(pwts = w_all, svydata = svy_data, 
  #                        stan_mod = mod_stan, stan_data = data_stan, 
  #                        par_stan = par_stan, u_pars = unc_par_hat)
  # }
  # Survey data frame for specifying survey design
  svy_data <- data.frame(x_mat = x_mat, w_all = w_all)
  # Get gradient of log-posterior: vector of length 272
  j_once <- grad_par(pwts = w_all, svydata = svy_data,
                     stan_mod = mod_stan, stan_data = data_stan,
                     par_stan = par_stan, u_pars = unc_par_hat)
  temp <- stats::var(j_once)

    # ### OLD CHANGED: CHECK THIS IS CORRECT for scale, rscales, combined.weights
    # # Create svrepdesign
    # ### I think this is wrong b/c replicate weights are not the same
    # svyrep <- survey::svrepdesign(design = svydes, repweights = wts_post, 
    #                               weights = res$data_vars$w_all, type = "other", 
    #                               scale = 1, rscales = 1)
    # # Get survey replicates
    # rep_temp <- survey::withReplicates(design = svyrep, theta = grad_par, 
    #                                    stan_mod = mod_stan, stan_data = data_stan, 
    #                                    par_stan = par_stan, u_pars = unc_par_hat)
    # # Create svrepdesign
    # svyrep <- survey::as.svrepdesign(design = svydes, type = "mrbbootstrap",
    #                                  replicates = num_reps)
    # # Get survey replicates
    # rep_temp <- survey::withReplicates(design = svyrep, theta = grad_par,
    #                                    stan_mod = mod_stan, stan_data = data_stan,
    #                                    par_stan = par_stan, u_pars = unc_par_hat)
    # J_hat <- stats::vcov(rep_temp)
  
  # Compute adjustment
  J_hat <- stats::var(j_l) + mean(j_l)^2
  temp <- stats::var(j_l)
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
  
  res$estimates_adjust <- estimates_adjust
  class(res) <- "wolca"
  
  # Save output
  if (save_res) {
    save(res, file = paste0(save_path, "_wolca_results.RData"))
  }
  
  return(res)
}



### Weight trimming function from Rafei thesis
trimw <- function(w, m='t1', c=10, max=Inf){
  if(m=='t1'){
    kn <- sqrt(c*mean(w^2, na.rm=T))
    i <- 0
    while(sum(w>kn & !is.na(w))>=1 & i<max){
      s <- sum(w, na.rm=T)
      w[w>kn & !is.na(w)] <- kn
      w[!is.na(w)] <- w[!is.na(w)]*s/sum(w, na.rm=T)
      kn <- sqrt(c*sum(w^2, na.rm=T)/sum(!is.na(w), na.rm=T))
      i <- i+1
    }
  } else if(m=='t2'){
    kn <- median(w, na.rm=T)+c*IQR(w, na.rm=T)
    i <- 0
    while(sum(w>kn, na.rm=T)>=1 & i<max){
      s <- sum(w, na.rm=T)
      w[w>kn & !is.na(w)] <- kn
      w[!is.na(w)] <- w[!is.na(w)]*s/sum(w, na.rm=T)
      kn <- median(w, na.rm=T)+c*IQR(w, na.rm=T)
      i <- i+1
    }
  } else if(m=='t3'){
    a <- 2+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T))
    b <- (sum(w, na.rm=T)-1)*(1+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T)))
    kn <- 1/(length(w[!is.na(w)])*qbeta(0.99, a, b))
    i <- 0
    while(sum(w>kn, na.rm=T)>=1 & i<max){
      s <- sum(w, na.rm=T)
      w[w>kn & !is.na(w)] <- kn
      w[!is.na(w)] <- w[!is.na(w)]*s/sum(w, na.rm=T)
      a <- 2+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T))
      b <- (sum(w, na.rm=T)-1)*(1+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T)))
      kn <- 1/(length(w[!is.na(w)])*qbeta(0.99, a, b))
      i <- i+1
    }
  }
  w	
}


get_summ_stats <- function(res, res_true = NULL) {
  if (!is.null(res$estimates_adjust)) {
    est <- res$estimates_adjust
  } else {
    est <- res$estimates
  }
  # This is a rough fix right now. Should be updated to proper version using theta
  est_pi <- sort(est$pi_med)
  #true_pi <- sort(res_true$estimates$pi_med)
  true_pi <- sort(prop.table(table(sim_pop$c_all)))
  ci_est <- apply(est$pi_red, 2, function(x) quantile(x, c(0.025, 0.975)))
  ci_est <- t(apply(ci_est, 1, sort))
  
  abs_bias <- mean(abs(est_pi - true_pi))
  ci_width <- mean(ci_est[2, ] - ci_est[1, ])
  cover <- mean(sapply(1:ncol(ci_est), function(x) 
    ci_est[1, x] <= true_pi[x] & true_pi[x] <= ci_est[2, x]))
  
  return(list(abs_bias = abs_bias, ci_width = ci_width, cover = cover))
}


#====================== OLD CODE ===============================================
# library(BART)
# library(baysc)
# wd <- "~/Documents/GitHub/npswolca/"  # working directory
# set.seed(1)
# 
# i <- 1
# load(paste0(wd, "sim_samp", i, "_B_wolca.RData"))
# load(paste0(wd, "sim_samp", i, "_R_wolca.RData"))
# 
# dat_B <- data.frame(sim_samp_B$covs)
# dat_R <- data.frame(sim_samp_R$covs)
# num_post <- 100
# pred_covs_R <- c("A1", "A2", "A1A2", "logA2")
# pred_covs_B <- c("A1", "A2", "A1A2", "logA2")
# frame_R <- 1
# frame_B <- 1
# 
# # Number of overlapping individuals
# length(intersect(sim_samp_R$ind_R, sim_samp_B$ind_B))
# 
# 
# # Create stacked sample
# n0 <- nrow(dat_R)
# n1 <- nrow(dat_B)
# samp_comb <- rbind(dat_B, dat_R)
# # Get known sampling probabilities for i in S_R
# pi_R <- sim_samp_R$pi_R
# logit_pi_R <- log(pi_R / (1 - pi_R))
# # Indicator for NPS given inclusion in stacked sample
# z <- rep(1:0, c(n1, n0))
# # Predict logit(pi_R) for i in S_B
# # Change so that it only predicts for those in S_B rather than all
# fit_pi_R <- wbart(x.train = samp_comb[z == 0, pred_covs_R, drop = FALSE], 
#                   y.train = logit_pi_R, 
#                   x.test = samp_comb[z == 1, pred_covs_R, drop = FALSE],
#                   ntree = 50L, nskip = 100L, ndpost = num_post) 
# hat_logit_pi_R <- fit_pi_R$yhat.test.mean
# # Convert to pi_R scale using expit
# hat_pi_R <- exp(hat_logit_pi_R) / (1 + exp(hat_logit_pi_R))
# # Check predicted pi_R against true pi_R 
# hat_pi_R[1:5]
# sim_samp_B$true_pi_R[1:5]
# # Get distribution of hat_pi_R
# hat_logit_pi_R_dist <- fit_pi_R$yhat.test
# hat_pi_R_dist <- exp(hat_logit_pi_R_dist) / (1 + exp(hat_logit_pi_R_dist))
# 
# # Predict pi_B for i in S_B
# fit_pi_B <- pbart(x.train = samp_comb[, pred_covs_B, drop = FALSE],
#                   y.train = z, ntree = 50L, nskip = 100L, ndpost = num_post)
# hat_pi_z <- fit_pi_B$prob.train.mean
# hat_pi_B <- hat_pi_z[z == 1] * hat_pi_R * frame_R / 
#   (frame_B * (1 - hat_pi_z[z == 1]))
# # Check predicted pi_B against true pi_B
# hat_pi_B[1:5]
# sim_samp_B$true_pi_B[1:5]
# # Get distribution of hat_pi_B, incorporating distribution of hat_pi_R
# hat_pi_z_dist <- fit_pi_B$prob.train
# hat_pi_B_dist <- matrix(NA, nrow = num_post, ncol = n1)
# # hat_pi_B_dist <- vector(mode = "list", length = num_post)
# for (m in 1:num_post) {
#   pi_z_m <- hat_pi_z_dist[m, ]
#   pi_B_m <- pi_z_m[z == 1] * hat_pi_R_dist[m, ] * frame_R / 
#     (frame_B * (1 - pi_z_m[z == 1]))
#   # Restrict to [0, 1]
#   hat_pi_B_dist[m, ] <- sapply(1:n1, function(x) max(min(pi_B_m[x], 1), 0))
# }
# # Example plot for individual 1
# hist(hat_pi_B_dist[, 1], breaks = 30)
# abline(v = hat_pi_B[1], col = "red", lwd = 2)
# abline(v = sim_samp_B$true_pi_B[1], col = "forestgreen", lwd = 2)
# 
# # Form pseudo-weights
# weights <- numeric(length = nrow(samp_comb))
# weights[z == 1] <- 1 / hat_pi_B
# weights[z == 0] <- 1 / pi_R
# wts <- weights[z == 1]
# 
# # Normalize weights to sum to N
# pi_B_m <- pi_B_m * N / sum(1 / pi_B_m)
# 
# ### Get distribution of weights for variance estimation
# wts_post <- t(1 / hat_pi_B_dist)
# # Normalize to sum to sample size
# wts_post <- apply(wts_post, 2, function(x) x / sum(x / n1))
# range(wts_post)
# range(res$data_vars$w_all)
# 
# ### Weight trimming
# wts <- trimw(wts, m = "t1", c = 10, max = 10)
# wts_post <- apply(wts_post, 2, function(x) trimw(x, m = "t1", c = 10, max = 10)) 
# 
# ### Run WOLCA model
# X_data <- sim_samp_B$X_data
# res <- wolca(x_mat = X_data, sampling_wt = wts, cluster_id = NULL, 
#              stratum_id = NULL, run_sampler = "both", K_max = 30, adapt_seed = 1,
#              class_cutoff = 0.05, n_runs = 20000, burn = 10000, thin = 5, 
#              update = 100, save_res = TRUE, save_path = paste0(wd, "test"))
# res$data_vars$cluster_id <- 1:425
# res_adj <- wolca_var_adjust(res = res, num_reps = 100, save_res = FALSE)
# 
# res_unwt <- wolca(x_mat = X_data, sampling_wt = rep(1, nrow(X_data)), cluster_id = NULL, 
#                   stratum_id = NULL, run_sampler = "both", K_max = 30, adapt_seed = 1,
#                   class_cutoff = 0.05, n_runs = 20000, burn = 10000, thin = 5, 
#                   update = 100, save_res = TRUE, save_path = paste0(wd, "test_unwt"))
# 
# res_true <- wolca(x_mat = X_data, sampling_wt = (1/sim_samp_B$true_pi_B), 
#                   cluster_id = NULL, 
#                   stratum_id = NULL, run_sampler = "both", K_max = 30, adapt_seed = 1,
#                   class_cutoff = 0.05, n_runs = 20000, burn = 10000, thin = 5, 
#                   update = 100, save_res = TRUE, save_path = paste0(wd, "test_true"))
