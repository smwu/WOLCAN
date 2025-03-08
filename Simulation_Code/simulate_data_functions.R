#=================================================
# Functions for simulating data
# Author: Stephanie Wu
# Date created: 2024/05/29
# Date updated: 2024/05/29
#=================================================


### Helper function to obtain offset to match desired sample size
root_n <- function(x, n) {
  uniroot(function(t){mean(exp(t + x) / (1 + exp(t + x))) - (n / N)}, 
          interval = c(-50, 50))$root
}

# Create and save simulated population
# Inherit parameters from baysc::simulate_pop
# n_B: Sample size for non-probability sample, used for selection probabilities offset
# n_R: Sample size for reference sample, used for selection probabilities offset
# rho: Correlation between covariates influencing selection
# high_overlap: Boolean indicating if the reference and non-probability samples 
# selection mechanism should have high overlap (default = TRUE) or low overlap.
sim_pop_wolcan <- function(N, J, K, R, rho, high_overlap = TRUE, n_B, n_R, 
                           formula_c, beta_mat_c, formula_x, beta_list_x, 
                           formula_y = NULL, xi_vec_y = NULL,
                           save_res = TRUE, save_path = NULL, pop_seed = 1) {
  
  # Set seed
  set.seed(pop_seed)
  
  ### Catch errors
  # Check save_res and save_path
  if (!is.logical(save_res)) {
    stop("save_res must be a boolean specifying if results should be saved")
  } else if (save_res) {
    if (is.null(save_path) | !is.character(save_path)) {
      stop("save_path must be a string specifying a path and file name, such as 
           '~/Documents/stratified_classes'")
    } else {
      last_slash_ind <- regexpr("\\/[^\\/]*$", save_path)
      if (!dir.exists(substr(save_path, start = 1, stop = last_slash_ind))) {
        stop("directory specified in save_path does not exist")
      }
      if (last_slash_ind == length(save_path)) {
        stop("please append the start of a file name to the end of save_path. 
             For example, '~/Documents/stratified_classes' will produce a saved 
             file named 'stratified_classes_sim_pop.RData'")
      }
    }
  }
  
  ### Generate covariates A1 and A2 in the population
  cov <- matrix(c(1, rho, rho, 1), 2, 2)
  mu <- c(0, 0)
  pop <- as.data.frame(mvrnorm(n = N, mu = mu, Sigma = cov))
  names(pop) <- c("A1", "A2")
  pop$A12 <- pop$A1^2
  pop$A1A2 <- pop$A1*pop$A2
  pop$logA2 <- log(abs(pop$A2))
  pop$sinA1A2 <- sin(pop$A1A2)
  pop$A3 <- rnorm(n = N, 2, 1)
  
  ### Generate selection probabilities
  ## High-overlap setting
  if (high_overlap) {
    # Offsets to obtain correct sample size
    offset_B <- root_n(x = (- 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
                            + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3),
                       n = n_B)
    offset_R <- root_n(x = (- 0.6 * pop$A1 + 0.4 * pop$A12 + 0.7 * pop$A2
                            + 0.1 * pop$logA2 - 0.05 * pop$sinA1A2 + 0.4 * pop$A3), n = n_R)
    # Selection probabilities for non-probability sample
    pop$pi_B <- exp(offset_B - 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
                    + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3) /
      (1 + exp(offset_B - 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
               + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3))
    # Selection probabilities for reference sample
    pop$pi_R <- exp(offset_R - 0.6 * pop$A1 + 0.4 * pop$A12 + 0.7 * pop$A2
                    + 0.1 * pop$logA2 - 0.05 * pop$sinA1A2 + 0.4 * pop$A3) /
      (1 + exp(offset_R - 0.6 * pop$A1 + 0.4 * pop$A12 + 0.7 * pop$A2
               + 0.1 * pop$logA2 - 0.05 * pop$sinA1A2 + 0.4 * pop$A3))
    ## Low-overlap setting
  } else {
    # Offsets to obtain correct sample size
    offset_B <- root_n(x = (- 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
                            + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3),
                       n = n_B)
    offset_R <- root_n(x = (0.7 * pop$A1 - 0.6 * pop$A2 + 0.1 * pop$logA2 
                            + 0.1 * pop$A1A2 - 0.1 * pop$A3), n = n_R)
    # Selection probabilities for non-probability sample
    pop$pi_B <- exp(offset_B - 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
                    + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3) /
      (1 + exp(offset_B - 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
               + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3))
    # Selection probabilities for reference sample
    pop$pi_R <- exp(offset_R + 0.7 * pop$A1 - 0.6 * pop$A2 + 0.1 * pop$logA2 
                    + 0.1 * pop$A1A2 - 0.1 * pop$A3) /
      (1 + exp(offset_R + 0.7 * pop$A1 - 0.6 * pop$A2 + 0.1 * pop$logA2 
               + 0.1 * pop$A1A2 - 0.1 * pop$A3))
  }
  
  ### Generate categorical latent class assignment C
  regr_vars <- labels(stats::terms(stats::as.formula(formula_c)))
  regr_vars <- regr_vars[!(stringr::str_detect(regr_vars, ":"))] 
  # Create design matrix
  design_mat_c <- stats::model.matrix(stats::as.formula(formula_c), data = pop)
  # Create latent class assignments
  out_vars_c <- create_categ_var(beta_mat = beta_mat_c, 
                                 design_mat = design_mat_c, split_dim = NULL,
                                 V = pop)
  c_all <- out_vars_c$categ_var
  prop.table(table(c_all))
  # Add latent class to data frame of variables
  pop$c_all <- as.factor(c_all)
  
  ### Generate observed manifest variables X
  regr_vars <- labels(terms(as.formula(formula_x)))
  regr_vars <- regr_vars[!(stringr::str_detect(regr_vars, ":"))] 
  # Get design matrix
  design_mat_x <- stats::model.matrix(stats::as.formula(formula_x), data = pop)
  # Create multivariate categorical exposure variable
  X_data <- matrix(NA, nrow = N, ncol = J)
  # Obtain underlying thetas for each item, class, and category
  true_global_thetas <- array(NA, dim = c(J, K, R))
  # Obtain underlying modal patterns for each item and class
  true_global_patterns <- matrix(NA, nrow = J, ncol = K)
  for (j in 1:J) {
    out_vars_x <- create_categ_var(beta_mat = beta_list_x[[j]], 
                                   design_mat = design_mat_x, 
                                   split_dim = "c_all", V = pop)
    X_data[, j] <- out_vars_x$categ_var
    true_global_thetas[j, , ] <- out_vars_x$pi_split
    true_global_patterns[j, ] <- apply(out_vars_x$pi_split, 1, which.max)
  }
  
  ### Generate outcome variable Y
  if (!is.null(formula_y)) {
    # Design matrix
    design_mat_y <- stats::model.matrix(stats::as.formula(formula_y), data = pop)
    # Get vector of individual linear predictors
    lin_pred <- design_mat_y %*% xi_vec_y
    # Get vector of individual underlying outcome probabilities
    true_Phi_under <- exp(lin_pred) / (1 + exp(lin_pred))
    # Outcome data
    Y_data <- stats::rbinom(n = N, size = 1, prob = true_Phi_under)
    
    # hist(true_Phi_under)
    # summary(true_Phi_under)
  } else {
    Y_data <- NULL
  }
  
  
  ### Add to population data and save
  sim_pop <- list(N = N, J = J, R = R, K = K, 
                  pop = pop[, c("A1", "A2", "A12", "A1A2", "logA2", "sinA1A2", "A3")], 
                  # A1 = pop$A1, A2 = pop$A2, A12 = pop$A12, A1A2 = pop$A1A2, 
                  # logA2 = pop$logA2, sinA1A2 = pop$sinA1A2, 
                  pi_R = pop$pi_R, pi_B = pop$pi_B, 
                  c_all = pop$c_all, X_data = X_data, Y_data = Y_data,
                  true_xi = xi_vec_y,
                  true_global_thetas = true_global_thetas, 
                  true_global_patterns = true_global_patterns)
  
  if (save_res) {
    save(sim_pop, file = paste0(save_path, "sim_pop_wolcan.RData"))
  }
  
  return(sim_pop)
}


# Create and save simulated population with informative sampling dependent on Y
# Inherit parameters from baysc::simulate_pop
# n_B: Sample size for non-probability sample, used for selection probabilities offset
# n_R: Sample size for reference sample, used for selection probabilities offset
# rho: Correlation between covariates influencing selection
# high_overlap: Boolean indicating if the reference and non-probability samples 
# selection mechanism should have high overlap (default = TRUE) or low overlap.
sim_pop_wolcan_inf <- function(N, J, K, R, rho, high_overlap = TRUE, n_B, n_R, 
                           formula_c, beta_mat_c, formula_x, beta_list_x, 
                           formula_y = NULL, xi_vec_y = NULL,
                           save_res = TRUE, save_path = NULL, pop_seed = 1) {
  
  # Set seed
  set.seed(pop_seed)
  
  ### Catch errors
  # Check save_res and save_path
  if (!is.logical(save_res)) {
    stop("save_res must be a boolean specifying if results should be saved")
  } else if (save_res) {
    if (is.null(save_path) | !is.character(save_path)) {
      stop("save_path must be a string specifying a path and file name, such as 
           '~/Documents/stratified_classes'")
    } else {
      last_slash_ind <- regexpr("\\/[^\\/]*$", save_path)
      if (!dir.exists(substr(save_path, start = 1, stop = last_slash_ind))) {
        stop("directory specified in save_path does not exist")
      }
      if (last_slash_ind == length(save_path)) {
        stop("please append the start of a file name to the end of save_path. 
             For example, '~/Documents/stratified_classes' will produce a saved 
             file named 'stratified_classes_sim_pop.RData'")
      }
    }
  }
  
  ### Generate covariates A1 and A2 in the population
  cov <- matrix(c(1, rho, rho, 1), 2, 2)
  mu <- c(0, 0)
  pop <- as.data.frame(mvrnorm(n = N, mu = mu, Sigma = cov))
  names(pop) <- c("A1", "A2")
  pop$A12 <- pop$A1^2
  pop$A1A2 <- pop$A1*pop$A2
  pop$logA2 <- log(abs(pop$A2))
  pop$sinA1A2 <- sin(pop$A1A2)
  pop$A3 <- rnorm(n = N, 2, 1)
  
  ### Generate categorical latent class assignment C
  regr_vars <- labels(stats::terms(stats::as.formula(formula_c)))
  regr_vars <- regr_vars[!(stringr::str_detect(regr_vars, ":"))] 
  # Create design matrix
  design_mat_c <- stats::model.matrix(stats::as.formula(formula_c), data = pop)
  # Create latent class assignments
  out_vars_c <- create_categ_var(beta_mat = beta_mat_c, 
                                 design_mat = design_mat_c, split_dim = NULL,
                                 V = pop)
  c_all <- out_vars_c$categ_var
  prop.table(table(c_all))
  # Add latent class to data frame of variables
  pop$c_all <- as.factor(c_all)
  
  
  ### Generate observed manifest variables X
  regr_vars <- labels(terms(as.formula(formula_x)))
  regr_vars <- regr_vars[!(stringr::str_detect(regr_vars, ":"))] 
  # Get design matrix
  design_mat_x <- stats::model.matrix(stats::as.formula(formula_x), data = pop)
  # Create multivariate categorical exposure variable
  X_data <- matrix(NA, nrow = N, ncol = J)
  # Obtain underlying thetas for each item, class, and category
  true_global_thetas <- array(NA, dim = c(J, K, R))
  # Obtain underlying modal patterns for each item and class
  true_global_patterns <- matrix(NA, nrow = J, ncol = K)
  for (j in 1:J) {
    out_vars_x <- create_categ_var(beta_mat = beta_list_x[[j]], 
                                   design_mat = design_mat_x, 
                                   split_dim = "c_all", V = pop)
    X_data[, j] <- out_vars_x$categ_var
    true_global_thetas[j, , ] <- out_vars_x$pi_split
    true_global_patterns[j, ] <- apply(out_vars_x$pi_split, 1, which.max)
  }
  
  ### Generate outcome variable Y
  # Design matrix
  design_mat_y <- stats::model.matrix(stats::as.formula(formula_y), data = pop)
  # Get vector of individual linear predictors
  lin_pred <- design_mat_y %*% xi_vec_y
  # Get vector of individual underlying outcome probabilities
  true_Phi_under <- exp(lin_pred) / (1 + exp(lin_pred))
  # Outcome data
  Y_data <- stats::rbinom(n = N, size = 1, prob = true_Phi_under)
  pop$y_all <- Y_data
  
  # hist(true_Phi_under)
  # summary(true_Phi_under)
  
  ### Generate selection probabilities with informative sampling
  ## High-overlap setting
  if (high_overlap) {
    # Offsets to obtain correct sample size
    offset_B <- root_n(x = (- 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
                            + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3
                            + 0.4 * pop$y_all),
                       n = n_B)
    offset_R <- root_n(x = (- 0.6 * pop$A1 + 0.4 * pop$A12 + 0.7 * pop$A2
                            + 0.1 * pop$logA2 - 0.05 * pop$sinA1A2 + 0.4 * pop$A3
                            + 0.3 * pop$y_all), n = n_R)
    # Selection probabilities for non-probability sample
    pop$pi_B <- exp(offset_B - 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
                    + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3
                    + 0.4 * pop$y_all) /
      (1 + exp(offset_B - 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
               + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3
               + 0.4 * pop$y_all))
    # Selection probabilities for reference sample
    pop$pi_R <- exp(offset_R - 0.6 * pop$A1 + 0.4 * pop$A12 + 0.7 * pop$A2
                    + 0.1 * pop$logA2 - 0.05 * pop$sinA1A2 + 0.4 * pop$A3
                    + 0.3 * pop$y_all) /
      (1 + exp(offset_R - 0.6 * pop$A1 + 0.4 * pop$A12 + 0.7 * pop$A2
               + 0.1 * pop$logA2 - 0.05 * pop$sinA1A2 + 0.4 * pop$A3
               + 0.3 * pop$y_all))
    ## Low-overlap setting
  } else {
    # Offsets to obtain correct sample size
    offset_B <- root_n(x = (- 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
                            + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3
                            + 0.4 * pop$y_all),
                       n = n_B)
    offset_R <- root_n(x = (0.7 * pop$A1 - 0.6 * pop$A2 + 0.1 * pop$logA2 
                            + 0.1 * pop$A1A2 - 0.1 * pop$A3
                            - 0.3 * pop$y_all), n = n_R)
    # Selection probabilities for non-probability sample
    pop$pi_B <- exp(offset_B - 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
                    + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3
                    + 0.4 * pop$y_all) /
      (1 + exp(offset_B - 0.9 * pop$A1 + 0.2 * pop$A12 + 0.8 * pop$A2
               + 0.2 * pop$logA2 - 0.1 * pop$sinA1A2 + 0.3 * pop$A3
               + 0.4 * pop$y_all))
    # Selection probabilities for reference sample
    pop$pi_R <- exp(offset_R + 0.7 * pop$A1 - 0.6 * pop$A2 + 0.1 * pop$logA2 
                    + 0.1 * pop$A1A2 - 0.1 * pop$A3
                    - 0.3 * pop$y_all) /
      (1 + exp(offset_R + 0.7 * pop$A1 - 0.6 * pop$A2 + 0.1 * pop$logA2 
               + 0.1 * pop$A1A2 - 0.1 * pop$A3
               - 0.3 * pop$y_all))
  }
  
  ### Add to population data and save
  sim_pop <- list(N = N, J = J, R = R, K = K, 
                  pop = pop[, c("A1", "A2", "A12", "A1A2", "logA2", "sinA1A2", "A3")], 
                  # A1 = pop$A1, A2 = pop$A2, A12 = pop$A12, A1A2 = pop$A1A2, 
                  # logA2 = pop$logA2, sinA1A2 = pop$sinA1A2, 
                  pi_R = pop$pi_R, pi_B = pop$pi_B, 
                  c_all = pop$c_all, X_data = X_data, Y_data = Y_data,
                  true_xi = xi_vec_y,
                  true_global_thetas = true_global_thetas, 
                  true_global_patterns = true_global_patterns)
  
  if (save_res) {
    save(sim_pop, file = paste0(save_path, "sim_pop_wolcan.RData"))
  }
  
  return(sim_pop)
}


### Create and save simulated sample data 
# num_samps: number of samples to select
# n_B: Sample size for non-probability sample
# n_R: Sample size for reference sample
# scenario: Data generating scenario
sim_samp_wolcan <- function(i, sim_pop, scenario, 
                            save_res = TRUE, save_path, samp_seed) {
  
  # Check save_res and save_path
  if (!is.logical(save_res)) {
    stop("save_res must be a boolean specifying if results should be saved")
  } else if (save_res) {
    if (is.null(save_path) | !is.character(save_path)) {
      stop("save_path must be a string specifying a path and file name, such as 
           '~/Documents/stratified_classes'")
    } else {
      last_slash_ind <- regexpr("\\/[^\\/]*$", save_path)
      if (!dir.exists(substr(save_path, start = 1, stop = last_slash_ind))) {
        stop("directory specified in save_path does not exist")
      }
      if (last_slash_ind == length(save_path)) {
        stop("please append the start of a file name to the end of save_path. 
             For example, '~/Documents/stratified_classes' will produce a saved 
             file named 'stratified_classes_sim_pop.RData'")
      }
    }
  }
  
  # Population size
  N <- sim_pop$N
  
  ### Generate samples using Poisson sampling with the selection probabilities
  set.seed(samp_seed)
  # Non-probability sample indicators
  delta_B <- rbinom(n = N, size = 1, prob = sim_pop$pi_B) 
  # Reference probability sample indicators
  delta_R <- rbinom(n = N, size = 1, prob = sim_pop$pi_R)
  table(delta_B)
  table(delta_R)
  # Indices of selection individuals
  ind_B <- which(delta_B == 1)
  ind_R <- which(delta_R == 1)
  
  # Subset to the non-probability sample
  sim_samp_B <- list(
    # indices, covariates, true (unknown) pi_B, true (unknown) pi_R, c_all, 
    # and X_data for selected individuals
    ind_B = ind_B, 
    covs = sim_pop$pop[ind_B, c("A1", "A2", "A12", "A1A2", "logA2", "sinA1A2", "A3")], 
    true_pi_B = sim_pop$pi_B[ind_B], 
    true_pi_R = sim_pop$pi_R[ind_B],
    c_all = sim_pop$c_all[ind_B], 
    X_data = sim_pop$X_data[ind_B, ], 
    Y_data = sim_pop$Y_data[ind_B],
    delta_B = delta_B)
  # Subset to the reference sample
  sim_samp_R <- list(
    # indices, covariates, known pi_R, c_all, and X_data for selected individuals
    ind_R = ind_R, 
    covs = sim_pop$pop[ind_R, c("A1", "A2", "A12", "A1A2", "logA2", "sinA1A2", "A3")], 
    pi_R = sim_pop$pi_R[ind_R], 
    c_all = sim_pop$c_all[ind_R], 
    X_data = sim_pop$X_data[ind_R, ], 
    Y_data = sim_pop$Y_data[ind_R],
    delta_R = delta_R)
  
  # Save and return results
  if (save_res) {
    save(sim_samp_B, file = paste0(save_path, "sim_samp_", i, "_B_wolcan.RData"))
    save(sim_samp_R, file = paste0(save_path, "sim_samp_", i, "_R_wolcan.RData"))
  }
  
  sim_samps <- list(sim_samp_B = sim_samp_B, sim_samp_R = sim_samp_R)
  return(sim_samps)
}


# Plot overlap
plot_overlap <- function(wd, data_dir, scenario, samp_i) {
  # Load convenience sample
  load(paste0(wd, data_dir,  "scen_", scenario, "/", "sim_samp_", samp_i, 
              "_B_wolcan.RData"))
  # Load reference sample
  load(paste0(wd, data_dir,  "scen_", scenario, "/", "sim_samp_", samp_i, 
              "_R_wolcan.RData"))
  # Load population
  load(paste0(wd, data_dir,  "scen_", scenario, "/", "sim_pop_wolcan.RData"))
  n_B <- nrow(sim_samp_B$covs)
  n_R <- nrow(sim_samp_R$covs)
  
  # Get inclusion probabilities
  comb_ind <- c(sim_samp_B$ind_B, sim_samp_R$ind_R)
  comb_pi_B <- sim_pop$pi_B[comb_ind]
  comb_pi_R <- sim_pop$pi_R[comb_ind]
  
  # Create df
  df <- data.frame(comb_ind, comb_pi_B, comb_pi_R)
  # Sort by increasion pi_B
  df <- df %>%
    arrange(comb_pi_B, comb_pi_R)
  # high_df$sample <- c(rep("Convenience", n_B),
  #                     rep("Reference", n_R))
  df$new_ind <- 1:length(comb_ind)
  df_long <- df %>%
    pivot_longer(
      cols = comb_pi_B:comb_pi_R,
      names_to = "Sample",
      values_to = "Incl_Prob"
    ) 
  df_long %>%
    ggplot(aes(x = new_ind, y = Incl_Prob, col = Sample)) + 
    geom_line() + theme_bw() + 
    facet_grid(~Sample) + 
    theme(strip.text = element_blank()) + 
    theme(legend.position = "right") +
    scale_color_manual(values = c("#66c2a5", "#fc8d62"), name = "", 
                       labels = expression(pi[B], pi[R])) + 
    xlab("Unit") + ylab("Inclusion Probability")
}
