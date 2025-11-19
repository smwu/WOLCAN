#=================================================
# Functions for simulating data
# Author: Stephanie Wu
# Date created: 2024/05/29
# Date updated: 2025/11/18
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


# Create and save simulated population for additional supplementary scenarios
# for more variables and more varied NPS weights
# Inherit parameters from baysc::simulate_pop
# N: population size
# J, K, R: latent class dimensions
# high_overlap: TRUE = more overlap between B and R selection; FALSE = low overlap
# n_B, n_R: target expected sample sizes for B (non-probability) and R (reference)
# scale_B: controls dispersion of non-probability selection (and thus weight variability)
# formula_c, beta_mat_c: latent class model
# formula_x, beta_list_x: manifest items model
# formula_y, xi_vec_y: optional outcome model (logistic)
# save_res, save_path, pop_seed: saving + reproducibility
sim_pop_wolcan_pr <- function(N, J, K, R, high_overlap = TRUE, n_B, n_R, 
                              scale_B = 0.8,
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
  
  #==================================
  ### Generate population covs: age, sex, educ, income, ethnicity, urban, disability 
  pop <- data.frame(id = 1:N)
  
  # Age: continuous, mean ~53 (truncate to 30–75)
  pop$age <- rnorm(N, mean = 53, sd = 8)
  pop$age <- pmin(pmax(pop$age, 30), 75)
  
  # Sex (0 = Male, 1 = Female); 52% female
  pop$sex <- rbinom(N, size = 1, prob = 0.52)
  pop$sex_f <- factor(pop$sex, levels = c(0, 1), labels = c("Male", "Female"))
  
  # Education (4 levels, fixed marginal proportions)
  # <HS 9%, HS 35%, College 47%, Graduate 9%
  edu_levels <- c("<HS", "HS", "College", "Graduate")
  edu_probs  <- c(0.09, 0.35, 0.47, 0.09)
  pop$educ <- factor(
    sample(edu_levels, N, replace = TRUE, prob = edu_probs),
    levels = edu_levels,
    ordered = TRUE
  )
  
  # Household income (3 levels, fixed marginal proportions)
  # Levels: <10k 19%, 10–20k 18%, >20k 63% 
  inc_levels <- c("<10k", "10-20k", ">20k")
  pop$hhinc <- NA_character_
  # Income probabilities by educationlevel. Higher education -> higher income
  inc_prob_by_educ <- list(
    "<HS"     = c(0.40, 0.30, 0.30),  # more low income
    "HS"      = c(0.25, 0.25, 0.50),
    "College" = c(0.10, 0.15, 0.75),
    "Graduate"= c(0.05, 0.10, 0.85)   # most >20k
  )
  # Sample income by probabilities for each education level
  for (e in levels(pop$educ)) {
    idx   <- which(pop$educ == e)
    probs <- inc_prob_by_educ[[e]]
    pop$hhinc[idx] <- sample(inc_levels, length(idx),
                             replace = TRUE, prob = probs)
  }
  pop$hhinc <- factor(pop$hhinc, levels = inc_levels, ordered = TRUE)
  
  # Ethnicity: 95% Puerto Rican, 5% not Puerto Rican
  eth_levels <- c("PR", "NotPR")
  eth_probs  <- c(0.95, 0.05)
  pop$ethnicity <- factor(
    sample(eth_levels, N, replace = TRUE, prob = eth_probs),
    levels = eth_levels
  )
  
  ## Unobserved drivers: urbanicity & disability
  # Urbanicity: 62% urban
  pop$urban <- rbinom(N, size = 1, prob = 0.62)
  pop$urban_f <- factor(pop$urban,
                        levels = c(0, 1),
                        labels = c("Non-urban", "Urban"))
  
  # Disability: more common among older & lower income
  dis_lp <- -2.5 +
    0.04 * (pop$age - 53) +
    -0.3 * as.numeric(pop$hhinc == ">20k") +
    rnorm(N, 0, 0.3)
  pop$disability <- rbinom(N, size = 1, prob = plogis(dis_lp))
  
  
  #====================================
  ### Generate selection probabilities
  # Convert variables to numeric 
  pop$age_cent <- pop$age - mean(pop$age)  # centered
  pop$sex_f_num   <- as.numeric(pop$sex_f == "Female")  # 1 = F
  pop$educ_num <- as.numeric(pop$educ) - 1          # 0, 1, 2, 3
  pop$hhinc_num  <- as.numeric(pop$hhinc) - 1         # 0, 1, 2
  pop$ethn_num   <- as.numeric(pop$ethnicity == "PR")  # 1 = PR
  
  ## True logit for NPS (S_B). Similar to PROSPECT
  ## - Strongly favors older, female, urban
  ## - Nonlinear in age, education, income
  ## - Depends on unobserved urban & disability (for misspecification)
  lin_B_raw <- 
    # Age: increasing quadratic function of centered age
    0.08 * pop$age_cent + 0.002 * (pop$age_cent^2) +
    # Sex: females more likely
    0.40 * pop$sex_f_num +
    # Education: increases with education, then slightly flattens at graduate
    (0.30 * pop$educ_num - 0.005 * pop$educ_num^2) +
    # Income: lowest income less likely, highest more likely
    (-0.40 * (pop$hhinc_num == 0)) + 
    (0.20 * (pop$hhinc_num == 2)) +
    # Ethnicity: PR more likely
    0.30 * pop$ethn_num +
    # Unobserved urban & disability: both strongly increase selection
    0.12 * pop$urban +
    0.15 * pop$disability +
    # Interactions: urban women and older disabled especially likely
    0.12 * pop$sex_f_num * pop$urban +
    0.05 * pop$age_cent * pop$disability
  
  # Center and scale to control dispersion and weight variability
  lin_B <- scale_B * (lin_B_raw - mean(lin_B_raw))
  # Offset to hit target sample sizes
  offset_B <- root_n(x = lin_B, n = n_B)
  # Selection probabilities for non-probability sample
  pop$pi_B <- plogis(offset_B + lin_B)
  # # Check
  plot(density(1/pop$pi_B))
  summary(1/pop$pi_B)
  
  ## Reference sample (S_R). Similar to PRCS
  # High overlap: similar to NPS, but milder, closer to "probability" sampling
  if (high_overlap) {
    lin_R <- 
      0.05 * pop$age_cent + 0.001 * (pop$age_cent^2) +
      0.30 * pop$sex_f_num +
      0.10 * pop$educ_num +
      (-0.05 * (pop$hhinc_num == 0)) + (0.05 * (pop$hhinc_num == 2)) +
      0.05 * pop$ethn_num +
      0.10 * pop$urban +
      0.10 * pop$disability +
      0.08 * pop$sex_f_num * pop$urban + 
      0.05 * pop$age_cent * pop$disability
    
    # Offsets to hit target sample sizes
    offset_R <- root_n(x = lin_R, n = n_R)
    # Selection probabilities for reference probability sample
    pop$pi_R <- plogis(offset_R + lin_R)
    
  } else {
    ## Low overlap: different from NPS. Reference sample is younger, more rural, etc.
    lin_R <- 
      0.25 * (1 - pop$sex_f_num) +  # male more likely
      -0.05 * pop$age_cent - 0.001 * (pop$age_cent^2) +  # younger more likely
      -0.15 * pop$educ_num +           # lower education more likely
      (0.05 * (pop$hhinc_num == 0)) -  # lower income more likely
      (0.10 * (pop$hhinc_num == 2)) +  
      -0.10 * pop$urban +              # non-urban more likely
      -0.05 * pop$ethn_num +           # non-PR more likely
      -0.05 * pop$urban +              # rural more likely
      -0.05 * pop$disability +         # disability less likely
      0.05 * (1 - pop$sex_f_num) * (1 - pop$urban)  # rural male more likely

    # Offsets to hit target sample sizes
    offset_R <- root_n(x = lin_R, n = n_R)
    # Selection probabilities for reference probability sample
    pop$pi_R <- plogis(offset_R + lin_R)
  }
  
  
  #============================================
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
  
  # # Check
  # # item 1 has association with hhinc_num (higher level 2)
  # t(sapply(1:3, function(k) prop.table(table(sim_pop$X_data[sim_pop$c_all == k, 1]))))
  # # item 30 has association with age (higher level 2)
  # t(sapply(1:3, function(k) prop.table(table(sim_pop$X_data[sim_pop$c_all == k, 30]))))
  # # item 5 has no associations
  # t(sapply(1:3, function(k) prop.table(table(sim_pop$X_data[sim_pop$c_all == k, 5]))))
  
  #=========================================
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
  
  #=========================================
  ### Add to population data and save
  sim_pop <- list(N = N, J = J, R = R, K = K, 
                  # Keep all the covariates so can choose what to use in models
                  pop = pop[, c("age", "sex", "sex_f_num",
                                "educ", "educ_num", "hhinc", "hhinc_num",
                                "ethnicity", "ethn_num", "urban", "disability")],                  
                  pi_R = pop$pi_R, 
                  pi_B = pop$pi_B, 
                  c_all = pop$c_all, 
                  X_data = X_data, 
                  Y_data = Y_data,
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


### Create and save simulated sample data for additional supplementary scenarios
### for more variables and more varied NPS weights
# num_samps: number of samples to select
# n_B: Sample size for non-probability sample
# n_R: Sample size for reference sample
# scenario: Data generating scenario
sim_samp_wolcan_pr <- function(i, sim_pop, scenario, 
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
    covs = sim_pop$pop[ind_B, c("age", "sex", "sex_f_num",
                                "educ", "educ_num", "hhinc", "hhinc_num",
                                "ethnicity", "ethn_num", "urban", "disability")], 
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
    covs = sim_pop$pop[ind_R, c("age", "sex", "sex_f_num",
                                "educ", "educ_num", "hhinc", "hhinc_num",
                                "ethnicity", "ethn_num", "urban", "disability")], 
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
