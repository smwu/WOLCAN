#=================================================
# Generate population and sample data for NPSWOLCA
# Author: Stephanie Wu
# Date created: 2024/04/16
# Date updated: 2024/05/29
#=================================================

rm(list = ls())

library(MASS)  # multivariate normal distribution
library(baysc) # source functions from baysc package
library(tidyverse) # data wrangling
# library(baysc)  # swolca and data generation
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/"
data_dir <- "Data/"                 # Data directory
res_dir <- "Results/"               # Results directory
code_dir <- "Simulation_Code/"      # Simulation code directory

# Source simulation functions
source(paste0(wd, code_dir, "simulate_data_functions.R"))

# # Source functions from baysc package
# file_list <- c("simulate_data.R")
# invisible(lapply(file_list, function(x) 
#   source(paste0(wd, "Model_Code/", "baysc_functions/", x))))

#=================== SCENARIO 1: BASELINE OLD ==================================
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size 2000/40,000 (5%), non-overlapping patterns, X~S, 
# high overlap for pi_R and pi_B
scenario <- 1

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 40000     # Population size
pop_seed <- 1  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A3"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A3")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                           modal_theta_prob = modal_theta_prob,
                           formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4))))
# Items 29-30 are affected as follows: level 2 probability decrease with A1
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 1, 0, 0), A3 = rep(0, 4))))

V_unique <- as.data.frame(expand.grid(c_all = as.factor(1:K),
                                      A1 = c(-4, 0, 4),
                                      A3 = c(-8, 0, 8)))
round(get_categ_probs(beta_mat = beta_list_x[[1]], formula = formula_x,
V_unique = V_unique), 3)

### Generate population
n_B <- 2000  # Sample size for non-probability sample
n_R <- 2000  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = TRUE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}

#=================== SCENARIO 11: X indep of S =================================
# X !~ S
scenario <- 11

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 40000     # Population size
pop_seed <- 1  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x <- get_betas_x(profiles = profiles, R = R,
                           modal_theta_prob = modal_theta_prob,
                           formula_x = formula_x, V_unique = V_unique)
# V_unique <- as.data.frame(expand.grid(c_all = as.factor(1:K), 
#                                       A2 = c(-8, -1.5, 0, 1.5, 8)))
# get_categ_probs(beta_mat = beta_list_x[[1]], formula = formula_x, V_unique = V_unique)

### Generate population
n_B <- 2000  # Sample size for non-probability sample
n_R <- 2000  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = TRUE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}

#=================== SCENARIO 0: BASELINE ======================================
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size 2000/40,000 (5%), non-overlapping patterns, X~S, 
# high overlap for pi_R and pi_B
scenario <- 0

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 40000     # Population size
pop_seed <- 1  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                                modal_theta_prob = modal_theta_prob,
                                formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3 + c_all:A1"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0), 
                                             `c_all2:A1` = rep(0, 4),
                                             `c_all3:A1` = rep(0, 4)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4), 
        `c_all2:A1` = rep(0, 4),
        `c_all3:A1` = rep(0, 4))))
# Items 29-30 are affected as follows: A1 associated with r=3 for k=1, 
# r=2 for k=2, and r=1 for k=3
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 0, 2, 0), A3 = rep(0, 4), 
        `c_all2:A1` = c(0, 2, -2, 0), `c_all3:A1` = c(0, -1, -2, -1))))

# V_unique <- as.data.frame(expand.grid(c_all = as.factor(1:K),
#                                       A1 = c(-4, 0, 4),
#                                       A3 = c(-8, 0, 8)))
# round(get_categ_probs(beta_mat = beta_list_x[[1]], formula = formula_x,
#                       V_unique = V_unique), 3)

### Generate population
n_B <- 2000  # Sample size for non-probability sample
n_R <- 2000  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = TRUE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}

#=================== SCENARIO 7: OVERLAPPING PATTERNS ==========================
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size 2000/40,000 (5%), non-overlapping patterns, X~S, 
# high overlap for pi_R and pi_B
scenario <- 7

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 40000     # Population size
pop_seed <- 7  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(1, times = 0.4 * J), 
                                        rep(4, times = 0.1 * J),
                                        rep(3, times = 0.5 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                                modal_theta_prob = modal_theta_prob,
                                formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3 + c_all:A1"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0), 
                                             `c_all2:A1` = rep(0, 4),
                                             `c_all3:A1` = rep(0, 4)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4), 
        `c_all2:A1` = rep(0, 4),
        `c_all3:A1` = rep(0, 4))))
# Items 29-30 are affected as follows: A1 associated with r=3 for k=1, 
# r=2 for k=2, and r=1 for k=3
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 0, 2, 0), A3 = rep(0, 4), 
        `c_all2:A1` = c(0, 2, -2, 0), `c_all3:A1` = c(0, -1, -2, -1))))

# V_unique <- as.data.frame(expand.grid(c_all = as.factor(1:K),
#                                       A1 = c(-4, 0, 4),
#                                       A3 = c(-8, 0, 8)))
# round(get_categ_probs(beta_mat = beta_list_x[[1]], formula = formula_x,
#                       V_unique = V_unique), 3)

### Generate population
n_B <- 2000  # Sample size for non-probability sample
n_R <- 2000  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = TRUE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}

#=================== SCENARIO 8: SAMPLE SIZE = 1% ==============================
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size 400/40,000 (1%), non-overlapping patterns, X~S, 
# high overlap for pi_R and pi_B
scenario <- 8

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 40000     # Population size
pop_seed <- 8  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                                modal_theta_prob = modal_theta_prob,
                                formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3 + c_all:A1"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0), 
                                             `c_all2:A1` = rep(0, 4),
                                             `c_all3:A1` = rep(0, 4)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4), 
        `c_all2:A1` = rep(0, 4),
        `c_all3:A1` = rep(0, 4))))
# Items 29-30 are affected as follows: A1 associated with r=3 for k=1, 
# r=2 for k=2, and r=1 for k=3
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 0, 2, 0), A3 = rep(0, 4), 
        `c_all2:A1` = c(0, 2, -2, 0), `c_all3:A1` = c(0, -1, -2, -1))))

### Generate population
n_B <- 400  # Sample size for non-probability sample
n_R <- 400  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = TRUE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}


#=================== SCENARIO 9: SAMPLE SIZE PROSPECT ==========================
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size similar to PROSPECT, non-overlapping patterns, X~S, 
# high overlap for pi_R and pi_B
scenario <- 9

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 3000000     # Population size
pop_seed <- 9  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                                modal_theta_prob = modal_theta_prob,
                                formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3 + c_all:A1"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0), 
                                             `c_all2:A1` = rep(0, 4),
                                             `c_all3:A1` = rep(0, 4)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4), 
        `c_all2:A1` = rep(0, 4),
        `c_all3:A1` = rep(0, 4))))
# Items 29-30 are affected as follows: A1 associated with r=3 for k=1, 
# r=2 for k=2, and r=1 for k=3
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 0, 2, 0), A3 = rep(0, 4), 
        `c_all2:A1` = c(0, 2, -2, 0), `c_all3:A1` = c(0, -1, -2, -1))))

### Generate population
n_B <- 2000  # Sample size for non-probability sample
n_R <- 100000  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = TRUE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}


#=================== SCENARIO 10: LOW OVERLAP ==================================
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size 2000/40,000 (5%), non-overlapping patterns, X~S, 
# low overlap for pi_R and pi_B
scenario <- 10

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 40000     # Population size
pop_seed <- 10  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                                modal_theta_prob = modal_theta_prob,
                                formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3 + c_all:A1"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0), 
                                             `c_all2:A1` = rep(0, 4),
                                             `c_all3:A1` = rep(0, 4)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4), 
        `c_all2:A1` = rep(0, 4),
        `c_all3:A1` = rep(0, 4))))
# Items 29-30 are affected as follows: A1 associated with r=3 for k=1, 
# r=2 for k=2, and r=1 for k=3
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 0, 2, 0), A3 = rep(0, 4), 
        `c_all2:A1` = c(0, 2, -2, 0), `c_all3:A1` = c(0, -1, -2, -1))))

### Generate population
n_B <- 2000  # Sample size for non-probability sample
n_R <- 2000  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = FALSE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}


#=================== SCENARIO 13: n_B 1%, n_R 5% ===============================
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size 1% for n_B, 5% for n_R, non-overlapping patterns, X~S, 
# high overlap for pi_R and pi_B
scenario <- 13

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 40000     # Population size
pop_seed <- 13  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                                modal_theta_prob = modal_theta_prob,
                                formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3 + c_all:A1"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0), 
                                             `c_all2:A1` = rep(0, 4),
                                             `c_all3:A1` = rep(0, 4)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4), 
        `c_all2:A1` = rep(0, 4),
        `c_all3:A1` = rep(0, 4))))
# Items 29-30 are affected as follows: A1 associated with r=3 for k=1, 
# r=2 for k=2, and r=1 for k=3
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 0, 2, 0), A3 = rep(0, 4), 
        `c_all2:A1` = c(0, 2, -2, 0), `c_all3:A1` = c(0, -1, -2, -1))))

### Generate population
n_B <- 400  # Sample size for non-probability sample
n_R <- 2000  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = TRUE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}


#=================== SCENARIO 14: n_B 5%, n_R 1% ===============================
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size 5% for n_B, 1% for n_R, non-overlapping patterns, X~S, 
# high overlap for pi_R and pi_B
scenario <- 14

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 40000     # Population size
pop_seed <- 14  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                                modal_theta_prob = modal_theta_prob,
                                formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3 + c_all:A1"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0), 
                                             `c_all2:A1` = rep(0, 4),
                                             `c_all3:A1` = rep(0, 4)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4), 
        `c_all2:A1` = rep(0, 4),
        `c_all3:A1` = rep(0, 4))))
# Items 29-30 are affected as follows: A1 associated with r=3 for k=1, 
# r=2 for k=2, and r=1 for k=3
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 0, 2, 0), A3 = rep(0, 4), 
        `c_all2:A1` = c(0, 2, -2, 0), `c_all3:A1` = c(0, -1, -2, -1))))

### Generate population
n_B <- 2000  # Sample size for non-probability sample
n_R <- 400  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = TRUE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}


#=================== SCENARIO 18: SAMPLE SIZE PROSPECT low overlap =============
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size similar to PROSPECT, non-overlapping patterns, X~S, 
# high overlap for pi_R and pi_B
scenario <- 18

# Create data folder if it doesn't exist
dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
if (!dir.exists(dir_path)) {
  dir.create(file.path(dir_path))
}

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 3000000     # Population size
pop_seed <- 9  # Set seed

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                                modal_theta_prob = modal_theta_prob,
                                formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3 + c_all:A1"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0), 
                                             `c_all2:A1` = rep(0, 4),
                                             `c_all3:A1` = rep(0, 4)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4), 
        `c_all2:A1` = rep(0, 4),
        `c_all3:A1` = rep(0, 4))))
# Items 29-30 are affected as follows: A1 associated with r=3 for k=1, 
# r=2 for k=2, and r=1 for k=3
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 0, 2, 0), A3 = rep(0, 4), 
        `c_all2:A1` = c(0, 2, -2, 0), `c_all3:A1` = c(0, -1, -2, -1))))

### Generate population
n_B <- 2000  # Sample size for non-probability sample
n_R <- 100000  # Sample size for reference sample
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                          n_R = n_R, high_overlap = FALSE, formula_c = formula_c, 
                          beta_mat_c = beta_mat_c, formula_x = formula_x, 
                          beta_list_x = beta_list_x, pop_seed = pop_seed, 
                          save_res = TRUE, save_path = dir_path) 

### Generate samples
num_samps <- 100  # Number of samples
for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}

#=================== SCENARIOS 19, 20, 21: Low overlap sample sizes ============
# Variance adjustment, pi_R unknown and predicted using continuous BART, 
# sample size 1% + low n_B + low n_R, non-overlapping patterns, X~S, 
# low overlap for pi_R and pi_B
scenario_list <- c(19, 20, 21)

### General parameters
rho <- 0.5     # Correlation between selection variables A1 and A2
N <- 40000     # Population size

### Parameters for generating categorical latent class assignment C
formula_c <- "~ A1 + A2 + A1A2"
beta_mat_c <- matrix(c(0, 0, 0, 0, 
                       0.4, -0.5, 0.75, 0.1,  
                       -0.2, -1, 1.2, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2", "A1A2")

### Parameters for generating observed manifest variables X
J <- 30; R <- 4; K <- 3
formula_x <- "~ c_all"
V_unique <- data.frame(c_all = as.factor(1:K))
profiles <- as.matrix(data.frame(C1 = c(rep(1, times = 0.5 * J),
                                        rep(3, times = 0.5 * J)),
                                 C2 = c(rep(4, times = 0.2 * J),
                                        rep(2, times = 0.8 * J)),
                                 C3 = c(rep(3, times = 0.3 * J),
                                        rep(4, times = 0.4 * J),
                                        rep(1, times = 0.3 * J))))
modal_theta_prob <- 0.85
beta_list_x_temp <- get_betas_x(profiles = profiles, R = R,
                                modal_theta_prob = modal_theta_prob,
                                formula_x = formula_x, V_unique = V_unique)
# Add in coefficients for A3, updating formula_x and beta_list_x
formula_x <- "~ c_all + A1 + A3 + c_all:A1"
# Items 1-2 are affected in the following manner: 
# level 2 probability increases as A3 increases
beta_list_x <- lapply(1:2, function(j) cbind(beta_list_x_temp[[j]],
                                             A1 = rep(0, 4),
                                             A3 = c(0, 0.5, 0, 0), 
                                             `c_all2:A1` = rep(0, 4),
                                             `c_all3:A1` = rep(0, 4)))
beta_list_x <- c(beta_list_x, lapply(3:(J-2), function(j) 
  cbind(beta_list_x_temp[[j]], A1 = rep(0, 4), A3 = rep(0, 4), 
        `c_all2:A1` = rep(0, 4),
        `c_all3:A1` = rep(0, 4))))
# Items 29-30 are affected as follows: A1 associated with r=3 for k=1, 
# r=2 for k=2, and r=1 for k=3
beta_list_x <- c(beta_list_x, lapply((J-1):J, function(j) 
  cbind(beta_list_x_temp[[j]], A1 = c(0, 0, 2, 0), A3 = rep(0, 4), 
        `c_all2:A1` = c(0, 2, -2, 0), `c_all3:A1` = c(0, -1, -2, -1))))


for (i in 1:length(scenario_list)) {
  scenario <- scenario_list[i]
  pop_seed <- scenario  # Set seed
  
  # Create data folder if it doesn't exist
  dir_path <- paste0(wd, data_dir, "scen_", scenario, "/")
  if (!dir.exists(dir_path)) {
    dir.create(file.path(dir_path))
  }
  
  ### Generate population
  n_B <- 400  # Sample size for non-probability sample
  n_R <- 400  # Sample size for reference sample
  if (scenario == 20) {
    n_R <- 2000
  } else if (scenario == 21) {
    n_B <- 2000
  }
  sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho, n_B = n_B, 
                            n_R = n_R, high_overlap = FALSE, formula_c = formula_c, 
                            beta_mat_c = beta_mat_c, formula_x = formula_x, 
                            beta_list_x = beta_list_x, pop_seed = pop_seed, 
                            save_res = TRUE, save_path = dir_path) 
  
  ### Generate samples
  num_samps <- 100  # Number of samples
  for (i in 1:num_samps) {
    sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, 
                                 scenario = scenario, save_res = TRUE, samp_seed = i,
                                 save_path = dir_path)
  }
}




#================= Check values for simulated population =======================
# Selection probabilities and weights
summary(sim_pop$pi_R)
summary(sim_pop$pi_B)
summary(1 / sim_pop$pi_R)
summary(1 / sim_pop$pi_B)
hist(sim_pop$pi_R, breaks = 30)
hist(sim_pop$pi_B, breaks = 30)
cor(sim_pop$pop$A1, sim_pop$pi_R)
cor(sim_pop$pop$A2, sim_pop$pi_R)
cor(sim_pop$pop$A3, sim_pop$pi_R)
cor(sim_pop$pop$A1, sim_pop$pi_B)
cor(sim_pop$pop$A2, sim_pop$pi_B)
cor(sim_pop$pop$A3, sim_pop$pi_B)

# Latent class assignment correlations with selection covariates
cor(as.numeric(sim_pop$c_all), sim_pop$pop$A1)
cor(as.numeric(sim_pop$c_all), sim_pop$pop$A2)
cor(as.numeric(sim_pop$c_all), sim_pop$pop$A3)

# Observed manifest variables correlations with selection covariates
cor(sim_pop$X_data, sim_pop$pop$A1)
cor(sim_pop$X_data, sim_pop$pop$A2)
cor(sim_pop$X_data, sim_pop$pop$A3)

# Observed thetas
# item 1 has association with A3 (higher level 2)
t(sapply(1:3, function(k) prop.table(table(sim_pop$X_data[sim_pop$c_all == k, 1]))))
# # item 4 has association with A2 (higher level 4)
# t(sapply(1:3, function(k) prop.table(table(sim_pop$X_data[sim_pop$c_all == k, 4]))))
# item 30 has association with A1 (higher level 2)
t(sapply(1:3, function(k) prop.table(table(sim_pop$X_data[sim_pop$c_all == k, 30]))))
# item 5 has no associations
t(sapply(1:3, function(k) prop.table(table(sim_pop$X_data[sim_pop$c_all == k, 5]))))


# Plot overlaps
plot_overlap(wd = wd, data_dir = data_dir, scenario = 0, samp_i = 1)
plot_overlap(wd = wd, data_dir = data_dir, scenario = 10, samp_i = 1)

# Obtain number of overlapping individuals
get_n_overlap <- function(wd, data_dir, scenario, num_samps) {
  num_overlap <- numeric(num_samps)
  prop_overlap <- numeric(num_samps)
  for (i in 1:num_samps) {
    samp_i <- i
    # Load convenience sample
    load(paste0(wd, data_dir,  "scen_", scenario, "/", "sim_samp_", samp_i, 
                "_B_wolcan.RData"))
    # Load reference sample
    load(paste0(wd, data_dir,  "scen_", scenario, "/", "sim_samp_", samp_i, 
                "_R_wolcan.RData"))
    ind_overlap <- intersect(sim_samp_B$ind_B, sim_samp_R$ind_R)
    num_overlap[i] <- length(ind_overlap)
    prop_overlap[i] <- length(ind_overlap) / 
      length(unique(c(sim_samp_B$ind_B, sim_samp_R$ind_R)))
  }
  res_overlap <- list(num_overlap = num_overlap, prop_overlap = prop_overlap)
  return(res_overlap)
}

scen_0 <- get_n_overlap(wd = wd, data_dir = data_dir, scenario = 0, 
                        num_samps = 100)
scen_10 <- get_n_overlap(wd = wd, data_dir = data_dir, scenario = 10, 
                         num_samps = 100)
scen <- scen_10
mean(scen$num_overlap)
summary(scen$num_overlap)
mean(scen$prop_overlap)
summary(scen$prop_overlap)

# Plot overlap
df_overlap <- data.frame(High = 100 * scen_0$prop_overlap, 
                         Low = 100 * scen_10$prop_overlap)
df_long <- df_overlap %>%
  pivot_longer(cols = everything(), names_to = "Scenario", values_to = "Percent")
df_long %>%
  ggplot(aes(x = Scenario, y = Percent, fill = Scenario)) + 
  geom_violin() + theme_bw() + ylab("Percent Overlap") + 
  theme(legend.position = "none") + ylim(0, 8)


