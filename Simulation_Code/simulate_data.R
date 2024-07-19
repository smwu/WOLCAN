#=================================================
# Generate population and sample data for NPSWOLCA
# Author: Stephanie Wu
# Date created: 2024/04/16
# Date updated: 2024/05/29
#=================================================

rm(list = ls())

library(MASS)  # multivariate normal distribution
library(baysc) # source functions from baysc package
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

#=================== Generate population: SCENARIO 1 ===========================
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


#=================== Generate population: SCENARIO 11 ==========================
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


#================= Check values for simulated population =======================
# Selection probabilities and weights
summary(sim_pop$pi_R)
summary(sim_pop$pi_B)
summary(1 / sim_pop$pi_R)
summary(1 / sim_pop$pi_B)
hist(sim_pop$pi_R)
hist(sim_pop$pi_B)
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



#================== Generate samples ===========================================
### Load population data 
scenario <- 1
load(paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData"))
num_samps <- 100  # Number of samples

### Parameters
n_B <- 2000  # Sample size for non-probability sample
n_R <- 2000  # Sample size for reference sample

for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, n_B = n_B, n_R = n_R, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = dir_path)
}




