#=================================================
# Generate population and sample data for NPSWOLCA
# Author: Stephanie Wu
# Date created: 2024/04/16
# Date updated: 2024/05/29
#=================================================

rm(list = ls())

library(MASS)  # multivariate normal distribution
library(baysc)  # swolca and data generation
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
data_dir <- "Data/"                 # Data directory
res_dir <- "Results/"               # Results directory
code_dir <- "Simulation_Code/"      # Simulation code directory

# Source simulation functions
source(paste0(wd, code_dir, "simulate_data_functions.R"))

#=================== Generate population =======================================

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
# formula_c <- "~ A1 + A12 + A2 + A1A2"
# beta_mat_c <- matrix(c(0, 0, 0, 0, 0, 
#                        0.4, -0.5, -0.275, 0.1, 0.175, 
#                        -0.2, -1, 0.25, 1.5, 0.15), nrow = 3, byrow = TRUE)
# colnames(beta_mat_c) <- c("Intercept", "A1", "A12", "A2", "A1A2")

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
# Add in coefficients for S to each j matrix, updating formula_x and beta_list_x
beta_list_x_temp <- beta_list_x
formula_x <- "~ c_all + A2"
# Items 1-3 are affected in the following manner: 
# level 4 probability increases as A2 increases
beta_list_x <- lapply(1:3, function(j) cbind(beta_list_x_temp[[j]],
                                             A2 = c(0, 0, 0, 5)))
beta_list_x <- c(beta_list_x, lapply(4:J, function(j) cbind(beta_list_x_temp[[j]],
                                                            A2 = rep(0, 4))))
# beta_list_x <- lapply(1:3, function(j) cbind(beta_list_x_temp[[j]],
#                                              s_all = c(0, 0.7, -0.7, -0.7)))
# beta_list_x <- c(beta_list_x, lapply(4:J, function(j) cbind(beta_list_x_temp[[j]],
#                                                             s_all = rep(0, 4))))
# V_unique <- as.data.frame(expand.grid(c_all = as.factor(1:K), 
#                                       A2 = c(-8, -1.5, 0, 1.5, 8)))
# get_categ_probs(beta_mat = beta_list_x[[1]], formula = formula_x, V_unique = V_unique)

### Generate population
sim_pop <- sim_pop_wolcan(N = N, J = J, K = K, R = R, rho = rho,
                      high_overlap = TRUE, formula_c = formula_c, 
                      beta_mat_c = beta_mat_c, formula_x = formula_x, 
                      beta_list_x = beta_list_x, pop_seed = pop_seed, 
                      save_res = TRUE, save_path = paste0(wd, data_dir, "scen_1/")) 

### Check values for simulated data
# Selection probabilities and weights
summary(sim_pop$pi_R)
summary(sim_pop$pi_B)
summary(1 / sim_pop$pi_R)
summary(1 / sim_pop$pi_B)
hist(sim_pop$pi_R)
hist(sim_pop$pi_B)
cor(sim_pop$A1, sim_pop$pi_R)
cor(sim_pop$A2, sim_pop$pi_R)
cor(sim_pop$A1, sim_pop$pi_B)
cor(sim_pop$A2, sim_pop$pi_B)

# Latent class assignment correlations with selection covariates
cor(sim_pop$c_all, sim_pop$A1)
cor(sim_pop$c_all, sim_pop$A2)

# Observed manifest variables correlations with selection covariates
cor(sim_pop$X_data, sim_pop$A1)
cor(sim_pop$X_data, sim_pop$A2)

# Observed thetas
# item 1 has association with A2
t(sapply(1:3, function(k) prop.table(table(sim_pop$X_data[sim_pop$c_all == k, 1]))))
# item 4 does not have association with A2
t(sapply(1:3, function(k) prop.table(table(sim_pop$X_data[sim_pop$c_all == k, 4]))))

#================== Generate samples ===========================================
### Load population data 
scenario <- 1
load(paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData"))
num_samps <- 10  # Number of samples

### Parameters
n_B <- 2000  # Sample size for non-probability sample
n_R <- 2000  # Sample size for reference sample

for (i in 1:num_samps) {
  sim_samps <- sim_samp_wolcan(i = i, sim_pop = sim_pop, n_B = n_B, n_R = n_R, 
                               scenario = scenario, save_res = TRUE, samp_seed = i,
                               save_path = paste0(wd, data_dir, "scen_1/"))
}




