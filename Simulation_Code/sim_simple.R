#=================================================
# Generate population and sample data for NPSWOLCA
# Simplified logistic regression setting
# Author: Stephanie Wu
# Date created: 2024/04/16
# Date updated: 2024/04/16
#=================================================

rm(list = ls())

library(MASS)  # multivariate normal distribution
library(baysc)  # swolca and data generation
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
data_dir <- "Data/"                 # Data directory
res_dir <- "Results/"               # Results directory
code_dir <- "Simulation_Code/"      # Simulation code directory

set.seed(1)

#=================== Generate population =======================================
# 10% sample for PS and NPS
rho <- 0.8
N <- 5000
n_B <- 500
n_R <- 500

### Helper function to obtain offset to match desired sample size
root_n <- function(x, n) {
  uniroot(function(t){mean(exp(t + x) / (1 + exp(t + x))) - (n / N)}, 
          interval = c(-50, 50))$root
}

### Generate covariates A1 and A2 in the population
cov <- matrix(c(1, rho, rho, 1), 2, 2)
mu <- c(0, 0)
pop <- as.data.frame(mvrnorm(n = N, mu = mu, Sigma = cov))
names(pop) <- c("A1", "A2")
pop$A12 <- pop$A1^2
pop$A1A2 <- pop$A1*pop$A2
pop$logA2 <- log(abs(pop$A2))

offset_R <- root_n(x = (- 0.5 * pop$A1 + pop$A2), n = n_R)
offset_B <- root_n(x = (- pop$A1 + 0.3 * pop$A2), n = n_B)

### Generate selection probabilities as simple logistic regression
pop$pi_R <- exp(offset_R - 0.5 * pop$A1 + pop$A2) / 
  (1 + exp(offset_R - 0.5 * pop$A1 + pop$A2))
pop$pi_B <- exp(offset_R - pop$A1 + 0.3 * pop$A2) / 
  (1 + exp(offset_R - pop$A1 + 0.3 * pop$A2))
summary(pop$pi_R)
summary(pop$pi_B)
summary(1 / pop$pi_R)
summary(1 / pop$pi_B)
hist(pop$pi_R)
hist(pop$pi_B)

### Generate categorical latent class assignment C
formula_c <- "~ A1 + A2"
beta_mat_c <- matrix(c(0, 0, 0, 
                       0.4, -0.5, -0.275,  
                       -0.2, -1, 0.25), nrow = 3, byrow = TRUE)
colnames(beta_mat_c) <- c("Intercept", "A1", "A2")
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
# Check correlations
cor(c_all, pop$A1)
cor(c_all, pop$A2)

### Generate observed manifest variables X
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
cor(X_data, pop$A1)
cor(X_data, pop$A2)
# Add to population data

sim_pop <- list(A1 = pop$A1, A2 = pop$A2, A12 = pop$A12, A1A2 = pop$A1A2, 
                logA2 = pop$logA2, pi_R = pop$pi_R, pi_B = pop$pi_B, 
                c_all = c_all, X_data = X_data, 
                true_global_thetas = true_global_thetas, 
                true_global_patterns = true_global_patterns)

save(sim_pop, file = paste0(wd, data_dir, "sim_pop_simple_wolca.RData"))

#================== Generate samples ===========================================
### Create samples
load(paste0(wd, data_dir, "sim_pop_simple_wolca.RData"))

num_samps <- 10

for (i in 1:num_samps) {
  set.seed(i)
  cat("Creating sample ", i, "\n")
  # Non-probability sample: n_B
  delta_B <- rbinom(n = N, size = 1, prob = sim_pop$pi_B) 
  # Reference probability sample: n_R
  delta_R <- rbinom(n = N, size = 1, prob = sim_pop$pi_R)
  cat("n_B:", sum(delta_B), "\n")
  cat("n_R", sum(delta_R), "\n")
  # Selected individuals for each sample
  ind_B <- which(delta_B == 1)
  ind_R <- which(delta_R == 1)
        # # Sampling probabilities for selected individuals in each sample
        # samp_pi_B <- sim_pop$pi_B[ind_B]
        # samp_pi_R <- sim_pop$pi_R[ind_R]
        
        # # Renormalize probabilities so sum(1/probs) = N
        # const_B <- sum(1 / samp_pi_B) / N
        # pi_B_renorm <- samp_pi_B * const_B
        # sum(1 / pi_B_renorm)
        # const_R <- sum(1 / samp_pi_R) / N
        # pi_R_renorm <- samp_pi_R * const_R
        # sum(1 / pi_R_renorm)
  
  # Subset to the sample data
  sim_samp_B <- list(ind_B = ind_B, 
                     covs = pop[ind_B, c("A1", "A2", "A12", "A1A2", "logA2")], 
                     true_pi_B = sim_pop$pi_B[ind_B],  # pi_B for those in S_B
                     true_pi_R = sim_pop$pi_R[ind_B],  # pi_R for those in S_B
                     c_all = sim_pop$c_all[ind_B], 
                     X_data = sim_pop$X_data[ind_B, ], delta_B = delta_B)
  sim_samp_R <- list(ind_R = ind_R, 
                     covs = pop[ind_R, c("A1", "A2", "A12", "A1A2", "logA2")], 
                     true_pi_B = sim_pop$pi_B[ind_R],  # pi_B for those in S_R
                     true_pi_R = sim_pop$pi_R[ind_R],  # pi_R for those in S_R
                     c_all = sim_pop$c_all[ind_R], 
                     X_data = sim_pop$X_data[ind_R, ], delta_R = delta_R)
  # Check sum of weights
  cat("wt_B", sum(1 / sim_samp_B$true_pi_B), "\n")
  cat("wt_R", sum(1 / sim_samp_R$true_pi_R), "\n")
  
  save(sim_samp_B, file = paste0(wd, data_dir, "sim_samp_simple", i, "_B_wolca.RData"))
  save(sim_samp_R, file = paste0(wd, data_dir, "sim_samp_simple", i, "_R_wolca.RData"))
}

