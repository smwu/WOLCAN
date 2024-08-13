#===========================================================
# Run and summarize scenarios for predicting pseudo-weights
# Author: Stephanie Wu
# Date created: 2024/08/10
# Date updated: 2024/08/10
#===========================================================

rm(list = ls())

library(tidyverse)
library(gridExtra)
library(knitr)
library(kableExtra)
library(BART)  # BART

# Set directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/"
data_dir <- "Data/"    # Data directory
res_dir <- "Results/"  # Results directory
code_dir <- "Summary_Code/"  # Model code directory
sum_dir <- "Summary_Results/"  # Summary results directory

# Source R model functions
source(paste0(wd, "Model_Code/model_functions.R"))

#============== Run and store weights ==========================================

### Default settings
hat_pi_R <- NULL  # RS selection probabilities for those in NPS
frame_B <- 1  # Coverage probability of NPS frame
frame_R <- 1  # Coverage probability of RS frame
trim_method = "t2"  # Trimming method using IQR
trim_c = 20         # Trimming constant

### Models and scenarios
# High overlap: 0, 8, 13, 14
# Low overlap: 10, 19, 20, 21
data_scen <- list(0, 8, 13, 14, 10, 19, 20, 21)
model_scen <- list("none", "logistic", "bart_500", "bart_1000", "bart_2000", 
                   "logistic_cov", "bart_1000_cov")

# Total number of scenarios
num_scen <- length(model_scen) * length(data_scen)
num_samps <- 100 # Number of samples

# Initialize results data frame
weights_res <- as.data.frame(matrix(NA, nrow = num_scen, ncol = num_samps))

# Counter to keep track of total scenarios
counter <- 1
for (i in 1:length(data_scen)) {  # Data-generating scenario
  scenario <- data_scen[[i]]
      # load(paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData"))
  
  for (j in 1:length(model_scen)) {  # Prediction model for weights
    model <- model_scen[[j]]
    # Covariates to include in prediction model
    if (model %in% c("bart_1000_cov", "logistic_cov")) {
      pred_covs_B <- pred_covs_R <- c("A1", "A2")
    } else {
      pred_covs_B <- pred_covs_R <- c("A1", "A2", "A1A2", "A3")
    }
    # Number of BART samples to use
    if (model == "bart_500") {
      num_post <- 500
    } else if (model == "bart_2000") {
      num_post <- 2000
    } else {
      num_post <- 1000
    }
    # Type of prediction model
    if (model %in% c("logistic", "logistic_cov")) {  # Logistic regression
      pred_model <- "glm"
    } else {  # BART model
      pred_model <- "bart"
    }
    
    for (k in 1:num_samps) {
      samp_i <- k
      load(paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", samp_i, "_B_wolcan.RData"))
      load(paste0(wd, data_dir, "scen_", scenario, "/sim_samp_", samp_i, "_R_wolcan.RData"))
      
      # Other variables
      x_mat <- sim_samp_B$X_data  # Multivariate categorical variables
      dat_B <- data.frame(sim_samp_B$covs)  # Covariates for NPS
      dat_R <- data.frame(sim_samp_R$covs)  # Covariates for reference
      pi_R <- sim_samp_R$pi_R  # RS selection probabilites for those in RS
      
      set.seed(samp_i)  # set seed
      
      if (model == "none") {  # no weights
        wts <- rep(1, nrow(x_mat))
        res <- mean(abs(wts - 1/(sim_samp_B$true_pi_B)))
        
      } else {
        invisible(capture.output(
          est_weights <- get_weights_bart(dat_B = dat_B, dat_R = dat_R, 
                                          pred_covs_B = pred_covs_B, 
                                          pred_covs_R = pred_covs_R, 
                                          num_post = num_post, pi_R = pi_R, 
                                          hat_pi_R = hat_pi_R,
                                          frame_B = frame_B, frame_R = frame_R,
                                          trim_method = trim_method, trim_c = trim_c,
                                          pred_model = pred_model)))
        
        # Mean absolute difference of weights
        res <- mean(abs(est_weights$wts - (1/sim_samp_B$true_pi_B)))
      }
      
      weights_res[counter, k] <- res
    }
    # Increment counter and print progress
    counter <- counter + 1
    print(paste0("Scenario ", scenario, " model ", model, " done!"))
    
    temp <- weights_res[counter, ]
    save(temp, file = paste0(wd, sum_dir, "scen_", scenario, "model_", model, 
                             ".RData"))
  }
}

rownames(weights_res) <- paste0(rep(data_scen, each = length(model_scen)), "_", 
                                rep(model_scen, times = length(data_scen)))

save(weights_res, file = paste0(wd, sum_dir, "weights_res.RData"))

# Reshape data for tables
weights_res_mean <- data.frame(mean_abs_bias = rowMeans(weights_res))
weights_res_mean$scenario <- rep(data_scen, each = length(model_scen))
weights_res_mean$model <- rep(model_scen, times = length(data_scen))
weights_res_mean <- weights_res_mean %>% 
  pivot_wider(names_from = model, values_from = mean_abs_bias)
weights_res_mean <- weights_res_mean %>%
  mutate(scenario = unlist(scenario))


# Sample size 5%: scenario 0 and 10
tb_1 <- as.data.frame(weights_res_mean) %>% 
  filter(scenario %in% c(0, 10)) %>%
  select(-scenario)
rownames(tb_1) <- c("High Overlap", "Low Overlap")
tb_1 %>%
  kbl(caption = "5% Sample Size", digits = 3) %>%
  kable_classic(full_width = F)

# Sample size 1%: scenario 8 and 19
tb_2 <- as.data.frame(weights_res_mean) %>% 
  filter(scenario %in% c(8, 19)) %>%
  select(-scenario)
rownames(tb_2) <- c("High Overlap", "Low Overlap")
tb_2 %>%
  kbl(caption = "1% Sample Size", digits = 3) %>%
  kable_classic(full_width = F)

# Sample size n_R 5% n_B 1%: scenario 13 and 20
tb_3 <- as.data.frame(weights_res_mean) %>% 
  filter(scenario %in% c(13, 20)) %>%
  select(-scenario)
rownames(tb_3) <- c("High Overlap", "Low Overlap")
tb_3 %>%
  kbl(caption = "n_R 5%, n_B 1%", digits = 3) %>%
  kable_classic(full_width = F)

# Sample size n_R 1% n_B 5%: scenario 14 and 21
tb_4 <- as.data.frame(weights_res_mean) %>% 
  filter(scenario %in% c(14, 21)) %>%
  select(-scenario)
rownames(tb_4) <- c("High Overlap", "Low Overlap")
tb_4 %>%
  kbl(caption = "n_R 1%, n_B 5%", digits = 3) %>%
  kable_classic(full_width = F)


# Plot comparison of BART number of posterior samples
weights_res_mean %>%
  pivot_longer(cols = bart_500:bart_2000, names_to = "num_post", 
               values_to = "mean_abs_bias") %>%
  mutate(`BART Samples` = factor(num_post, levels = c("bart_500", "bart_1000", "bart_2000"),
                           labels = c("500", "1000", "2000"))) %>%
  ggplot(aes(x = as.factor(scenario), y = mean_abs_bias, fill = `BART Samples`)) + 
  theme_bw() + geom_bar(stat = "identity", position = position_dodge()) + 
  xlab("Scenario") + ylab("Mean Absolute Bias for Weights")
