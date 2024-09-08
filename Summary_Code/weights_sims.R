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
pi_res <- as.data.frame(matrix(NA, nrow = num_scen, ncol = num_samps))

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
        pi <- mean(abs(1/wts - sim_samp_B$true_pi_B))
        
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
        # Mean absolute different of inclusion probabilities
        pi <- mean(abs(est_weights$hat_pi_B - sim_samp_B$true_pi_B))
      }
      
      weights_res[counter, k] <- res
      pi_res[counter, k] <- pi
    }
    # Increment counter and print progress
    counter <- counter + 1
    print(paste0("Scenario ", scenario, " model ", model, " done!"))
    
    temp <- weights_res[counter, ]
    save(temp, file = paste0(wd, sum_dir, "scen_", scenario, "/model_", model, 
                             ".RData"))
    temp2 <- pi_res[counter, ]
    save(temp2, file = paste0(wd, sum_dir, "scen_", scenario, "/model_", model, 
                              "_pi.RData"))
  }
}

rownames(weights_res) <- paste0(rep(data_scen, each = length(model_scen)), "_", 
                                rep(model_scen, times = length(data_scen)))
rownames(pi_res) <- paste0(rep(data_scen, each = length(model_scen)), "_", 
                                rep(model_scen, times = length(data_scen)))

save(weights_res, file = paste0(wd, sum_dir, "weights_res.RData"))
save(pi_res, file = paste0(wd, sum_dir, "pi_res.RData"))

#load(paste0(wd, sum_dir, "weights_res.RData"))

# Reshape data for tables
weights_res_mean <- data.frame(mean_abs_bias = rowMeans(weights_res))
weights_res_mean$scenario <- rep(data_scen, each = length(model_scen))
weights_res_mean$model <- rep(model_scen, times = length(data_scen))
weights_res_mean <- weights_res_mean %>% 
  pivot_wider(names_from = model, values_from = mean_abs_bias)
weights_res_mean <- weights_res_mean %>%
  mutate(scenario = unlist(scenario))

pi_res_mean <- data.frame(mean_abs_bias = rowMeans(pi_res))
pi_res_mean$scenario <- rep(data_scen, each = length(model_scen))
pi_res_mean$model <- rep(model_scen, times = length(data_scen))
pi_res_mean <- pi_res_mean %>% 
  pivot_wider(names_from = model, values_from = mean_abs_bias)
pi_res_mean <- pi_res_mean %>%
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


# For pi
# Sample size 5%: scenario 0 and 10
tb_1 <- as.data.frame(pi_res_mean) %>% 
  filter(scenario %in% c(0, 10)) %>%
  select(-scenario)
rownames(tb_1) <- c("High Overlap", "Low Overlap")
tb_1 %>%
  kbl(caption = "5% Sample Size", digits = 3) %>%
  kable_classic(full_width = F)

# Sample size 1%: scenario 8 and 19
tb_2 <- as.data.frame(pi_res_mean) %>% 
  filter(scenario %in% c(8, 19)) %>%
  select(-scenario)
rownames(tb_2) <- c("High Overlap", "Low Overlap")
tb_2 %>%
  kbl(caption = "1% Sample Size", digits = 3) %>%
  kable_classic(full_width = F)

# Sample size n_R 5% n_B 1%: scenario 13 and 20
tb_3 <- as.data.frame(pi_res_mean) %>% 
  filter(scenario %in% c(13, 20)) %>%
  select(-scenario)
rownames(tb_3) <- c("High Overlap", "Low Overlap")
tb_3 %>%
  kbl(caption = "n_R 5%, n_B 1%", digits = 3) %>%
  kable_classic(full_width = F)

# Sample size n_R 1% n_B 5%: scenario 14 and 21
tb_4 <- as.data.frame(pi_res_mean) %>% 
  filter(scenario %in% c(14, 21)) %>%
  select(-scenario)
rownames(tb_4) <- c("High Overlap", "Low Overlap")
tb_4 %>%
  kbl(caption = "n_R 1%, n_B 5%", digits = 3) %>%
  kable_classic(full_width = F)

# Rearrange and print full weights df
displ <- weights_res_mean %>%
  slice(match(c(0, 10, 14, 21, 13, 20, 8, 19), scenario)) %>%
  mutate(Overlap = rep(c("High", "Low"), times = 4)) %>%
  select(Overlap, none, logistic, logistic_cov, bart_500:bart_2000, bart_1000_cov)
displ <- as.data.frame(displ)
colnames(displ) <- c("Overlap", "No Model", "LogReg", "LogRegMiss", "BART500", 
                     "BART1000", "BART2000", "BART1000Miss")
displ$`Sample Size` <- rep(c("$n_B 5\\%, n_R 5\\%$", "$n_B 5\\%, n_R 1\\%$", 
                      "$n_B 1\\%, n_R 5\\%$", "$n_B 1\\%, n_R 1\\%$"), each = 2)
displ <- displ %>% select(`Sample Size`, Overlap:`BART1000Miss`)
displ %>%
  kbl(digits = 3, format = "html", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE)

#================ Plot comparison for high and low overlap
plot_df <- displ 
plot_df %>%
  pivot_longer(cols = `No Model`:BART1000Miss, names_to = "Model", 
               values_to = "mean_abs_bias") %>%
  mutate(Model = factor(Model, 
    levels = c("No Model", "LogReg", "BART500",  "BART1000", "BART2000", 
               "LogRegMiss","BART1000Miss"))) %>%
  ggplot(aes(x = `Sample Size`, y = mean_abs_bias, fill = Model)) + 
  theme_bw() + geom_bar(stat = "identity", position = position_dodge()) + 
  facet_grid(.~Overlap, labeller = as_labeller(c("High" = "High Overlap",
                                                 "Low" = "Low Overlap"))) + 
  scale_fill_brewer(palette = "Set2") + 
  ylab("Mean Absolute Bias for Pseudo-Weights") + 
  labs(fill = "Model") + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="aliceblue")) +
  guides(fill = guide_legend(nrow = 1)) + 
  scale_x_discrete(name = "Sample Size",
                   limits = c("$n_B 5\\%, n_R 5\\%$", 
                              "$n_B 5\\%, n_R 1\\%$",
                              "$n_B 1\\%, n_R 5\\%$",
                              "$n_B 1\\%, n_R 1\\%$"),
                   labels = c(expression(atop(n[B]*" 5%", n[R]*" 5%")),
                              expression(atop(n[B]*" 5%", n[R]*" 1%")),
                              expression(atop(n[B]*" 1%", n[R]*" 5%")),
                              expression(atop(n[B]*" 1%", n[R]*" 1%"))))
ggsave(paste0(wd, "Figures/weights_sims.png"), width = 8.5, height = 5, 
       dpi = 700, device = png)

#================ Plot comparison of BART number of posterior samples
weights_res_mean %>%
  pivot_longer(cols = bart_500:bart_2000, names_to = "num_post", 
               values_to = "mean_abs_bias") %>%
  mutate(`BART Samples` = factor(num_post, levels = c("bart_500", "bart_1000", "bart_2000"),
                           labels = c("500", "1000", "2000"))) %>%
  ggplot(aes(x = as.factor(scenario), y = mean_abs_bias, fill = `BART Samples`)) + 
  theme_bw() + geom_bar(stat = "identity", position = position_dodge()) + 
  xlab("Scenario") + ylab("Mean Absolute Bias for Weights")

#================ Plot single realization of high overlap vs low overlap
# High overlap
load(paste0(wd, data_dir, "scen_0/", "sim_pop_wolcan.RData"))
load(paste0(wd, data_dir, "scen_0/", "sim_samp_1_B_wolcan.RData"))
load(paste0(wd, data_dir, "scen_0/", "sim_samp_1_R_wolcan.RData"))
n_B <- length(sim_samp_B$ind_B)
n_R <- length(sim_samp_R$ind_R)
comb_df <- data.frame(unit_sample = c(rep("S[B]", n_B), rep("S[R]", n_R)))
comb_df$pi_B <- c(sim_samp_B$true_pi_B, sim_pop$pi_B[sim_samp_R$ind_R])
comb_df$pi_R <- c(sim_samp_B$true_pi_R, sim_pop$pi_R[sim_samp_R$ind_R])
# Sort by ind_B, then ind_R, and then sort within by pi_B
comb_df <- comb_df %>%
  arrange(unit_sample, pi_B) %>%
  mutate(Unit = 1:nrow(comb_df))
# Convert to long format
comb_df_long <- comb_df %>%
  pivot_longer(cols = pi_B:pi_R, names_to = "Sample", 
               values_to = "Inclusion Probability")
comb_df_long %>% 
  ggplot(aes(x = Unit, y = `Inclusion Probability`)) + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2", 
                     labels = c(expression(pi[i]^B), expression(pi[i]^R))) + 
  facet_grid(. ~ unit_sample, scales = "free_x", 
             labeller = label_parsed) + 
  # facet_grid(Sample ~ unit_sample, scales = "free_x") + 
  geom_line(aes(size = Sample, alpha = Sample, col = Sample)) + 
  scale_size_manual(values = c(0.6, 0.3), 
                    labels = c(expression(pi[i]^B), expression(pi[i]^R))) + 
  scale_alpha_manual(values = c(1, 0.6), 
                     labels = c(expression(pi[i]^B), expression(pi[i]^R))) + 
  theme(strip.background = element_rect(fill="aliceblue"))

ggsave(paste0(wd, "Figures/high_overlap.png"), width = 8.5, height = 3, 
       dpi = 700, device = png)

# Low overlap
load(paste0(wd, data_dir, "scen_10/", "sim_pop_wolcan.RData"))
load(paste0(wd, data_dir, "scen_10/", "sim_samp_1_B_wolcan.RData"))
load(paste0(wd, data_dir, "scen_10/", "sim_samp_1_R_wolcan.RData"))
n_B <- length(sim_samp_B$ind_B)
n_R <- length(sim_samp_R$ind_R)
comb_df <- data.frame(unit_sample = c(rep("S[B]", n_B), rep("S[R]", n_R)))
comb_df$pi_B <- c(sim_samp_B$true_pi_B, sim_pop$pi_B[sim_samp_R$ind_R])
comb_df$pi_R <- c(sim_samp_B$true_pi_R, sim_pop$pi_R[sim_samp_R$ind_R])
# Sort by ind_B, then ind_R, and then sort within by pi_B
comb_df <- comb_df %>%
  arrange(unit_sample, pi_B) %>%
  mutate(Unit = 1:nrow(comb_df))
# Convert to long format
comb_df_long <- comb_df %>%
  pivot_longer(cols = pi_B:pi_R, names_to = "Sample", 
               values_to = "Inclusion Probability")
comb_df_long %>% 
  ggplot(aes(x = Unit, y = `Inclusion Probability`)) + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2", 
                     labels = c(expression(pi[i]^B), expression(pi[i]^R))) + 
  facet_grid(. ~ unit_sample, scales = "free_x", 
             labeller = label_parsed) + 
  # facet_grid(Sample ~ unit_sample, scales = "free_x") + 
  geom_line(aes(size = Sample, alpha = Sample, col = Sample)) + 
  scale_size_manual(values = c(0.6, 0.3), 
                    labels = c(expression(pi[i]^B), expression(pi[i]^R))) + 
  scale_alpha_manual(values = c(1, 0.6), 
                     labels = c(expression(pi[i]^B), expression(pi[i]^R))) + 
  theme(strip.background = element_rect(fill="aliceblue"))

ggsave(paste0(wd, "Figures/low_overlap.png"), width = 8.5, height = 3, 
       dpi = 700, device = png)
