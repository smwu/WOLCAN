#==================================
## Summarize and Plot Model Results
## Programmer: SM Wu   
## Last Updated: 2023/05/20
#==================================

rm(list = ls())

library(baysc)
library(abind)  # binding arrays
library(flextable)
library(dplyr)
library(bayesplot)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(knitr)
library(kableExtra)
library(gt)
library(parallel)
library(stringr)
library(Hmisc)  # plotting
library(ggbeeswarm)  # beeswarm plot

# Set directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/"
data_dir <- "Data/"    # Data directory
res_dir <- "Results/"  # Results directory
code_dir <- "Summary_Code/"  # Model code directory
sum_dir <- "Summary_Results/"  # Summary results directory

# Source summary functions
source(paste0(wd, code_dir, "summary_functions.R"))

#================ Summarize and save results ===================================
# scenario <- 9
scenarios <- c(0, 7:10)
scenarios <- c(0, 10, 8, 19, 9, 18, 7)
scenarios <- c(6, 15)
for (i in 1:length(scenarios)) {
  scenario <- scenarios[i]
  samp_i_seq <- 1:100
  # Define path to save results
  save_path <- paste0(wd, sum_dir, "scen_", scenario, "/")
  
  # Create scenario results folder if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(file.path(save_path))
  }
  
  # Get metrics for all samples
  save_scen_metrics(scenario = scenario, samp_i_seq = samp_i_seq, WOLCAN = TRUE, 
                    WOLCA = TRUE, save_path = save_path, wd = wd, 
                    data_dir = data_dir, res_dir = res_dir, subset = FALSE, 
                    dist_type = "mean_abs", parallel = FALSE)
}


# bias comparison
metrics_all$metrics_wolca$theta_bias
metrics_all$metrics_wolcan$theta_bias
metrics_all$metrics_wolca$theta_mode_bias
metrics_all$metrics_wolcan$theta_mode_bias
# samp_i_seq <- 1:10
# save_path <- paste0(wd, sum_dir, "scen_", scenario, "/", "top10_")
# save_scen_metrics(scenario = scenario, samp_i_seq = samp_i_seq, WOLCAN = TRUE, 
#                   WOLCA = TRUE, save_path = save_path, wd = wd, 
#                   data_dir = data_dir, res_dir = res_dir, subset = FALSE, 
#                   dist_type = "mean_abs", parallel = TRUE)

#============== Create Tables ==================================================

scenarios <- c(2, 3, 4, 5)
scenarios <- c(0, 7, 8, 10, 9)
scenarios <- c(6, 15)
# scen_names <- c("Baseline: n = 2000 (5%)")
scen_names <- paste0("Scen ", scenarios)
# scen_names <- c("WS all", "WS mean", "MI no var adj", "no wts adj")
save_names <- rep("metrics_scen", length(scenarios))
save_paths <- paste0(wd, sum_dir, "scen_", scenarios, "/")
create_app_tables_wolcan(save_paths = save_paths, scenarios = scenarios, 
                         scen_names = scen_names, overall_name = "Scenario",
                         format = "html", WOLCA = TRUE, WOLCAN = TRUE)

# Fill in scenario 6 and 15 unweighted with baseline
load(paste0(wd, sum_dir, "scen_0/summary.RData"))
scen_0_metrics <- metrics_all
load(paste0(wd, sum_dir, "scen_6/summary.RData"))
metrics_all$metrics_wolca <- scen_0_metrics$metrics_wolca
save(metrics_all, file = paste0(wd, sum_dir, "scen_6/summary.RData"))
load(paste0(wd, sum_dir, "scen_15/summary.RData"))
metrics_all$metrics_wolca <- scen_0_metrics$metrics_wolca
save(metrics_all, file = paste0(wd, sum_dir, "scen_15/summary.RData"))

# Table of all scenarios together                         
scenarios <- c(0, 10, 8, 19, 9, 18, 7, 6, 15)
scen_names <- c("Sample size 5% high overlap (baseline)", 
                "Sample size 5% low overlap", "Sample size 1% high overlap", 
                "Sample size 1% low overlap", 
                "Sample size PROSPECT high overlap", 
                "Sample size PROSPECT low overlap",
                "Baseline with non-disjoint patterns", 
                "Baseline with D = 10",
                "Baseline with missing selection covariates")
save_names <- rep("metrics_scen", length(scenarios))
save_paths <- paste0(wd, sum_dir, "scen_", scenarios, "/")
create_app_tables_wolcan(save_paths = save_paths, scenarios = scenarios, 
                         scen_names = scen_names, overall_name = "Scenario",
                         format = "html", WOLCA = TRUE, WOLCAN = TRUE)

# Table of all scenarios using modal theta for bias
create_app_tables_wolcan(save_paths = save_paths, scenarios = scenarios, 
                         scen_names = scen_names, overall_name = "Scenario",
                         format = "html", WOLCA = TRUE, WOLCAN = TRUE, 
                         modal = TRUE)

#================= Grouped param boxplots across scenarios =====================
# Use modal theta for theta bias
# Order of rows: For each scenario, L iterations of wolca and wolcan
L <- 100
scenarios <- c(0, 10, 8, 19, 9, 18)
num_scen <- length(scenarios)
num_rows <- num_scen*2*L
plot_df <- as.data.frame(matrix(NA, nrow = num_rows, ncol = 6))
colnames(plot_df) <- c("pi_bias", "pi_cov", "pi_var", "theta_bias", "theta_cov", "theta_var")
plot_df$scen <- rep(scenarios, each = 2*L)
plot_df$`Sample Size` <- rep(c("5%", "1%", "PROSPECT"), each = 2*2*L)
plot_df$Overlap <- rep(c("High", "Low"), each = 2*L, times = num_scen/2)
plot_df$Model <- rep(c("Unweighted", "WOLCAN"), each = L, times = 2)

counter <- 0  # help propagate the correct row
for (i in 1:num_scen) {
  load(paste0(wd, sum_dir, "scen_", scenarios[i], "/summary.RData"))
  print(paste0(metrics_all$metrics_wolca$pi_bias))  # sanity check
  print(counter)
  
  plot_df[(counter + 1):(counter + 2*L), "pi_bias"] <- 
    c(metrics_all$metrics_wolca$pi_dist, 
      metrics_all$metrics_wolcan$pi_dist)
  plot_df[(counter + 1):(counter + 2*L), "pi_cov"] <- 
    c(metrics_all$metrics_wolca$pi_cover_all, 
      metrics_all$metrics_wolcan$pi_cover_all)
  plot_df[(counter + 1):(counter + 2*L), "pi_var"] <- 
    c(metrics_all$metrics_wolca$pi_var_all, 
      metrics_all$metrics_wolcan$pi_var_all)
  plot_df[(counter + 1):(counter + 2*L), "theta_bias"] <- 
    c(metrics_all$metrics_wolca$theta_mode_dist, 
      metrics_all$metrics_wolcan$theta_mode_dist)
  plot_df[(counter + 1):(counter + 2*L), "theta_cov"] <- 
    c(metrics_all$metrics_wolca$theta_cover_all, 
      metrics_all$metrics_wolcan$theta_cover_all)
  plot_df[(counter + 1):(counter + 2*L), "theta_var"] <- 
    c(metrics_all$metrics_wolca$theta_var_all, 
      metrics_all$metrics_wolcan$theta_var_all)
  counter <- counter + 2*L  # Increment counter by 2*L for next scenario
} 

# Sanity check of values
mean(plot_df[plot_df$scen == 0 & plot_df$Model == "Unweighted", "pi_bias"])
mean(plot_df[plot_df$scen == 9 & plot_df$Model == "Unweighted", "pi_bias"])


plot_df_format <- plot_df %>% 
  pivot_longer(cols = pi_bias:theta_var, names_to = "metric_all", values_to = "Value")
temp <- nrow(plot_df_format)
plot_df_format <- plot_df_format %>%
  mutate(Parameter = str_split_i(metric_all, "_", i = 1), 
         Metric = str_split_i(metric_all, "_", i = 2))
dodge <- position_dodge(width = 0.7)
# plot_df_format %>% 
#   filter(Overlap == "High") %>%
#   ggplot(aes(x =`Sample Size`, y = Value, col = Model)) + 
#   geom_violin(width = 0.8, trim = FALSE, position = dodge) +
#   stat_summary(fun=mean, geom="point", size=1, position = dodge) + 
#   #stat_summary(fun.data=mean_sdl, geom="pointrange", position = dodge, size = 0.5) + 
#   #geom_boxplot(width = 0.1, position = dodge, alpha = 0.2) + 
#   theme_bw() + 
#   facet_grid(Parameter ~ Metric, scales = "free")

plot_df_format2 <- plot_df_format %>%
  mutate(Scenario = paste(`Sample Size`, Overlap),
         Scenario = factor(Scenario, levels = c("5% High", "5% Low", 
                                                 "1% High", "1% Low", 
                                                 "PROSPECT High", "PROSPECT Low"),
                           labels = c("5%\nHigh", "5%\nLow", 
                                      "1%\nHigh", "1%\nLow", 
                                      "PROSP\nHigh", "PROSP\nLow")),
         Parameter = factor(Parameter, levels = c("pi", "theta"), 
                            labels = c("Pi", "Theta")),
         Metric = factor(Metric, levels = c("bias", "var", "cov"),
                         labels = c("Mean Absolute Bias", 
                                    "Posterior Interval Width", "Coverage")))
plot_df_format2 %>% 
  ggplot(aes(x = Scenario, y = Value, fill = Model, col = Model)) + 
  geom_violin(width = 1, trim = TRUE, position = dodge, alpha = 0.1,
              linewidth = 0.1) +
  #draw_quantiles = c(0.25, 0.5, 0.75)
  #geom_beeswarm(corral = "gutter", dodge.width = 0.9, cex=0.1) + 
  geom_quasirandom(dodge.width = 0.7, cex=0.15, alpha=0.5, method = "tukeyDense",
                   width = 0.25) +
  scale_fill_brewer(palette = "Set2") + 
  scale_color_brewer(palette = "Set2") + 
  #geom_boxplot(width = 0.8, position = dodge) + 
  stat_summary(fun=mean, geom="point", size=1, position = dodge, col = "black") + 
  theme_bw() + 
  facet_grid(Metric ~ Parameter, scales = "free_y") + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="aliceblue"))



#================= Plot theta patterns =========================================

# Baseline
scenario <- 0
save_path <- paste0(wd, sum_dir, "scen_", scenario, "/")
p_theta <- plot_theta_patterns_wolcan(wd = wd, data_dir = data_dir, 
                                            scenario = scenario, 
                                      save_path = save_path) 
annotate_figure(p_theta, 
                top = text_grob("Modal pattern per latent class, averaged across simulations"))


### Investigation of WOLCA extra classes
# Load summary
load(paste0(save_path, "summary.RData"))
summ_wolca <- metrics_all$metrics_wolca
# Obtain WOLCA iterations with extra class
samp_extra <- which(summ_wolca$K_all != 3)
length(samp_extra) # 37 iterations
# Get theta_mode for all iterations with extra class
theta_mode_extra <- array(NA, dim = c(length(samp_extra), dim(summ_wolca$theta_mode)))
for (i in 1:length(samp_extra)) {
  temp <- get_metrics_wolcan(wd = wd, data_dir = data_dir, 
                             res_dir = res_dir, scenario = scenario, 
                             model = "wolca", samp_i_seq = samp_extra[i],
                             subset = subset, dist_type = dist_type,
                             parallel = FALSE, 
                             save_path = save_path)
  theta_mode_extra[i, , ] <- temp$theta_mode
}
table(theta_mode_extra[, 29, 4])
table(theta_mode_extra[, 30, 4])
table(theta_mode_extra[, 29, 4], theta_mode_extra[, 30, 4])
# Most common is 29=3, 30=4 (e.g., samp_i = 2)

# Plot comparison with WOLCA samp_i = 2
# Process WOLCA samp_i = 2
load(paste0(wd, sum_dir, "scen_0/samp_2_wolca.RData"))
wolca_mode <- summ_i$theta_mode
# Load summary
load(paste0(save_path, "summary.RData"))
# Load true population data
pop_data_path <- paste0(wd, data_dir, "scen_", scenario, "/sim_pop_wolcan.RData")
load(pop_data_path)
# Obtain true observed population parameters
true_params <- get_true_params_wolcan(sim_pop = sim_pop)
p_true <- theta_mode_plot_wolcan(sim_pop$true_global_patterns, "True Classes") + 
  guides(fill = guide_legend(reverse = FALSE)) +
  labs(fill = "Modal θ Level")
p_WOLCAN <- theta_mode_plot_wolcan(metrics_all$metrics_wolcan$theta_mode, "WOLCAN Classes")
p_unwt <- theta_mode_plot_wolcan(wolca_mode, "Unweighted Classes")
#p_unwt <- theta_mode_plot(metrics_all$metrics_wolca$th, "Unweighted Classes")
ggarrange(p_true, 
          p_WOLCAN + theme(axis.title.y = element_blank()), 
          p_unwt + theme(axis.title.y = element_blank()), 
          nrow = 1, common.legend = TRUE, legend = "right", widths = c(1, 1, 1.2))

#================== Plot pi ====================================================

plot_pi_patterns_wolcan(wd = wd, data_dir = data_dir, scenario = scenario, 
                        samp_i_seq = samp_i_seq, save_path = save_path,
                        y_lim = c(0,1))

#================== Run simulations and summary for outcome model ==============
library(survey)
pop$Y <- Y_data

test <- glm(as.formula(Y ~ c_all + A1 + A2 + c_all:A1), data = pop, 
            family = binomial(link = "logit"))
round(summary(test)$coefficients[, 1], 3)
xi_vec_y

samp_data <- data.frame(c_all = sim_samp_B$c_all, Y = sim_samp_B$Y_data, 
                        sim_samp_B$covs, wts = 1/sim_samp_B$true_pi_B, 
                        clus = 1:nrow(sim_samp_B$covs))
svy_des <- survey::svydesign(ids = ~clus, weights = ~wts, data = samp_data)
svy_test <- survey::svyglm(formula = as.formula(Y ~ c_all + A1 + A2 + c_all:A1),
                           design = svy_des, 
                           family = stats::binomial(link = "logit"))
round(summary(svy_test)$coefficients[, 1], 3)

# True parameters
xi_vec_y

scenarios <- c(0, 6:12)
xi_est <- as.data.frame(matrix(NA, nrow = length(scenarios), 
                               ncol = length(xi_vec_y)))
for (scenario in scenarios) {
  # Load data
  
  # Load model output
  
  # Run survey-weighted regression
  
  # Save results
  
  ### Summary metrics
  # Get best order
  
  # Calculate bias
  
  # Calculate variance
  
  # Calculate coverage
}
