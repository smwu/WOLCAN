#==================================
## Summarize and Plot Model Results
## Programmer: SM Wu   
## Last Updated: 2023/05/20
#==================================

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
scenario <- 1
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

# samp_i_seq <- 1:10
# save_path <- paste0(wd, sum_dir, "scen_", scenario, "/", "top10_")
# save_scen_metrics(scenario = scenario, samp_i_seq = samp_i_seq, WOLCAN = TRUE, 
#                   WOLCA = TRUE, save_path = save_path, wd = wd, 
#                   data_dir = data_dir, res_dir = res_dir, subset = FALSE, 
#                   dist_type = "mean_abs", parallel = TRUE)

#============== Create Tables ==================================================

scenarios <- c(0, 1)
# scen_names <- c("Baseline: n = 2000 (5%)")
scen_names <- paste0("Scen ", scenarios)
save_names <- rep("metrics_scen", 3)
save_paths <- paste0(wd, sum_dir, "scen_", scenarios, "/")
create_app_tables_wolcan(save_paths = save_paths, scenarios = scenarios, 
                         scen_names = scen_names, overall_name = "Scenario",
                         format = "html")


#================= Plot theta patterns =========================================

# Baseline
p_theta <- plot_theta_patterns_wolcan(wd = wd, data_dir = data_dir, 
                                            scenario = scenario, 
                                      save_path = save_path) 
annotate_figure(p_theta, 
                top = text_grob("Modal pattern per latent class, averaged across simulations"))


#================== Plot pi ====================================================

plot_pi_patterns_wolcan(wd = wd, data_dir = data_dir, scenario = scenario, 
                        samp_i_seq = samp_i_seq, save_path = save_path,
                        y_lim = c(0,1))


