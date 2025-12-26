#===========================================================
# Summarize scenarios for predicting pseudo-weights
# Author: Stephanie Wu
# Date created: 2025/11/24
# Date updated: 2025/11/24
#===========================================================

rm(list = ls())

library(tidyverse)
library(gridExtra)
library(knitr)
library(kableExtra)
library(BART)  # BART

# Set directories
wd <- "~/Documents/GitHub/WOLCAN/"  # Working directory
wd <- "/n/netscratch/stephenson_lab/Lab/stephwu18/WOLCAN/"
data_dir <- "Data/"    # Data directory
res_dir <- "Results/"  # Results directory
code_dir <- "Summary_Code/"  # Model code directory
sum_dir <- "Summary_Results/"  # Summary results directory

#=============== Read in data ==================================================

# Data scenarios:
# High overlap: 0, 8, 13, 14, 24
# Low overlap: 10, 19, 20, 21, 25
scenarios <- c(0, 8, 13, 14, 24, 10, 19, 20, 21, 25)

model_scen <- list("none", "logistic_true", "logistic", 
                   "bart_500", "bart_1000", "bart_2000", 
                   "logistic_cov", "bart_1000_cov")

# Initialize results list for all models
# MAB = mean absolute bias, averaged across individuals in each sample
# SE and CI width are averaged across individuals in each sample
wts_res_all <- as.data.frame(matrix(NA, nrow = 1, ncol = 8))
colnames(wts_res_all) <- c("mab_weights", "mab_pi_B", "se_w_B_mean", 
                           "ci_width_w_B_mean", "se_pi_B_mean", 
                           "ci_width_pi_B_mean", "model", "scenario")

for (i in 1:length(scenarios)) {
  scenario <- scenarios[i]
  for (j in 1:length(model_scen)) {
    model <- model_scen[j]
    # Read in results
    load(paste0(wd, sum_dir, "scen_", scenario, "/model_", model, 
                "_wts_res.RData"))
    # Transpose 
    wts_res <- as.data.frame(t(wts_res))
    # Average results across independent samples
    wts_res_avg <- as.data.frame(t(colMeans(wts_res, na.rm = TRUE)))
    
    # Rename and add variables to stack with other results
    wts_res_avg <- wts_res_avg %>%
      rename(mab_weights = weights,
             mab_pi_B = pi_B)
    wts_res_avg$model <- model
    wts_res_avg$scenario <- scenario
    
    # Stack with other results
    wts_res_all <- rbind(wts_res_all, wts_res_avg)
  }
}

# Remove empty first row
wts_res_all <- wts_res_all[-1, ]


# Reshape data for tables
# Mean absolute bias
weights_res_mean <- wts_res_all %>%
  select(mab_weights, model, scenario) %>%
  rename(mean_abs_bias = mab_weights)
weights_res_mean <- weights_res_mean %>%
  pivot_wider(names_from = model, values_from = mean_abs_bias)
weights_res_mean <- weights_res_mean %>%
  mutate(scenario = unlist(scenario))

# Mean SE/posterior SD
se_w_B_mean <- wts_res_all %>%
  select(se_w_B_mean, model, scenario)
se_w_B_mean <- se_w_B_mean %>%
  pivot_wider(names_from = model, values_from = se_w_B_mean)
se_w_B_mean <- se_w_B_mean %>%
  mutate(scenario = unlist(scenario))

pi_res_mean <- wts_res_all %>%
  select(mab_pi_B, model, scenario) %>%
  rename(mean_abs_bias = mab_pi_B)
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

# Application-based variable weights: scenario 24 and 25
tb_5 <- as.data.frame(weights_res_mean) %>% 
  filter(scenario %in% c(24, 25)) %>%
  select(-scenario)
rownames(tb_5) <- c("High Overlap", "Low Overlap")
tb_5 %>%
  kbl(caption = "Variable Weights", digits = 3) %>%
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

# Sample size n_B 1% n_R 5%: scenario 13 and 20
tb_3 <- as.data.frame(pi_res_mean) %>% 
  filter(scenario %in% c(13, 20)) %>%
  select(-scenario)
rownames(tb_3) <- c("High Overlap", "Low Overlap")
tb_3 %>%
  kbl(caption = "n_B 1%, n_R 5%", digits = 3) %>%
  kable_classic(full_width = F)

# Sample size n_B 5% n_R 1%: scenario 14 and 21
tb_4 <- as.data.frame(pi_res_mean) %>% 
  filter(scenario %in% c(14, 21)) %>%
  select(-scenario)
rownames(tb_4) <- c("High Overlap", "Low Overlap")
tb_4 %>%
  kbl(caption = "n_B 5%, n_R 1%", digits = 3) %>%
  kable_classic(full_width = F)

# Application-based variable weights
tb_5 <- as.data.frame(pi_res_mean) %>% 
  filter(scenario %in% c(24, 25)) %>%
  select(-scenario)
rownames(tb_5) <- c("High Overlap", "Low Overlap")
tb_5 %>%
  kbl(caption = "Application-Based", digits = 3) %>%
  kable_classic(full_width = F)


# Rearrange and print full weights df
# Mean absolute bias
displ <- weights_res_mean %>%
  slice(match(c(0, 10, 14, 21, 13, 20, 8, 19, 24, 25), scenario)) %>%
  mutate(Overlap = rep(c("High", "Low"), times = 5)) %>%
  select(Overlap, none, logistic_true, logistic, logistic_cov, 
         bart_500:bart_2000, bart_1000_cov)
displ <- as.data.frame(displ)
colnames(displ) <- c("Overlap", "No Model", "LogRegTrue", "LogReg", "LogRegMiss", 
                     "BART500", "BART1000", "BART2000", "BART1000Miss")
displ$`Sample Size` <- rep(c("$n_B 5\\%, n_R 5\\%$", "$n_B 5\\%, n_R 1\\%$", 
                             "$n_B 1\\%, n_R 5\\%$", "$n_B 1\\%, n_R 1\\%$",
                             "Application-Based"), each = 2)
displ <- displ %>% select(`Sample Size`, Overlap:`BART1000Miss`)
displ %>%
  kbl(digits = 3, format = "html", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE)


# Uncertainty
displ_se <- se_w_B_mean %>%
  slice(match(c(0, 10, 14, 21, 13, 20, 8, 19, 24, 25), scenario)) %>%
  mutate(Overlap = rep(c("High", "Low"), times = 5)) %>%
  select(Overlap, none, logistic_true, logistic, logistic_cov, 
         bart_500:bart_2000, bart_1000_cov)
displ_se <- as.data.frame(displ_se)
colnames(displ_se) <- c("Overlap", "No Model", "LogRegTrue", "LogReg", "LogRegMiss", 
                     "BART500", "BART1000", "BART2000", "BART1000Miss")
displ_se$`Sample Size` <- rep(c("$n_B 5\\%, n_R 5\\%$", "$n_B 5\\%, n_R 1\\%$", 
                             "$n_B 1\\%, n_R 5\\%$", "$n_B 1\\%, n_R 1\\%$",
                             "Application-Based"), each = 2)
displ_se <- displ_se %>% select(`Sample Size`, Overlap:`BART1000Miss`)
displ_se %>%
  kbl(digits = 3, format = "html", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE)

### Weights and SE in same table
displ_comb <- displ
# Combine mean absolute bias and SE into "mean (SE)" with rounding to nearest tenth
method_cols <- c("No Model", "LogRegTrue", "LogReg", "LogRegMiss",
                 "BART500", "BART1000", "BART2000", "BART1000Miss")

for (col in method_cols) {
  displ_comb[[col]] <- sprintf("%.1f (%.1f)",
                          round(displ[[col]], 1),
                          round(displ_se[[col]], 1))
}

# Final table
displ_comb %>%
  kbl(format = "html", booktabs = TRUE) %>%
  kable_classic(full_width = FALSE)

#================ Plot comparison for high and low overlap
plot_df <- displ 
plot_df %>%
  pivot_longer(cols = `No Model`:BART1000Miss, names_to = "Model", 
               values_to = "mean_abs_bias") %>%
  mutate(Model = factor(Model, 
                        levels = c("No Model", "LogRegTrue", "LogReg", 
                                   "BART500",  "BART1000", "BART2000", 
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
                              "$n_B 1\\%, n_R 1\\%$",
                              "Application-Based"),
                   labels = c(expression(atop(n[B]*" 5%", n[R]*" 5%")),
                              expression(atop(n[B]*" 5%", n[R]*" 1%")),
                              expression(atop(n[B]*" 1%", n[R]*" 5%")),
                              expression(atop(n[B]*" 1%", n[R]*" 1%")),
                              "App-Based"))
ggsave(paste0(wd, "Figures/weights_sims_new.png"), width = 8.5, height = 5, 
       dpi = 700, device = png)

### Previous plot
### Plot without "none" model
plot_df_old <- displ %>% 
  select(-c(`LogRegTrue`)) %>%
  filter(`Sample Size` != "Application-Based")
plot_df_old %>%
  pivot_longer(cols = `No Model`:BART1000Miss, names_to = "Model", 
               values_to = "mean_abs_bias") %>%
  mutate(Model = factor(Model, 
                        levels = c("No Model","LogReg", 
                                   "BART500",  "BART1000", "BART2000", 
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

### Plot without "none" model
plot_df_excl <- displ %>% select(-c(`No Model`))
plot_df_excl %>%
  pivot_longer(cols = `LogRegTrue`:BART1000Miss, names_to = "Model", 
               values_to = "mean_abs_bias") %>%
  mutate(Model = factor(Model, 
                        levels = c("LogRegTrue", "LogReg", 
                                   "BART500",  "BART1000", "BART2000", 
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
                              "$n_B 1\\%, n_R 1\\%$",
                              "Application-Based"),
                   labels = c(expression(atop(n[B]*" 5%", n[R]*" 5%")),
                              expression(atop(n[B]*" 5%", n[R]*" 1%")),
                              expression(atop(n[B]*" 1%", n[R]*" 5%")),
                              expression(atop(n[B]*" 1%", n[R]*" 1%")),
                              "App-Based"))
ggsave(paste0(wd, "Figures/weights_sims_excl_no_model.png"), width = 8.5, height = 5, 
       dpi = 700, device = png)


### Facet by scenario
plot_df %>%
  pivot_longer(
    cols = `No Model`:BART1000Miss,
    names_to = "Model", 
    values_to = "mean_abs_bias"
  ) %>%
  mutate(
    Model = factor(
      Model, 
      levels = c("No Model", "LogRegTrue", "LogReg", 
                 "BART500", "BART1000", "BART2000", 
                 "LogRegMiss", "BART1000Miss")
    ),
    scenario = factor(scenario)  # optional but usually helpful for faceting
  ) %>%
  ggplot(aes(x = `Sample Size`, y = mean_abs_bias, fill = Model)) + 
  theme_bw() +
  geom_bar(stat = "identity", position = position_dodge()) + 
  
  # facet by scenario (rows) and Overlap (columns)
  facet_grid(
    scenario ~ Overlap,
    labeller = labeller(
      Overlap = c("High" = "High Overlap",
                  "Low"  = "Low Overlap")
    )
  ) + 
  scale_fill_brewer(palette = "Set2") + 
  ylab("Mean Absolute Bias for Pseudo-Weights") + 
  labs(fill = "Model") + 
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "aliceblue")
  ) +
  guides(fill = guide_legend(nrow = 1)) + 
  scale_x_discrete(
    name = "Sample Size",
    limits = c("$n_B 5\\%, n_R 5\\%$", 
               "$n_B 5\\%, n_R 1\\%$",
               "$n_B 1\\%, n_R 5\\%$",
               "$n_B 1\\%, n_R 1\\%$",
               "Application-Based"),
    labels = c(
      expression(atop(n[B]*" 5%", n[R]*" 5%")),
      expression(atop(n[B]*" 5%", n[R]*" 1%")),
      expression(atop(n[B]*" 1%", n[R]*" 5%")),
      expression(atop(n[B]*" 1%", n[R]*" 1%")),
      "App-Based"
    )
  ) + 
  scale_y_continuous(trans='log10')

ggsave(paste0(wd, "Figures/weights_sims_log10.png"), width = 8.5, height = 5, 
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


#================ Plot weights density for single realization ==================
# Pick one sample
samp_i <- 1

# Data scenarios:
# High overlap: 0, 8, 13, 14, 24
# Low overlap: 10, 19, 20, 21, 25
scenarios <- c(0, 8, 13, 14, 24, 10, 19, 20, 21, 25)

scen_names <- c(expression(n[B]*" 5%, "*n[R]*" 5%, High Overlap"),
                expression(n[B]*" 1%, "*n[R]*" 1%, High Overlap"),
                expression(n[B]*" 1%, "*n[R]*" 5%, High Overlap"),
                expression(n[B]*" 5%, "*n[R]*" 1%, High Overlap"),
                "App-Based, High Overlap",
                expression(n[B]*" 5%, "*n[R]*" 5%, Low Overlap"),
                expression(n[B]*" 1%, "*n[R]*" 1%, Low Overlap"),
                expression(n[B]*" 1%, "*n[R]*" 5%, Low Overlap"),
                expression(n[B]*" 5%, "*n[R]*" 1%, Low Overlap"),
                "App-Based, Low Overlap")

# Open PNG device
png(paste0(wd, "Figures/weights_sim_densities.png"),
    width = 10, height = 6, units = "in", res = 300)  # 300 dpi for nice quality

# Open a multi-panel layout: 2 rows x 5 columns
par(mfrow = c(2, 5), mar = c(5, 4, 2, 1))  # adjust margins if needed

for (i in seq_along(scenarios)) {
  scen_i <- scenarios[i]
  
  # NOTE: use scen_i here, not 'scenario'
  load(paste0(wd, data_dir, "scen_", scen_i, "/sim_samp_", samp_i, "_B_wolcan.RData"))
  
  dens <- density(1 / sim_samp_B$true_pi_B)
  plot(
    dens,
    main = scen_names[i],
    xlab = "NPS Weight",
    ylab = "Density",
    cex.main = 1.1
  )
}

dev.off()  # important: closes the PNG device and writes the file

# (Optional) reset layout afterward
par(mfrow = c(1, 1))

