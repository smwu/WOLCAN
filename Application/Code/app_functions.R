#=================================================
# PROSPECT application functions
# Author: Stephanie Wu
# Date created: 2024/09/29
# Date updated: 2024/09/29
#=================================================

### Survey-weighted logistic regression
wtd_logreg <- function(res, data, formula_y) {
  svy_data <- data.frame(dbh_class = as.factor(res$estimates_adjust$c_all),
                         restrict_data, 
                         wts = res$data_vars$sampling_wt)
  svy_des <- survey::svydesign(ids = ~1, weights = ~wts, data = svy_data)
  log_reg <- survey::svyglm(formula = as.formula(formula_y), 
                            design = svy_des, 
                            family = stats::quasibinomial(link = "logit"))
  return(summary(log_reg))
}

# Display plots of results
res_plots <- function(res) {
  print(baysc::plot_pattern_profiles(res, item_labels = item_labels, 
                                     class_title = class_title, 
                                     categ_title = categ_title,
                                     categ_labels = categ_labels))
  print(baysc::plot_pattern_probs(res, item_labels = item_labels, 
                                  class_title = class_title, 
                                  categ_title = categ_title,
                                  categ_labels = categ_labels))
  print(baysc::plot_class_dist(res))
  # Traceplot for pi
  plot(res$estimates_adjust$pi_red[, 1], type = "l", ylim = c(0, 1))
  hist(res$data_vars$sampling_wt, breaks = 30)
}

# Function to run weighted logistic regression, incorporating variability from 
# estimated weights
# Inputs:
#   wts_draws: nxD matrix of the D draws from the weights posterior distribution
#   subset: String specifying whether to subset to those unaware ("unaware") to 
# measure impact of diet behavior pattern on disease, or to subset to those with 
# the disease ("disease") to measure impact of diagnosis on diet behavior. Must 
# be one of "unaware" or "disease".
#   condition: String specifying the condition to focus on. Must be one of "htn"
# (hypertension), "t2d" (type 2 diabetes), or "chol" (high cholesterol)
# Add variability of the weights
wtd_logreg_wolcan <- function(wts_draws, subset, condition, save_res = TRUE,
                              save_path) {
  # Check input arguments
  if (subset == "unaware") {
    filter_categ <- 3
  } else if (subset == "disease") {
    filter_categ <- 1
  } else {
    stop("Input argument 'subset' must be one of 'unaware' or 'disease'")
  }
  
  if (condition == "htn") {
    cond_categ <- "htn_categ"
    cond <- "hypertension"
  } else if (condition == "t2d") {
    cond_categ <- "t2d_categ"
    cond <- "diabetes"
  } else if (condition == "chol") {
    cond_categ <- "chol_categ"
    cond <- "cholesterol"
  } else {
    stop("Input argument 'condition' must be one of 'htn', 't2d', or 'chol'")
  }
  
  # Number of draws
  D <- ncol(wts_draws)  
  # Initialize lists
  adj_parms_draws <- list()
  sampled_parms_draws <- list()
  stan_fit_draws <- list()
  
  # Define brms formula
  model_formula <- as.formula(paste0(cond, " | weights(wts_d) ~ dbh_class ",
                                     "+ Age + Sex + Educ+ Inc_hh + Urban",
                                     " + Physical_activity + Smoking_status ",
                                     "+ Food_security")) 
  brms_mod <- brms::brmsformula(model_formula, center = TRUE)
  
  # For each draw, run the variance-adjusted survey-weighted regression and 
  # store results
  for (d in 1:D) {
    print("d: ", d)
    set.seed(d)
    svy_data_d <- svy_data %>%
      mutate(wts_d = wts_draws[, d])
    svy_data_subset_d <- svy_data_d %>%
      filter(cond_categ != filter_categ) %>% 
      drop_na(cond, Age, Sex, Educ, Inc_hh, Urban, Physical_activity, 
              Smoking_status, Food_security)
    svy_des_subset_d <- survey::svydesign(ids = ~1, weights = ~wts_d,
                                          data = svy_data_subset_d)
    print(weights(svy_des_subset_d)[1:10])
    fit_d <- csSampling::cs_sampling_brms(svydes = svy_des_subset_d, 
                                          brmsmod = brms_mod, 
                                          data = svy_data_subset_d, 
                                          family = bernoulli(link = "logit"),
                                          ctrl_stan = list(chains = 3, iter = 2000, 
                                                           warmup = 1000, thin = 5))
    sampled_parms_draws[[d]] <- fit_d$sampled_parms
    adj_parms_draws[[d]] <- fit_d$adjusted_parms
    stan_fit_draws[[d]] <- fit_d$stan_fit
  }
  
  all_adj_parms <- as.data.frame(do.call(rbind, adj_parms_draws))
  
  if (save_res) {
    wtd_logreg_res <- list(all_adj_parms = all_adj_parms, 
                           sampled_parms_draws = sampled_parms_draws,
                           stan_fit_draws = stan_fit_draws)
    save(wtd_logreg_res, file = paste0(save_path, "logreg.RData"))
  }
}