#=================================================
# PROSPECT application functions
# Author: Stephanie Wu
# Date created: 2024/09/29
# Date updated: 2024/09/29
#=================================================

### Survey-weighted logistic regression
wtd_logreg <- function(res, data, formula_y, filter_inds) {
  svy_data <- data.frame(dbh_class = as.factor(res$estimates_adjust$c_all),
                         data, 
                         wts = res$data_vars$sampling_wt)
  svy_data <- svy_data[filter_inds, ]
  svy_des <- survey::svydesign(ids = ~1, weights = ~wts, data = svy_data)
  log_reg <- survey::svyglm(formula = as.formula(formula_y), 
                            design = svy_des, 
                            family = stats::quasibinomial(link = "logit"))
  return(log_reg)
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
#   model_covs: String specifying whether the model will be "marginal", with 
# "core" covariates (default), or with "full" covariate adjustment
# Add variability of the weights
wtd_logreg_wolcan <- function(wts_draws, subset, condition, save_res = TRUE,
                              save_path, model_covs = "core", svy_data) {
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
    cond <- "high_cholesterol"
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
  if (model_covs == "marginal") {
    add_covs <- ""
  } else if (model_covs == "core") {
    if (subset == "unaware") {
      add_covs <- paste0("+ Age + Sex + Educ + Inc_hh + Ethnicity + Urban ", 
                         "+ Physical_activity + Smoking_status + Drinking_status", 
                         "+ hypertension + high_cholesterol")
    }
  } else if (model_covs == "full") {
    add_covs <- paste0("+ Age + Sex + Educ + Inc_hh + Ethnicity + Urban ",
                       "+ Physical_activity + Smoking_status + Drinking_status", 
                       "+ hypertension + high_cholesterol + Food_security",
                       "+ WIC_SNAP + Social_support + Perceived_stress",
                       "+ Depression + Anxiety")
  } else if (model_covs == "core_int") {
    if (subset == "unaware") {
      add_covs <- paste0("* (Age + Sex + Inc_hh + Urban + Physical_activity)",
                         "+ Educ + Smoking_status + Drinking_status",
                         "+ hypertension + high_cholesterol")
    }
  } else if (model_covs == "full_int") {
    add_covs <- paste0("* (Age + Sex + Educ + Inc_hh + Urban + Physical_activity",
                       "+ Smoking_status + Food_security",
                       "+ Ethnicity + Drinking_status + WIC_SNAP + Social_support",
                       "+ Perceived_stress + Depression + Anxiety)")
  } else {
    stop("Invalid model_covs input specified.")
  }
  if (subset == "unaware") {
    model_formula <- as.formula(paste0(cond, " | weights(wts_d) ~ dbh_class ",
                                       add_covs)) 
  } else if (subset == "disease") {
    model_formula <- as.formula(paste0("dbh_class | weights(wts_d) ~ ", cond_categ,
                                       add_covs)) 
  } else {
    stop("Invalid subset input specified.")
  }
  
  brms_mod <- brms::brmsformula(model_formula, center = TRUE)
  
  # For each draw, run the variance-adjusted survey-weighted regression and 
  # store results
  for (d in 1:D) {
    print(paste0("d: ", d))
    set.seed(d)
    svy_data_d <- svy_data %>%
      mutate(wts_d = wts_draws[, d]) 
    svy_data_subset_d <- svy_data_d %>%
      filter(!!rlang::sym(cond_categ) != filter_categ)
    # Remove NAs
    svy_data_subset_d <- svy_data_subset_d %>%
      drop_na(hypertension, diabetes, high_cholesterol, 
              Age, Sex, Educ, Inc_hh, Ethnicity, Urban, Physical_activity, 
              Smoking_status, Drinking_status, Food_security, WIC_SNAP, 
              Social_support, Perceived_stress, Depression, Anxiety)
    # Renormalize weights
    svy_data_subset_d <- svy_data_subset_d %>%
      mutate(wts_d = wts_d * nrow(svy_data_subset_d) / sum(wts_d))
      
    # # Remove NAs
    # if (model_covs == "marginal") {
    #   svy_data_subset_d <- svy_data_subset_d %>%
    #     drop_na(!!rlang::sym(cond), !!rlang::sym(cond_categ))
    # } else if (model_covs %in% c("core", "core_int")) {
    #   svy_data_subset_d <- svy_data_subset_d %>%
    #     drop_na(!!rlang::sym(cond), !!rlang::sym(cond_categ), 
    #             Age, Sex, Educ, Inc_hh, Urban, Physical_activity, 
    #             Smoking_status, Food_security)
    # } else if (model_covs %in% c("full", "full_int")) {
    #   svy_data_subset_d <- svy_data_subset_d %>%
    #     drop_na(!!rlang::sym(cond), !!rlang::sym(cond_categ), 
    #             Age, Sex, Educ, Inc_hh, Urban, Physical_activity, 
    #             Smoking_status, Food_security, Ethnicity, Drinking_status, 
    #             WIC_SNAP, Social_support, Perceived_stress, Depression, Anxiety)
    # } else {
    #   stop("Invalid model_covs input specified.")
    # }
      
    svy_des_subset_d <- survey::svydesign(ids = ~1, weights = ~wts_d,
                                          data = svy_data_subset_d)
    # print(weights(svy_des_subset_d)[1:10])
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
  
  # Output
  wtd_logreg_res <- list(all_adj_parms = all_adj_parms, 
                         sampled_parms_draws = sampled_parms_draws,
                         stan_fit_draws = stan_fit_draws,
                         brms_mod = brms_mod, data = svy_data_subset_d)
  
  if (save_res) {
    save(wtd_logreg_res, file = paste0(save_path, "logreg.RData"))
  }
  
  return(wtd_logreg_res)
}


# fit_obj: Fit object resulting from a call to `wtd_logreg_wolcan()`
summarize_parms <- function(fit_obj, quant_lb = 0.025, quant_ub = 0.975, 
                            round_digits = 3, parm_names = TRUE) {
  all_adj_parms <- fit_obj$all_adj_parms
  summary_adj_parms <- t(apply(all_adj_parms, 2, function(x) c(mean(x, na.rm = TRUE), 
                                                               quantile(x, c(0.025, 0.975), na.rm = TRUE),
                                                               mean(x > 0, na.rm = TRUE),
                                                               mean(x < 0, na.rm = TRUE))))
  summary_adj_parms <- round(summary_adj_parms, round_digits)
  colnames(summary_adj_parms) <- c("Mean", "2.5%", "97.5%", "P(xi>0)", "P(xi<0)")
  if (parm_names) {
    mod_mat <- brms::make_standata(formula = fit_obj$brms_mod$formula, 
                                   data = fit_obj$data, family = "bernoulli")
    parm_names <- colnames(mod_mat$X)
    parm_names_cs <- parm_names[-1]  # Remove intercept
    rownames(summary_adj_parms)[1:length(parm_names_cs)] <- parm_names_cs
  }
  
  return(summary_adj_parms)
}


# Table of estimated weights by covariates
get_mean_props_wolcan <- function(var_vec, wts, data, digits = 1) {
  
  out_all <- list()
  
  for (i in 1:length(var_vec)) {
    var_name <- var_vec[i]
    x <- data[, var_name]
    
    if (is.factor(x)) {
      num_levels <- length(levels(x))
      lev_lab <- levels(x)
    } else {
      num_levels <- 1
      lev_lab <- ""
    }
    
    out_df <- as.data.frame(matrix(NA, nrow = num_levels, ncol = 3))
    colnames(out_df) <- c("Variable", "Level", "Mean")
    out_df[, 1] <- rep(var_name, times = num_levels)
    out_df[, 2] <- lev_lab
    
    # Calculate the proportion
    svy_des <- survey::svydesign(ids = ~1, weights = ~wts, data = data)
    prop_temp <- as.data.frame(
      survey::svymean(stats::as.formula(paste0("~", var_name)), 
                      svy_des, na.rm = TRUE))
    
    if (is.factor(x)) {
      prop <- round(prop_temp[, 1] * 100, digits = digits)
      prop <- paste0(prop, "%")
    } else {
      prop <- round(prop_temp[, 1], digits = digits)
    }
    out_df[, 3] <- prop
    
    out_all[[i]] <- out_df
  }
  # Row-bind all variables
  out_all_df <- do.call("rbind", out_all)
  return(out_all_df)
}


# Function to plot the pattern probabilities horizontally
plot_pattern_probs_wolcan <- function (res, item_labels = NULL, 
                                       categ_labels = NULL, categ_title = "Risk Level", 
                                       class_labels = NULL, 
                                       class_title = "Dietary Behavior Pattern", 
                                       x_title = "Risk Level Probability", 
                                       num_rows = 12,
                                       ...) {
  if (is(res, "wolcan")) {
    est_item_probs <- res$estimates_adjust$theta_med
  } else if (is(res, "wolca")) {
    est_item_probs <- res$estimates$theta_med
  }
  
  K <- dim(est_item_probs)[2]
  # item_labels <- 1:res$data_vars$J
  class_labels <- 1:K
  # categ_labels <- 1:res$data_vars$R
  
  dimnames(est_item_probs)[[1]] <- item_labels
  dimnames(est_item_probs)[[2]] <- class_labels
  theta_plot <- data.frame(expand.grid(lapply(dim(est_item_probs), seq_len)), 
                           value = as.vector(est_item_probs))
  Item <- Class <- Probability <- Level <- NULL
  colnames(theta_plot) <- c("Item", "Class", "Level", "Probability")
  theta_plot %>% 
    ggplot2::ggplot(ggplot2::aes(y = factor(Class, labels = class_labels), 
                                 x = Probability, 
                                 fill = factor(Level, levels = c(3, 2, 1)))) + 
    ggplot2::geom_bar(stat = "identity", position = "stack") + 
    # ggplot2::facet_grid(factor(Item, labels = item_labels) ~ .) + 
    ggplot2::facet_wrap(. ~ factor(Item, labels = item_labels), nrow = num_rows, 
                        dir = "v") +
    ggplot2::scale_fill_brewer(type = "seq", palette = "RdYlBu", direction = 1, 
                               name = categ_title, labels = c("High", "Med", "Low")) + 
    ggplot2::guides(fill = guide_legend(reverse = TRUE)) + 
    ggplot2::theme_bw() + 
    ggplot2::labs(y = class_title, x = x_title) + 
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), 
                   # axis.text.x = ggplot2::element_text(size = 11, color = "black"), 
                   # axis.text.y = ggplot2::element_text(size = 11, color = "black"), 
                   # axis.title.x = ggplot2::element_text(size = 13, color = "black", 
                   #                                      face = "bold"), 
                   # axis.title.y = ggplot2::element_text(size = 13, color = "black", 
                   #                                      face = "bold"), 
                   # legend.title = ggplot2::element_text(size = 13, color = "black", 
                   #                                      face = "bold"), 
                   # legend.text = ggplot2::element_text(size = 11, color = "black"), 
                   legend.position = "right", 
                   strip.text = ggplot2::element_text(size = 9, 
                                                      margin = margin(0.05,0,0.05,0, "cm")),
                   strip.background = ggplot2::element_rect(fill = "gray90"))
}


vars_across_class_wolcan <- function(c_all, cov_df, sampling_wt, res, stratum_id = NULL, 
                              cluster_id = NULL, digits = 1, col_props = TRUE) {
  if (!is.factor(c_all)) {
    stop("c_all must be a factor")
  }
  # Check object class
  if (!(inherits(res, c("swolca", "wolca")))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  # Covariate variable names
  var_names <- colnames(cov_df)
  # Covariate variable levels
  var_levels_list <- sapply(var_names, function(x) get_levels(df = cov_df, var = x))
  # Continuous variables
  cont_vars <- which(sapply(cov_df, function(x) !is.factor(x)))
  # Get number of rows and columns for the dataframe
  num_rows <- length(unlist(sapply(cov_df, levels))) + length(cont_vars) + 2
  num_cols <- length(levels(c_all)) + 3
  
  # Set survey design
  if (is.null(cluster_id)) {
    cluster_id = 1:length(c_all)
  }
  if (is.null(stratum_id)) {
    svy_design <- survey::svydesign(ids = ~cluster_id,
                                    weights = ~sampling_wt,
                                    data = data.frame(cov_df, Class = c_all))
  } else {
    svy_design <- survey::svydesign(ids = ~cluster_id,
                                    weights = ~sampling_wt,
                                    strata= ~stratum_id,
                                    data = data.frame(cov_df, Class = c_all))
  }
  
  # Initialize dataframe
  output_df <- as.data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  # Column names, including latent class names
  colnames(output_df) <- c("Variable", "Level", levels(c_all), "Overall")
  # Variable names
  output_df[, 1] <- c("N: % (posterior SD %)", "n: %", 
                      unlist(sapply(1:length(var_names), function(x) 
                        if(x %in% cont_vars){  # continuous variables
                          c(paste0(var_names[x], ": Mean"))
                        } else {  # factor variables
                          c(paste0(var_names[x], ": %"), rep("", length(var_levels_list[[x]]) - 1))
                        })))
  # Variable levels
  output_df[, 2] <- c("", "", unlist(var_levels_list))
  # Get estimates for population and sample sizes
  output_df[1, -c(1, 2)] <- unlist(get_cov_props(svy_design = svy_design, 
                                                 cov = "population", var_levels = "", 
                                                 res = res, digits = digits))
  output_df[2, -c(1, 2)] <- unlist(get_cov_props(svy_design = svy_design, 
                                                 cov = "sample", var_levels = "", 
                                                 digits = digits))
  # Get estimates for covariates
  row <- 2
  for (i in 1:length(var_names)) {
    var_levels <- var_levels_list[[i]]
    output_df[row + 1:length(var_levels), -c(1, 2)] <- 
      unlist(get_cov_props(svy_design = svy_design, cov = var_names[i], 
                           var_levels = var_levels, digits = digits, 
                           col_props = col_props))
    row <- row + length(var_levels)
  }
  # Output table
  return(output_df)
}
