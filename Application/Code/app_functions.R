#=================================================
# PROSPECT application functions
# Author: Stephanie Wu
# Date created: 2024/09/29
# Date updated: 2024/09/29
#=================================================

### Frequentist survey-weighted logistic regression
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
wtd_logreg_wolcan <- function(wts_draws, exposure, filter_categ, condition, 
                              save_res = TRUE,
                              save_path, model_covs = "core", svy_data) {
  
  # Print filter category
  if (filter_categ == 1) {
    print("Subsetting to exclude those without the disease outcome")
  } else if (filter_categ == 3) {
    print("Subsetting to exclude those with self-reported diagnosis")
  } else if (filter_categ == 2) {
    print("Subsetting to exclude those with non-reported outcome")
  } else {
    print("No subsetting performed. Full sample is used.")
  }
  
  # Specify the outcome of interest 
  if (condition == "htn") {
    cond_categ <- "htn_categ"
    cond_diag <- "htn_aware"
    cond <- "hypertension"
    other_y <- c("diabetes + high_cholesterol")  # comorbidities
    other_y_diag <- c("t2d_aware + chol_aware")  # comorbid diagnoses
  } else if (condition == "t2d") {
    cond_categ <- "t2d_categ"
    cond_diag <- "t2d_aware"
    cond <- "diabetes"
    other_y <- c("hypertension + high_cholesterol") # comorbidities
    other_y_diag <- c("htn_aware + chol_aware") # comorbid diagnoses
  } else if (condition == "chol") {
    cond_categ <- "chol_categ"
    cond_diag <- "chol_aware"
    cond <- "high_cholesterol"
    other_y <- c("hypertension + diabetes") # comorbidities
    other_y_diag <- c("t2d_aware + htn_aware") # comorbid diagnoses
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
  } else if (model_covs == "core") {  # Core covariates
    if (exposure == "dbh_class") {  # DBH -> disease
      add_covs <- paste0("+ Age_cent + Sex + Educ + Inc_hh + Ethnicity + Urban ", 
                         "+ Physical_activity + Smoking_status + Drinking_status +", 
                         other_y)
    } else if (exposure == "diagnosis") {  # Diagnosis -> DBH
      add_covs <- paste0("+ Age_cent + Sex + Educ + Inc_hh + Ethnicity + Urban ", 
                         "+ Physical_activity + Smoking_status + Drinking_status +", 
                         other_y_diag)
    }
    
  } else if (model_covs == "full") {  # All covariates
    if (exposure == "dbh_class") { # DBH -> disease
      add_covs <- paste0("+ Age_cent + Sex + Educ + Inc_hh + Ethnicity + Urban ",
                         "+ Physical_activity + Smoking_status + Drinking_status +", 
                         other_y, " + Food_security",
                         "+ WIC_SNAP + Social_support_cent + Perceived_stress_cent",
                         "+ Depression + Anxiety")
    } else if (exposure == "diagnosis") { # Diagnosis -> DBH
      add_covs <- paste0("+ Age_cent + Sex + Educ + Inc_hh + Ethnicity + Urban ",
                         "+ Physical_activity + Smoking_status + Drinking_status +", 
                         other_y_diag, " + Food_security",
                         "+ WIC_SNAP + Social_support_cent + Perceived_stress_cent",
                         "+ Depression + Anxiety")
    }
    
  } else if (model_covs == "core_int") {  # Core covariates, with interactions
    if (exposure == "dbh_class") { # DBH -> disease
      add_covs <- paste0("* (Age_cent + Sex) + Physical_activity + Inc_hh + Urban",
                         "+ Educ + Smoking_status + Drinking_status +",
                         other_y)
    } else if (exposure == "diagnosis") { # Diagnosis -> DBH
      add_covs <- paste0("* (Age_cent + Sex + Physical_activity) + Inc_hh + Educ",
                         " + Urban + Smoking_status + Drinking_status +",
                         other_y_diag)
    }
  } else if (model_covs == "full_int") { # All covariates, with interactions
    add_covs <- paste0("* (Age_cent + Sex + Educ + Inc_hh + Urban + Physical_activity",
                       "+ Smoking_status + Food_security",
                       "+ Ethnicity + Drinking_status + WIC_SNAP + Social_support_cent",
                       "+ Perceived_stress_cent + Depression + Anxiety)")
  } else {
    stop("Invalid model_covs input specified.")
  }
  
  # Create brms model formula
  if (exposure == "dbh_class") {
    model_formula <- as.formula(paste0(cond, " | weights(wts_d) ~ dbh_class ",
                                       add_covs)) 
  } else if (exposure == "diagnosis") {
    model_formula <- as.formula(paste0("dbh_class | weights(wts_d) ~ ", cond_diag,
                                       add_covs)) 
  } else {
    stop("Invalid exposure input specified.")
  }
  brms_mod <- brms::brmsformula(model_formula, center = TRUE)
  
  
  # For each draw, run the variance-adjusted survey-weighted regression and 
  # store results
  for (d in 1:D) {
    print(paste0("d: ", d))
    set.seed(d)
    # Specify weights for the draw
    svy_data_d <- svy_data %>%
      mutate(wts_d = wts_draws[, d]) 
    # Subset sample if specified
    svy_data_subset_d <- svy_data_d %>%
      filter(!!rlang::sym(cond_categ) != filter_categ)
    # Remove NAs and center continuous variables
    if (exposure == "dbh_class") {
      svy_data_subset_d <- svy_data_subset_d %>%
        drop_na(hypertension, diabetes, high_cholesterol, 
                Age, Sex, Educ, Inc_hh, Ethnicity, Urban, Physical_activity, 
                Smoking_status, Drinking_status, Food_security, WIC_SNAP, 
                Social_support, Perceived_stress, Depression, Anxiety)%>%
        # Center continuous variables 
        mutate(Age_cent = Age - mean(Age, na.rm = TRUE), 
               Social_support_cent = Social_support - 
                 mean(Social_support, na.rm = TRUE),
               Perceived_stress_cent = Perceived_stress - 
                 mean(Perceived_stress, na.rm = TRUE))
    } else if (exposure == "diagnosis") {
      svy_data_subset_d <- svy_data_subset_d %>%
        drop_na(htn_aware, t2d_aware, chol_aware, 
                Age, Sex, Educ, Inc_hh, Ethnicity, Urban, Physical_activity, 
                Smoking_status, Drinking_status, Food_security, WIC_SNAP, 
                Social_support, Perceived_stress, Depression, Anxiety) %>%
        # Center continuous variables 
        mutate(Age_cent = Age - mean(Age, na.rm = TRUE), 
               Social_support_cent = Social_support - 
                 mean(Social_support, na.rm = TRUE),
               Perceived_stress_cent = Perceived_stress - 
                 mean(Perceived_stress, na.rm = TRUE))
    }
    
    # Renormalize weights
    svy_data_subset_d <- svy_data_subset_d %>%
      mutate(wts_d = wts_d * nrow(svy_data_subset_d) / sum(wts_d))
    
    # Specify survey design  
    svy_des_subset_d <- survey::svydesign(ids = ~1, weights = ~wts_d,
                                          data = svy_data_subset_d)
    
    # Fit weighted logistic regression for the draw
    if (exposure == "dbh_class") { # DBH -> disease
      fit_d <- csSampling::cs_sampling_brms(svydes = svy_des_subset_d, 
                                            brmsmod = brms_mod, 
                                            data = svy_data_subset_d, 
                                            family = bernoulli(link = "logit"),
                                            ctrl_stan = list(chains = 3, iter = 2000, 
                                                             warmup = 1000, thin = 1))
    } else if (exposure == "diagnosis") { # diagnosis -> DBH
      fit_d <- csSampling::cs_sampling_brms(svydes = svy_des_subset_d, 
                                            brmsmod = brms_mod, 
                                            data = svy_data_subset_d, 
                                            family = categorical(link = "logit"),
                                            ctrl_stan = list(chains = 3, iter = 2000, 
                                                             warmup = 1000, thin = 1))
    }
    
    # Store results for this draw
    sampled_parms_draws[[d]] <- fit_d$sampled_parms
    adj_parms_draws[[d]] <- fit_d$adjusted_parms
    stan_fit_draws[[d]] <- fit_d$stan_fit
  }
  
  # Row-bind together results from all draws
  all_adj_parms <- as.data.frame(do.call(rbind, adj_parms_draws))
  
  # Save and return output for all draws
  # Note: 'b_Intercept' is intercept on original scale. 'Intercept' is after 
  # internal centering and should not be used
  wtd_logreg_res <- list(all_adj_parms = all_adj_parms,  # adjusted parameters
                         sampled_parms_draws = sampled_parms_draws, # Original parameters
                         stan_fit_draws = stan_fit_draws, 
                         brms_mod = brms_mod, 
                         data = svy_data_subset_d)
  
  if (save_res) {
    save(wtd_logreg_res, file = paste0(save_path, "logreg.RData"))
  }
  
  return(wtd_logreg_res)
}

# Function to summarize output from `wtd_logreg_wolcan()` 
# Inputs:
#   fit_obj: Fit object resulting from a call to `wtd_logreg_wolcan()`
#   quant_lb: CI lower bound quantile. Default is `0.025`.
#   quant_ub: CI upper bound quantile. Default is `0.975`.
#   round_digits: Number of digits to round to. Default is `3`. 
#   parm_names: Boolean specifying if coefficient names should be outputted. 
# Default is `TRUE`.
#   exponentiate: Boolean specifying if results should be exponentiated to the 
# odds scale (`TRUE`; default) or left on the log-odds scale (`FALSE`).
#   diag: Boolean specifying if the exposure is disease diagnosis, which would 
# result in a categorical outcome of dietary behavior class. Default is `FALSE`, 
# indicating dietary behavior class as the exposure and disease status as the 
# binary outcome. 
# Output: Table of summarized regression coefficients
summarize_parms <- function(fit_obj, quant_lb = 0.025, quant_ub = 0.975, 
                            round_digits = 3, parm_names = TRUE, 
                            exponentiate = TRUE, diag = FALSE) {
  all_adj_parms <- fit_obj$all_adj_parms
  summary_adj_parms <- t(apply(all_adj_parms, 2, function(x) 
    c(mean(x, na.rm = TRUE), 
      quantile(x, c(quant_lb, quant_ub), na.rm = TRUE),
      mean(x > 0, na.rm = TRUE),
      mean(x < 0, na.rm = TRUE))))
  if (exponentiate) {
    summary_adj_parms[, 1:3] <- exp(summary_adj_parms[, 1:3])
  }
  summary_adj_parms <- as.data.frame(round(summary_adj_parms, round_digits))
  colnames(summary_adj_parms) <- c("Mean", paste0(quant_lb*100, "%"), 
                                   paste0(quant_ub*100, "%"), "P(xi>0)", "P(xi<0)")
  if (parm_names) {
    if (diag) {  # categorical outcome
      mod_mat <- brms::make_standata(formula = fit_obj$brms_mod$formula, 
                                     data = fit_obj$data, family = "categorical")
      parm_names <- colnames(mod_mat$X_mu2)
      categs <- levels(fit_obj$data$dbh_class)
      n_categs <- length(categs)
      parm_names <- c(parm_names[-1], parm_names[1])  # Intercept at end
      parm_names_cs <- rep(parm_names, times = (n_categs-1))
      summary_adj_parms$dbh_class <- NA
      summary_adj_parms$dbh_class[1:length(parm_names_cs)] <-
        rep(categs[-1], each = length(parm_names))
    } else {  # binary outcome
      mod_mat <- brms::make_standata(formula = fit_obj$brms_mod$formula, 
                                     data = fit_obj$data, family = "bernoulli")
      parm_names <- colnames(mod_mat$X)
      parm_names_cs <- parm_names[-1]  # Remove intercept
    }
    summary_adj_parms$Variable <- NA
    summary_adj_parms$Variable[1:length(parm_names_cs)] <- parm_names_cs
  }
  
  return(summary_adj_parms)
}


# Table of estimated weights by covariates
# Inputs:
#   var_vec: String vector indicating variables over which to summarize.
#   wts: nx1 vector of weights (or pseudo-weights) for all individuals in the sample
#   data: Matrix or dataframe of sample data
#   digits: Number of digits to round to. Default is `1`.
# Output: Returns a dataframe where each row is the survey-weighted mean of a 
# variable in `var_vec`. Proportions for each level are used if the variable is 
# categorical. 
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
# Similar to baysc::plot_pattern_probs() function, but the barplot is horizontal
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


# Create table of distribution of classes across covariate variables
# Similar to baysc::vars_across_class() function
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



# Create forest plot of outcome regression odds ratios comparing the dietary 
# behavior patterns
# Inputs:
#   summ_mod_list: List of models summarized using `summarize_parms()`
# Output: Forest plot of odds ratios for type 2 diabetes, hypertension, and 
# high cholesterol
create_mod_plot <- function(summ_mod_list) {
  n_mods <- length(summ_mod_list)
  plot_df <- as.data.frame(matrix(NA, nrow = 3*n_mods, ncol = 5))
  colnames(plot_df) <- c("Outcome", "Pattern", "OR", "Lower", "Upper")
  plot_df[, 1] <- c(rep("Type 2 Diabetes", 3), rep("Hypertension", 3), 
                    rep("High Cholesterol", 3))
  plot_df[, 2] <- rep(c("DBP2", "DBP3", "DBP4"), times = 3)
  counter <- 0
  for (i in 1:length(summ_mod_list)) {
    plot_df[1:3 + (counter * 3), 3:5] <- summ_mod_list[[i]][1:3, 1:3]
    counter <- counter + 1
  }
  
  plot_df %>% 
    ggplot(aes(x=Outcome, y=OR, ymin=Lower, ymax=Upper,col=fct_rev(Pattern),
               fill=fct_rev(Pattern))) + 
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    theme_bw() + 
    #specify position here
    geom_linerange(size=5,position=position_dodge(width = 0.5)) +
    geom_hline(yintercept=1, lty=2) +
    #specify position here too
    geom_point(size=3, shape=21, colour="white", stroke = 0.5,
               position=position_dodge(width = 0.5)) +
    xlab("Outcome") + ylab("Odds Ratio Scale") + 
    geom_hline(yintercept = 1, linetype = "dashed", col = "gray") + 
    guides(fill = guide_legend(reverse = TRUE),
           color = guide_legend(reverse = TRUE)) + 
    labs(fill = "Pattern", col = "Pattern") +
    coord_flip() 
}

# Create forest plot of regression odds ratios with faceting for multiple plots
# comparing subsetted and non-subsetted models
# Inputs:
#   summ_mod_list_multiple: List of lists, where the inner list is models 
# summarized using `summarize_parms()`, and the outer list is the non-subsetted 
# vs. subsetted models.
#   mod_labels: String vector specifying the model labels for the non-subsetted
# vs. subsetted models.
#   legend_labels: Labels for the dietary behavior classes. Default is `NULL`, 
# which results in the number of the pattern appendd to "DBP" as the label.
# Output: Facetted forest plots of odds ratios for type 2 diabetes, hypertension, 
# and high cholesterol, for the non-subsetted and subsetted models.
create_mod_plot_multiple <- function(summ_mod_list_multiple, mod_labels,
                                     legend_labels = NULL) {
  if (is.null(legend_labels)) {
    legend_labels <- c("DBP2", "DBP3", "DBP4")
  }
  n_mult <- length(summ_mod_list_multiple)
  plot_df_list <- list()
  for (i in 1:n_mult) {
    summ_mod_list <- summ_mod_list_multiple[[i]]
    n_mods <- length(summ_mod_list)
    plot_df_i <- as.data.frame(matrix(NA, nrow = 3*n_mods, ncol = 5))
    colnames(plot_df_i) <- c("Outcome", "Pattern", "OR", "Lower", "Upper")
    plot_df_i[, 1] <- c(rep("Type 2 Diabetes", 3), rep("Hypertension", 3), 
                      rep("High Cholesterol", 3))
    plot_df_i[, 2] <- rep(legend_labels, times = 3)
    counter <- 0
    for (j in 1:length(summ_mod_list)) {
      plot_df_i[1:3 + (counter * 3), 3:5] <- summ_mod_list[[j]][1:3, 1:3]
      counter <- counter + 1
    }
    n_rows <- nrow(plot_df_i)
    plot_df_list[[i]] <- plot_df_i
  }
  
  plot_df <- as.data.frame(do.call("rbind", plot_df_list))
  plot_df$Subset <- factor(rep(mod_labels, each = n_rows), levels = mod_labels)
  
  plot_df %>% 
    ggplot(aes(x=Outcome, y=OR, ymin=Lower, ymax=Upper,col=Pattern,
               fill=Pattern)) + 
    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
    theme_bw() + 
    #specify position here
    geom_linerange(size=4.5,position=position_dodge(width = 0.5)) +
    geom_hline(yintercept=1, lty=2) +
    #specify position here too
    geom_point(size=2.5, shape=21, colour="white", stroke = 0.5,
               position=position_dodge(width = 0.5)) +
    xlab("Outcome") + ylab("Odds Ratio Scale") + 
    geom_hline(yintercept = 1, linetype = "dashed", col = "gray") + 
    guides(fill = guide_legend(reverse = TRUE),
           color = guide_legend(reverse = TRUE)) + 
    labs(fill = "Pattern", col = "Pattern") +
    coord_flip() + 
    # Facet by whether or not subsetting was done
    facet_grid(Subset ~ .) + 
    theme(strip.background = element_rect(fill = "aliceblue"))
}


# Basic plots of results for troubleshooting
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
