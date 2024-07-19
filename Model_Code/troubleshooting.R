# Checking how many draws have different K
for (samp_i in 1:100) {
  res_path <- paste0(wd, res_dir, "scen_", scenario, "/samp_", samp_i, "_wolcan_results.RData")
  load(res_path)
  if (res$estimates_adjust$K != 3) {
    print(paste0("i: ", samp_i, ", K_fixed:", res$estimates_adjust$K))
  }
  if (res$estimates_adjust$K_red != 3) {
    print(paste0("i: ", samp_i, ", K_red:", res$estimates_adjust$K_red))
  }
}

# Troubleshooting thetas for each draw
temp1 <- apply(theta_red_stack[201:400, , , ], c(2, 3, 4), stats::median)
temp2 <- apply(temp1, c(1, 2), which.max)



### Simulated data summaries
summary(1/pop$pi_B)
IQR(1/pop$pi_B)

# Get sample theta
sim_samp <- sim_samp_B
samp_theta <- array(NA, dim=c(J, K, R))
for (j in 1:J) {
  for (k in 1:K) {
    for (r in 1:R) {
      samp_theta[j,k,r] <- sum((sim_samp$X_data[,j]==r) & 
                                 (sim_samp$c_all==k)) / sum(sim_samp$c_all==k) 
    }
  }
}

# Get population theta
pop_theta <- array(NA, dim=c(J, K, R))
for (j in 1:J) {
  for (k in 1:K) {
    for (r in 1:R) {
      pop_theta[j,k,r] <- sum((sim_pop$X_data[,j]==r) & 
                                 (sim_pop$c_all==k)) / sum(sim_pop$c_all==k) 
    }
  }
}

samp_theta[1,,]
pop_theta[1,,]
samp_theta[3,,]
pop_theta[3,,]
samp_theta[30,,]
pop_theta[30,,]
hist(sim_pop$pop$A3, freq = FALSE, col = "darkgreen")
hist(sim_samp$covs$A3, freq = FALSE, col=rgb(1,0,0,0.5), add = TRUE)



# Quantiles of weights
load("/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/Results/scen_11/samp_1_wolcan_weights.RData")

quantile(colMaxs(est_weights$wts_post), seq(0.1, 1, by = 0.1))

# Plot quantiles
d <- data.frame(x = colMaxs(est_weights$wts_post))

breaks <- seq(min(d$x), max(d$x), length.out = 50)
quantiles <- quantile(d$x, seq(0, 1, 0.1))
quantiles2 <- sapply(quantiles, function(x) breaks[which.min(abs(x - breaks))])

d$bar <- as.numeric(as.character(cut(d$x, breaks, na.omit((breaks + dplyr::lag(breaks)) / 2))))
d$fill <- cut(d$x, quantiles2, na.omit((quantiles2 + dplyr::lag(quantiles2)) / 2))

ggplot(d, aes(bar, y = 1, fill = fill)) +
  geom_col(position = 'stack', col = 1, show.legend = FALSE, width = diff(breaks)[1], size = 0.3) +
  scale_fill_brewer(type = 'qual', palette = 3) +
  theme_classic() +
  coord_fixed(diff(breaks)[1], expand = FALSE) + # makes square blocks
  labs(x = 'max(weight)', y = 'count')

### Why is unweighted performing well?
load("/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/Results/scen_0/samp_1_wolca_results.RData")
res_unwt <- res
load("/n/holyscratch01/stephenson_lab/Users/stephwu18/WOLCAN/Results/scen_0/samp_1_wolcan_results.RData")
round(res_unwt$estimates$theta_med[1,,], 3)
round(res$estimates_adjust$theta_med[1,,], 3)
round(samp_theta[1,,], 3)
round(pop_theta[1,,], 3)
round(res_unwt$estimates$theta_med[3,,], 3)
round(res$estimates_adjust$theta_med[3,,], 3)
round(samp_theta[3,,], 3)
round(pop_theta[3,,], 3)
round(res_unwt$estimates$theta_med[30,,], 3)
round(res$estimates_adjust$theta_med[30,,], 3)
round(samp_theta[30,,], 3)
round(pop_theta[30,,], 3)


# Theta modes
est_theta_modal <- apply(estimates$theta_med[, order_sub_est, ], c(1,2), max)
true_theta_modal <- apply(true_params$true_theta[ , order_sub_true, ], c(1,2), max) 
dist_modal_wolcan <- abs(est_theta_modal - true_theta_modal)
sum(dist_modal_wolcan)  
mean(dist_modal_wolcan)
## Have to re-run summary_function code first for wolca
est_theta_modal <- apply(estimates$theta_med[, order_sub_est, ], c(1,2), max)
true_theta_modal <- apply(true_params$true_theta[ , order_sub_true, ], c(1,2), max) 
dist_modal_wolca <- abs(est_theta_modal - true_theta_modal)
sum(dist_modal_wolca)  
mean(dist_modal_wolca)
  