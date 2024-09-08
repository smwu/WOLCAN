/* Weighted logistic regression*/

data {
    int<lower=1> n; // number of observations
	  int<lower=1> q; // number of linear predictors
    int<lower=0, upper = 1> y[n]; // Response variable
    matrix[n, q] V; // coefficient matrix
    vector<lower=0>[n] weights;  // individual-level weights
}

parameters{
  vector[q] xi; /* regression coefficients from linear predictor */
}

transformed parameters{
  vector[n] mu;
  mu = V * xi; /* linear predictor */
} /* end transformed parameters block */

model{
  
  /* recommended normal informative prior rather than uniform improper prior in 
  (-inf, inf), following https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations */
  // xi ~ normal(0, 2.5);
  
  /* directly update the log-probability for sampling */
  // for (i in 1:n) {
  //   target += weights[i] * bernoulli_logit_lpmf(y[i] | mu[i]);
  // }
  target += bernoulli_logit_lpmf(y | mu);
  
} /* end model{} block */
