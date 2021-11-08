
// 04_twopart_lognormal_meas_err_alpha_param.stan
// These are the parameters for the lognormal model with measurement error and
// and estimated alpha param.
// Meant to include in 04_twopart_lognormal_meas_err_alpha_model and
// 04_twopart_lognormal_meas_err_alpha_generate
parameters {
  vector[Kc] b_y;  // population-level effects for observed
  real temp_y_inter;  // temporary intercept for centered predictors
  real<lower=0> sigma_y;  // residual SD
  vector[Kc] b_obs;  // population-level effects for hurdle
  real temp_obs_inter;  // temporary intercept for centered predictors
  vector<lower=0>[N_obs] temp_y_latent; // latent variable for y with measurement error
  // Later R code will ignore anything with the "temp_" prefix
  real<lower=-1, upper=0> temp_inv_alpha; // cost param (1 / alpha)
}
transformed parameters {
}
