
// 03_twopart_lognormal_meas_err_param.stan
// These are the parameters for the lognormal model
// Meant to include in 03_twopart_lognormal_meas_err_model and
// 03_twopart_lognormal_meas_err_generate
parameters {
  vector[Kc] b_y;  // population-level effects for observed
  real temp_y_inter;  // temporary intercept for centered predictors
  real<lower=0> sigma_y;  // residual SD
  vector[Kc] b_obs;  // population-level effects for hurdle
  real temp_obs_inter;  // temporary intercept for centered predictors
  vector<lower=0>[N_obs] temp_y_latent; // latent variable for y with measurement error
  // Later R code will ignore anything with the "temp_" prefix
}
transformed parameters {
}
