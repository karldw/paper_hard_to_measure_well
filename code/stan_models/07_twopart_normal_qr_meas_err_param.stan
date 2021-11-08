
// 07_twopart_normal_qr_meas_err_param.stan
// These are the parameters for the normal-expected value model
// Meant to include in 07_twopart_normal_qr_meas_err_model
parameters {
  // Coefficients for the measured (lognormal) outcome
  vector[Kc] b_y;  // population-level effects
  real temp_y_inter;  // temporary intercept for centered predictors
  real<lower=0> sigma_y;  // residual SD
  // Coefficients for the observed-or-not indicator
  vector[K_qr] b_obs_qr;  // population-level effects
  real temp_obs_inter;  // temporary intercept for centered predictors
  vector<lower=0>[N_obs] temp_y_latent; // latent variable for y with measurement error
  // Later R code will ignore anything with the "temp_" prefix
}
transformed parameters {
}
