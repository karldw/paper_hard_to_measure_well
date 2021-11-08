
// 06_twopart_normal_qr_alpha_param.stan
// These are the parameters for the normal-expected value model with cost param
// Meant to include in 06_twopart_normal_qr_alpha_model and
// 06_twopart_normal_qr_alpha_generate
parameters {
  // Coefficients for the measured (lognormal) outcome
  vector[Kc] b_y;  // population-level effects
  real temp_y_inter;  // temporary intercept for centered predictors
  real<lower=0> sigma_y;  // residual SD
  // Coefficients for the observed-or-not indicator
  vector[K_qr] b_obs_qr;  // population-level effects
  real temp_obs_inter;  // temporary intercept for centered predictors
  real<upper=0> cost_alpha; // common cost param
  real<lower=0> sigma_alpha; // residual SD of cost param
}
transformed parameters {
}
