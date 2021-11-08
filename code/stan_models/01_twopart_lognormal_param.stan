
// 01_twopart_lognormal_param.stan
// These are the parameters for the lognormal model
// Meant to include in 01_twopart_lognormal_model and 01_twopart_lognormal_generate
parameters {
  vector[Kc] b_y;  // population-level effects for observed
  real temp_y_inter;  // temporary intercept for centered predictors
  real<lower=0> sigma_y;  // residual SD
  vector[Kc] b_obs;  // population-level effects for hurdle
  real temp_obs_inter;  // temporary intercept for centered predictors
}
transformed parameters {
}
