
// 08_twopart_lognormal_heterog_alpha_param.stan
// These are the parameters for the lognormal model with cost param alpha
// Meant to include in 08_twopart_lognormal_heterog_alpha_model and
// 08_twopart_lognormal_heterog_alpha_generate
parameters {
  vector[Kc] b_y;  // population-level effects for observed
  real temp_y_inter;  // temporary intercept for centered predictors
  real<lower=0> sigma_y;  // residual SD
  vector[Kc] b_obs;  // population-level effects for hurdle
  real temp_obs_inter;  // temporary intercept for centered predictors
  // real<lower=-1, upper=0> temp_inv_alpha; // cost param (1 / alpha)

  real temp_alpha_inter;
  vector[Kc] b_alpha;
}
transformed parameters {
}
