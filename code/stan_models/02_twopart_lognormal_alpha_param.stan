
// 02_twopart_lognormal_alpha_param.stan
// These are the parameters for the lognormal model with cost param alpha
// Meant to include in 02_twopart_lognormal_alpha_model and
// 02_twopart_lognormal_alpha_generate
parameters {
  vector[Kc] b_y;  // population-level effects for observed
  real temp_y_inter;  // temporary intercept for centered predictors
  real<lower=0> sigma_y;  // residual SD

  // common_A needs to be lower than min(e * p) for the fraction to be well defined
  real<lower=0, upper=shift_amount_times_price_min_times_T> common_A;
  real<lower=-1, upper=0> temp_inv_alpha; // cost param (1 / alpha)
}
transformed parameters {
}
