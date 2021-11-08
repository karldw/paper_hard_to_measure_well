
// Import the shared code:
#include 00_twopart_data_shared.stan
#include 04_twopart_lognormal_meas_err_alpha_param.stan
model {
  // Priors:
  b_y ~ student_t(3, 0, 3);
  temp_y_inter ~ student_t(3, 0, 3);
  sigma_y ~  student_t(3, 0, 3);
  // These priors for the logit were chosen so the predicted probability was
  // roughly very flat (even priors as wide as normal(0, 1) end up putting a
  // large fraction of the probability quite close to 0 or 1)
  // Further discussion: https://mikedecr.github.io/post/nonflat-implications/
  // NOTE: these will change if I standardize the variables or QR decompose Xc.
  b_obs ~ normal(0, 0.75);
  temp_obs_inter ~ normal(0, 0.75);

  if (!prior_only) {

    // linear predictor for outcome:
    vector[N_obs] mu_y = temp_y_inter + Xc_obs * b_y;

    // See notes on Duan smearing in 02_twopart_lognormal_alpha_model.stan
    // (Except here we're using log(temp_y_latent) instead of y_shift_log)
    real smear_y = mean(exp(log(temp_y_latent) - mu_y));
    // Calculate the expectation of e * p * T (private value of a leak when
    // leaking), re-adding the shift_amount
    vector[N] leak_private_value = time_period_hr * price .*
      (exp(temp_y_inter + Xc * b_y) * smear_y + shift_amount);
    vector[N] A_i = inv_logit(temp_obs_inter + Xc * b_obs) * shift_amount_times_price_min_times_T;
    vector[N] prob_leak = 1 - (leak_private_value ./ A_i) ^ temp_inv_alpha;

    // Note: do we want to think about a special case for no measurement error?

    // Observed data is distributed normally according to latent data and
    // reported std err.
    y_shift ~ normal(temp_y_latent, y_se);
    // Then model the latent data as a function of X and sigma.
    // Note that obs_dummy is currently unaffected.
    temp_y_latent ~ lognormal(mu_y, sigma_y);
    obs_dummy ~ bernoulli(prob_leak);
  }
}
generated quantities {
  // Undo the centering from before to calculate the intercepts intercept
  real b_y_intercept = temp_y_inter - dot_product(means_X, b_y);
  real b_obs_intercept = temp_obs_inter - dot_product(means_X, b_obs);
  real alpha = inv(temp_inv_alpha);
}
