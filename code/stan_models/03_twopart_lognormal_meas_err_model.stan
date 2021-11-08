
// Import the shared code:
#include 00_twopart_data_shared.stan
#include 03_twopart_lognormal_meas_err_param.stan
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
    // Note: do we want to think about a special case for no measurement error?

    // Observed data is distributed normally according to latent data and
    // reported std err.
    y_shift ~ normal(temp_y_latent, y_se);
    // Then model the latent data as a function of X and sigma.
    // Note that obs_dummy is currently unaffected.
    temp_y_latent ~ lognormal(mu_y, sigma_y);
    obs_dummy ~ bernoulli_logit_glm(Xc, temp_obs_inter, b_obs);
  }
}
generated quantities {
  // Undo the centering from before to calculate the intercepts intercept
  real b_y_intercept = temp_y_inter - dot_product(means_X, b_y);
  real b_obs_intercept = temp_obs_inter - dot_product(means_X, b_obs);
}
