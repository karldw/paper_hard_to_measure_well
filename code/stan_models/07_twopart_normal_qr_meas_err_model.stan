
// Import the shared code:
#include 00_twopart_data_shared.stan
#include 07_twopart_normal_qr_meas_err_param.stan
model {
  // These priors were chosen to very roughly fit the range of the data
  // (note that the X are centered, but not all unit-scale)
  // Priors:
  b_y ~ student_t(3, 0, 3);
  temp_y_inter ~ student_t(3, 0, 3);
  sigma_y ~  student_t(3, 0, 3);
  // These priors for the logit were chosen so the predicted probability was
  // roughly very flat (even priors as wide as normal(0, 1) end up putting a
  // large fraction of the probability quite close to 0 or 1)
  // Further discussion: https://mikedecr.github.io/post/nonflat-implications/
  b_obs_qr ~ normal(0, 0.75);
  temp_obs_inter ~ normal(0, 0.75);

  if (!prior_only) {
    // Make the matrix and its QR decomp.
    // See comments in 05_twopart_normal_qr_model.stan
    vector[N] log_emiss = temp_y_inter + Xc * b_y + (sigma_y ^ 2 / 2);
    matrix[N, Kc + 1] X_combo = append_col(Xc, exp(log_emiss));
    matrix[N, K_qr] Q = qr_thin_Q(X_combo) * sqrt(N - 1);

    // Observed data is distributed normally according to latent data and
    // reported std err.
    y_shift ~ normal(temp_y_latent, y_se);
    // Then model the latent data as a function of X and sigma.
    // Note that obs_dummy is only affected via b_y.
    log(temp_y_latent) ~ normal_id_glm(Xc_obs, temp_y_inter, b_y, sigma_y);

    obs_dummy ~ bernoulli_logit_glm(Q, temp_obs_inter, b_obs_qr);
  }
}
generated quantities {
  // Declare these variables here so they get written. We'll assign to them in
  // the block below.
  vector[K_qr] b_obs;

  // Optional code to simulate outcomes:
  // (variable definitions can't be inside an `if` block if we want them printed)
  // real y_pred[N] = normal_rng(temp_y_inter + Xc * b_y, sigma_y);
  // int obs_pred[N];

  // Reverse the QR reparameterization to get the b_obs back
  // Local block here to avoid writing out the whole X_combo matrix
  // See https://github.com/stan-dev/cmdstan/issues/553
  {
    matrix[N, Kc + 1] X_combo = append_col(
      Xc,
      exp(temp_y_inter + Xc * b_y + (sigma_y ^ 2 / 2))
    );
    matrix[N, K_qr] R_inverse = inverse(qr_thin_R(X_combo) / sqrt(N - 1));
    b_obs = R_inverse * b_obs_qr;
    // See note above about predicted outcomes
    // obs_pred = bernoulli_logit_rng(temp_obs_inter + X_combo * b_obs);
  }
  // actual population-level intercepts
  real b_y_intercept = temp_y_inter - dot_product(means_X, b_y);
  real b_obs_intercept = temp_obs_inter - dot_product(means_X, b_obs[:Kc]);
}
