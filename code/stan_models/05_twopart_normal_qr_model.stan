
// Import the shared code:
#include 00_twopart_data_shared.stan
#include 05_twopart_normal_qr_param.stan
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
    // Note: we could calculate the smearing with:
    // real smear_y = mean(exp(y_shift_log - (temp_y_inter + Xc_obs * b_y)));
    // And then use `exp(temp_y_inter + Xc * b_y) * smear_y` instead of
    // `exp(temp_y_inter + Xc * b_y + (sigma_y ^ 2 / 2))`.
    // Note that this doesn't matter for predictions -- we'd just end up with
    // a different coef -- but including sigma_y makes it easier to calculate
    // predicted values later.
    // Note that I'm not using the generated quantities block for predicted
    // values because it's a half-baked feature and I'm tired of fighting it.

    // Make the matrix and its QR decomp. Manual recommends scaling by sqrt(N-1)
    // Order matters here - Xc first, then anything else - because of b_obs_intercept
    vector[N] log_emiss = temp_y_inter + Xc * b_y + (sigma_y ^ 2 / 2);
    matrix[N, Kc + 1] X_combo = append_col(Xc, exp(log_emiss));
    matrix[N, K_qr] Q = qr_thin_Q(X_combo) * sqrt(N - 1);

    y_shift_log ~ normal_id_glm(Xc_obs, temp_y_inter, b_y, sigma_y);
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
