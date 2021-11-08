
// Import the shared code:
#include 00_twopart_data_shared.stan
#include 05_twopart_normal_qr_param.stan
generated quantities {
  // Note: using the lognormal_rng here is a higher-variance choice than
  // using the expected value. Shouldn't matter for point values, but might
  // matter for CI calc.
  // Note: could also do duan-smoothed version for leak_size_expect
  vector[N] leak_size_expect = exp(temp_y_inter + Xc * b_y + (sigma_y ^ 2 / 2)) + shift_amount;
  vector[N] leak_size_draw = to_vector(lognormal_rng(temp_y_inter + Xc * b_y, sigma_y)) + shift_amount;

  // Note: this prob_leak re-builds b_obs from the generated quantities block of
  // 05_twopart_normal_qr_model.stan (need to reverse the QR transform).
  vector[N] prob_leak;
  // Local block here to avoid writing out the whole X_combo matrix
  // See https://github.com/stan-dev/cmdstan/issues/553
  {
    // Note: there are other ways we could reasonably include the leak value
    // here: (a) Duan smearing, (b) pred_size directly, (c) current approach
    matrix[N, Kc + 1] X_combo = append_col(Xc, leak_size_expect);
    // matrix[N, K_qr] R_inverse = inverse(qr_thin_R(X_combo) / sqrt(N - 1));
    matrix[N, K_qr] Q = qr_thin_Q(X_combo) * sqrt(N - 1);

    // Other option here for pred_size would be `exp(temp_y_inter + Xc * b_y) *
    // smear_y` if I really wanted to commit to the smearing transform.
    prob_leak = inv_logit(temp_obs_inter + Q * b_obs_qr);
  }
}
