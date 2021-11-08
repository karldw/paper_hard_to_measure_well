
// Import the shared code:
#include 00_twopart_data_shared.stan
#include 01_twopart_lognormal_param.stan
generated quantities {
  // Note: using the lognormal_rng here is a higher-variance choice than
  // using the expected value. Shouldn't matter for point values, but might
  // matter for CI calc.
  // Note: could also do duan-smoothed version for leak_size_expect
  vector[N] leak_size_expect = exp(temp_y_inter + Xc * b_y + (sigma_y ^ 2 / 2)) + shift_amount;
  vector[N] leak_size_draw = to_vector(lognormal_rng(temp_y_inter + Xc * b_y, sigma_y)) + shift_amount;
  vector[N] prob_leak = inv_logit(temp_obs_inter + Xc * b_obs);

  vector[N] prob_size_above_threshold;
  // Note: `detect_threshold` is a single real number. For different
  // detect_threshold values, create different variables.
  for (i in 1:N) {
    prob_size_above_threshold[i] = 1 - lognormal_cdf(
      detect_threshold_shifted | temp_y_inter + Xc[i] * b_y, sigma_y
    );
  }
}
