
// Import the shared code:
#include 00_twopart_data_shared.stan
#include 04_twopart_lognormal_meas_err_alpha_param.stan
generated quantities {
  // Note: using the lognormal_rng here is a higher-variance choice than
  // using the expected value. Shouldn't matter for point values, but might
  // matter for CI calc.

  vector[N] leak_size_expect;
  vector[N] leak_size_draw;
  vector[N] prob_leak;
  vector[N] prob_size_above_threshold;

  // These are the A_i and alpha in the stan model, but adding more descriptive
  // names here for saving.
  vector[N] cost_param_A = inv_logit(temp_obs_inter + Xc * b_obs) * shift_amount_times_price_min_times_T;
  vector[N] cost_param_alpha = rep_vector(inv(temp_inv_alpha), N);

  {
    // Inner block here to avoid printing temporary vars.
    // Most of this math is the same as in 02_twopart_lognormal_alpha_model.stan
    // See that file for details.
    vector[N_obs] mu_y_obs = temp_y_inter + Xc_obs * b_y;
    vector[N] mu_y_all = temp_y_inter + Xc * b_y;
    leak_size_draw = to_vector(lognormal_rng(mu_y_all, sigma_y)) + shift_amount;
    // This is all the same as 02_twopart_lognormal_alpha_generate.stan, except
    // here we're using log(temp_y_latent) instead of y_shift_log
    real smear_y = mean(exp(log(temp_y_latent) - mu_y_obs));
    // Calculate the expectation of e * p, re-adding the shift_amount
    leak_size_expect = exp(mu_y_all) * smear_y + shift_amount;
    prob_leak = (leak_size_expect .* price * time_period_hr ./ cost_param_A) ^ temp_inv_alpha;

    // Note: `detect_threshold_shifted` is a single real number. For different
    // detect_threshold_shifted values, create different variables.
    for (i in 1:N) {
      prob_size_above_threshold[i] = 1 - lognormal_cdf(
        detect_threshold_shifted | mu_y_all[i], sigma_y
      );
    }
  }
}
// I thik this is still true:
// This generate code currently doesn't work with bootstrap, in almost all cases having mismatches in the temp_y_latent param.
// I think what's going on here is different random number generation, so the bootstrap samples don't end up with the same number of latent variables.
// This could be solved by:
// 1. Allowing stan to ignore differences in parameters
// 2. Allowing stan to not write some parameters
// 3. Post-editing the fits file to harmonize or delete the param that are causing trouble.
// 4. Making the random number generation line up exactly.
