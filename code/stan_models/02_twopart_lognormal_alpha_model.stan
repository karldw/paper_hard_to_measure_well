
// Import the shared code:
#include 00_twopart_data_shared.stan
#include 02_twopart_lognormal_alpha_param.stan
model {
  // Priors:
  // Scale the prior std dev on leak size params by `time_period_hr` (different
  // `time_period_hr` values give widely ranging `leak_private_value`)
  b_y ~ student_t(3, 0, 3 * time_period_hr);
  temp_y_inter ~ student_t(3, 0, 3 * time_period_hr);
  sigma_y ~  student_t(3, 0, 3 * time_period_hr);
  // These priors for the logit were chosen so the predicted probability was
  // roughly flat (even priors as wide as normal(0, 0.75) end up putting a
  // large fraction of the probability quite close to 0 or 1)
  // Further discussion: https://mikedecr.github.io/post/nonflat-implications/
  // NOTE: these will change if I standardize the variables or QR decompose Xc.

  // Implicit uniform prior on common_A
  // Remember this is the prior on 1/alpha, and is combined with the [-1, 0]
  // bounds on temp_inv_alpha to be a pretty diffuse folded-normal prior.
  temp_inv_alpha ~ normal(-1, 0.7);

  if (!prior_only) {
    // linear predictor for outcome:
    vector[N_obs] mu_y = temp_y_inter + Xc_obs * b_y;

    // Calculate Duan's smearing (robust to non-normal log(y) residuals, but
    // still assumes conditional homoskedasticity)
    // smear_y_log replaces (sigma_y ^ 2 / 2) in the expectation of e
    // y_shift_log and mu_y are the values and mean on the log scale, so
    // (y_shift_log - mu_y) are the residuals on the log scale.
    // We want the mean of the exp(log-scale residual) to as part of our expression
    // for the expectation of e * p.
    real smear_y = mean(exp(y_shift_log - mu_y));
    // Calculate the expectation of e * p * T (private value of a leak when
    // leaking), re-adding the shift_amount
    vector[N] leak_private_value = time_period_hr * price .*
      (exp(temp_y_inter + Xc * b_y) * smear_y + shift_amount);
    vector[N] prob_leak = (leak_private_value ./ common_A) ^ temp_inv_alpha;

    y_shift ~ lognormal(mu_y, sigma_y);
    obs_dummy ~ bernoulli(prob_leak);
  }
}
generated quantities {
  // Undo the centering from before to calculate the intercepts intercept
  real b_y_intercept = temp_y_inter - dot_product(means_X, b_y);
  real alpha = inv(temp_inv_alpha);
}
