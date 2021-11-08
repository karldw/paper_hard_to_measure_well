
// Import the shared code:
#include 00_twopart_data_shared.stan
#include 08_twopart_lognormal_heterog_alpha_param.stan
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
  b_obs ~ normal(0, 0.5);
  temp_obs_inter ~ normal(0, 0.5);
  // Remember this is the prior on 1/alpha, and is combined with the [-1, 0]
  // bounds on temp_inv_alpha to be a pretty diffuse folded-normal prior.
  // temp_inv_alpha ~ normal(-1, 0.7);
  temp_alpha_inter ~ normal(0, 0.75);
  b_alpha ~ normal(0, 0.75);

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
    // We need temp_obs_inter + Xc * b_obs to be in (0, e * p) for every i, which
    // is infeasible for almost all b_obs. But we can do an ad-hoc logit scale to
    // force it to be true, then multiply by shift_amount * min(price)
    // This approach forces A_i to be in a specific range, which is a little
    // unusual, but seems fine.
    // Doesn't work to use min(e * p) here because of identifiability challenges
    vector[N] A_i = inv_logit(temp_obs_inter + Xc * b_obs) * shift_amount_times_price_min_times_T;

    vector[N] inv_alpha = -inv_logit(temp_alpha_inter + Xc * b_alpha);
    vector[N] prob_leak = (leak_private_value ./ A_i) .^ inv_alpha;

    y_shift ~ lognormal(mu_y, sigma_y);
    obs_dummy ~ bernoulli(prob_leak);
  }
}
generated quantities {
  // Undo the centering from before to calculate the intercepts intercept
  real b_y_intercept = temp_y_inter - dot_product(means_X, b_y);
  real b_obs_intercept = temp_obs_inter - dot_product(means_X, b_obs);
  real b_alpha_intercept = temp_alpha_inter - dot_product(means_X, b_alpha);
}
