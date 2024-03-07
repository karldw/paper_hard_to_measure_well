
functions {

  // Which elements of x are != 0?
  array[] int which(array[] int x) {
    int N = size(x);
    array[N] int is_nonzero;
    for (i in 1:N) {
      is_nonzero[i] = x[i] != 0;
    }
    int n_nonzero = sum(is_nonzero);
    array[n_nonzero] int result = zeros_int_array(n_nonzero);
    int j = 0;
    for (i in 1:N) {
      if (is_nonzero[i]) {
        j += 1;
        result[j] = i;
      }
    }
    return result;
  }
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  vector<lower=0>[N] noise; // reported std err of response (not always used)
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=0, upper=1> prior_only;  // should the likelihood be ignored?
  real<lower=0> shift_amount;
  int<lower=0, upper=1> bootstrap; // should we use the original data or bootstrap?
  vector<lower=0, upper=10000>[N] price; // commodity_price
  real<lower=0.01> time_period_hr; // assumed leak duration (only useful for cost_coef_models)
}
transformed data {
  array[N] int<lower=1, upper=N> boot_idx;
  simplex[N] uniform = rep_vector(1.0 / N, N);
  for (i in 1:N) {
    // boot_idx is either 1:N if bootstrap == 0, or a random draw if bootstrap == 1
    boot_idx[i] = bootstrap ? categorical_rng(uniform) : i;
  }
  // Can't write to Y, so make a new Y_reindex
  // NB: later code *should not* use the original unshuffled X or Y.
  vector[N] Y_reindex = Y[boot_idx];
  matrix[N, K] X_reindex = X[boot_idx, ];

  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    // Xc and means_X are the shuffled, centered variables.
    means_X[i - 1] = mean(X_reindex[, i]);
    Xc[, i - 1] = X_reindex[, i] - means_X[i - 1];
  }

  // Construct the y_shift variable
  array[N] int obs_dummy;
  for (i in 1:N) {
    obs_dummy[i] = Y_reindex[i] > shift_amount;
  }
  int N_obs = sum(obs_dummy); // number of observed y, after shifting

  if (N_obs < 2 || N_obs > (N - 2)) {
    reject("Can't estimate model");
  }
  array[N_obs] int<lower=1> obs_idx = which(obs_dummy); // index of observed y
  // Models will use either y_shift or y_shift_log, but that's fine.
  vector[N_obs] y_shift = Y_reindex[obs_idx] - shift_amount;
  vector[N_obs] y_shift_log = log(y_shift);
  matrix[N_obs, Kc] Xc_obs = Xc[obs_idx, ];
  // Note: might want y_se for the log-scale, but need to think more.
  vector[N_obs] y_se = noise[boot_idx][obs_idx];

  vector[N] log_price = log(price);
  real shift_amount_times_price_min_times_T = shift_amount * min(price) * time_period_hr;

  // Note: hard-coding detect_threshold for now. Needs to match detect_threshold
  // in the outcomes_analysis.py program.
  real<lower=0.01> detect_threshold_shifted = 100 - shift_amount;
}

parameters {
  vector[Kc] b_y;  // population-level effects for observed
  real temp_y_inter;  // temporary intercept for centered predictors
  real<lower=0> sigma_y;  // residual SD
  vector[Kc] b_obs;  // population-level effects for hurdle
  real temp_obs_inter;  // temporary intercept for centered predictors

  real temp_alpha_inter;
  vector[Kc] b_alpha;
}
transformed parameters {
}

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
