
functions {

  // Which elements of x are != 0?
  array[] int which(array[] int x) {
    int N = size(x);
    array[N] int is_nonzero;
    for (i in 1:N) {
      is_nonzero[i] = x[i] != 0;
    }
    array[sum(is_nonzero)] int result;
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
    // boot_idx is either 1:N if bootstrap == 0
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

  // Number of cols in the QR reparameterization
  // (could add more than 1 if we wanted additional covariates as fn of b_y)
  int K_qr = Kc + 1;

  vector[N] log_price = log(price);
  vector[N] price_times_T = price * time_period_hr;

  // Additional params, sometimes useful in the generated block:
  // empty arrays and vec so we can match function signatures exactly.
  array[0] real empty_real;
  array[0] int empty_int;
  vector[0] empty_vec;

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
  // real<lower=-1, upper=0> temp_inv_alpha; // cost param (1 / alpha)

  real temp_alpha_inter;
  vector[Kc] b_alpha;
}
transformed parameters {
}

generated quantities {
  // Note: using the lognormal_rng here is a higher-variance choice than
  // using the expected value. Shouldn't matter for point values, but might
  // matter for CI calc.

  vector[N] leak_size_expect;
  vector[N] leak_size_draw;

  // Note that the prob_leak output is not actually used for this model in
  // outcomes_analysis.py -- it was for a previous model -- but keeping it around
  // is easier than changing all the code.
  vector[N] prob_leak;
  vector[N] prob_size_above_threshold;

  // These are the A_i and alpha in the stan model, but adding more descriptive
  // names here for saving.
  vector[N] cost_param_A;
  vector[N] cost_param_alpha;

  {
    // Inner block here to avoid printing temporary vars.
    // Most of this math is the same as in 02_twopart_lognormal_alpha_model.stan
    // See that file for details.
    vector[N_obs] mu_y_obs = temp_y_inter + Xc_obs * b_y;
    vector[N] mu_y_all = temp_y_inter + Xc * b_y;

    leak_size_draw = to_vector(lognormal_rng(mu_y_all, sigma_y)) + shift_amount;

    real smear_y = mean(exp(y_shift_log - mu_y_obs));
    // Calculate the expectation of e * p, re-adding the shift_amount
    leak_size_expect = exp(mu_y_all) * smear_y + shift_amount;
    vector[N] leak_private_value_expect = price_times_T .* leak_size_expect;
    cost_param_A = inv_logit(temp_obs_inter + Xc * b_obs) .* leak_private_value_expect;

    // There's a hassle where, if leak_size_draw is much smaller than
    // leak_size_expect, the probability will be malformed (because A_i will be
    // outside the acceptable range). This happens quite rarely (about 0.00058%).
    int count_replaced = 0;
    int replaced_draw_bool = 0;
    for (i in 1:N) {
      replaced_draw_bool = 0;
      while (cost_param_A[i] > leak_size_draw[i] * price_times_T[i]) {
        replaced_draw_bool = 1;
        leak_size_draw[i] = lognormal_rng(mu_y_all[i], sigma_y) + shift_amount;
      }
      count_replaced = count_replaced + replaced_draw_bool;
    }
    if (count_replaced > 0) {
      // Important: we thin out the MCMC draws after this, so the total number
      // of values printed can be greater than the value that actually matter.
      print("  INFO: replaced leak_size_draw well pad count: ", count_replaced);
    }

    vector[N] inv_alpha = -inv_logit(temp_alpha_inter + Xc * b_alpha);
    cost_param_alpha = inv(inv_alpha);
    // See note above -- prob_leak is not used.
    prob_leak = (leak_private_value_expect ./ cost_param_A) ^ inv_alpha;

    // Note: `detect_threshold_shifted` is a single real number. For different
    // detect_threshold_shifted values, create different variables.
    for (i in 1:N) {
      prob_size_above_threshold[i] = 1 - lognormal_cdf(
        detect_threshold_shifted | mu_y_all[i], sigma_y
      );
    }
  }
}
