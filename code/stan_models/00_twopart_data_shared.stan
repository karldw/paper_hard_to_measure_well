
// 00_twopart_data_shared.stan
// data and transformed data blocks, meant to be included in other programs.
// (this file can be compiled, but won't run)
functions {
#include 00_shared_functions.stan
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
  int<lower=1, upper=N> boot_idx[N];
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
  int obs_dummy[N];
  for (i in 1:N) {
    obs_dummy[i] = Y_reindex[i] > shift_amount;
  }
  int N_obs = sum(obs_dummy); // number of observed y, after shifting

  if (N_obs < 2 || N_obs > (N - 2)) {
    reject("Can't estimate model");
  }
  int<lower=1> obs_idx[N_obs] = which(obs_dummy); // index of observed y
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
  real shift_amount_times_price_min_times_T = shift_amount * min(price) * time_period_hr;

  // Additional params, sometimes useful in the generated block:
  // empty arrays and vec so we can match function signatures exactly.
  real empty_real[0];
  int empty_int[0];
  vector[0] empty_vec;

  // Note: hard-coding detect_threshold for now. Needs to match detect_threshold
  // in the outcomes_analysis.py program.
  real<lower=0.01> detect_threshold_shifted = 100 - shift_amount;
}
