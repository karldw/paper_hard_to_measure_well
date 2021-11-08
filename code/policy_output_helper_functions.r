# Helper functions for number formatting:
signif2 <- function(x) {
  signif(zapsmall(x, 5), 2)
}

signif3 <- function(x) {
  signif(zapsmall(x, 6), 3)
}

fill_redundant_audit_rules <- function(df, time_H) {
  # For speed, we try not to estimate redundant audit rules (e.g. only estimate
  # the 'none' rule once, since the param don't matter), but it makes things
  # easier later if those results exist for all combinations (e.g. we might
  # want to compare outcomes for low-tau, but only have outcomes for high tau)
  stopifnot(noNAs(df))
  id_vars <- c("audit_frac", "audit_rule", "detect_threshold", "tau_T", "audit_cost")
  est_vars <- setdiff(names(df), id_vars)
  maybe_not_est_vars <- !grepl("_(mean|conf_low|conf_high)$", est_vars, perl=TRUE)
  if (any(maybe_not_est_vars)) {
    warning(
      "These variables were assumed to be estimates, but they have unexpected names: ",
      paste(est_vars[maybe_not_est_vars], collapse=", ")
    )
  }
  # This is basically the reverse of the paring-down done in make_audit_possibilities
  full_options <- tidyr::expand(df, !!!rlang::syms(id_vars)) %>%
    # Drop cases we don't ever want:
    dplyr::filter(
      audit_cost  > 0 | audit_frac  > 0 | audit_rule %in% c("none", "remote"),
      audit_cost == 0 | audit_frac == 0,
    )
  # Note: we do inner joins here. If values for some combinations are expected
  # but aren't estimated in outcomes_analysis.py, this function won't raise an
  # error.

  # For audit_rule == none, fill all possibilities
  df_none <- full_options %>%
    dplyr::filter(audit_rule == "none") %>%
    safejoin::safe_inner_join(
      dplyr::select(df, -audit_frac, -audit_cost, -tau_T, -detect_threshold),
      by="audit_rule",
      check="B C L T"
    )

  # For uniform and target_x:
  # Within groups of audit_rule, audit_frac, audit_cost, and tau_T, fill in
  # the estimates across values of detect_threshold
  # It won't be the cases that all rows match in full_options, or in df
  df_uniform_targetx <- full_options %>%
    dplyr::filter(audit_rule %in% c("uniform", "target_x")) %>%
    safejoin::safe_inner_join(
      dplyr::select(df, -detect_threshold),
      by=c("audit_rule", "audit_frac", "audit_cost", "tau_T"),
      check="B C L T"
    )
  # For remote:
  # Within groups of detect_threshold and tau_T, fill in the estimates across
  # values of audit_frac and audit_cost.
  df_remote <- full_options %>%
    dplyr::filter(audit_rule == "remote") %>%
    safejoin::safe_inner_join(
      dplyr::select(df, -audit_frac, -audit_cost),
      by=c("audit_rule","detect_threshold", "tau_T"),
      check="B C L T"
    )
  # For target_e:
  # No filling -- no values can be filled in across estimates.
  df_target_e <- dplyr::filter(df, audit_rule == "target_e")
  df_filled <- dplyr::bind_rows(df_none, df_uniform_targetx, df_remote, df_target_e)

  stopifnot(noNAs(df_filled))
  df_filled
}


annualize_estimates <- function(df) {
  time_H <- unique(df$time_H)
  stopifnot(length(time_H) == 1, time_H > 0)
  hours_per_year <- 8760 # Ignoring leapyears
  periods_per_year <- hours_per_year / time_H
  stopifnot(periods_per_year >= 0.99) # accomodate rounding

  # Here we annualize the columns that make sense to think of an annual total.
  # We don't annualize the fee_per_kg_* variables, shadow_price_* variables,
  # any of the *_rel_* variables (we could, but the annual version is identical
  # to the non annual version, since you'd need to multiply and divide by
  # `periods_per_year`).
  # Note that audit_frac_annual may be > 1, indicating the average well gets
  # audited more than once a year.
  # The annualized versions are created as new variables starting with 'annual_'
  out <- dplyr::mutate(df, dplyr::across(
    .cols = c(
    "audit_frac",
    "dwl_mean_mean",
    "dwl_mean_conf_low",
    "dwl_mean_conf_high",
    "dwl_tot_mean",
    "dwl_tot_conf_low",
    "dwl_tot_conf_high",
    "emis_mean_mean",
    "emis_mean_conf_low",
    "emis_mean_conf_high",
    "emis_reduce_mean_mean",
    "emis_reduce_mean_conf_low",
    "emis_reduce_mean_conf_high",
    "emis_tot_mean",
    "emis_tot_conf_low",
    "emis_tot_conf_high"),
    .fns = list(annual = ~ . * periods_per_year),
    .names = "{.fn}_{.col}"
  ))
  out
}


make_tau_T <- function() {
  # read_constants() is in shared_functions.r, reading constants.json
  const <- read_constants()
  df <- expand.grid(tau=const$TAU_LEVELS, time_T=const$T_LEVELS, KEEP.OUT.ATTRS=TRUE) %>%
    dplyr::mutate(tau_T = tau * time_T)
  stopifnot(rlang::is_named(df$tau), rlang::is_named(df$time_T))
  tau_T_str <- paste0(names(df$tau), "-", names(df$time_T))
  df <- dplyr::mutate(df, tau_T_str = tau_T_str)
  if (anyDuplicated(df$tau_T) != 0) {
    dups <- dplyr::group_by(df, tau_T) %>% dplyr::filter(dplyr::n() > 1) %>% dplyr::ungroup()
    print(dups)
    stop("Fix duplicates above:")
  }
  # Note: we check that all values match later.
  stopifnot(nrow(df) > 2)
  df
}


read_one_summary <- function(filename, fill_redundant) {
  stopifnot(length(filename) == 1)
  # 1. Parse filename for length
  # 2. Read parquet file.
  # 3. fill_redundant_audit_rules

  time_H <- filename_to_time_period_hr(filename)
  model_name <- dirname(filename) %>%
    file.path("model_fit.rds") %>%
    filename_to_model_name()
  df <- arrow::read_parquet(filename)
  if (fill_redundant) {
    # This fill_redundant_audit_rules was a pain to get right, so run it before
    # adding other variables.
    df %<>% fill_redundant_audit_rules(time_H)
  }
  imputed_N <- mean(df$dwl_tot_mean / df$dwl_mean_mean)
  if (!isTRUE(all.equal(imputed_N, read_constants()$N_WELL_PADS))) {
    stop(glue::glue("Number of well pads seems wrong {imputed_N} vs {N_WELL_PADS}"))
  }
  tau_T_df <- make_tau_T()
  df %<>% dplyr::mutate(time_H = !!time_H, model_name = !!model_name) %>%
    annualize_estimates() %>%
    # Merge with the tau_T info and check that it's a m:1 join with all LHS rows matching
    safejoin::safe_left_join(tau_T_df, by="tau_T", check="B C V L T")
  df
}

read_policy_summaries <- function(filenames, fill_redundant) {
  stopifnot(length(filenames) >= 1)
  purrr::map_dfr(filenames, read_one_summary, fill_redundant)
}


format_policy_details <- function(audit_rule, detect_threshold, latex=TRUE) {
  combo_labels <- c(
    none_low = "None",
    none_high = "None",
    uniform_low = "Uniform",
    uniform_high = "Uniform",
    target_x_low = "Target covariates",
    target_x_high = "Target covariates",
    target_e_low = "Target leaks, low threshold",
    target_e_high = "Target leaks, high threshold",
    remote_low = "Remote, low threshold",
    remote_high = "Remote, high threshold"
  )
  if (!latex) {
    combo_labels %<>% tidy_gsub("$", "", fixed=TRUE)
  }

  # Now bundle all of that into a factor so things sort correctly.
  combo <- paste(
    audit_rule,
    ifelse(detect_threshold > 0, "high", "low"),
    sep="_"
  )
  combo_levels <- names(combo_labels)
  combo_labels <- unname(combo_labels)
  if (!all(combo %in% combo_levels)) {
    print(combo)
  }
  stopifnot(
    length(audit_rule) == length(detect_threshold),
    combo %in% combo_levels,
    length(combo_levels) == length(combo_labels)
  )
  factor(combo, levels = combo_levels, labels = combo_labels)
}


audit_amount_to_number <- function(audit_amount_str, type=c("fixed", "optimal")) {
  type <- match.arg(type)
  stopifnot(length(audit_amount_str) > 0)
  if (type == "fixed") {
    number_part <- tidy_gsub(audit_amount_str, "pct$", "")
  } else {
    number_part <- tidy_gsub(audit_amount_str, "^optimal-([0-9+])usd$", "\\1")
  }
  num <- as.numeric(number_part)
  stopifnot(!anyNA(num))
  num
}
