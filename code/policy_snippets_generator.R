suppressMessages(
  here::i_am("code/policy_snippets_generator.R", uuid="380db8d7-171a-4058-983a-81e0df94a7f0")
)

source(here::here("code/shared_functions.r"))
source(snakemake@input[["policy_output_helper_functions"]])
options(scipen = 99, mc.cores=snakemake@threads)
# Note: memory limiting doesn't play well, even though actual memory
# consumption isn't that high. Maybe something with arrow memory allocation
if (FALSE) {
  memory_limit(snakemake@resources[["mem_mb"]])
}

results_to_snippets <- function(results_df, snippet_files, wildcards) {
  stopifnot(length(snippet_files) >= 1)
  output_strings <- purrr::map_chr(snippet_files, make_one_snippet, df=results_df, wildcards=wildcards)
  purrr::walk2(output_strings, snippet_files, writeLines)
  invisible(output_strings)
}


parse_inputs <- function(filename, wildcards) {
  filename <- basename(filename)
  stopifnot(
    length(filename) == 1,
    startsWith(filename, "intext_OUTCOME"),
    endsWith(filename, ".tex"),
    nrow(df) >= 1
  )
  message(filename)

  match_outcome <- wildcards$outcome %||% stop("missing outcome")
  match_rule <- wildcards$audit_rule %||% stop("missing audit_rule")
  # tauT can be in two formats: a string like "low-1week" or a number like "0.3"
  # Support both here.
  match_tauT <- wildcards$audit_tauT %||% stop("missing audit_tauT")
  number_regex <- "^\\d+\\.?\\d*$" # numbers (basic, and decimals require 0)
  if (isTRUE(grepl(number_regex, match_tauT, perl = TRUE))) {
    match_tauT_num <- as.numeric(match_tauT)
  } else {
    match_tauT_num <- make_tau_T() %>%
      dplyr::filter(tau_T_str == !!match_tauT) %>%
      purrr::chuck("tau_T")
  }
  stopifnot(length(match_tauT_num) == 1, !anyNA(match_tauT_num))

  if (endsWith(match_rule, "_high")) {
    match_threshold <- const$POLICY_DETECT_THRESHOLD_HIGH
  } else {
    match_threshold <- const$POLICY_DETECT_THRESHOLD_LOW
  }

  outcome_compare_across_rules <- c(
    "welfare_pct_of_target_e_low",
    "welfare_gain_pct_vs_uniform"
  )
  input_param <- list(
    filename = filename,
    match_outcome = match_outcome,
    audit_rule = tidy_gsub(match_rule, "(_low|_high)$", ""),
    match_threshold = match_threshold,
    match_tauT_num = match_tauT_num,
    compare_across_rules = match_outcome %in% outcome_compare_across_rules,
    # Note: optimal-* for could be supported here, but currently isn't.
    match_pct = audit_amount_to_number(wildcards$audit_amount, type="fixed"),
    match_CI = wildcards$CI %||% stop("missing CI")
  )
  input_param
}



make_one_snippet <- function(filename, wildcards, df) {
  param <- parse_inputs(filename, wildcards)
   # This is a kludge, but we want to have a couple of snippets that compare
   # across audit rules.
  if (param$compare_across_rules) {
    return(make_one_snippet_compare_across_rules(param, df))
  }
  match_outcome <- param$match_outcome
  audit_rule <- param$audit_rule
  df_in_case_of_error <- dplyr::select(df, audit_rule, audit_frac, tau_T, detect_threshold, audit_cost)

  df %<>% dplyr::filter(
      audit_rule == !!audit_rule,
      approx_equal(audit_frac * 100, !!param$match_pct),
      approx_equal(tau_T, !!param$match_tauT_num),
      detect_threshold == !!param$match_threshold,
    ) %>%
    dplyr::mutate(dplyr::across(
      .cols=c(
        "annual_emis_mean_mean",
        "annual_emis_mean_conf_low",
        "annual_emis_mean_conf_high",
        "annual_emis_reduce_mean_mean",
        "annual_emis_reduce_mean_conf_low",
        "annual_emis_reduce_mean_conf_high"
      ),
      .fns = ~ . * METHANE_GWP / 1000,
      .names = "{.col}_tons_co2e"
    ))
  if (nrow(df) > 1) {
    print(df)
    stop("Found multiple matching result rows. See above.\n", filename)
  }
  if (nrow(df) == 0) {
    print(wildcards)
    print(df_in_case_of_error)
    stop("Failed to find any audit results to match. See above.\n", filename)
  }

  if (match_outcome == "welfare_gain_pct") {
    df %<>% dplyr::transmute(
      estimate = dwl_tot_rel_pct_mean,
      conf_low = dwl_tot_rel_pct_conf_low,
      conf_high = dwl_tot_rel_pct_conf_high,
      unit = "",
    )
  } else if (match_outcome == "emission_reduce_tonCO2e") {
    df %<>% dplyr::transmute(
      estimate = annual_emis_reduce_mean_mean_tons_co2e,
      conf_low = annual_emis_reduce_mean_conf_low_tons_co2e,
      conf_high = annual_emis_reduce_mean_conf_high_tons_co2e,
      unit = "~tons",
    )
  } else if (match_outcome == "emission_tonCO2e") {
    df %<>% dplyr::transmute(
      estimate = annual_emis_mean_mean_tons_co2e,
      conf_low = annual_emis_mean_conf_low_tons_co2e,
      conf_high = annual_emis_mean_conf_high_tons_co2e,
      unit = "~tons",
    )
  } else if (match_outcome == "emission_reduce_pct") {
    df %<>% dplyr::transmute(
      # Note that this is a shortcut -- would be more correct to calculate the
      # ratio for every bootstrap iter, then take the mean, rather than taking
      # the ratio of the means. Shouldn't matter too much for the means, would
      # definitely matter for the CI, but I don't actually use those at all.
      # Set them to zero here so the rest of the code works.
      # Also note that emissions pct change =
      # emissions reduced / (new emission total + emissions reduced)
      estimate = 100 * annual_emis_reduce_mean_mean_tons_co2e / (annual_emis_reduce_mean_mean_tons_co2e + annual_emis_mean_mean_tons_co2e),
      conf_low = 0,
      conf_high = 0,
      unit = "",
    )
  } else if (match_outcome == "fee_mean") {
    # Note: it's possible to add other fee measures by editing the outcomes_analysis.py file
    # Note: All of these measures also have CI available.
    # $2 / kg CH4 = $2000 / ton CH4 = $51.55 / ton CO2e
    df %<>% dplyr::transmute(
      estimate = fee_per_kg_mean_mean * 1000 / METHANE_GWP,
      conf_low = fee_per_kg_mean_conf_low * 1000 / METHANE_GWP,
      conf_high = fee_per_kg_mean_conf_high * 1000 / METHANE_GWP,
      unit = "",
    )
  } else if (grepl("^fee_p[0-9]+$", match_outcome, perl=TRUE)) {
    # Create names like fee_per_kg_p10_conf_high or fee_per_kg_p50_mean without
    # making a case for every one.
    pct_str <- tidy_gsub(match_outcome, "^fee_(p[0-9]+)$", "\\1")
    mean_var <- rlang::sym(glue::glue("fee_per_kg_{pct_str}_mean"))
    conf_low_var <- rlang::sym(glue::glue("fee_per_kg_{pct_str}_conf_low"))
    conf_high_var <- rlang::sym(glue::glue("fee_per_kg_{pct_str}_conf_high"))
    df %<>% dplyr::transmute(
      estimate  = !!mean_var      * 1000 / METHANE_GWP,
      conf_low  = !!conf_low_var  * 1000 / METHANE_GWP,
      conf_high = !!conf_high_var * 1000 / METHANE_GWP,
      unit = "",
    )
  } else if (match_outcome == "net_private_cost_per_mcf_pct_price") {
    df %<>% dplyr::transmute(
      estimate = net_private_cost_per_mcf_pct_price_mean,
      conf_low = net_private_cost_per_mcf_pct_price_conf_low,
      conf_high = net_private_cost_per_mcf_pct_price_conf_high,
      unit = "",
    )
  } else {
    stop("Unknown outcome type: ", match_outcome)
  }
  glue_df(df, param$match_CI)
}


glue_df <- function(df, match_CI) {
  df %<>% dplyr::mutate(
    estimate = signif3(estimate),
    conf_low = signif3(conf_low),
    conf_high = signif3(conf_high),
  )

  if (match_CI != "") {
    pattern <- "{estimate}{unit} \\mbox{{[{conf_low}, {conf_high}]}}%"
  } else {
    pattern <- "{estimate}{unit}%"
  }

  glue::glue(pattern, .envir = df)
}


# Load the data summary that corresponds to this outcome, but a different
# audit rule.
# NB: that this assumes the file exists, and is from the same model run.
load_counterpart <- function(snakemake, alt_audit_rule) {
  outcome_summary_filename <- snakemake@input[["outcome_summaries"]]
  stopifnot(length(outcome_summary_filename) == 1)
  alt_audit_rule <- match.arg(alt_audit_rule, c("none", "uniform", "target_x", "target_e_low", "target_e_high", "remote_low", "remote_high"))
  orig_audit_rule_wildcard <- snakemake@wildcards$audit_rule %||% stop("missing audit_rule")
  new_basename <- tidy_gsub(
    basename(outcome_summary_filename),
    paste0("rule=", orig_audit_rule_wildcard),
    paste0("rule=", alt_audit_rule),
    fixed=TRUE
  )
  filename <- file.path(dirname(outcome_summary_filename), new_basename)
  alt_df <- read_policy_summaries(filename, fill_redundant=TRUE)
  alt_df
}


make_one_snippet_compare_across_rules <- function(param, df) {
  # parse_inputs already validated the inputs above
  match_outcome <- param$match_outcome

  # welfare_gain_pct_vs_uniform reports the percentage difference from the
  # uniform policy DWL reduction. This code could be generalized to break apart
  # DWL and uniform, or to report percentage point differences, but we don't
  # need that.
  if (match_outcome == "welfare_gain_pct_vs_uniform") {
    df_compare <- load_counterpart(snakemake, "uniform")
    stopifnot(nrow(df) == 1, nrow(df_compare) == 1)
    dwl_compare <- df_compare$dwl_tot_rel_pct_mean
    df %<>% dplyr::transmute(
      estimate = 100 * (dwl_tot_rel_pct_mean - !!dwl_compare) / !!dwl_compare,
      conf_low = NA_real_,
      conf_high = NA_real_,
      unit = "",
    )
  } else if (match_outcome == "welfare_pct_of_target_e_low") {
    df_compare <- load_counterpart(snakemake, "target_e_low")
    stopifnot(nrow(df) == 1, nrow(df_compare) == 1)
    dwl_compare <- df_compare$dwl_tot_rel_pct_mean
    # Note the slightly different construction - we're asking how much worse
    # this policy is than target_e_low
    df %<>% dplyr::transmute(
      estimate = 100 * dwl_tot_rel_pct_mean / (!!dwl_compare),
      conf_low = NA_real_,
      conf_high = NA_real_,
      unit = "",
    )
    # Note: if you add cases here, also edit outcome_compare_across_rules in parse_inputs
  } else {
    # Should not be here
    stop("programming error")
  }
  glue_df(df, param$match_CI)
}


const <- read_constants()
METHANE_GWP <- const$METHANE_GWP


all_res <- read_policy_summaries(
  snakemake@input[["outcome_summaries"]],
  fill_redundant=TRUE
)
results_to_snippets(all_res, snakemake@output$snippets, wildcards=snakemake@wildcards)
