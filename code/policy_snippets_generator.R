
source(here::here("code/shared_functions.r"))
source(snakemake@input[["policy_output_helper_functions"]])
options(scipen = 99, mc.cores=snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])

results_to_snippets <- function(results_df, snippet_files, wildcards) {
  stopifnot(length(snippet_files) >= 1)
  output_strings <- purrr::map_chr(snippet_files, make_one_snippet, df=results_df, wildcards=wildcards)
  purrr::walk2(output_strings, snippet_files, writeLines)
  invisible(output_strings)
}


make_one_snippet <- function(filename, wildcards, df) {
  filename <- basename(filename)
  stopifnot(
    length(filename) == 1,
    startsWith(filename, "OUTCOME"),
    endsWith(filename, ".tex"),
    nrow(df) >= 1
  )
  message(filename)
  df_orig <- df
  match_outcome <- wildcards$outcome %||% stop("missing outcome")
  match_rule <- wildcards$audit_rule %||% stop("missing audit_rule")
  match_tauT <- wildcards$audit_tauT %||% stop("missing audit_tauT")
  match_tauT_num <- make_tau_T() %>% dplyr::filter(tau_T_str == !!match_tauT) %>%
    purrr::chuck("tau_T")
  stopifnot(length(match_tauT_num) == 1)
  # Note: optimal-* for could be supported here, but currently isn't.
  match_pct <- audit_amount_to_number(wildcards$audit_amount, type="fixed")
  match_CI <- wildcards$CI %||% stop("missing CI")

  const <- read_constants()
  if (endsWith(match_rule, "_high") ) {
    match_threshold <- const$POLICY_DETECT_THRESHOLD_HIGH
  } else {
    match_threshold <- const$POLICY_DETECT_THRESHOLD_LOW
  }
  audit_rule <- tidy_gsub(match_rule, "(_low|_high)$", "")
  df %<>% dplyr::filter(
      audit_rule == !!audit_rule,
      approx_equal(audit_frac * 100, !!match_pct),
      approx_equal(tau_T, !!match_tauT_num),
      detect_threshold == !!match_threshold,
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
    print(dplyr::select(df_orig, audit_rule, audit_frac, tau_T, detect_threshold, audit_cost))
    stop("Failed to find any audit results to match. See above.\n", filename)
  }

  if (match_outcome == "welfare_gain_pct") {
    df %<>% dplyr::transmute(
      estimate = dwl_tot_rel_pct_mean,
      conf_low = dwl_tot_rel_pct_conf_low,
      conf_high = dwl_tot_rel_pct_conf_high,
      unit = "~percent",
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
  } else if (match_outcome == "fee_p10") {
    df %<>% dplyr::transmute(
      estimate = fee_per_kg_p10_mean * 1000 / METHANE_GWP,
      conf_low = fee_per_kg_p10_conf_low * 1000 / METHANE_GWP,
      conf_high = fee_per_kg_p10_conf_high * 1000 / METHANE_GWP,
      unit = "",
    )
  } else if (match_outcome == "fee_p90") {
    df %<>% dplyr::transmute(
      estimate = fee_per_kg_p90_mean * 1000 / METHANE_GWP,
      conf_low = fee_per_kg_p90_conf_low * 1000 / METHANE_GWP,
      conf_high = fee_per_kg_p90_conf_high * 1000 / METHANE_GWP,
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

  glue::glue(pattern, .envir=df)
}

const <- read_constants()
TAU_LEVELS <- const$TAU_LEVELS
T_LEVELS <- const$T_LEVELS
METHANE_GWP <- const$METHANE_GWP


all_res <- read_policy_summaries(
  snakemake@input[["outcome_summaries"]],
  fill_redundant=TRUE
)
results_to_snippets(all_res, snakemake@output$snippets, wildcards=snakemake@wildcards)
