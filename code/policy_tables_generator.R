
source(here::here("code/shared_functions.r"))
source(snakemake@input[["policy_output_helper_functions"]])
options(scipen = 99, mc.cores=snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])


fee_table_shared <- function(df, outfile) {
  df %<>% dplyr::mutate(
      term = format_policy_details(.data$audit_rule, .data$detect_threshold),
    ) %>%
    dplyr::arrange(term) %>% # use factor levels for sorting
    dplyr::mutate(term = as.character(term))
  stopifnot(nrow(df) > 0, !anyDuplicated(df$term))
  audit_res_table <- list(
    # Rename each and process individually because the format_estimate_above_interval
    # code expects columns named term, estimate, conf_low, and conf_high
    dplyr::transmute(df,
      term = term,
      estimate  = fee_per_kg_mean_mean,
      conf_low  = fee_per_kg_mean_conf_low,
      conf_high = fee_per_kg_mean_conf_high
    ),
    dplyr::transmute(df,
      term = term,
      estimate  = fee_per_kg_med_mean,
      conf_low  = fee_per_kg_med_conf_low,
      conf_high = fee_per_kg_med_conf_high
    ),
    dplyr::transmute(df, term = term,
      estimate  = fee_per_kg_p10_mean,
      conf_low  = fee_per_kg_p10_conf_low,
      conf_high = fee_per_kg_p10_conf_high
    ),
    dplyr::transmute(df, term = term,
      estimate  = fee_per_kg_p90_mean,
      conf_low  = fee_per_kg_p90_conf_low,
      conf_high = fee_per_kg_p90_conf_high
    )
  ) %>%
    purrr::map(~dplyr::mutate(.,
      estimate  = signif3(100 * estimate / SCM_PER_KG),
      conf_low  = signif2(100 * conf_low / SCM_PER_KG),
      conf_high = signif2(100 * conf_high / SCM_PER_KG),
    )) %>%
    purrr::map(format_estimate_above_interval, align="@") %>%
    merge_estimates_df()
  # These are added in the latex comments to help catch errors.
  comment_names <- "fee mean, fee median, fee p10, fee p90"
  audit_res_table %>%
    make_table_fragment(escaped = FALSE, add_comments=comment_names) %>%
    writeLines(outfile)
}


make_one_table_fee_fixed_pct <- function(results_df, wildcards, outfile) {
  # Assemble the expected fee table, with a name like
  # "expected_fee_10pct_low_tau.tex". Read the desired percentage from the
  # filename.
  # Columns are mean, median, p10 and p90 audit frac, with CI for each.
  # The point estimates are the mean across draws of the statistic, so the
  # variable name for the mean is `fee_per_kg_mean_mean`, with CI
  # `fee_per_kg_mean_conf_low` and `fee_per_kg_mean_conf_high`, and so on.
  message(outfile)
  desired_pct <- audit_amount_to_number(wildcards$audit_amount, type="fixed")
  tau_level <- tidy_gsub(wildcards$audit_tau, "\\-$", "")
  T_level <- wildcards$audit_T

  results_df %>% dplyr::filter(
      approx_equal(audit_frac * 100, desired_pct),
      !audit_rule %in% c("remote", "none"),
      tau == TAU_LEVELS[[tau_level]],
      time_T == T_LEVELS[[T_level]],
      audit_rule == "target_e" | detect_threshold == 0,
    ) %>%
    fee_table_shared(outfile)
  outfile
}


make_one_table_fee_optim_pct <- function(results_df, wildcards, outfile) {
  message(outfile)
  dollar_amt <- audit_amount_to_number(wildcards$audit_amount, type="optimal")
  tau_level <- tidy_gsub(wildcards$audit_tau, "\\-$", "")
  T_level <- wildcards$audit_T

  df <- results_df %>% dplyr::filter(
      approx_equal(audit_cost, dollar_amt),
      !audit_rule %in% c("remote", "none"),
      tau == TAU_LEVELS[[tau_level]],
      time_T == T_LEVELS[[T_level]],
      audit_rule == "target_e" | detect_threshold == 0,
    )
  if (nrow(df) == 0) {
    msg <- glue::glue("No results with this combination of values: tau={tau_level}, T={T_level}, audit_cost={dollar_amt}")
    stop(msg)
  }
  fee_table_shared(df, outfile)
  outfile
}


make_one_table_dwl_fixed <- function(results_df, wildcards, outfile) {
  # Note that both this case and the optim case already include audit costs
  desired_pct <- audit_amount_to_number(wildcards$audit_amount, type="fixed")
  tau_level <- tidy_gsub(wildcards$audit_tau, "\\-$", "")
  T_level <- wildcards$audit_T
  stopifnot(tau_level == "all")
  df <- dplyr::filter(results_df,
      approx_equal(audit_frac * 100, desired_pct),
    ) %>%
    dplyr::rename(
      estimate = dwl_tot_rel_pct_mean,
      conf_low = dwl_tot_rel_pct_conf_low,
      conf_high = dwl_tot_rel_pct_conf_high,
    )
  make_one_table_dwl_emis_shared(df, outfile, T_level)
}


make_one_table_dwl_optim <- function(results_df, wildcards, outfile) {
  dollar_amt <- audit_amount_to_number(wildcards$audit_amount, type="optimal")
  tau_level <- tidy_gsub(wildcards$audit_tau, "\\-$", "")
  T_level <- wildcards$audit_T
  stopifnot(tau_level == "all")
  df <- dplyr::filter(results_df,
      approx_equal(audit_cost, dollar_amt)
    ) %>%
    dplyr::rename(
      estimate = dwl_tot_rel_pct_mean,
      conf_low = dwl_tot_rel_pct_conf_low,
      conf_high = dwl_tot_rel_pct_conf_high,
    )
  make_one_table_dwl_emis_shared(df, outfile, T_level)
}


make_one_table_emis_fixed <- function(results_df, wildcards, outfile) {
  desired_pct <- audit_amount_to_number(wildcards$audit_amount, type="fixed")
  tau_level <- tidy_gsub(wildcards$audit_tau, "\\-$", "")
  T_level <- wildcards$audit_T
  stopifnot(tau_level == "all")
  df <- dplyr::filter(results_df,
      approx_equal(audit_frac * 100, desired_pct),
    ) %>%
    dplyr::rename(
      estimate = emis_tot_rel_pct_mean,
      conf_low = emis_tot_rel_pct_conf_low,
      conf_high = emis_tot_rel_pct_conf_high,
    )
  make_one_table_dwl_emis_shared(df, outfile, T_level)
}


make_one_table_emis_optim <- function(results_df, wildcards, outfile) {
  dollar_amt <- audit_amount_to_number(wildcards$audit_amount, type="optimal")
  tau_level <- tidy_gsub(wildcards$audit_tau, "\\-$", "")
  T_level <- wildcards$audit_T
  stopifnot(tau_level == "all")
  df <- dplyr::filter(results_df,
      approx_equal(audit_cost, dollar_amt)
    ) %>%
    dplyr::rename(
      estimate = emis_tot_rel_pct_mean,
      conf_low = emis_tot_rel_pct_conf_low,
      conf_high = emis_tot_rel_pct_conf_high,
    )
  make_one_table_dwl_emis_shared(df, outfile, T_level)
}


make_one_table_dwl_emis_shared <- function(df, outfile, T_level) {
  message(outfile)
  to_table <- df %>%
    dplyr::filter(
      audit_rule != "none",
      audit_rule %in% c("remote", "target_e") | detect_threshold == 0,
      time_T == T_LEVELS[[T_level]],
    ) %>%
    dplyr::mutate(term = format_policy_details(audit_rule, detect_threshold)) %>%
    # Sort by factor level
    dplyr::arrange(term) %>%
    dplyr::mutate(term = as.character(term))
  stopifnot(nrow(df) > 3, nrow(to_table) > 3, all(TAU_LEVELS %in% df$tau))


  # Columns (in order below):
  # tau high, tau med, tau low
  # Actual outcome variable depends on what's in df -- see the functions that
  # call this one.
  tab <- list(
    dplyr::select(
      dplyr::filter(to_table, tau == TAU_LEVELS[["high"]]),
      term, estimate, conf_low, conf_high,
    ),
    dplyr::select(
      dplyr::filter(to_table, tau == TAU_LEVELS[["med"]]),
      term, estimate, conf_low, conf_high,
    ),
    dplyr::select(
      dplyr::filter(to_table, tau == TAU_LEVELS[["low"]]),
      term, estimate, conf_low, conf_high,
    )
  ) %>%
  purrr::map(~dplyr::mutate(.,
    estimate  = signif3(estimate),
    conf_low  = signif3(conf_low),
    conf_high = signif3(conf_high),
  )) %>%
  purrr::map(format_estimate_above_interval, align="@") %>%
  merge_estimates_df() %>%
  make_table_fragment(escaped = FALSE,
    add_comments=glue::glue("high tau, med tau, low tau (all using T={T_level})")
  )

  writeLines(tab, outfile)
  outfile
}


results_to_table <- function(results_df, wildcards, table_file) {
  stopifnot(length(table_file) == 1)
  table_type <- wildcards$table_type %||% stop("missing table_type")
  amount_type <- ifelse(endsWith(wildcards$audit_amount, "pct"), "fixed", "optimal")
  if (table_type == "expected_fee" && amount_type == "fixed") {
    make_one_table_fee_fixed_pct(results_df, wildcards, table_file)
  } else if (table_type == "expected_fee" && amount_type == "optimal") {
    make_one_table_fee_optim_pct(results_df, wildcards, table_file)
  } else if (table_type == "outcome_dwl" && amount_type == "fixed") {
    make_one_table_dwl_fixed(results_df, wildcards, table_file)
  } else if (table_type == "outcome_dwl" && amount_type == "optimal") {
    make_one_table_dwl_optim(results_df, wildcards, table_file)
  } else if (table_type == "outcome_emis" && amount_type == "fixed") {
    make_one_table_emis_fixed(results_df, wildcards, table_file)
  } else if (table_type == "outcome_emis" && amount_type == "optimal") {
    make_one_table_emis_optim(results_df, wildcards, table_file)
  } else {
    print(wildcards)
    stop("No table match.")
  }
}

const <- read_constants()
TAU_LEVELS <- const$TAU_LEVELS
T_LEVELS <- const$T_LEVELS
SCM_PER_KG <- const$SOCIAL_COST_METHANE_PER_KG

all_res <- read_policy_summaries(
  snakemake@input[["outcome_summaries"]],
  fill_redundant=TRUE
)
results_to_table(all_res, snakemake@wildcards, snakemake@output$tables)
