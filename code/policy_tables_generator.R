
suppressMessages(
  here::i_am("code/policy_tables_generator.R", uuid="7f35bda8-6c07-4a60-9636-5b37caa6d6c0")
)

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
    dplyr::select(df,
      term = term,
      estimate  = fee_per_kg_mean_mean,
      conf_low  = fee_per_kg_mean_conf_low,
      conf_high = fee_per_kg_mean_conf_high
    ),
    dplyr::select(df,
      term = term,
      estimate  = fee_per_kg_med_mean,
      conf_low  = fee_per_kg_med_conf_low,
      conf_high = fee_per_kg_med_conf_high
    ),
    dplyr::select(df,
      term = term,
      estimate  = fee_per_kg_p25_mean,
      conf_low  = fee_per_kg_p25_conf_low,
      conf_high = fee_per_kg_p25_conf_high
    ),
    dplyr::select(df,
      term = term,
      estimate  = fee_per_kg_p75_mean,
      conf_low  = fee_per_kg_p75_conf_low,
      conf_high = fee_per_kg_p75_conf_high
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
  comment_names <- "fee mean, fee median, fee p25, fee p75"
  audit_res_table %>%
    make_table_fragment(escaped = FALSE, add_comments=comment_names) %>%
    writeLines(outfile)
}


make_one_table_fee_fixed_pct <- function(results_df, wildcards, outfile) {
  # Assemble the expected fee table, with a name like
  # "tableA04_expected_fee_10pct_low_tau.tex". Read the desired percentage from
  # the filename.
  # Columns are mean, median, p25 and p75 audit frac, with CI for each.
  # The point estimates are the mean across draws of the statistic, so the
  # variable name for the mean is `fee_per_kg_mean_mean`, with CI
  # `fee_per_kg_mean_conf_low` and `fee_per_kg_mean_conf_high`, and so on.
  message(outfile)
  desired_pct <- audit_amount_to_number(wildcards$audit_amount, type="fixed")
  tau_level <- tidy_gsub(wildcards$audit_tau, "\\-$", "")
  T_level <- wildcards$audit_T
  results_df %>% dplyr::filter(
      approx_equal(audit_frac * 100, desired_pct) | audit_rule %in% c("none", "remote"),
      !audit_rule %in% c("none"),
      tau_str == tau_level,
      T_str == T_level,
      audit_rule %in% c("target_e", "remote") | detect_threshold == 0,
    ) %>%
    fee_table_shared(outfile)
  invisible(outfile)
}

make_one_table_dwl_emis_fixed <- function(results_df, wildcards, outfile) {
  desired_pct <- audit_amount_to_number(wildcards$audit_amount, type="fixed")
  tau_level <- tidy_gsub(wildcards$audit_tau, "\\-$", "")
  T_level <- wildcards$audit_T

  df <- results_df %>%
    dplyr::filter(
      approx_equal(audit_frac * 100, desired_pct) | audit_rule %in% c("none", "remote"),
      !audit_rule %in% c("none"),
      tau_str == tau_level,
      T_str == T_level,
      audit_rule %in% c("target_e", "remote") | detect_threshold == 0,
    ) %>%
    dplyr::mutate(
      term = format_policy_details(.data$audit_rule, .data$detect_threshold),
    ) %>%
    dplyr::arrange(term) %>% # use factor levels for sorting
    dplyr::mutate(term = as.character(term))
  lst <- list(
    dplyr::transmute(df,
      term = term,
      estimate = signif3(dwl_tot_rel_pct_mean),
      conf_low = signif2(dwl_tot_rel_pct_conf_low),
      conf_high = signif2(dwl_tot_rel_pct_conf_high),
    ),
    dplyr::transmute(df,
      term = term,
      estimate = signif3(emis_tot_rel_pct_mean),
      conf_low = signif2(emis_tot_rel_pct_conf_low),
      conf_high = signif2(emis_tot_rel_pct_conf_high),
    )
  )

  comment_names <- "dwl, emis"
  tab <- lst %>%
    purrr::map(format_estimate_above_interval, align="@") %>%
    merge_estimates_df()  %>%
    make_table_fragment(escaped = FALSE, add_comments=comment_names)

  writeLines(tab, outfile)
  invisible(outfile)
}


results_to_table <- function(results_df, wildcards, table_file) {
  stopifnot(length(table_file) == 1)
  outcome_type <- wildcards$outcome_type %||% stop("missing outcome_type")
  amount_type <- ifelse(endsWith(wildcards$audit_amount, "pct"), "fixed", "optimal")
  if (amount_type == "optimal") {
    stop("no longer used")
  }
  if (outcome_type == "expected_fee" && amount_type == "fixed") {
    make_one_table_fee_fixed_pct(results_df, wildcards, table_file)
  } else if (outcome_type == "outcome_dwl_emis" && amount_type == "fixed") {
    make_one_table_dwl_emis_fixed(results_df, wildcards, table_file)
  } else {
    print(wildcards)
    stop("No table match.")
  }
}

SCM_PER_KG <- read_constants()$SOCIAL_COST_METHANE_PER_KG

all_res <- read_policy_summaries(
  snakemake@input[["outcome_summaries"]],
  fill_redundant=TRUE
)
results_to_table(all_res, snakemake@wildcards, snakemake@output$tables)
