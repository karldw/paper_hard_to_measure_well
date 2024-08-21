suppressMessages(
  here::i_am("code/policy_audit_gains_rel_plot.R", uuid="2547096c-b358-4e82-a172-12fb98006241")
)

source(here::here("code/shared_functions.r"))


if (!exists("snakemake")) {
  snakemake <- SnakemakePlaceholder(
    input = list(
       outcome_summaries = here::here(
          "data/generated/policy_outcomes/main_spec/08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours",
          c(
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=high-1week.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=low-1week.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=med-1week.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=med-3month.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=2500.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=7500.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=high-1week.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=low-1week.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=med-1week.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=med-3month.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=2500.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=7500.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=high-1week.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=low-1week.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=med-1week.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=med-3month.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=2500.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=7500.parquet"
          )
        ),
      r_lib = "scratch/snakemake_flags/setup_r_library",
      script = "code/policy_audit_gains_rel_plot.R",
      constants = "code/constants.json",
      policy_output_helper_functions = "code/policy_output_helper_functions.r"
    ),
    output = list(
      "graphics/audit_gains_rel_plot_dwl_frac=1pct.pdf"
    ),
    wildcards = list(
      audit_amount = "1pct",
      out_measure = "dwl",
      simple = ""
    ),
    threads = 1
  )
}

source(snakemake@input[["policy_output_helper_functions"]])
options(scipen = 99, mc.cores=snakemake@threads)


plot_rel_outcome <- function(df, outcome_var, outfile, fee_labels = c("detailed", "simple")) {
  stopifnot(
    length(outfile) == 1,
    length(outcome_var) == 1,
    nrow(df) > 0,
    # Labels below assume a year
    df$time_H == 8760
  )
  df %<>% dplyr::ungroup()
  df %<>% dplyr::filter(detect_threshold > 0 | !audit_rule %in% c("target_e", "remote"))
  should_be_empty <- dplyr::count(df, audit_rule, tau_T_str) %>% dplyr::filter(n > 1)
  if (nrow(should_be_empty) > 0) {
    print(should_be_empty)
    stop("Unexpected combos above")
  }
  fee_labels <- match.arg(fee_labels)
  if (outcome_var == "dwl") {
    df %<>% dplyr::rename(
      estimate  = dwl_tot_rel_pct_mean,
      conf_low  = dwl_tot_rel_pct_conf_low,
      conf_high = dwl_tot_rel_pct_conf_high,
    )
    subtitle <- "Percent DWL improvement"
    ylims <- c(-0.01, 100.01)
  } else if (outcome_var == "emis") {
    df %<>% dplyr::rename(
      estimate  = emis_tot_rel_pct_mean,
      conf_low  = emis_tot_rel_pct_conf_low,
      conf_high = emis_tot_rel_pct_conf_high
    )
    subtitle <- "Percent emissions improvement"
    ylims <- c(-0.01, 100.01)
  } else if (outcome_var == "fee_per_kg_p90") {
    # Plotting median here; could do other moments too.
    df %<>% dplyr::mutate(
      estimate  = fee_per_kg_p90_mean / SCM_PER_KG * 100,
      conf_low  = fee_per_kg_p90_conf_low / SCM_PER_KG * 100,
      conf_high = fee_per_kg_p90_conf_high / SCM_PER_KG * 100,
    )
    subtitle <- "90th percentile expected fee (% of external cost)"
    ylims <- c(-0.01, 100.01)
  } else if (outcome_var == "fee_per_kg_med") {
    # Plotting median here; could do other moments too.
    df %<>% dplyr::mutate(
      estimate  = fee_per_kg_med_mean / SCM_PER_KG * 100,
      conf_low  = fee_per_kg_med_conf_low / SCM_PER_KG * 100,
      conf_high = fee_per_kg_med_conf_high / SCM_PER_KG * 100,
    )
    subtitle <- "Median expected fee (% of external cost)"
    ylims <- c(-0.01, 100.01)
  } else if (outcome_var == "fee_per_kg_mean") {
    # Plotting median here; could do other moments too.
    df %<>% dplyr::mutate(
      estimate  = fee_per_kg_mean_mean / SCM_PER_KG * 100,
      conf_low  = fee_per_kg_mean_conf_low / SCM_PER_KG * 100,
      conf_high = fee_per_kg_mean_conf_high / SCM_PER_KG * 100,
    )
    subtitle <- "Mean expected fee (% of external cost)"
    ylims <- c(-0.01, 100.01)
  } else {
    stop("Add a case for ", outcome_var)
  }

  df %<>% dplyr::mutate(
    policy = factor(audit_rule,
      levels = c("none", "uniform", "target_x", "target_e", "remote"),
      labels = c("None", "No targeting", "Target on well\ncharacteristics", "Target on\nremote sensing", "Remote only")
    )
  )
  count_tau_T_levels <- dplyr::n_distinct(df$tau_T)
  stopifnot(count_tau_T_levels > 0)
  relative_label_list <- list(
    "",
    c("Low", "High"),
    c("Low", "Med", "High"),
    c("Low", "Med-low", "Med-high", "High"),
    c("Low", "Med-low", "Med", "Med-high", "High"),
    c("V Low", "Low", "Med-low", "Med", "Med-high", "High"),
    c("V Low", "Low", "Med-low", "Med", "Med-high", "High", "V High")
  )
  if (count_tau_T_levels > length(relative_label_list)) {
    stop("Sorry, don't have labels for ", count_tau_T_levels, " different levels")
  }
  relative_labels <- relative_label_list[[count_tau_T_levels]]
  # Make a small df here to make the labeling easier.
  fee_label_df <- dplyr::distinct(df, tau_T_str, tau_T_str_fct, tau_T) %>%
    dplyr::arrange(tau_T) %>%
    dplyr::mutate(
      relative_label = relative_label_list[[count_tau_T_levels]],
    )

  if (fee_labels == "detailed") {
    fee_label_df %<>% dplyr::mutate(tau_T_fct = tau_T_str_fct)
  } else {
    fee_label_df %<>% dplyr::mutate(
      tau_T_fct = factor(tau_T, levels=sort(unique(tau_T)), labels=relative_label)
    )
  }
  fee_label_df %<>% dplyr::select(tau_T_str, tau_T_fct)
  # Check that all rows match
  df %<>% powerjoin::power_inner_join(
    fee_label_df, by="tau_T_str",
    check=merge_specs(duplicate_keys_left = "ignore")
  )
  fee_label <- "Fee"

  if (!all(df$audit_frac %in% c(0, 0.01))) {
    stop("The label in the graph below will be wrong")
  }
  plt_bar <- df %>%
    ggplot2::ggplot(ggplot2::aes(policy, estimate, ymin=conf_low, ymax=conf_high, fill=tau_T_fct)) +
    ggplot2::geom_col(position="dodge") +
    ggplot2::scale_fill_brewer(palette="Blues") +
    ggplot2::theme_bw() +
    ggplot2::labs(x="", y="", fill=fee_label,
      subtitle=paste(subtitle, "auditing 1% of wells per year")
    ) +
    ggplot2::ylim(ylims[1], ylims[2])
  if (fee_label == "detailed") {
    plt <- plt + ggplot2::geom_errorbar(width=0.2, position=ggplot2::position_dodge(width=0.9))
  }
  save_plot(plt_bar, outfile, scale_mult=1.2)
}


make_rel_outcome_barplots <- function(df, plot_file_list) {
  stopifnot(length(plot_file_list) >= 1)
  # This could be improved by doing it on demand, but for now just make all the plots
  for (out_measure in c("dwl", "emis", "fee_per_kg_p90", "fee_per_kg_med")) {
    local_path <- glue::glue("graphics/audit_gains_rel_plot_{out_measure}_frac=1pct.pdf")
    filename <- here::here(local_path)
    message(local_path)
    to_plot <- dplyr::filter(df,
      !audit_rule %in% c("none", "remote"),
      tau_T_str %in% c(
        "low-1week",
        "med-1week", "med-3month",
        "high-1week", "high-3month",
        "2500", "7500"
      ),
      audit_frac == 0.01 | audit_rule == "remote",
    )
    plot_rel_outcome(to_plot, out_measure, filename)
  }
  # For the blog:

  filename <- here::here("graphics/audit_gains_rel_plot_emis_frac=1pct_simple.pdf")
  message(local_path)
  # This is a kludge, but for this plot only, rescale the relative outcome so
  # it's the percent reduction below BAU emissions, rather than the percent of
  # the way to first best.
  # Do this by multiplying by ~0.81, the emissions reduction at first best.
  # (The standard errors will be slightly wrong.)
  emis_bau <- dplyr::filter(df, audit_rule == "none")$emis_tot_mean[1]

  to_plot <- dplyr::filter(df,
    !audit_rule %in% c("none", "remote"),
    # small set of intermediate fee levels:
    tau_T_str %in% c(
      "med-3month",
      "high-3month",
      "med-1week"
    ),
    audit_frac == 0.01 | audit_rule == "remote",
  ) %>% dplyr::mutate(
    emis_tot_rel_pct_mean = 100 * emis_reduce_tot_mean / !!emis_bau,
    emis_tot_rel_pct_conf_low = 100 * emis_reduce_tot_conf_low / !!emis_bau,
    emis_tot_rel_pct_conf_high = 100 * emis_reduce_tot_conf_high / !!emis_bau,
  )
  message("  fee values: tau * T = {δ * 1wk, δ * 3mo, δ * 6mo}")
  plot_rel_outcome(to_plot, "emis", filename, fee_labels="simple")

  invisible(NULL)
}


SCM_PER_KG <- read_constants()$SOCIAL_COST_METHANE_PER_KG

if (isTRUE(snakemake@wildcards$simple == "_simple")) {
  stop("Simple version of the graph needs to be repaired before it can be used.")
}

all_res_no_dup <- read_policy_summaries(
  snakemake@input[["outcome_summaries"]],
  fill_redundant=FALSE
)
make_rel_outcome_barplots(all_res_no_dup, snakemake@output)
