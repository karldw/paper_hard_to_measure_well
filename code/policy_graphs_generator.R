suppressMessages(
  here::i_am("code/policy_graphs_generator.R", uuid="4bc31bae-8143-4ee3-93a1-d55c925caa19")
)
source(here::here("code/shared_functions.r"))

if (!exists("snakemake")) {
  snakemake <- SnakemakePlaceholder(
      input = list(
        outcome_summaries = here::here(
          "data/generated/policy_outcomes/main_spec/08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours",
          c(
            "audit_outcome_summary_rule=none_frac=0pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=remote_high_frac=0pct_tauT=high-1week.parquet",
            "audit_outcome_summary_rule=remote_high_frac=0pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=remote_high_frac=0pct_tauT=low-1week.parquet",
            "audit_outcome_summary_rule=remote_high_frac=0pct_tauT=low-3month.parquet",
            "audit_outcome_summary_rule=remote_high_frac=0pct_tauT=med-1week.parquet",
            "audit_outcome_summary_rule=remote_high_frac=0pct_tauT=med-3month.parquet",
            "audit_outcome_summary_rule=remote_high_frac=0pct_tauT=2500.parquet",
            "audit_outcome_summary_rule=remote_high_frac=0pct_tauT=7500.parquet",
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT=high-1week.parquet",
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT=low-1week.parquet",
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT=low-3month.parquet",
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT=med-1week.parquet",
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT=med-3month.parquet",
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT=2500.parquet",
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT=7500.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=high-1week.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=low-1week.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=low-3month.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=med-1week.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=med-3month.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=2500.parquet",
            "audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=7500.parquet",
            "audit_outcome_summary_rule=target_e_low_frac=1pct_tauT=high-1week.parquet",
            "audit_outcome_summary_rule=target_e_low_frac=1pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=target_e_low_frac=1pct_tauT=low-1week.parquet",
            "audit_outcome_summary_rule=target_e_low_frac=1pct_tauT=low-3month.parquet",
            "audit_outcome_summary_rule=target_e_low_frac=1pct_tauT=med-1week.parquet",
            "audit_outcome_summary_rule=target_e_low_frac=1pct_tauT=med-3month.parquet",
            "audit_outcome_summary_rule=target_e_low_frac=1pct_tauT=2500.parquet",
            "audit_outcome_summary_rule=target_e_low_frac=1pct_tauT=7500.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=high-1week.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=low-1week.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=low-3month.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=med-1week.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=med-3month.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=2500.parquet",
            "audit_outcome_summary_rule=target_x_frac=1pct_tauT=7500.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=high-1week.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=high-3month.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=low-1week.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=low-3month.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=med-1week.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=med-3month.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=2500.parquet",
            "audit_outcome_summary_rule=uniform_frac=1pct_tauT=7500.parquet"
          )
        ),
        r_lib = "scratch/snakemake_flags/setup_r_library",
        script = "code/policy_graphs_generator.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r"
      ),
      output = list(
        plot_fee = "graphics/outcomes_fee_frac=1pct.pdf",
        plot_dwl_emis_ordinal = "graphics/outcomes_dwl_emis_frac=1pct.pdf",
        plot_dwl_emis_cardinal = "graphics/outcomes_dwl_emis_frac=1pct_cardinal.pdf"
      ),
      wildcards = list(
        audit_amount = "1pct"
      ),
      threads = 1
    )
}

source(snakemake@input[["policy_output_helper_functions"]])
options(scipen = 99, mc.cores = snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])


format_audit_rules_variables <- function(results_df) {
  stopifnot(
    is.character(GRAPH_POLICY_ORDERING),
    length(GRAPH_POLICY_ORDERING) >= 1
  )
  check_n_unique_audit_frac(results_df)

  df <- results_df %>%
    dplyr::mutate(
      policy_details = factor(
        format_policy_details(audit_rule, detect_threshold),
        levels = !!GRAPH_POLICY_ORDERING
      ),
      uses_remote_info = factor(
        audit_rule %in% c("target_e", "remote"),
        levels = c(FALSE, TRUE),
        labels = c("Uniform and covariates", "Uses remote info")
      ),
    ) %>%
    dplyr::filter(.data$policy_details %in% !!GRAPH_POLICY_ORDERING)
  df
}

add_shared_plot_elements <- function(plt) {
  out <- plt +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_bw() +
    #ggplot2::scale_color_viridis_d(begin=0.2, end=0.95) +
    ggplot2::scale_color_brewer(palette = "Paired") +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::scale_shape_manual(values = c(0, 15, 5, 18, 1, 19)) +
    ggplot2::labs(
      shape = "",
      color = "",
      fill = "",
      size = "",
      x = "Fee when found leaking (τ × T)",
    )
  out
}


make_fee_graph <- function(results_df, wildcards, plot_file) {
  df_fee <- format_audit_rules_variables(results_df) %>%
    # These tau_T values were added for the DWL/emis plot, but are less helpful here.
    dplyr::filter(!tau_T_str %in% c("2500", "7500")) %>%
    dplyr::mutate(
      fee_interquartile = fee_per_kg_p75_mean - fee_per_kg_p25_mean,
      # Note: a more Tukey-standard way to plot this would be:
      # whisker_max = pmin(fee_per_kg_p99_mean, fee_per_kg_p75_mean + 1.5 * fee_interquartile),
      # whisker_min = pmax(fee_per_kg_p01_mean, fee_per_kg_p25_mean - 1.5 * fee_interquartile),

      # But that's not what people really expect, so instead plotting the p1 and p99
      whisker_max = fee_per_kg_p99_mean,
      whisker_min = fee_per_kg_p01_mean,
    )

  # For the fees, we're not showing confidence intervals (though we could).
  # Instead, we're showing summary stats on the wells.
  # Values are the mean across MCMC draws.
  plt <- df_fee %>%
    ggplot2::ggplot(ggplot2::aes(
      x = tau_T_str_fct,
      lower = fee_per_kg_p25_mean / !! SCM_PER_KG,
      middle = fee_per_kg_med_mean / !! SCM_PER_KG,
      upper = fee_per_kg_p75_mean / !! SCM_PER_KG,
      ymin = whisker_min / !! SCM_PER_KG,
      ymax = whisker_max / !! SCM_PER_KG,
      fill = policy_details,
      color = policy_details,
      shape = policy_details,
    )) %>%
    add_shared_plot_elements() +
    ggplot2::geom_point(
      ggplot2::aes(y=fee_per_kg_med_mean / !! SCM_PER_KG, size=policy_details),
      position = ggplot2::position_dodge2(width=0.9),
    ) +
    ggplot2::scale_size_manual(values = c(2.5, 2.5, 2.5, 3.5, 2.5, 2.5)) +
    ggplot2::geom_boxplot(
      stat="identity",
      # Note: for ggplot2 >= 3.4.0, would use `linewidth` instead of `size`
      size=0.3,
      alpha=0.7,
    ) +
    ggplot2::facet_grid(
      rows = dplyr::vars(uses_remote_info),
      scales = "free_y"
    ) +
    ggplot2::labs(
      y = STATISTIC_LABELS["fee"],
    ) +
    ggplot2::theme(
      legend.position = c(0.17, 0.83)
    )
  save_plot(plot = plt, filename = plot_file, scale_mult = 1.2)

}

make_emission_dwl_graph <- function(results_df, wildcards, plot_file, x_axis_form = c("ordinal", "cardinal")) {
  # Layout is somewhat different than fee graph, because we're not emphasizing the differences between the fees charged to different wells.
  # y-axis is % change, x-axis is same as fee axis
  # Facet cols: one for DWL, one for emiss
  # Facet rows: same as fee (divide up remote sensing and not)
  x_axis_form <- match.arg(x_axis_form)
  df_tmp <- format_audit_rules_variables(results_df)
  df_long <- dplyr::bind_rows(
      dplyr::transmute(df_tmp,
        uses_remote_info, policy_details,
        tau_T, tau_T_str, tau_T_str_fct,
        point_val = dwl_tot_rel_pct_mean / 100,
        bound_low = dwl_tot_rel_pct_conf_low / 100,
        bound_high = dwl_tot_rel_pct_conf_high / 100,
        statistic_type = STATISTIC_LABELS["dwl"],
      ),
      dplyr::transmute(df_tmp,
        uses_remote_info, policy_details,
        tau_T, tau_T_str, tau_T_str_fct,
        point_val = emis_tot_rel_pct_mean / 100,
        bound_low = emis_tot_rel_pct_conf_low / 100,
        bound_high = emis_tot_rel_pct_conf_high / 100,
        statistic_type = STATISTIC_LABELS["emis"],
      ),
    ) %>%
    dplyr::mutate(
      statistic_type = factor(statistic_type, levels=unname(STATISTIC_LABELS)),
    )

  stopifnot(noNAs(df_long), is.factor(df_long$tau_T_str_fct))

  pointsize <- 1.5
  dodge_width <- pointsize / 3
  if (x_axis_form == "ordinal") {
    x_var <- rlang::sym("tau_T_str_fct")
    # These values are included to let the reader see the curvature in the
    # cardinal plot, but are not relevant for the ordinal plot.
    df_long %<>% dplyr::filter(!tau_T_str %in% c("2500", "7500"))
  } else {
    x_var <- rlang::sym("tau_T")
    df_long %<>% dplyr::filter(!tau_T_str %in% c("low-3month"))
  }
  plt <- df_long %>%
    ggplot2::ggplot(ggplot2::aes(
        x = !!x_var,
        y = point_val,
        ymin = bound_low,
        ymax = bound_high,
        # 'group' is necessary to make ggplot connect dots with factor x-axis
        # group = policy_details,
        color = policy_details,
        shape = policy_details,
      )
    ) %>%
    add_shared_plot_elements() +
    # Use geom_point and geom_errorbar(..., width=0) rather than geom_pointrange
    # so that I can adjust the line width of the error bar and the point separately
    ggplot2::geom_point(
      ggplot2::aes(size=policy_details),
      position=ggplot2::position_dodge(width=dodge_width)
    ) +
    ggplot2::geom_errorbar(
      position=ggplot2::position_dodge(width=dodge_width),
      width = 0,
      # Note: for ggplot2 >= 3.4.0, would use `linewidth` instead of `size`
      size=pointsize / 3
    ) +
    ggplot2::scale_size_manual(values = c(1, 1, 1, 1.2, 1, 1)) +
    # ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size=2))) +
    ggplot2::facet_grid(
      rows = dplyr::vars(statistic_type),
      # scales = "free_y"
    ) +
    ggplot2::labs(
      y = "Percent reduction",
    )
  if (x_axis_form == "ordinal") {
    plt <- plt +
      ggplot2::theme(legend.position = c(0.17, 0.83))
  } else {
    plt <- plt +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 2L)) +
      ggplot2::theme(legend.position = "bottom")
  }

  label_tau_T <- x_axis_form != "ordinal"
  if (label_tau_T) {
    plt <- plt + ggplot2::geom_text(
      data = dplyr::filter(df_long, policy_details == "Remote, low threshold", !as.character(tau_T_str_fct) %in% c("2500", "7500")),
      mapping = ggplot2::aes(label = tau_T_str_fct),
      color="gray",
      nudge_y = 0.08,
      size=2,
      lineheight = 0.8,
    )
  }
  if (x_axis_form == "ordinal") {
    save_plot(plot = plt, filename = plot_file, scale_mult = 1.2)
  } else {
    save_plot(plot = plt, filename = plot_file, scale_mult = c(1, 1.5))
  }
}

SCM_PER_KG <- read_constants()$SOCIAL_COST_METHANE_PER_KG

STATISTIC_LABELS <- c(
  fee = "Expected fee % of med × 1 year",
  dwl = "DWL % reduced from BAU",
  emis = "Emis. % reduced from BAU"
)

# Note: the order here matters, and will be used in the "Paired" color palette
GRAPH_POLICY_ORDERING <- c(
  # "None",
  "Uniform",
  "Target covariates",
  "Target leaks, high threshold",
  "Remote, high threshold",
  "Target leaks, low threshold",
  "Remote, low threshold"
)

all_res <- read_policy_summaries(
    snakemake@input[["outcome_summaries"]],
    fill_redundant = TRUE
)

make_fee_graph(
  all_res, snakemake@wildcards,
  snakemake@output[["plot_fee"]]
)
make_emission_dwl_graph(
  all_res, snakemake@wildcards,
  snakemake@output[["plot_dwl_emis_ordinal"]],
  x_axis_form = "ordinal"
)
make_emission_dwl_graph(
  all_res, snakemake@wildcards,
  snakemake@output[["plot_dwl_emis_cardinal"]],
  x_axis_form = "cardinal"
)
