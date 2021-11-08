
source(here::here("code/shared_functions.r"))
source(snakemake@input[["policy_output_helper_functions"]])

options(scipen = 99, mc.cores=snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])

plot_shadow_prices <- function(results_df, plot_file) {
  constants <- read_constants()
  price_comparison <- constants$AUDIT_COST_TO_COMPARE_FIXED_VS_OPTIM_DWL
  plt <- results_df %>%
    dplyr::filter(
      !audit_rule %in% c("none", "remote"),
      tau == constants$TAU_LEVELS[["med"]],
      time_T == constants$T_LEVELS[["3month"]],
      (audit_rule == "target_e") + (detect_threshold == 0) == 1,
      audit_frac > 0,
    ) %>%
    dplyr::mutate(
      policy = format_policy_details(audit_rule, detect_threshold, latex=FALSE),
      policy = factor(policy, levels=levels(policy), labels=tidy_gsub(levels(policy), ", ", "\n", fixed=TRUE)),
    ) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = audit_frac,
      y = -shadow_price_mean,
      ymin = -shadow_price_conf_low,
      ymax = -shadow_price_conf_high,
      color = policy,
      fill = policy,
    )) +
    ggplot2::geom_hline(yintercept=price_comparison, linetype="dashed", alpha=0.6) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha=0.2, linetype="blank") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Audit budget (fraction of wells per years)",
      y = "Shadow price",
      color="Audit rule"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::scale_fill_brewer(palette = "Dark2", guide="none") +
    ggplot2::theme(legend.position = c(0.85, 0.8))

  save_plot(plt, plot_file, reproducible=TRUE)
}

# Note: with changing design, this read_summary / fill_redundant_audit_rules
# isn't the cleanest, but it's easier to reuse code.

all_res <- read_policy_summaries(
  snakemake@input[["outcome_summaries"]],
  fill_redundant=TRUE
)

plot_shadow_prices(all_res, snakemake@output$shadow_price_plot)
