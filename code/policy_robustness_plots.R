
# Note: the graphs below are not binary-reproducible run to run, so will always
# show as different.

suppressMessages(
  here::i_am("code/policy_robustness_plots.R", uuid="f1c51536-d8eb-4022-be69-f7c1b61e4af1")
)
source(here::here("code/shared_functions.r"))

library(ggplot2)

if (!exists("snakemake")) {
  snakemake <- SnakemakePlaceholder(
    input = list(
      outcome_summaries = c(
        glue::glue("data/generated/policy_outcomes/robustness{1:52}/08_twopart_lognormal_heterog_alpha-period_8760_hours/audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=med-3month.parquet"),
        "data/generated/policy_outcomes/main_spec/08_twopart_lognormal_heterog_alpha-period_8760_hours/audit_outcome_summary_rule=target_e_high_frac=1pct_tauT=med-3month.parquet"
      ),
      policy_output_helper_functions = "code/policy_output_helper_functions.r"
    ),
    output = list(
      plot_file_dwl = "graphics/robustness_audit_result_dwl_rule=target_e_high_frac=1pct_tauT=med-3month.pdf",
      plot_file_emiss = "graphics/robustness_audit_result_emission_rule=target_e_high_frac=1pct_tauT=med-3month.pdf"
    ),
    resources = list(mem_mb = 4000),
    threads = 1,
    wildcards = list(
      audit_rule = ""
    )
  )
}


source(snakemake@input[["policy_output_helper_functions"]])


options(scipen = 5, mc.cores=snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])


df = snakemake@input[["outcome_summaries"]] %>%
  purrr::map_dfr(arrow::read_parquet, .id="file_idx") %>%
  dplyr::mutate(
    source_filename = snakemake@input[["outcome_summaries"]][as.integer(file_idx)],
    spec_type = dplyr::case_when(
      grepl("/main_spec/", source_filename, fixed = TRUE) ~ "Main",
      grepl("/robustness", source_filename, fixed = TRUE) ~ "Robustness",
      TRUE ~ NA_character_
    )
  )
stopifnot(!anyNA(df))

plotting_shared <- function(df) {
  plt = df %>%
    dplyr::arrange(coef_value) %>%
    dplyr::mutate(idx = dplyr::row_number()) %>%
    ggplot(aes(
      x = coef_value,
      y = idx,
      xmin = coef_conf_low,
      xmax = coef_conf_high,
      color = spec_type,
    )) +
    scale_color_manual(values = c("black", "gray")) +
    geom_point() +
    geom_errorbarh() +
    labs(
      y = "",
      color = ""
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.length.y.left = unit(0, "pt"),
      legend.position = c(0.15, 0.75)
    ) +
    expand_limits(x = c(0.40, 0.75))
  plt
}

plt_emis <- df %>%
  dplyr::mutate(
    coef_value = emis_mean_rel_pct_mean / 100,
    coef_conf_low = emis_mean_rel_pct_conf_low / 100,
    coef_conf_high = emis_mean_rel_pct_conf_high / 100,
  ) %>%
  plotting_shared() +
  labs(
    x = "Emis. % reduced from BAU"
  )
save_plot(
  plt_emis,
  snakemake@output[["plot_file_emiss"]]
)

plt_dwl <- df %>%
  dplyr::mutate(
    coef_value = dwl_tot_rel_pct_mean / 100,
    coef_conf_low = dwl_tot_rel_pct_conf_low / 100,
    coef_conf_high = dwl_tot_rel_pct_conf_high / 100,
  ) %>%
  plotting_shared() +
  labs(
    x = "DWL % reduced from BAU"
  )
save_plot(
  plt_dwl,
  snakemake@output[["plot_file_dwl"]]
)
