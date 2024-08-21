suppressMessages(
  here::i_am("code/policy_aggregate_abatement_plot.R", uuid="7d491c25-7a1e-433f-b850-f73660b5d341")
)

source(here::here("code/shared_functions.r"))


if (!exists("snakemake")) {
  snakemake <- SnakemakePlaceholder(
    input = list(
      outcome_summaries = fs::dir_ls(
        here::here("data/generated/policy_outcomes/main_spec/08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours"),
        regexp = "audit_outcome_summary_rule=remote_low_frac\\=0pct_tauT=[0-9]+\\.?[0-9]*.parquet$",
      ),
      r_lib = "scratch/snakemake_flags/setup_r_library",
      script = "code/policy_aggregate_abatement_plot.R",
      constants = "code/constants.json",
      policy_output_helper_functions = "code/policy_output_helper_functions.r"
    ),
    output = list(
      aggregate_abatement_plot = "graphics/aggregate_abatement_curve.pdf"
    ),
    threads = 1
  )
}


source(snakemake@input[["policy_output_helper_functions"]])
options(scipen = 99, mc.cores=snakemake@threads)

const <- read_constants()
METHANE_GWP <- const$METHANE_GWP
SCM_PER_KG <- const$SCM_PER_KG

plot_aggregate_abatement_curve <- function(results_df, plot_file) {
  stopifnot(length(plot_file) == 1, nrow(results_df) > 1)
  df <- results_df %>%
    dplyr::filter(
      audit_rule == "remote",
      audit_frac == 0,
      detect_threshold == 0,
    ) %>%
    dplyr::transmute(
      # Convert units: cost was in $/kg CH4, now in $/ton CO2e
      # emiss was in kg CH4, now in millions of metric tons CO2e
      # tot_cost_per_kg_mean includes both private cost (commodity value) and
      # the fee amount
      tot_cost_per_ton_co2e = tot_cost_per_kg_mean_mean * 1000 / METHANE_GWP,
      emis_reduce_mn_ton_co2e = emis_reduce_tot_mean * METHANE_GWP / 1e9,
      emis_reduce_mn_ton_co2e_low = emis_reduce_tot_conf_low * METHANE_GWP / 1e9,
      emis_reduce_mn_ton_co2e_high = emis_reduce_tot_conf_high * METHANE_GWP / 1e9,
    )
  line_private <- min(df$tot_cost_per_ton_co2e)
  line_optimal <- line_private + (SCM_PER_KG * 1000 / METHANE_GWP)

  plt <- ggplot2::ggplot(df, ggplot2::aes(
      x=tot_cost_per_ton_co2e,
      y=emis_reduce_mn_ton_co2e,
      ymin=emis_reduce_mn_ton_co2e_low,
      ymax=emis_reduce_mn_ton_co2e_high
    )) +
    ggplot2::geom_vline(xintercept=line_private, alpha=0.3, linetype="dashed") +
    ggplot2::geom_vline(xintercept=line_optimal, alpha=0.3, linetype="dashed") +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha=0.3) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip(xlim = c(0, NA)) +
    # because of the coord_flip, x and y labelds are flipped
    ggplot2::labs(
      y = "Annual reductions in methane emissions (millions of tons CO2e)", # x axis
      x = "Marginal cost ($ per ton CO2e)",
      title="Aggregate abatement supply curve"
    )

  save_plot(plt, plot_file)
}

# Note: here we don't need the extra processing that happens in read_policy_summaries
all_res_no_dup <- purrr::map_dfr(snakemake@input[["outcome_summaries"]], arrow::read_parquet)
plot_aggregate_abatement_curve(all_res_no_dup, snakemake@output$aggregate_abatement_plot)
