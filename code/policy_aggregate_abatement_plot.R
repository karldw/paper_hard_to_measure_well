
source(here::here("code/shared_functions.r"))
source(snakemake@input[["policy_output_helper_functions"]])
options(scipen = 99, mc.cores=snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])

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

  save_plot(plt, plot_file, reproducible=TRUE)
}

all_res_no_dup <- read_policy_summaries(
  snakemake@input[["outcome_summaries"]],
  fill_redundant=FALSE
)
plot_aggregate_abatement_curve(all_res_no_dup, snakemake@output$aggregate_abatement_plot)
