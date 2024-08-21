
suppressMessages(
  here::i_am("code/plot_well_operator_stats.R", uuid="69e06e85-ee2b-4c59-8bdb-1afd55468de6")
)
source(here::here("code/shared_functions.r"))
source(here::here("code/model_data_prep.r"))
options(scipen=99)

plot_well_pads_per_operator_cdf <- function(well_pad_crosswalk, output_file) {
  stopifnot(nrow(well_pad_crosswalk) > 1,  length(output_file) == 1)

  operator_well_pads <- well_pad_crosswalk %>%
    dplyr::distinct(well_pad_id, operator_company_name)

  should_be_empty <- dplyr::group_by(operator_well_pads, well_pad_id) %>%
    dplyr::count(operator_company_name) %>%
    dplyr::filter(n > 1)
  stopifnot(nrow(should_be_empty) == 0)

  pad_count_df <- operator_well_pads %>%
    dplyr::ungroup() %>%  # should already be ungrouped, but be sure
    dplyr::count(operator_company_name)

  plt <- ggplot2::ggplot(pad_count_df, ggplot2::aes(n)) +
    ggplot2::stat_ecdf(geom = "step") +
    ggplot2::scale_x_log10() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x="Well pads per operator",
      y="CDF"
    )
  save_plot(plt, output_file)
}

plot_wells_per_pad_histogram <- function(well_pad_crosswalk, output_file) {
  stopifnot(nrow(well_pad_crosswalk) > 1, length(output_file) == 1)
  pad_count_df <- well_pad_crosswalk %>%
    dplyr::ungroup() %>%  # should already be ungrouped, but be sure
    dplyr::count(well_pad_id) %>%
    dplyr::mutate(n = winsorize(n, c(0, 0.01)))
  plt <- ggplot2::ggplot(pad_count_df, ggplot2::aes(n)) +
    ggplot2::geom_histogram(binwidth=1) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x="Wells per well pad",
      y="Count"
    )
  save_plot(plt, output_file)
}


load_crosswalk_for_relevant_wells <- function(measurement_file, crosswalk_file) {
  stopifnot(length(measurement_file) == 1, length(crosswalk_file) == 1)
  # Load the well pad IDs for the wells that were flown over.
  well_pads <- prep_measurement_data(measurement_file) %>%
    dplyr::select(well_pad_id)
  # Subset the crosswalk for those wells, checking that the well pads are unique
  # and all rows in the `well_pads` data match
  well_pad_crosswalk <- arrow::read_parquet(crosswalk_file) %>%
    powerjoin::power_inner_join(well_pads, by="well_pad_id",
      check=merge_specs(
        duplicate_keys_left = "ignore",
        unmatched_keys_left = "ignore",
      )
    )

  well_pad_crosswalk
}

main <- function(snakemake) {
  well_pad_crosswalk <- load_crosswalk_for_relevant_wells(
    snakemake@input$matched_leakage,
    snakemake@input$well_pad_crosswalk
  )

  plot_wells_per_pad_histogram(
    well_pad_crosswalk,
    snakemake@output$wells_per_pad_histogram
  )
  plot_well_pads_per_operator_cdf(
    well_pad_crosswalk,
    output_file=snakemake@output$well_pads_per_operator_cdf
  )
}


if (! exists("snakemake")) {
  message("Using placeholder snakemake")
  snakemake <- SnakemakePlaceholder(
    input = list(
      headers = "data/generated/production/well_headers/",
      # prod = "data/generated/production/monthly_production/",
      well_pad_crosswalk = "data/generated/production/well_pad_crosswalk_1970-2018.parquet",
      matched_leakage = "data/generated/methane_measures/matched_wells_all.parquet"
    ),
    output = list(
      well_pads_per_operator_cdf  = "graphics/well_pads_per_operator_cdf.pdf",
      wells_per_pad_histogram = "graphics/wells_per_pad_histogram.pdf"
    ),
    resources = list(
      mem_mb = 2000
    )
  )
}

try(memory_limit(snakemake@resources$mem_mb))

main(snakemake)
rm(snakemake) # easier debugging
