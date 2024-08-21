
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
      plot_file_dwl = "graphics/figureA11_robustness_audit_result_dwl_rule=target_e_high_frac=1pct_tauT=med-3month.pdf",
      plot_file_emiss = "graphics/figureA11_robustness_audit_result_emission_rule=target_e_high_frac=1pct_tauT=med-3month.pdf"
    ),
    resources = list(mem_mb = 6000),
    threads = 2,
    wildcards = list(
      audit_rule = ""
    )
  )
}


source(snakemake@input[["policy_output_helper_functions"]])


options(scipen = 5, mc.cores=snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])

PLOT_COLORS = c(
  "Main" = "black",
  "Robustness: add geo" = "plum2",
  "Robustness: subsets" = "gray"
)


df = snakemake@input[["outcome_summaries"]] |>
  purrr::map_dfr(arrow::read_parquet, .id="file_idx") |>
  dplyr::mutate(
    source_filename = snakemake@input[["outcome_summaries"]][as.integer(file_idx)],
    # Note: this hard-coding of robustness1, robustnes2, and robustness3 relies
    # on the ordering used in model_data_prep_robustness.r, where those
    # robustness specs were assigned #1-3. If adding more robustness specs,
    # need to edit here too.
    spec_type = dplyr::case_when(
      grepl("/main_spec/", source_filename, fixed = TRUE) ~ names(PLOT_COLORS)[1],
      grepl("/robustness(1|2|3)/", source_filename, perl = TRUE) ~ names(PLOT_COLORS)[2],
      grepl("/robustness", source_filename, fixed = TRUE) ~ names(PLOT_COLORS)[3],
      TRUE ~ NA_character_
    )
  )
stopifnot(!anyNA(df))

plotting_shared <- function(df) {
  plt = df |>
    dplyr::arrange(coef_value) |>
    dplyr::mutate(idx = dplyr::row_number()) |>
    ggplot(aes(
      x = coef_value,
      y = idx,
      xmin = coef_conf_low,
      xmax = coef_conf_high,
      shape = spec_type,
      color = spec_type,
    )) +
    scale_shape_manual(values = c(19, 2, 19)) +
    scale_color_manual(values = PLOT_COLORS) +
    geom_point() +
    geom_errorbarh() +
    labs(
      y = "",
      shape = "",
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


# This is copied with minimal modification from group_into_well_pads.R
add_count_wells_nearby_same_operator <- function(well_df) {
  # Work from most populous first, just so the parallel processing doesn't
  # wait too long on single straglers.
  operators <- well_df |>
    dplyr::count(operator_company_name) |>
    dplyr::arrange(dplyr::desc(n)) |>
    purrr::chuck("operator_company_name")

  counts_df <- furrr::future_map_dfr(
    operators, count_wells_nearby_one_operator,
    well_df=well_df,
    .options=furrr::furrr_options(seed=TRUE)
  )
  stopifnot(anyDuplicated(counts_df$well_pad_id) == 0)
  out <- dplyr::inner_join(well_df, counts_df, by="well_pad_id")
  out
}


count_wells_nearby_one_operator <- function(operator_company_name, well_df) {
  individual_well_geom <- well_df |>
    dplyr::filter(.data$operator_company_name == !!operator_company_name) |>
    dplyr::select(entity_id, surface_longitude_wgs84, surface_latitude_wgs84) |>
    lonlat_to_projected(c("surface_longitude_wgs84", "surface_latitude_wgs84"))

  well_pad_df <- well_df |>
    dplyr::filter(.data$operator_company_name == !!operator_company_name) |>
    dplyr::select(well_pad_id, well_pad_lon, well_pad_lat) |>
    dplyr::distinct()

  nearby_dist_km <- 10

  stopifnot(
    length(operator_company_name) == 1L,
    anyDuplicated(well_df[["entity_id"]]) == 0L,
    anyDuplicated(well_df[["well_pad_id"]]) != 0L,
    nrow(individual_well_geom) > 0,
    nrow(well_pad_df) > 0,
    nearby_dist_km > 0
  )

  well_pad_buffered <- well_pad_df |>
    lonlat_to_projected(c("well_pad_lon", "well_pad_lat")) |>
    sf::st_buffer(dist = nearby_dist_km * 1000) # default to 10km
  intersect_count <- sf::st_intersects(well_pad_buffered, individual_well_geom) |>
    purrr::map_int(length)

  # Note: this count includes wells in the well pad.
  out <- well_pad_df |>
    dplyr::select(well_pad_id) |>
    dplyr::mutate(same_operator_count_wells_within_10km = !!intersect_count)

  out
}


plot_well_operator_density_vs_well_density <- function() {
  well_pads_matched <- arrow::read_parquet(here::here("data/generated/methane_measures/matched_wells_all.parquet"))
  well_pad_crosswalk <- arrow::read_parquet(here::here("data/generated/production/well_pad_crosswalk_1970-2018.parquet"))

  # Note: of the flown-over data 833 well pads rows (of 14,197), had multiple
  # operator_company_name values. Here we take the modal one.
  well_df <- well_pad_crosswalk |>
    dplyr::group_by(well_pad_id) |>
    dplyr::mutate(operator_company_name = Mode(operator_company_name)) |>
    dplyr::ungroup()

  # Note: this takes a while - you might want to increase snakemake@threads
  future::plan("multicore", workers = snakemake@threads)
  to_plot <- add_count_wells_nearby_same_operator(well_df) |>
    dplyr::distinct(well_pad_id, tot_count_wells_within_10km, same_operator_count_wells_within_10km) |>
    dplyr::inner_join(
      dplyr::select(well_pads_matched, well_pad_id, basin),
      by="well_pad_id"
    )
  # Report the corelation of the asinh-transformed variables
  # (because that's the transformation used in the regression)
  cor_asinh <- cor(asinh(to_plot$tot_count_wells_within_10km), asinh(to_plot$same_operator_count_wells_within_10km))
  # Print out the (pearson) correlation:
  message(glue::glue("Correlation: {signif(cor_asinh, 2)}"))

  plt <- ggplot2::ggplot(to_plot, ggplot2::aes(tot_count_wells_within_10km, same_operator_count_wells_within_10km, color=basin)) +
    ggplot2::geom_point(alpha=0.1, size=0.01) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(transform = "log10") +
    ggplot2::scale_y_continuous(transform = "log10") +
    ggplot2::geom_abline(slope=1, intercept = 0, linetype="dashed", alpha=0.8, linewidth=1) +
    ggplot2::labs(
      x = "Total count of wells within 10 km",
      y = "Same-operator count of wells within 10 km",
      color = "Basin",
    ) +
    ggplot2::theme(
      legend.position = c(0.2, 0.8)
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=1, alpha=0.9)))
    save_plot(plt, here::here("graphics/robustness_well_count_vs_same_operator_count.pdf"))
  invisible(NULL)
}

plt_emis <- df |>
  dplyr::mutate(
    coef_value = emis_mean_rel_pct_mean / 100,
    coef_conf_low = emis_mean_rel_pct_conf_low / 100,
    coef_conf_high = emis_mean_rel_pct_conf_high / 100,
  ) |>
  plotting_shared() +
  labs(
    x = "Emis. % reduced from BAU"
  )
save_plot(
  plt_emis,
  snakemake@output[["plot_file_emiss"]]
)

plt_dwl <- df |>
  dplyr::mutate(
    coef_value = dwl_tot_rel_pct_mean / 100,
    coef_conf_low = dwl_tot_rel_pct_conf_low / 100,
    coef_conf_high = dwl_tot_rel_pct_conf_high / 100,
  ) |>
  plotting_shared() +
  labs(
    x = "DWL % reduced from BAU"
  )
save_plot(
  plt_dwl,
  snakemake@output[["plot_file_dwl"]]
)

# plot_well_operator_density_vs_well_density()

