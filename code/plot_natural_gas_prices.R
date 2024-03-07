suppressMessages(
  here::i_am("code/plot_natural_gas_prices.R", uuid="b3318ae8-e132-4a01-b327-c4e61a84c08d")
)

source(here::here("code/shared_functions.r"))

if (!exists("snakemake")) {
  warning("using placeholder snakemake")

  snakemake <- SnakemakePlaceholder(
    input = list(
      nat_gas_prices = here::here("data/generated/nat_gas_prices_by_basin.parquet"),
      matched_leakage = here::here("data/generated/methane_measures/matched_wells_all.parquet"),
      r_lib = here::here("scratch/snakemake_flags/setup_r_library")
    ),
    output = list(
      nat_gas_price_timeseries = here::here("graphics/nat_gas_price_timeseries.pdf"),
      nat_gas_price_histogram = here::here("graphics/nat_gas_price_histogram.pdf")
    ),
    threads = 1L
  )
}

plot_price <- function(price_file, outfile) {
  stopifnot(length(price_file) == 1, length(outfile) == 1)

  to_plot <- arrow::read_parquet(price_file, col_select=c("date", "basin", "price_real")) %>%
    dplyr::filter(
      date >= as.Date("2014-01-01"), date <= as.Date("2020-01-01"),
      basin %in% c("San Juan", "San Joaquin"),
    ) %>%
    dplyr::mutate(region_fct = factor(basin, levels = c("San Joaquin", "San Juan"), labels=c("CA", "NM / CO")))
  stopifnot(noNAs(to_plot))
  plt <- ggplot2::ggplot(to_plot, ggplot2::aes(date, price_real, color=region_fct)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::theme_bw() +
    ggplot2::labs(x="", y="Price ($2019/mcf)", color="") +
    ggplot2::ylim(0, NA) +
    ggplot2::theme(legend.position = c(0.1, 0.2))
  save_plot(plt, outfile, reproducible=TRUE)
  invisible(plt)
}

plot_histogram <- function(matched_leaks_file, outfile) {
  to_plot <- matched_leaks_file %>%
    arrow::read_parquet(col_select=c("basin", "emiss_kg_hr", "gas_price_per_mcf")) %>%
    dplyr::mutate(
      detect_emiss = factor(!is.na(emiss_kg_hr), levels=c(TRUE, FALSE), labels=c("Leak detected", "No leak detected")),
      region_fct = factor(dplyr::case_when(
        basin %in% c("San Joaquin", "Other California") ~ "California",
        basin == "San Juan" ~ "New Mexico / Colorado",
      )),
      basin_fct = factor(basin, levels = c("San Joaquin", "Other California", "San Juan")),
    )
  # Note: could do colors by basin, and the code is set up for it, but the graph
  # becomes kind of busy without adding much insight.
  plt <- ggplot2::ggplot(
      to_plot,
      ggplot2::aes(
        x = gas_price_per_mcf,
        # fill = basin_fct
      )
    ) +
    ggplot2::geom_histogram(bins=20) +
    ggplot2::facet_grid(
      cols=ggplot2::vars(region_fct),
      rows=ggplot2::vars(detect_emiss),
      scales="free_y"
    ) +
    ggplot2::scale_x_continuous(labels=scales::label_dollar()) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x="Gas price ($2019/mcf)",
      y="Count",
      fill="Basin FE",
    ) +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::theme(legend.position = c(0.15, 0.8))
  save_plot(plt, outfile, reproducible=TRUE)
  invisible(plt)
}


plot_ecdf <- function(matched_leaks_file, outfile) {
  to_plot <- matched_leaks_file %>%
    arrow::read_parquet(col_select=c("basin", "emiss_kg_hr", "gas_price_per_mcf")) %>%
    dplyr::mutate(
      detect_emiss = factor(!is.na(emiss_kg_hr), levels=c(TRUE, FALSE), labels=c("Leak detected", "No leak detected")),
      region_fct = factor(dplyr::case_when(
        basin %in% c("San Joaquin", "Other California") ~ "California",
        basin == "San Juan" ~ "New Mexico / Colorado",
      ))
    )
  plt <- ggplot2::ggplot(to_plot, ggplot2::aes(gas_price_per_mcf, color=detect_emiss)) +
    ggplot2::stat_ecdf(geom="step") +
    ggplot2::facet_grid(
      cols=ggplot2::vars(region_fct),
      scales="free_y"
    ) +
    ggplot2::scale_x_continuous(labels=scales::label_dollar()) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x="Gas price ($2019/mcf)",
      y="Count",
      color=""
    ) +
    ggplot2::theme(legend.position = c(0.15, 0.85))
  save_plot(plt, outfile, reproducible=TRUE)
  invisible(plt)
}

plot_price(
  snakemake@input$nat_gas_prices,
  snakemake@output$nat_gas_price_timeseries
)

plot_histogram(
  snakemake@input$matched_leakage,
  snakemake@output$nat_gas_price_histogram
)

# Note: if using, need to define nat_gas_price_ecdf filename
# plot_ecdf(
#   snakemake@input$matched_leakage,
#   snakemake@output$nat_gas_price_ecdf
# )