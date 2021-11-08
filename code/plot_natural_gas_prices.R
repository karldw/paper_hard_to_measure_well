
source(here::here("code/shared_functions.r"))

if (!exists("snakemake")) {
  stop("This program is meant to by run by snakemake")
}

plot_price <- function(price_file, outfile) {
  stopifnot(length(price_file) == 1, length(outfile) == 1)

  to_plot <- arrow::read_parquet(price_file, col_select=c("date", "basin", "price_real")) %>%
    dplyr::filter(
      date >= as.Date("2014-01-01"), date <= as.Date("2020-01-01"),
      basin %in% c("San Juan", "San Joaquin"),
    ) %>%
    dplyr::mutate(basin_fct = factor(basin, levels = c("San Joaquin", "San Juan"), labels=c("CA", "NM / CO")))
  stopifnot(noNAs(to_plot))
  plt <- ggplot2::ggplot(to_plot, ggplot2::aes(date, price_real, color=basin_fct)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::theme_bw() +
    ggplot2::labs(x="", y="Price ($2019/mcf)", color="") +
    ggplot2::ylim(0, NA) +
    ggplot2::theme(legend.position = c(0.1, 0.2))
  save_plot(plt, outfile, reproducible=TRUE)
  invisible(plt)
}

plot_price(
  snakemake@input$nat_gas_prices,
  snakemake@output$nat_gas_price_timeseries
)
