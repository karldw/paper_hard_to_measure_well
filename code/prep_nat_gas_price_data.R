source(here::here("code/shared_functions.r"))

read_sheet <- function(filename) {
  stopifnot(length(filename) == 1)
  # Read as text instead of numeric because read_excel thinks they're dates and
  # gives a bunch of warnings.
  df <- readxl::read_excel(
      filename,
      range="Data!A16:N136",
      na=c("", "NA"),
      col_types = c("date", rep("text", 13)),
    ) %>%
    dplyr::rename_all(make_better_names) %>%
    dplyr::mutate(date = lubridate::date(date)) %>%
    dplyr::mutate(dplyr::across(-date, as.numeric))
  names(df) <- tidy_gsub(names(df), "snl_natural_gas_(.+)_bidweek_index_final_price_value", "\\1")

  df %<>% tidyr::pivot_longer(-date, names_to = "price_series", values_to="price_nominal_per_mmbtu")
  df
}


match_basin <- function(df, basin_mapping_csv) {
  stopifnot(length(basin_mapping_csv) == 1)
  basin_mapping <- data.table::fread(basin_mapping_csv, skip=3, data.table=FALSE)
  # Note this is a m:m (cartesian) match -- a single basin can have multiple
  # series, and a series can be matched to multiple basins. We collapse by
  # taking the (unweighted) mean.
  out <- safejoin::safe_inner_join(df, basin_mapping, by="price_series", check="B C M N L T") %>%
    dplyr::group_by(date, basin) %>%
    dplyr::summarize(
      price_nominal = mean_(price_nominal_per_mmbtu),
      price_real = mean_(price_real),
      .groups="drop"
    )
  out
}


adjust_for_inflation <- function(df, price_index_file) {
  stopifnot(length(price_index_file) == 1)
  price_index <- load_price_index(price_index_file, base_year=2019)
  out <- safejoin::safe_inner_join(df, price_index, by="date", check="B C V L T") %>%
    dplyr::mutate(price_real = price_nominal_per_mmbtu / price_index)
  out
}


write_parquet <- function(x, sink, ...) {
  stopifnot(length(sink) == 1)
  arrow::write_parquet(x, sink, ..., version="2.0")
}


if (!exists("snakemake")) {
  stop("This program is meant to by run by snakemake")
}

df <- read_sheet(snakemake@input$snl_price_file) %>%
  adjust_for_inflation(snakemake@input$price_index_file) %>%
  match_basin(snakemake@input$basin_mapping) %>%
  write_parquet(snakemake@output$nat_gas_prices)
