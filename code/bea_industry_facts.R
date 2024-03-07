suppressMessages(
  here::i_am("code/bea_industry_facts.R", uuid="c4bb8916-4aea-47c6-99ab-97909e20f68d")
)
source(here::here("code/shared_functions.r"))
memory_limit(snakemake@resources[["mem_mb"]])
suppressWarnings(loadNamespace("lubridate")) # avoid https://github.com/tidyverse/lubridate/issues/965
options(warn = 2, mc.cores=1)

df_gdp <- readxl::read_excel(snakemake@input$gdp_file, range="A6:W907")
years <- 1997:2017
stopifnot(names(df_gdp) == c("Line", "...2", as.character(years)))
df_gdp %<>% dplyr::rename(line = Line, item = ...2) %>%
  dplyr::mutate(line = as.integer(line), item = trimws(item))


assert_item <- function(df, item) {
  stopifnot(setequal(df$item, item))
  df
}

df_deprec <- readxl::read_excel(snakemake@input$deprec_file, range="A6:BW85")
deprec_years <- 1947:2019
stopifnot(names(df_deprec) == c("Line", "...2", as.character(deprec_years)))
df_deprec %<>% dplyr::rename(line = Line, item = ...2) %>%
  dplyr::filter(!is.na(line)) %>%
  dplyr::mutate(line = as.integer(line), item = trimws(item)) %>%
  dplyr::filter(line == 6) %>%
  dplyr::select(-line) %>%
  data.table::transpose(make.names=TRUE) %>%
  dplyr::mutate(year = !!deprec_years) %>%
  dplyr::rename(deprec_current_cost = "Oil and gas extraction")


prices = load_price_index(snakemake@input$price_index_file, base_year=2019) %>%
  dplyr::mutate(year = as.integer(lubridate::year(date))) %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(price_index = mean(price_index), .groups="drop")

oil_gas <- dplyr::filter(df_gdp, line %in% 56:57) %>%
  assert_item(c("Value added", "Compensation of employees")) %>%
  dplyr::select(-line) %>%
  data.table::transpose(make.names=TRUE) %>%
  dplyr::mutate(year = !!years) %>%
  dplyr::mutate_if(is.character, as.numeric) %>%
  powerjoin::power_inner_join(df_deprec, by="year", check=
    merge_specs(
      unmatched_keys_left = "ignore",
      unmatched_keys_right = "ignore",
    )
  ) %>%
  powerjoin::power_inner_join(prices, by="year", check=
    merge_specs(
      unmatched_keys_left = "ignore",
      unmatched_keys_right = "ignore",
    )
  ) %>%
  dplyr::rename_all(make_better_names) %>%
  dplyr::mutate(
    net_value_nominal_bn = value_added - compensation_of_employees - deprec_current_cost,
    net_value_real_bn = net_value_nominal_bn / price_index,
  ) %>%
  dplyr::as_tibble()


stopifnot(!anyNA(oil_gas))

# 2017 is lastest easily available. 2013 is arbitrary, but 5 years seems nice.
mean_net_val <- dplyr::filter(oil_gas, dplyr::between(year, 2013, 2017))$net_value_real_bn %>% mean() %>%
  signif(2)
to_write <- c(
  "% Def: value added - employee compensation - current-cost depreciation, expressed in $2019 billions",
  "% Data from BEA",
  paste0(mean_net_val, "%") # "%" to avoid extra space in latex
)

writeLines(to_write, snakemake@output$net_value)
