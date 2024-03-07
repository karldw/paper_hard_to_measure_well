suppressMessages(
  here::i_am("code/summary_stats.R", uuid="d35c83a9-af54-4dc8-a035-d5ded6531bb1")
)

library(ggplot2)
suppressWarnings(loadNamespace("lubridate")) # avoid https://github.com/tidyverse/lubridate/issues/965

source(here::here("code/shared_functions.r"))
source(here::here("code/model_data_prep.r"))

loadNamespace("sf")
options(scipen=99)


count_gas_completions_by_size <- function(header_files) {
  # Note: cut requires size cuts to cover the whole domain of the variable being
  # cut. Therefore, have 0 and Inf.
  size_cuts <- c(0, 10, 100, 1000, Inf)
  used_cols <- c("completion_date", "first_prod_date", "production_type", "prac_ip_gas_daily")
  # Filter to wells that start producing within 1 year of completion. Anything
  # else seems like a data error or weird well. Because completion is daily but
  # production is monthly, also check the month.
  df <- load_headers_by_filename(header_files, used_cols) %>%
    dplyr::filter(
      !is.na(completion_date),
      production_type == "Gas",
      prac_ip_gas_daily > 0,
    ) %>%
    dplyr::mutate(
      completion_month = lubridate::make_date(
        lubridate::year(completion_date),
        lubridate::month(completion_date),
        1L),
      size_class = cut(prac_ip_gas_daily, size_cuts, right=FALSE, dig.lab=50)
    ) %>%
    dplyr::filter(
      is.na(first_prod_date) |
      dplyr::between(as.numeric(first_prod_date - completion_date), 0, 365) |
      first_prod_date == completion_month
    )
  out <- df %>%
    dplyr::group_by(size_class, completion_month) %>%
    dplyr::summarize(N = dplyr::n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(size_class, completion_month, fill=list(N = 0L))
  out
}


get_pretty_varnames <- function(orig_varnames) {
  # data might have other variables, but we only want to display these.
  # (variables in the list below but not in the data will not raise an error)
  all_pretty_vars <- c(
    age_yr = "Age (yr)",
    gas_avg_mcfd = "Gas (mcfd)",
    oil_avg_bbld = "Oil (bbld)",
    detect_emiss_pct = "Detect leak (\\%)",
    emiss_kg_hr = "Leak size (kg/hr)",
    gas_price_per_mcf = "Gas price (\\textdollar/mcf)",
    distance_to_nearest_pad_km = "Dist. to next pad (km)",
    tot_count_wells_within_10km = "N wells within 10 km"
  )
  all_pretty_vars[intersect(names(all_pretty_vars), orig_varnames)]
}

add_emiss_vars <- function(df) {
  if ("emiss_kg_hr" %in% names(df)) {
    df %<>% dplyr::mutate(
        emiss_kg_hr = dplyr::na_if(emiss_kg_hr, 0),
        detect_emiss = not_na(emiss_kg_hr),
      )
  }
  # Lyon et al. 2016 has detect_emiss but not emiss_kg_hr
  if ("detect_emiss" %in% names(df)) {
    df$detect_emiss_pct <- df$detect_emiss * 100
  }
  df
}


convert_km_to_nearest_well_pad <- function(df) {
  if ("distance_to_nearest_pad_m" %in% names(df)) {
    df <- dplyr::mutate(df, distance_to_nearest_pad_km = distance_to_nearest_pad_m / 1000)
  }
  df
}

summary_stats_wells <- function(df, output_file) {
  summary_stats <- list(
    mean = mean_,
    sd = sd_,
    p10 = function(x) quantile_(x, 0.1),
    p90 = function(x) quantile_(x, 0.9)
    # min = min_,
    # max = max_
    # n = function(x) sum(!is.na(x))
  )

  # As currently written, these names can only letters a-z, A-Z, and 0-9, not "_"
  # Change the regex in names_pattern if you want more flexibility.
  stopifnot(all(grepl("^[a-z0-9A-Z]+$", names(summary_stats), perl=TRUE)))
  cols_to_format <- names(summary_stats) %>% setdiff("n")

  df %<>% add_emiss_vars() %>%
    convert_km_to_nearest_well_pad()

  pretty_vars <- get_pretty_varnames(names(df))

  df_sum <- df %>%
    dplyr::summarize_at(
      names(pretty_vars),
      summary_stats,
    ) %>%
    tidyr::pivot_longer(
      everything(),
      names_pattern = "(.+)_([a-zA-Z0-9]+)$",
      names_to = c("var", "stat")
    ) %>%
    tidyr::pivot_wider(
      id_cols = c("var"),
      names_from="stat"
    ) %>%
    dplyr::mutate(var_pretty = pretty_vars[var])
  if (any(df_sum$sd == 0)) {
    varnames <- dplyr::filter(df_sum, sd == 0)$var_pretty %>% paste(collapse="', ''")
    message("\nDropping '", varnames, "'\nfrom the table because of zero variance")
    df_sum %<>% dplyr::filter(sd != 0)
  }
  x <- df_sum %>%
    gt::gt() %>%
    gt::fmt_percent(
      columns = !!cols_to_format,
      rows = var == "detect_emiss",
      locale="en_US",
    ) %>%
    gt::fmt_number(
      columns = !!cols_to_format,
      rows = var %in% names(pretty_vars),
      decimals = 1,
      sep_mark = ""
    ) %>%
    gt::cols_hide("var") %>%
    gt::cols_move_to_start(c("var_pretty"))
  # Run gt::as_latex, but remove all the headers and footers
  make_table_fragment(x, escaped=FALSE) %>%
    # also remove all the $ they put around every number
    tidy_gsub("$", "", fixed=TRUE) %>%
    writeLines(output_file)
  invisible(df)
}

write_n_obs <- function(df_list, file_list) {
  stopifnot(length(df_list) == length(file_list))
  lengths <- purrr::map_int(df_list, nrow)
  lengths <- ifelse(lengths >= 10000, prettyNum(lengths, big.mark=","), prettyNum(lengths)) %>%
    paste0("%") # prevent extra space in latex
  purrr::walk2(lengths, file_list, ~writeLines(.x, .y, sep=""))
}

load_all_wells <- function(prod_dir, header_dir, well_pad_crosswalk) {
  # Steps:
  # 1. Load prod and headers
  # 2. Match to well pads
  # 3. Match to prices

  as_of_date <- as.Date("2018-06-01")
  as_of_year <- lubridate::year(as_of_date)
  as_of_month <- lubridate::month(as_of_date)
  as_of_days_per_month <- as.vector(lubridate::days_in_month(as_of_date)) # 30

  monthly_prod <- arrow::open_dataset(prod_dir) %>%
    dplyr::filter(year == !!as_of_year, month == !!as_of_month, daily_avg_gas > 0) %>%
    dplyr::select(entity_id, daily_avg_gas, daily_avg_oil) %>%
    dplyr::collect()

  # entity_id is assigned by drillinginfo
  # well_pad_id is assigned by us, and could change, so don't hardcode values
  well_pad_mapping <- well_pad_crosswalk %>%
    arrow::read_parquet(col_select=c(
      "entity_id", # assigned by drillinginfo
      "well_pad_id", # assigned in group_into_well_pads.R
      "aapg_geologic_province",
      "distance_to_nearest_pad_m",  # calculated in group_into_well_pads.R
      "tot_count_wells_within_10km" # calculated in group_into_well_pads.R
    )) %>%
    dplyr::rename(basin = aapg_geologic_province)

  well_pad_info <- arrow::open_dataset(header_dir) %>%
    dplyr::select(entity_id, first_prod_date, state) %>%
    dplyr::collect() %>%
    dplyr::filter(!is.na(first_prod_date)) %>%
    dplyr::mutate(age_yr = as.numeric(!!as_of_date - first_prod_date) / 365.25) %>%
    dplyr::select(-first_prod_date) %>%
    # 1:1 join, no column conflicts, type of "by" must be the same, and all
    # rows of right side must have a match
    powerjoin::power_inner_join(monthly_prod, by="entity_id", check=merge_specs(unmatched_keys_left="ignore")) %>%
    # 1:M join (U), no column conflicts (C), type of "by" must be the same (T)
    # We don't have all rows matching because:
    # - we filtered some wells by production above
    # - wells older than 1970 shouldn't show up in the well pad
    # - a small number of wells newer than 1970 aren't showing up
    # https://github.com/karldw/methane_abatement/issues/51
    powerjoin::power_inner_join(well_pad_mapping, by="entity_id",
      check=merge_specs(
        duplicate_keys_right = "ignore",
        unmatched_keys_left = "ignore",
        unmatched_keys_right = "ignore",
      )
    ) %>%
    # Harmonize basin name after the merge because the function needs a `state`
    # column
    harmonize_basin_name(group_small_CA_basins=TRUE) %>%
    # aggregate to well pads first:
    dplyr::group_by(well_pad_id)
  # Do a quick check here that we're matching the states we expect:
  expected_states <- c("CA", "CO", "NM")
  actual_states <- sort(unique(well_pad_info$state))
  if (! setequal(expected_states, actual_states)) {
    stop(
      "Expected state matches to be ", paste(expected_states, collapse=", "),
      " but they were actually ", paste(actual_states, collapse=", ")
    )
  }

  # This is vastly faster when not part of a pipeline. I don't know why.
  dplyr::summarize(well_pad_info,
    # varnames here need to match other data sources
    oil_avg_bbld = sum_(daily_avg_oil),
    gas_avg_mcfd = sum_(daily_avg_gas),
    age_yr = mean_(age_yr),
    basin = min(basin),
    distance_to_nearest_pad_m = mean_(distance_to_nearest_pad_m),
    tot_count_wells_within_10km = mean_(tot_count_wells_within_10km),
    .groups="drop"
  )
  # Note: could do histograms / sparklines of results here.
  # well_pad_info
}

count_gas_vs_nongas_wells <- function(prod_dir, well_pad_crosswalk, outfile) {
  as_of_date <- as.Date("2018-06-01")
  as_of_year <- lubridate::year(as_of_date)
  as_of_month <- lubridate::month(as_of_date)

  well_pad_mapping <- well_pad_crosswalk %>%
    arrow::read_parquet(col_select=c("entity_id", "well_pad_id"))

  well_prod <- arrow::open_dataset(prod_dir) %>%
    dplyr::filter(year == !!as_of_year, month == !!as_of_month) %>%
    dplyr::select(entity_id, gas) %>%
    dplyr::collect() %>%
    dplyr::filter(!is.na(gas)) %>%
    # 1:1 join, no column conflicts, type of "by" must be the same
    powerjoin::power_inner_join(well_pad_mapping, by="entity_id",
      check=merge_specs(
        unmatched_keys_left = "ignore",
        unmatched_keys_right = "ignore",
      )
    ) %>%
    dplyr::group_by(well_pad_id) %>%
    dplyr::summarize(gas = sum(gas), .groups="drop") %>%
    dplyr::mutate(any_gas = gas >= 0.01)

  stopifnot(noNAs(well_prod))
  # End in "%" to prevent extra whitespace
  mean_any_gas_pct <- paste0(signif(100 * mean(well_prod$any_gas), 3), "\\%%")
  writeLines(mean_any_gas_pct, outfile)

  invisible(well_prod)
}

calc_summaries_for_balance_tables <- function(df) {
  stopifnot(!is.null(df), nrow(df) > 10)
  df %<>% convert_km_to_nearest_well_pad()
  pretty_vars <- get_pretty_varnames(names(df))
  summary_stats <- list(
    mean = function(x) signif(mean_(x), 3),
    p10 = function(x) signif(quantile_(x, 0.1), 2),
    p90 = function(x) signif(quantile_(x, 0.9), 2)
  )
  df %>%
    dplyr::summarize_at(
      names(pretty_vars),
      summary_stats,
    )  %>%
    tidyr::pivot_longer(
      everything(),
      names_pattern = "(.+)_([a-zA-Z0-9]+)$",
      names_to = c("var", "stat")
    ) %>%
    tidyr::pivot_wider(
      id_cols = c("var"),
      names_from="stat"
    ) %>%
    # finally, rename so format_estimate_above_interval finds the right
    # variables. Note: these aren't confidence intervals, but it's easier
    # to set these names and reuse code.
    dplyr::transmute(
      term = pretty_vars[var],
      estimate = mean,
      conf_low = p10,
      conf_high = p90,
    )
}

balance_table_by_leak_status <- function(df_list, output_file) {
  stopifnot(length(output_file) == 1)
  df_list <- df_list[c("jpl_wells_all", "four_corners_all_wells", "lyon")] %>%
    purrr::map(add_emiss_vars) %>%
    purrr::map(~dplyr::select(., -any_of("emiss_kg_hr")))

  stopifnot(purrr::map_lgl(df_list, ~noNAs(.$detect_emiss)))

  # Six columns, with and without leaks, for three studies
  summary_lst <- list(
    # These names are written to the tex file as a comment to help avoid mistakes
    "California with emiss"   = dplyr::filter(df_list$jpl_wells_all,           detect_emiss),
    "California no emiss"     = dplyr::filter(df_list$jpl_wells_all,          !detect_emiss),
    "Four Corners with emiss" = dplyr::filter(df_list$four_corners_all_wells,  detect_emiss),
    "Four Corners no emiss"   = dplyr::filter(df_list$four_corners_all_wells, !detect_emiss),
    "Lyon with emiss"         = dplyr::filter(df_list$lyon,                    detect_emiss == 1),
    "Lyon no emiss"           = dplyr::filter(df_list$lyon,                    detect_emiss == 0)
  ) %>%
  purrr::map(calc_summaries_for_balance_tables)

  tab <- summary_lst %>%
    purrr::map(format_estimate_above_interval, align="@") %>%
    merge_estimates_df() %>%
    make_table_fragment(escaped=FALSE, add_comments = names(summary_lst)) %>%
    writeLines(output_file)
  invisible(summary_lst)
}

find_within_basin_non_match <- function(df, all_possible) {
  # Don't have to do this separately for each basin, since each well belongs to
  # exactly one basin.
  acceptable_basins <- dplyr::distinct(df, basin)
  stopifnot(noNAs(acceptable_basins), nrow(all_possible) > 0)
  out <- all_possible %>%
    dplyr::inner_join(acceptable_basins, by="basin") %>%
    dplyr::anti_join(df, by="well_pad_id")
  out
}

balance_table_by_flyover_status <- function(df_list, output_file) {
  stopifnot(length(output_file) == 1)
  non_matches <- df_list[c("jpl_wells_all", "four_corners_all_wells")] %>%
    purrr::map(find_within_basin_non_match, all_possible=df_list$all_possible_wells)
  summary_lst <- list(
      "California flown over"       = df_list$jpl_wells_all,
      "California not flown over"   = non_matches$jpl_wells_all,
      "Four Corners flown over"     = df_list$four_corners_all_wells,
      "Four Corners not flown over" = non_matches$four_corners_all_wells
    ) %>%
    purrr::map(calc_summaries_for_balance_tables)

  tab <- summary_lst %>%
    purrr::map(format_estimate_above_interval, align="@") %>%
    merge_estimates_df() %>%
    make_table_fragment(escaped=FALSE, add_comments = names(summary_lst))
  writeLines(tab, output_file)
  invisible(summary_lst)
}


if (! exists("snakemake")) {
  message("Using placeholder snakemake")
  TEX_FRAGMENTS <- fs::fs_path(here::here("output/tex_fragments"))
  snakemake <- SnakemakePlaceholder(
    input = list(
      # natural_gas_prices="data/eia/NG_PRI_FUT_S1_D.xls",
      # price_index = "data/fred/CPILFENS.csv",
      headers = "data/generated/production/well_headers/",
      prod = "data/generated/production/monthly_production/",
      well_pad_crosswalk = "data/generated/production/well_pad_crosswalk_1970-2018.parquet",
      # study_files = c(
      #   "alvarez_2018" = "data/studies/alvarez_etal_2018/aar7204_Database_S1.xlsx",
      #   "omara_2018" = "data/studies/omara_etal_2018/Omara_etal_SI_tables.csv",
      #   "duren_2019" = "data/studies/duren_etal_2019/Plume_list_20191031.csv"
      # ),
      matched_leakage = "data/generated/methane_measures/matched_wells_all.parquet",
      ground_studies = "data/generated/methane_measures/ground_studies.parquet",
      lyon_etal_2016 = "data/generated/methane_measures/lyon_etal_2016.parquet"
    ),
    output = list(
      # Order matters here (we'll index by position later)
      well_summary_stats = c(
        TEX_FRAGMENTS / "well_summary_stats_aviris_matched.tex", # 1
        TEX_FRAGMENTS / "well_summary_stats_lyon.tex",           # 2
        TEX_FRAGMENTS / "well_summary_stats_all_2018.tex"        # 3
      ),
      covariate_balance_by_leak = TEX_FRAGMENTS / "well_covariate_balance.tex",
      covariate_balance_by_flyover = TEX_FRAGMENTS / "well_covariate_balance_by_flyover.tex",
      # Order matters here (we'll index by position later)
      obs_counts = c(
        TEX_FRAGMENTS / "num_rows_aviris_matched.tex",
        TEX_FRAGMENTS / "num_rows_aviris_with_leak.tex",
        TEX_FRAGMENTS / "num_rows_lyon.tex",
        TEX_FRAGMENTS / "num_rows_all_2018.tex"
      ),
      percent_wells_with_gas = TEX_FRAGMENTS / "summ_percent_wells_nonzero_gas_in_us.tex"
    )
  )
}

header_dir <- snakemake@input[["headers"]]
prod_dir <- snakemake@input[["prod"]]


if (!exists("df_list")) {
  # for debugging speed.
  df_list <- prep_measurement_data_extra(
    filename_aviris=snakemake@input[["matched_leakage"]],
    filename_ground_studies=snakemake@input[["ground_studies"]]
  )
  df_list$lyon <- arrow::read_parquet(snakemake@input[["lyon_etal_2016"]])
  # Important note:
  # load_all_wells will load all wells *for the states that have been processed*
  # Currently we only process CO, NM, and CA, since those are the states we're
  # using.
  # Note: we might want to keep same-basin wells, rather than same-state ones.
  df_list$all_possible_wells <- load_all_wells(prod_dir, header_dir, snakemake@input[["well_pad_crosswalk"]])
}
purrr::walk2(
  # NOTE: this order has to be the same as the well_summary_stats list in snakemake
  df_list[c("aviris_all", "lyon", "all_possible_wells")],
  snakemake@output[["well_summary_stats"]],
  summary_stats_wells
)

write_n_obs(
  # NOTE: this order has to be the same as the obs_counts list in snakemake
  df_list[c("aviris_all", "aviris_trunc", "lyon", "all_possible_wells")],
  snakemake@output[["obs_counts"]]
)

balance_table_by_leak_status(df_list, snakemake@output$covariate_balance_by_leak)
balance_table_by_flyover_status(df_list, snakemake@output$covariate_balance_by_flyover)
count_gas_vs_nongas_wells(prod_dir, snakemake@input$well_pad_crosswalk,
  snakemake@output$percent_wells_with_gas
)

rm(snakemake) # easier debugging

