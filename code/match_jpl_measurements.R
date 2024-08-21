suppressMessages(
  here::i_am("code/match_jpl_measurements.R", uuid="bb6782d2-63d8-4013-a394-b10ba4a71e55")
)

source(here::here("code/shared_functions.r"))

options(scipen=5, warn=1)
set.seed(6350) # CRS as seed? Sure.

# Note on variable names:
# mcfd is mcf per day (thousands of cubic feet per day), but is reported monthly
# bbld is barrels per day, but is measured monthly.
# There are other variables that report the average on days when production was
# happening, but I don't use those (maybe I should?).

LON_LAT_CRS <- sf::st_crs(CRS_LONGLAT_int) # https://epsg.io/4326
OUTPUT_CRS  <- sf::st_crs(CRS_PROJECT_int) # https://epsg.io/6350 (CONUS Albers)
MAX_ACCEPTABLE_MATCH_DIST_METERS <- 500
# These are the source types we'll keep. There are others, like landfills and
# dairies, that we don't care about. This is a list of oil and gas types.
# I wasn't sure about gas compressors, but they seem to be far from wells.
SOURCE_TYPES <- c(
  # "gas compressor",
  # "gas distribution line",
  # "gas LNG station",
  # "gas processing plant",
  # "gas storage facility",
  "oil/gas compressor",
  "oil/gas drill rig",
  "oil/gas gathering line",
  "oil/gas possible plugged well",
  "oil/gas pumpjack",
  "oil/gas stack",
  "oil/gas tank",
  "oil/gas unknown infrastucture",
  "oil/gas waste lagoon"
)

read_jpl_plumes <- function(input_files) {
  df <- data.table::fread(input_files[["duren_2019_plumes"]], data.table=FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::rename_all(make_better_names) %>%
    dplyr::filter(source_type_best_estimate %in% !!SOURCE_TYPES) %>%
    dplyr::mutate(
      datetime_of_detection = lubridate::mdy_hm(paste(date_of_detection, time_of_detection_utc), tz="UTC"),
      date_of_detection = lubridate::mdy(date_of_detection),
      flight_name = stringr::str_match(candidate_identifier, "^(ang\\d{8}t\\d{6})")[, 2L, drop=TRUE]
    ) %>%
    dplyr::rename(
      emiss_kg_hr = qplume_kg_hr_plume_emissions,
      emiss_se_kg_hr = sigma_qplume_kg_hr_uncertainty_for_plume_emissions,
    ) %>%
    dplyr::select(-sectors_ipcc, -time_of_detection_utc) %>%
    ensure_id_vars(source_identifier, datetime_of_detection) %>%
    sf::st_as_sf(
      crs=LON_LAT_CRS,
      coords=c("plume_longitude_deg", "plume_latitude_deg"),
      agr=c(
        source_identifier = "identity", candidate_identifier = "identity",
        date_of_detection = "constant", source_type_best_estimate = "constant",
        emiss_kg_hr = "constant", emiss_se_kg_hr = "constant",
        datetime_of_detection = "constant", flight_name = "constant"
      )
    ) %>%
    sf::st_transform(OUTPUT_CRS)
  stopifnot(!anyNA(df$flight_name))
  df
}

read_jpl_sites <- function(input_files) {
  # Sites are different than plumes because they revisited
  # This excel file is translated from a PDF in the Durent et al. (2019)
  # supplementary materials.
  df <- readxl::read_excel(input_files[["duren_2019_sites"]]) %>%
    dplyr::rename_all(make_better_names) %>%
    dplyr::filter(source_type %in% !!SOURCE_TYPES) %>%
    dplyr::rename(
      emiss_kg_hr = qsource_kg_hr,
      persistence_frac = source_persistence_f,
    ) %>%
    # Run this rename twice because the sigma character gets assigned a different name on windows.
    rename_cols(c("emiss_se_kg_hr" = "x_q_kg_hr"), strict=FALSE) %>%
    rename_cols(c("emiss_se_kg_hr" = "x_u_f073_q_kg_hr"), strict=FALSE) %>%
    dplyr::select(-ipcc_sector, -confidence_in_persistence) %>%
    ensure_id_vars(source_identifier) %>%
    sf::st_as_sf(
      crs=LON_LAT_CRS,
      coords=c("source_longitude_deg", "source_latitude_deg"),
      agr=c(
        source_identifier = "identity", source_type = "constant",
        n_overflights = "constant", persistence_frac = "constant",
        emiss_kg_hr = "constant", emiss_se_kg_hr = "constant"
      )
    ) %>%
    sf::st_transform(OUTPUT_CRS)
  df
}


read_headers <- function(years, states) {
  header_dir <- here::here("data/generated/production/well_headers/")
  stopifnot(dir.exists(header_dir), !anyNA(years), !anyNA(states))
  date_min <- lubridate::make_date(min(years), 1, 1)
  date_max <- lubridate::make_date(max(years) + 1, 1, 1)

  well_headers <- arrow::open_dataset(header_dir) %>%
    dplyr::filter(
      state %in% !!states,
      # !(is.na(first_60_oil) & is.na(first_60_gas)),
      # !is.na(completion_date)
    ) %>%
    dplyr::select(
      county, state, production_type, drill_type, aapg_geologic_province,
      surface_latitude_wgs84, surface_longitude_wgs84, first_prod_date, last_prod_date,
      first_60_oil, first_60_gas, completion_date, spud_date, months_produced, entity_id,
    ) %>%
    dplyr::collect() %>%
    dplyr::filter(
      !is.na(surface_latitude_wgs84),
      !is.na(surface_longitude_wgs84),
      !is.na(first_prod_date),
      !is.na(last_prod_date),
      # Drop wells that started producing after the end of observation or ended
      # before the beginning (just for efficiency of not doing spatial operations
      # on a bunch of irrelevant wells)
      first_prod_date < !!date_max,
      last_prod_date > !!date_min,
      !production_type %in% c("WATER", "STEAMFLOOD"),
    ) %>%
    dplyr::rename(basin = aapg_geologic_province) %>%
    dplyr::mutate(production_type = dplyr::case_when(
      production_type %in% c("GAS", "OIL") ~ production_type,
      production_type == "GAS STORE" ~ "GAS",
      production_type == "OIL (CYCLIC STEAM)" ~ "OIL",
      TRUE ~ NA_character_
    )) %>%
    ensure_id_vars(entity_id) %>%
    sf::st_as_sf(
      crs=LON_LAT_CRS,
      coords=c("surface_longitude_wgs84", "surface_latitude_wgs84"),
      agr=c(
        entity_id = "identity", county = "constant", production_type = "constant",
        drill_type = "constant", basin = "constant", months_produced = "constant",
        surface_latitude_wgs84 = "constant", surface_longitude_wgs84 = "constant",
        first_prod_date = "constant", last_prod_date = "constant",
        first_60_oil = "constant", first_60_gas = "constant",
        completion_date = "constant", spud_date = "constant"
      )
    ) %>%
    sf::st_transform(OUTPUT_CRS)
  well_headers
}

match_wells_to_flight_day <- function(well_headers, flight_paths) {
  wells_flown_over <- well_headers %>%
    ensure_id_vars(entity_id) %>%
    # Only consider wells that are in the flight paths *and* are producing when
    # the plane flys over. This join generates duplicate rows wells because
    # the there are multiple flights over the same areas. Keep rows that had
    # any overpass while active, then drop duplicates.
    # (doing it this way because distinct doesn't work quite the same on sf objects)
    sf::st_join(dplyr::select(flight_paths, flight_date), join=sf::st_intersects, left=FALSE) %>%
    dplyr::filter(flight_date >= first_prod_date, flight_date <= last_prod_date) %>%
    sf::st_drop_geometry() %>%
    dplyr::distinct(entity_id, flight_date)

  dplyr::inner_join(well_headers, wells_flown_over, by="entity_id")
}

st_nn <- function(...){
  # Suppress nngeo's narration
  suppressMessages(nngeo::st_nn(...))
}

load_flight_paths <- function(plume_measurements) {
  # Flight paths:
  flight_df <- data.table::fread(
      here::here("data/studies/duren_etal_2019/AVIRIS-NG Flight Lines - AVIRIS-NG Flight Lines.csv"),
      data.table=FALSE
    ) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(flight_name = Name) %>%
    dplyr::rename_all(make_better_names) %>%
    dplyr::mutate(
      flight_date = lubridate::make_date(.data$year, .data$month, .data$day),
      utc_hour = dplyr::if_else(dplyr::between(utc_hour, 0, 23), utc_hour, NA_integer_),
      utc_minute = dplyr::if_else(dplyr::between(utc_minute, 0, 23), utc_minute, NA_integer_),
      flight_datetime_utc = lubridate::make_datetime(
        .data$year, .data$month, .data$day, .data$utc_hour, .data$utc_minute
      ),
    ) %>%
    dplyr::semi_join(plume_measurements, by="flight_name") %>% # keep matching flights
    dplyr::select(flight_name, flight_date, flight_datetime_utc, number_of_lines, investigator,
      dplyr::starts_with("lat"), dplyr::starts_with("lon")
    ) %>%
    dplyr::distinct() %>%
    ensure_id_vars(flight_name)
  stopifnot(setequal(flight_df$flight_name, plume_measurements$flight_name))

  flight_polygons <- flight_df %>%
    dplyr::select(flight_name, dplyr::starts_with("lon"), dplyr::starts_with("lat")) %>%
    # This is terrible.
    tidyr::pivot_longer(dplyr::starts_with("lon"), names_to="idx1", names_prefix="lon", values_to="lon") %>%
    tidyr::pivot_longer(dplyr::starts_with("lat"), names_to="idx2", names_prefix="lat", values_to="lat") %>%
    dplyr::filter(idx1 == idx2) %>%
    dplyr::filter(!is.na(lat), !is.na(lon)) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = LON_LAT_CRS) %>%
    dplyr::group_by(flight_name) %>%
    dplyr::summarize(
      geometry = sf::st_cast(sf::st_combine(geometry), "POLYGON"),
      .groups="drop"
    ) %>%
    sf::st_transform(OUTPUT_CRS)

  # Join back the variables (flight_polygons has to be first to get the s3 dispatch right)
  out <- dplyr::inner_join(flight_polygons,
    dplyr::select(flight_df, flight_name, flight_date, flight_datetime_utc, number_of_lines, investigator),
    by="flight_name"
  )
  if ("date_of_detection" %in% colnames(plume_measurements)) {
    # flight dates are detection dates:
    stopifnot(setequal(out$flight_date, plume_measurements$date_of_detection))
  }
  stopifnot(
    # detected plumes are inside flight paths:
    all(sf::st_intersects(plume_measurements, sf::st_union(out), sparse=FALSE))
  )
  out
}

load_monthly_prod <- function(wells_by_flight_day) {
  if (inherits(wells_by_flight_day, "sf")) {
    wells_by_flight_day %<>% sf::st_drop_geometry()
  }
  prod_dir <- here::here("data/generated/production/monthly_production/")
  stopifnot(
    dir.exists(prod_dir),
    is_id(wells_by_flight_day, entity_id, flight_date)
  )
  well_to_match <- wells_by_flight_day %>%
    dplyr::transmute(
      flight_date = flight_date,
      year = as.integer(lubridate::year(flight_date)),
      month = as.integer(lubridate::month(flight_date)),
      entity_id = entity_id
  )
  # These are just for speed, so we don't read anything we're positive we don't want
  desired_entity_id <- unique(well_to_match$entity_id)
  desired_year <- unique(well_to_match$year)
  desired_month <- unique(well_to_match$month)

  # Only diff from the stored schema is year, which currently registers as string
  # To select another column, add it to this schema.
  schema <- arrow::schema(
    year = arrow::int32(), month = arrow::int32(), entity_id = arrow::int32(),
    daily_avg_gas = arrow::float(), daily_avg_oil = arrow::float(), well_count = arrow::int32()
  )
  prod <- arrow::open_dataset(prod_dir, schema=schema) %>%
    dplyr::filter(
      daily_avg_gas > 0,
      year %in% !!desired_year,
      month %in% !!desired_month,
      entity_id %in% !!desired_entity_id,
    ) %>%
    dplyr::collect() %>%
    powerjoin::power_inner_join(
      well_to_match,
      by=c("entity_id", "year", "month"),
      check=merge_specs(
        duplicate_keys_right = "ignore",
        unmatched_keys_left = "ignore",
        unmatched_keys_right = "ignore",
      )
    )
  prod
}


load_well_records <- function(flight_paths, states, nat_gas_price) {
  years <- flight_paths$flight_date %>% lubridate::year() %>% unique()
  # First, find wells that were flown over and had production start before the
  # flight date and end after the flight date.
  # wells_by_flight_day has one row for each well each day it was flown over.
  wells_by_flight_day <- read_headers(years, states) %>%
    match_wells_to_flight_day(flight_paths)

  # Note: If the well was flown over in different months, this is the average
  # across months.
  monthly_prod_when_flown_over <- load_monthly_prod(wells_by_flight_day) %>%
    dplyr::group_by(entity_id) %>%
    dplyr::summarize(
      oil_avg_bbld = mean_(daily_avg_oil),
      gas_avg_mcfd = mean_(daily_avg_gas),
      .groups="drop"
    )
  # Recall wells_by_flight_day has one row for each well each day it was flown over.
  # Drop down to one row per well.
  # Do both well age and price here because it's easier than handling prices later
  well_age_and_price <- wells_by_flight_day %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(age_yr = as.numeric(flight_date - first_prod_date) / 365.25) %>%
    match_full_state_names() %>%
    dplyr::rename(state = state_full) %>%
    harmonize_basin_name() %>%
    match_commodity_prices(nat_gas_price) %>%
    dplyr::group_by(entity_id) %>%
    dplyr::summarize(
      age_yr = mean(age_yr),
      gas_price_per_mcf = mean(gas_price_per_mcf),
      gas_frac_methane = mean(gas_frac_methane),
      .groups="drop"
    )
  stopifnot(noNAs(well_age_and_price))

  wells <- wells_by_flight_day %>%
    dplyr::group_by(entity_id) %>%
    dplyr::select(-flight_date) %>%
    dplyr::distinct(.keep_all=TRUE) %>% # could speed this up here by not doing a spatial distinct
    ensure_id_vars(entity_id) %>%
    dplyr::inner_join(well_age_and_price, by="entity_id") %>% # 1:1 join
    harmonize_basin_name(group_small_CA_basins=TRUE)
  dplyr::inner_join(wells, monthly_prod_when_flown_over, by="entity_id") # 1:1 join
}

match_with_wells <- function(observed_sites, wells) {
  sites_matched <- sf::st_join(
      observed_sites, wells,
      join=st_nn, k=1, maxdist=MAX_ACCEPTABLE_MATCH_DIST_METERS,
      left=FALSE,
      progress=FALSE
    ) %>%
    # drop doubly-matched sites.
    dplyr::group_by(entity_id) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup() %>%
    ensure_id_vars(entity_id)
  if (nrow(sites_matched) != nrow(observed_sites)) {
    warning(paste0(
      nrow(observed_sites) - nrow(sites_matched),
      " leaks didn't have a matching well"
    ))
  }
  sites_matched
}


aggregate_to_well_pads <- function(wells_all, well_pad_mapping_file) {
  stopifnot(
    anyDuplicated(wells_all$entity_id) == 0,
    length(well_pad_mapping_file) == 1
  )
  well_pad_mapping <- arrow::read_parquet(
    well_pad_mapping_file,
    col_select=c("entity_id", "well_pad_id", "well_pad_lon", "well_pad_lat", "tot_count_wells_within_10km", "distance_to_nearest_pad_m")
  )
  # We drop study-specific columns here, but might want to bring them back.
  # n_overflights, persistence_frac, ...
  # NOTE: we're doing an inner join with well_pad_id here. That means wells that
  # aren't part of the well pad mapping will not be part of the output.
  well_pad_df <- wells_all %>%
    # Note: not currently true that all LHS are matched -- inform instead of erroring.
      powerjoin::power_inner_join(
        well_pad_mapping,  by="entity_id",
        check=merge_specs(
          duplicate_keys_left = "ignore",
          unmatched_keys_left = "inform",
          unmatched_keys_right = "ignore",
        )
    ) %>%
    dplyr::group_by(well_pad_id) %>%
    dplyr::summarize(
      # These aggregating functions are defined in shared_functions.r
      emiss_kg_hr = mean_(emiss_kg_hr),
      emiss_se_kg_hr = mean_(emiss_se_kg_hr),
      production_type = Mode(production_type),
      drill_type = Mode(drill_type),
      first_60_oil = sum_(first_60_oil),
      first_60_gas = sum_(first_60_gas),
      oil_avg_bbld = sum_(oil_avg_bbld),
      gas_avg_mcfd = sum_(gas_avg_mcfd),
      months_produced = mean_(months_produced),
      county = Mode(county),
      # Like Lyon et al. 2016, we'll define pad age by age of the most recently
      # drilled well.
      age_yr = min(age_yr),
      gas_price_per_mcf = mean(gas_price_per_mcf), # probably unique by well pad
      # already unique by well_pad_id; the particular aggregation doesn't matter
      well_pad_lon = min(well_pad_lon),
      well_pad_lat = min(well_pad_lat),
      basin = min(basin),
      gas_frac_methane = min(gas_frac_methane),
      tot_count_wells_within_10km = mean_(tot_count_wells_within_10km),
      distance_to_nearest_pad_m = mean_(distance_to_nearest_pad_m),
      .groups="drop"
    )

  stopifnot(!anyNA(well_pad_df$well_pad_id))

  # NOTES on well pads:
  # - There's still a vast number of well pads with no measurements
  # (currently 112 with measurements and 17006 without. This could be improved
  # slightly by thinking about deduplication, but it's minimal.)
  well_pad_df
}

load_ground_studies <- function(input_files) {
  # Drop studies that can't be used for our application. See the notes in
  # notes/measurement_paper_notes.csv
  studies_that_dont_work_for_us <- c("Rella et al.")

  df_alvarez <- readxl::read_excel(input_files[["alvarez_2018"]], sheet="inputs_meas_sites") %>%
    dplyr::rename_all(make_better_names) %>%
    dplyr::select(basin, study, methane_emissions_kgh, gas_production_mcfd) %>%
    dplyr::rename(emiss_kg_hr = methane_emissions_kgh) %>%
    dplyr::mutate(
      study = dplyr::case_when(
        study == "Omara et al." ~ "Omara et al. (2016)",
        study == "Robertson et al." ~ "Robertson et al. (2017)",
        TRUE ~ study
      ),
      basin = dplyr::case_when(
        basin == "SWPA" ~ "SW Pennsylvania",
        basin == "Weld County" ~ "Denver Julesburg",
        TRUE ~ basin
      )
    )

  df_omara_2018 <- data.table::fread(input_files[["omara_2018"]], data.table=FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::transmute(
      basin = dplyr::case_when(
        site == "Uinta Basin (Uintah County, UT)" ~ "Uinta", # match with Alvarez et al.
        site == "Denver Julesburg Basin (Weld County, CO)" ~ "Denver Julesburg",
        site == "NE PA (Bradford, Susquehanna, Wyoming, Sullivan Counties)" ~ "NE Pennsylvania",
        TRUE ~ site
      ),
      study = "Omara et al. (2018)",
      gas_production_mcfd = tot_gas_prod_mcfd,
      emiss_kg_hr = avg_emissions_kg_per_h,
    )

  df_zavala_araiza_2018 <- readxl::read_excel(input_files[["zavala_araiza_2018"]]) %>%
    dplyr::rename_all(make_better_names) %>%
    dplyr::filter(method == "tracer flux") %>% # different sampling strategies
    dplyr::transmute(
      emiss_kg_hr = ch4_emission_rate_kgh,
      gas_production_mcfd = gas_production_mcfd,
      study = "Zavala-Araiza et al. (2018)",
      basin = "Alberta",
      # other variables that could be interesting:
      # ch4_lb, ch4_ub, oil_production_bbld, age_yr, wells_per_site,
      # reported_emissions_kgh, gas_composition_c1_percent
    )
  out <- dplyr::bind_rows(df_alvarez, df_omara_2018, df_zavala_araiza_2018)
  stopifnot(all(studies_that_dont_work_for_us %in% out$study)) # guard against typos
  out <- dplyr::filter(out, !study %in% !!studies_that_dont_work_for_us)
  out
}

compare_ground_studies_with_jpl <- function(wells_all, ground_studies) {
  # Table:
  # - N non-zero
  # - N zero
  # - mean, median leak rate
  # - correlation with size, among non-zero leaks
  # - age?
  wells_num_zero <- dplyr::filter(wells_all, is.na(emiss_kg_hr) | emiss_kg_hr < 5) %>% nrow()
  wells_all_stats <- dplyr::filter(wells_all, emiss_kg_hr >= 5, !is.na(gas_avg_mcfd)) %>%
    dplyr::summarize(
      source = "JPL flights",
      n_positive = dplyr::n(),
      mean_if_positive = mean_(emiss_kg_hr),
      corr_with_size = stats::cor(emiss_kg_hr, gas_avg_mcfd),
      .groups="drop"
    ) %>%
    dplyr::mutate(n_zero = !!wells_num_zero)
  ground_studies_num_zero <- dplyr::filter(ground_studies, is.na(emiss_kg_hr) | emiss_kg_hr == 0) %>% nrow()
  ground_studies_stats <- dplyr::filter(ground_studies, emiss_kg_hr > 0) %>%
    dplyr::summarize(
      source = "Ground studies",
      n_positive = dplyr::n(),
      mean_if_positive = mean_(emiss_kg_hr),
      # corr_with_size = stats::cor(emiss_kg_hr, gas_avg_mcfd)
      corr_with_size = 0,
      .groups="drop"
    ) %>%
    dplyr::mutate(n_zero = !!ground_studies_num_zero)
  tab <- cbind(t(wells_all_stats), t(ground_studies_stats))
  return(tab)
}

read_lyon_etal_2016 <- function(input_files) {
  locations <- readxl::read_excel(input_files[["lyon_etal_2016_locations"]]) %>%
    dplyr::as_tibble() %>%
    dplyr::rename_all(make_better_names) %>%
    ensure_id_vars(pad_id) %>%
    dplyr::select(-basin) %>%
    dplyr::mutate(longitude = as.double(stringr::str_trim(longitude))) # excel issues

  df <- readxl::read_excel(input_files[["lyon_etal_2016_measures"]], sheet=1) %>%
    dplyr::as_tibble() %>%
    dplyr::rename_all(make_better_names) %>%
    dplyr::select(-video_id) %>%
    ensure_id_vars(pad_id) %>%
    dplyr::left_join(locations, by="pad_id") %>%
    # Rename to harmonize with other datasets
    dplyr::rename(
      gas_avg_mcfd = gas_production_mcf_pad_day,
      oil_avg_bbld = oil_production_bbl_pad_day,
      detect_emiss = emissions_detected_0_no_1_yes,
    ) %>%
    dplyr::mutate(
      # This is approx the same as our definition
      age_yr = well_age_months_since_initial_production_of_newest_well / 12,
      flight_date = as.Date("2014-08-01"),
      state = basin_to_state(basin),
    )

  df
}

basin_to_state <- function(basin) {
  stopifnot(is.character(basin))
  # These aren't perfect, but for basins that cross states, we pick the majority
  # state. (But if you're worried about this, you should worry about using
  # state citygate prices instead.)
  conversions = c(
    "Bakken" = "North Dakota",
    "Barnett" = "Texas",
    "EagleFord" = "Texas",
    "Fayetteville" = "Arkansas",
    "Marcellus" = "Pennsylvania",
    "PowderRiver" = "Wyoming",
    "Uintah" = "Utah",
    "San Juan" = "New Mexico",
    "San Joaquin" = "California"
  )
  missing_basins <- setdiff(unique(basin), names(conversions))
  if (length(missing_basins) > 0) {
    stop("Missing states for these basins: ", paste(missing_basins, collapse=", "))
  }
  conversions[basin]
}

read_frankenberg_etal_2016 <- function(input_files) {
  measures <- readxl::read_excel(input_files[["frankenberg_etal_2016_measures"]]) %>%
    dplyr::transmute(
      # Do the rename and select in one step. Drop "Thorpe ID" and
      # "Rank of size of flux for 245 sources (Frankenberg et al., 2016)."
      plume_id = `Thorpe ID`,
      latitude = Latitude,
      longitude = Longitude,
      source_type = `Source designation (Thorpe, performed after Frankenberg et al., 2016, total 273 sources with some repeat observations)`,
      emiss_kg_hr = `Estimated flux for 245 sources (Frankenberg et al., 2016). Units kg CH4/hr`
    ) %>%
    dplyr::slice(-1) %>% # omit spacer row
    dplyr::filter(
      not_na(emiss_kg_hr), # sources that had a plume, but it wasn't quantified.
      source_type %in% c("Tank", "Tanks", "Unknown", "Unknown facility",
        "Unknown infrastructure", "Well completion", "Wellpad infrastructure"
      ),
      emiss_kg_hr == 0 | emiss_kg_hr > 0.001 # one *very* small value seems like a mistake
    )
  # Source types in the data:
  # source_type                n
  # Coal mine vent shaft       1
  # Gas processing plant      10
  # Natural                    3
  # Pipeline                   3
  # Tank                      61
  # Tanks                      1
  # Unknown                   10
  # Unknown facility          17
  # Unknown infrastructure     2
  # Well completion            1
  # Wellpad infrastructure   135

  sources <- readxl::read_excel(input_files[["frankenberg_etal_2016_sources"]]) %>%
    dplyr::rename_all(make_better_names)
  multiple_divider_row <- which(sources$x == "Multiple Overpasses")
  total_rows <- nrow(sources)
  stopifnot(identical(multiple_divider_row, 179L), total_rows > 179)
  sources <- dplyr::slice(sources, -multiple_divider_row) %>%
    dplyr::select(longitude, latitude, file) %>%
    dplyr::mutate(
      flight_name =
        stringr::str_match(.data$file, "^(ang\\d{8}t\\d{6})_ch4_v1e_img$")[, 2, drop=TRUE],
      flight_date = lubridate::ymd(
        stringr::str_match(.data$file, "^ang(\\d{8})t\\d{6}_ch4_v1e_img$")[, 2, drop=TRUE]
      )
    ) %>%
    ensure_id_vars(longitude, latitude)
  stopifnot(
    !anyNA(sources$flight_name),
    nrow(dplyr::anti_join(measures, sources, by=c("longitude", "latitude"))) == 0
  )
  # For these two, longitude and latitude match exactly because they're from the
  # same underyling measurements. Some sources were detected but don't have
  # quantified measurements.
  out <- dplyr::left_join(measures, sources, by=c("longitude", "latitude")) %>%
    # Drop 4 obs that don't have associated flight paths in the online data.
    # (They're just east of the flight main zone; not sure if they were recorded the same)
    dplyr::filter(! flight_name %in% c("ang20150421t160633", "ang20150423t150648")) %>%
    dplyr::select(-file) %>%
    sf::st_as_sf(
      crs=LON_LAT_CRS,
      coords=c("longitude", "latitude"),
      agr=c(
        plume_id = "identity", source_type = "constant", emiss_kg_hr = "aggregate",
        flight_name = "constant", flight_date = "constant"
      )
    ) %>%
    sf::st_transform(OUTPUT_CRS) %>%
    # Add all-NA variables to make it easier to share some functions
    dplyr::mutate(
      emiss_se_kg_hr = NA_real_,
    )
  out
}

write_datasets <- function(data_lst, output_lst) {
  data_names <- names(data_lst)
  stopifnot(
    length(data_names) == length(data_lst),
    all(data_names %in% names(output_lst))
  )
  for (nm in data_names) {
    df <- data_lst[[nm]]
    if (inherits(df, "sf")) {
      df <- geometry_to_lonlat(df)
    }
    arrow::write_parquet(df, output_lst[[nm]])
  }
  NULL
}

standardize_columns_for_plotting <- function(df_list, censor_threshold=5) {
  # NOTE: in the output of emiss_kg_hr, zero means emissions could have been
  # detected and quantified but were not. NA means emissions could not have been
  # quantified.
  df_list$jpl_wells_all %<>%
    dplyr::transmute(
      emiss_kg_hr = dplyr::if_else(is.na(emiss_kg_hr), 0, emiss_kg_hr),
      detect_emiss = emiss_kg_hr > 0,
      src = "California",
    )
  df_list$four_corners_all_wells %<>%
    dplyr::transmute(
      emiss_kg_hr = dplyr::if_else(is.na(emiss_kg_hr), 0, emiss_kg_hr),
      detect_emiss = emiss_kg_hr > 0,
      src = "Four Corners",
    )
  df_list$lyon %<>%
    dplyr::transmute(
      emiss_kg_hr = NA_real_,
      detect_emiss = detect_emiss == 1,
      src = "Lyon et al.",
    )
  df_list$ground_studies_censored_5kgh <- df_list$ground_studies %>%
    dplyr::transmute(
      emiss_kg_hr = dplyr::if_else(.data$emiss_kg_hr > 5, .data$emiss_kg_hr, 0),
      detect_emiss = .data$emiss_kg_hr > 0,
      src = "Ground studies censored at 5 kg/hr",
    )
  df_list$ground_studies_censored_10kgh <- df_list$ground_studies %>%
    dplyr::transmute(
      emiss_kg_hr = dplyr::if_else(.data$emiss_kg_hr > 10, .data$emiss_kg_hr, 0),
      detect_emiss = .data$emiss_kg_hr > 0,
      src = "Ground studies censored at 10 kg/hr",
    )
  df_list$ground_studies %<>%
    # dplyr::filter(emiss_kg_hr > 0) %>%
    dplyr::transmute(
      emiss_kg_hr = .data$emiss_kg_hr,
      detect_emiss = .data$emiss_kg_hr > 0,
      src = "Ground studies",
    )

  df_list
}

filter_within_distance <- function(x, y, max_dist=units::as_units(1000, "m")) {
  # Keep elements of x that are within 5000 m of any element of y
  y_buffer <- sf::st_geometry(y) %>%
    sf::st_union() %>%
    sf::st_buffer(dist=max_dist)
  stopifnot(length(y_buffer) == 1)
  sf::st_intersection(x, y_buffer)
}

data_to_check_matches <- function(snakemake) {
  outfile <- snakemake@output[["data_to_check_matches"]] %||% stop("Need outfile")
  # lyon <- read_lyon_etal_2016(snakemake@input)
  well_pad_mapping <- snakemake@input[["well_pad_crosswalk"]] %>%
    arrow::read_parquet(col_select=c("entity_id", "well_pad_id")) %>%
    ensure_id_vars(entity_id)

  jpl_sites <- read_jpl_sites(snakemake@input) %>%
    dplyr::rename(source_id = source_identifier)
  wells_ca <- read_headers(2016:2017, "CA") %>%
    filter_within_distance(jpl_sites) %>%
    dplyr::inner_join(well_pad_mapping, by="entity_id")

  four_corners_df <- read_frankenberg_etal_2016(snakemake@input) %>%
    dplyr::rename(source_id = plume_id)

  wells_nm_co <- read_headers(2015, c("CO", "NM")) %>%
    filter_within_distance(four_corners_df) %>%
    dplyr::inner_join(well_pad_mapping, by="entity_id")

  df_list <- readRDS(snakemake@output[["cleaned_matched_obs"]])
  jpl_matches <- df_list$jpl_wells_all %>%
    dplyr::filter(!is.na(emiss_kg_hr)) %>%
    dplyr::select(source_identifier, entity_id) %>%
    dplyr::rename(source_id = source_identifier)
  four_corners_matches <- df_list$four_corners_all_wells %>%
    dplyr::filter(!is.na(emiss_kg_hr)) %>%
    dplyr::select(plume_id, entity_id) %>%
    dplyr::rename(source_id = plume_id)
  stopifnot(
    anyDuplicated(jpl_matches$source_id) == 0,
    anyDuplicated(jpl_matches$entity_id) == 0,
    anyDuplicated(four_corners_matches$source_id) == 0,
    anyDuplicated(four_corners_matches$entity_id) == 0
  )
  jpl_sites %<>% dplyr::left_join(jpl_matches, by="source_id")
  four_corners_df %<>% dplyr::left_join(four_corners_matches, by="source_id")
  out <- list(
    ca    = list(measures = jpl_sites,       wells = wells_ca),
    co_nm = list(measures = four_corners_df, wells = wells_nm_co)
  )
  saveRDS(out, outfile)
  invisible(out)
}

match_full_state_names <- function(df) {
  state_names <- tibble::tibble(state_abb = state.abb, state_full = state.name)
  # 1:m merge, require a match of all rows in df
  powerjoin::power_inner_join(state_names, df,
    by=c("state_abb"="state"),
    check=merge_specs(
      duplicate_keys_right = "ignore",
      unmatched_keys_left = "ignore",
    )
  )
}

match_commodity_prices <- function(df, price_file) {
  stopifnot("basin" %in% names(df), length(price_file) == 1)
  prices <- arrow::read_parquet(price_file, col_select=c("date", "basin", "price_real")) %>%
    dplyr::rename(gas_price_per_mcf = price_real)
  out <- df %>%
    dplyr::mutate(date_monthly = first_of_month(flight_date)) %>%
    powerjoin::power_inner_join(
      # Check that the join is m:1, all rows of left have a match, and there are
      # no conflicts in columns or types.
      prices,
      by=c("basin"="basin", "date_monthly"="date"),
      check=merge_specs(
        duplicate_keys_left = "ignore",
        unmatched_keys_right = "ignore",
      )
    ) %>%
    dplyr::select(-date_monthly) %>%
    dplyr::mutate(
      gas_frac_methane = 0.95, # APPROX!
    )
  out
}

write_match_percent <- function(count_total, count_match, outfile) {
  stopifnot(
    count_total > 100, count_total < 10000, count_match > 0,
    count_match <= count_total, length(outfile) == 1
  )
  count_unmatch <- count_total - count_match
  pct_drop <- signif(100 * count_unmatch / count_total, 2)
  writeLines(paste0(pct_drop, "\\%%"), outfile)
  invisible(NULL)
}

match_jpl_california <- function(input_files) {
  # Note:
  # Some wells are flown over multiple times. Ideally, we would use that
  # information to estimate leak probabilities for these wells. However, for the
  # short run, we're going to count wells that weren't leaking on their first
  # overpass as not leaking, becuase we definitely want to avoid overstating the
  # leak rate (as using ever-leak status would do)
  #
  # The steps we follow are:
  # - Match each flight polygon with wells.
  # - For each matched well, only keep one randomly drawn flight
  # - Within groups of flight ID, match wells to detected plumes.
  # - Report counts of all plumes, matched plumes, and plumes that would have
  #   matched if we weren't dropping the later overpasses.
  plumes <- read_jpl_plumes(input_files) %>%
    dplyr::filter(!is.na(emiss_kg_hr), !is.na(emiss_se_kg_hr))
  flight_paths <- load_flight_paths(plumes)

  # Load the well covariates.
  # Note: this is could be improved a bit, because for repeat flyovers we want the
  # specific draw we have, but load_well_records returns the average for the
  # well across all visits.
  well_info <- load_well_records(flight_paths, "CA", input_files$nat_gas_prices)


  # This is like match_wells_to_flight_day, but different because we pick out
  # individual flights instead of flight days (some days have multiple flights)
  # All flights are in CA
  wells_flight_match <- read_headers(
      years=unique(lubridate::year(flight_paths$flight_date)),
      states="CA"
    ) %>%
    # Keep only wells that we kept in load_well_records (dropping gas == 0)
    dplyr::filter(.data$entity_id %in% well_info$entity_id) %>%
    # Only consider wells that are in the flight paths *and* are producing when
    # the plane flys over. This join generates duplicate rows wells because
    # the there are multiple flights over the same areas. Keep rows that had
    # any overpass while active, then drop duplicates.
    # (doing it this way because distinct doesn't work quite the same on sf objects)
    sf::st_join(
      dplyr::select(flight_paths, flight_date, flight_name),
      join=sf::st_intersects,
      left=FALSE
    ) %>%
    dplyr::filter(flight_date >= first_prod_date, flight_date <= last_prod_date) %>%
    dplyr::select(entity_id, flight_date, flight_name) %>%
    ensure_id_vars(entity_id, flight_date, flight_name)

  flight_names <- unique(wells_flight_match$flight_name)
  stopifnot(noNAs(wells_flight_match), length(flight_names) > 10)

  wells_rand_overpass <- wells_flight_match %>%
    dplyr::group_by(entity_id) %>%
    dplyr::slice_sample(n = 1) %>%
    dplyr::ungroup()

  .match_one_overpass <- function(flight_name, plumes, wells) {
    plumes %<>% dplyr::filter(.data$flight_name == !!flight_name) %>%
      dplyr::select(-flight_name)
    wells  %<>% dplyr::filter(.data$flight_name == !!flight_name)
    if (nrow(plumes) == 0 || nrow(wells) == 0) {
      return(NULL)
    }
    matched <- sf::st_join(plumes, wells,
        join=st_nn, k=1, progress=FALSE, left=FALSE
      ) %>%
      sf::st_drop_geometry() %>%
      # Now, if there are multiple plumes observed and matched to a well
      # _from the same flight_, then average them.
      dplyr::group_by(entity_id) %>%
      dplyr::summarize(
        emiss_kg_hr = mean_(emiss_kg_hr),
        emiss_se_kg_hr = mean_(emiss_se_kg_hr),
        .groups="drop"
      )
    all_wells_flown_over <- wells %>%
      sf::st_drop_geometry() %>%
      powerjoin::power_left_join(
        matched, by="entity_id",
        check=merge_specs(unmatched_keys_left = "ignore")
      )
    all_wells_flown_over
  }
  # Loop over flight_name and match wells and plumes for each flight.
  # Note: there's absolutely a better way to do this.
  all_wells_flown_over <- purrr::map_dfr(
    flight_names, .match_one_overpass,
    plumes=plumes, wells=wells_rand_overpass
  )
  stopifnot(
    anyDuplicated(all_wells_flown_over$entity_id) == 0,
    all(all_wells_flown_over$entity_id %in% well_info$entity_id)
  )
  all_wells_flown_over %<>% dplyr::inner_join(well_info, by="entity_id")
  matched_wells <- dplyr::filter(all_wells_flown_over, !is.na(emiss_kg_hr))

  observed_well_pads <- aggregate_to_well_pads(all_wells_flown_over, input_files$well_pad_crosswalk)

  could_have_been_matches <- sf::st_join(
    plumes, wells_flight_match,
    join=st_nn, k=1,
    maxdist=MAX_ACCEPTABLE_MATCH_DIST_METERS,
    left=FALSE,
    progress=FALSE
  )

  list(
    observed_well_pads = observed_well_pads,
    matched_wells = matched_wells,
    count_total = nrow(plumes),
    count_matched = sum(!is.na(all_wells_flown_over$emiss_kg_hr)),
    count_would_have_matched = nrow(could_have_been_matches)
  )
}


if (!exists("snakemake")) {
  message("This script is meant to be run with snakemake. Using a placeholder.")
  snakemake <- SnakemakePlaceholder(
    input = list(
      well_pad_crosswalk = "data/generated/production/well_pad_crosswalk_1970-2018.parquet",
      headers = glue::glue("data/generated/production/well_headers/first_prod_year={year}/file.parquet", year=1990:2018),
      prod = glue::glue("data/generated/production/monthly_production/year={year}/file.parquet", year=1990:2018),
      alvarez_2018 = "data/studies/alvarez_etal_2018/aar7204_Database_S1.xlsx",
      omara_2018 = "data/studies/omara_etal_2018/Omara_etal_SI_tables.csv",
      duren_2019_plumes = "data/studies/duren_etal_2019/Plume_list_20191031.csv",
      duren_2019_sites = "data/studies/duren_etal_2019/41586_2019_1720_MOESM3_ESM.xlsx",
      lyon_etal_2016_locations = "data/studies/lyon_etal_2016/es6b00705_si_005.xlsx",
      lyon_etal_2016_measures = "data/studies/lyon_etal_2016/es6b00705_si_004.xlsx",
      frankenberg_etal_2016_sources = "data/studies/frankenberg_etal_2016/AVNG_sources_all2.xlsx",
      frankenberg_etal_2016_measures = "data/studies/frankenberg_etal_2016/FourCorners_AV_NG_detections_Werner.xlsx",
      zavala_araiza_2018 = "data/studies/zavala-araiza_etal_2018/elementa-6-284-s1.xlsx",
      nat_gas_prices = "data/generated/nat_gas_prices_by_basin.parquet"
    ),
    output = list(
      # plot_obs_count_jpl = "graphics/observation_count_jpl_flights.pdf",
      # plot_jpl_flights_qq = "graphics/jpl_flights_qq_plot.pdf",
      lyon_etal_2016 = "data/generated/methane_measures/lyon_etal_2016.parquet",
      ground_studies = "data/generated/methane_measures/ground_studies.parquet",
      cleaned_matched_obs = "data/generated/methane_measures/matched_wells_all.parquet",
      aviris_match_fraction_dropped = "output/tex_fragments/intext_aviris_match_fraction_dropped.tex"
    ),
    threads = 4,
    resources = list(mem_mb = 7000),
    rule = ""
  )
}

main <- function(snakemake) {
  ground_studies <- load_ground_studies(snakemake@input) %>%
    arrow::write_parquet(snakemake@output[["ground_studies"]])

  lyon <- read_lyon_etal_2016(snakemake@input) %>%
    match_commodity_prices(snakemake@input$nat_gas_prices) %>%
    arrow::write_parquet(snakemake@output[["lyon_etal_2016"]])

  four_corners_df <- read_frankenberg_etal_2016(snakemake@input)
  jpl_ca_list <- match_jpl_california(snakemake@input)
  jpl_wells_all <- jpl_ca_list$observed_well_pads
  jpl_sites_matched <- jpl_ca_list$matched_wells

  # Keep track of the total number of sites and how many were matched.
  # (Note that read_jpl_sites already drops some non-O&G sites)
  # Note: these numbers are used to determine geographic match quality, so I'm
  # including counts that would have matched but for the random-overpass filter
  # (see details in match_jpl_california)
  count_aviris_total <- jpl_ca_list$count_total
  count_aviris_match <- jpl_ca_list$count_would_have_matched

  wells_nm_co <- load_flight_paths(four_corners_df) %>%
    load_well_records(c("CO", "NM"), snakemake@input$nat_gas_prices)

  four_corners_matched <- match_with_wells(four_corners_df, wells_nm_co)
  # Keep track of the total number of sites and how many were matched.
  count_aviris_total <- count_aviris_total + nrow(four_corners_df)
  count_aviris_match <- count_aviris_match + nrow(four_corners_matched)

  four_corners_vars_to_keep <- c(setdiff(colnames(four_corners_df), "geometry"), "entity_id")
  four_corners_all_wells <- four_corners_matched %>%
    sf::st_drop_geometry() %>%
    dplyr::select(!!four_corners_vars_to_keep) %>%
    dplyr::right_join(wells_nm_co, by="entity_id") %>%
    aggregate_to_well_pads(snakemake@input[["well_pad_crosswalk"]])

  main_df <- dplyr::bind_rows(jpl_wells_all, four_corners_all_wells)
  arrow::write_parquet(main_df, snakemake@output[["cleaned_matched_obs"]])

  write_match_percent(count_aviris_total, count_aviris_match,
    snakemake@output[["aviris_match_fraction_dropped"]]
  )
}

arrow::set_cpu_count(snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])


main(snakemake)
