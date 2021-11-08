
# Note: see bottom of this file for some globals that get defined.


harmonize_one_variable <- function(df, decoder, varname) {
  stopifnot(
    setequal(colnames(decoder), c("raw_value", "harmonized_value")),
    anyDuplicated(decoder$raw_value) == 0,
    length(varname) == 1
  )
  by_vars <- setNames("raw_value", varname)
  non_matches <- dplyr::anti_join(df, decoder, by=by_vars)
  if (nrow(non_matches) > 0) {
    print(dplyr::count(df, !!varname))
    stop("Failed to match decoder status for the above")
  }
  # First rename is e.g. "production_status" to "production_status_raw"
  # Second rename is "harmonized_value" to "production_status"
  rename1 <- setNames(varname, paste0(varname, "_raw"))
  rename2 <- setNames("harmonized_value", varname)
  out <- dplyr::left_join(df, decoder, by=by_vars) %>%
    rename_cols(rename1) %>%
    rename_cols(rename2)
  out
}

harmonize_production_status <- function(df) {
  # Pulling some status codes from http://www.ogb.state.ms.us/well_status_list.php
  decoder <- tibble::tribble(
    ~raw_value, ~harmonized_value,
    "AB", "ABANDONED",
    "ABANDONED", "ABANDONED",
    "ABANDONED, DRY HOLE", "ABANDONED",
    "AC", "ACTIVE",
    "ACTIVE", "ACTIVE",
    "ACTIVE GAS PROD", "ACTIVE",
    "ACTIVE INJ", "INJECTION",
    "ACTIVE PERMIT", "ACTIVE",
    "APPLICATION FOR PERMIT TO DRILL", "PRE-PRODUCTION",
    "CA", "CLOSED",
    "CANCEL", "CLOSED",
    "CANCELLED", "CLOSED",
    "CM", "OTHER", # ?? small number of these
    "CO", "OTHER", # ?? small number of these
    "COMPLETED", "PRE-PRODUCTION",
    "CONSTRUCTION", "PRE-PRODUCTION",
    "CONVERTED", "PRE-PRODUCTION", # ?? small number of these
    "DC", "OTHER", # ?? small number of these
    "DD", "OTHER", # ?? small number of these
    "DG", "PRE-PRODUCTION", # drilling
    "DH", "OTHER", # ?? small number of these
    "DISPOSAL", "INJECTION",
    "DRILLING", "PRE-PRODUCTION",
    "DRILLING ACTIVE", "PRE-PRODUCTION",
    "DRILLING SUSPENDED", "CLOSED",
    "DRL", "PRE-PRODUCTION",
    "DRY", "CLOSED",
    "DRY HOLE", "CLOSED",
    "EX", "CLOSED", # small number of these
    "EXPIRED", "CLOSED",
    "EXPIRED PERMIT", "CLOSED",
    "FUTURE USE", "PRE-PRODUCTION",
    "IN", "CLOSED",
    "INACTIVE", "CLOSED",
    "INJ", "INJECTION",
    "INJECTION", "INJECTION",
    "ISSUED", "PRE-PRODUCTION",
    "NEW DRILL", "PRE-PRODUCTION",
    "NOT DRILLED", "OTHER",
    "NOTICE OF INTENT TO ABANDON", "CLOSED",
    "NR", "OTHER", # no report    "ORND", "CLOSED",
    "ORPHAN", "CLOSED",
    "ORND", "CLOSED",
    "OTHER", "OTHER",
    "P & A", "CLOSED", # "plugged and abandoned"
    "P&A", "CLOSED",
    "PA", "CLOSED",
    "PB", "OTHER",
    "PERM", "PRE-PRODUCTION",
    "PERMIT", "PRE-PRODUCTION",
    "PERMIT ISSUED", "PRE-PRODUCTION",
    "PERMITTED", "PRE-PRODUCTION",
    "PM", "PRE-PRODUCTION",
    "PLUGGED", "CLOSED",
    "PLUGGED BACK", "CLOSED",
    "PLUGGING/PLU", "CLOSED",
    "PR", "OTHER",
    "RC", "OTHER", # ?? small number of these
    "RECLASSIFIED TO ANOTHER TYPE OF ENTITY", "OTHER",
    "RELEASED", "CLOSED",
    "SHUT IN", "CLOSED",
    "SHUT-IN", "CLOSED",
    "SI", "CLOSED",
    "SPUDDED", "PRE-PRODUCTION",
    "SUSP", "CLOSED",
    "SUSPENDED", "CLOSED",
    "TA", "CLOSED", # temporarily abandoned
    "TR", "OTHER", # ?? small number of these
    "UL", "OTHER", # ?? small number of these
    "UM", "OTHER", # ?? small number of these
    "UN", "OTHER", # unknown
    "VENTING", "ACTIVE",
    "VP", "OTHER", # ?? small number of these
    "WORKED OVER TO ANOTHER ZONE", "CLOSED",
    "WSW", "OTHER", # ?? small number of these
    NA_character_, "OTHER", # small number of these
    )

  harmonize_one_variable(df, decoder, "production_status")
}

help_message_production_cols_available <- function() {
  paste(
    "You must specify which columns you want to load. Here are the available columns:",
    '  "entity_id",',  # int64
    '  "api",',  # string
    '  "api_list",',  # string
    '  "oil",',  # float
    '  "gas",',  # float
    '  "water",',  # float
    '  "well_count",',  # int32
    '  "days",',  # int32
    '  "daily_avg_oil",',  # float
    '  "daily_avg_gas",',  # float
    '  "daily_avg_water",',  # float
    '  "reservoir",',  # string
    '  "well_name",',  # string
    '  "well_number",',  # string
    '  "operator_alias",',  # string
    '  "production_type",',  # string
    '  "production_status",',  # string
    '  "entity_type",',  # string
    '  "month" (specify "month" to load "prod_date")',  # int32
  )
}

load_production_by_filename <- function(files, ...) {
  years <- stringr::str_match(files, "/year=(\\d{4})/")[, 2, drop=TRUE] %>% as.integer()
  stopifnot(!anyNA(years))
  load_production_by_year(years, ...)
}

load_production_by_year <- function(years, col_select, production_status = NULL, filter_fn = NULL) {
  if (missing(col_select)) {
    stop(help_message_production_cols_available())
  }
  stopifnot(length(years) > 0)
  years <- as.integer(years)
  if (!is.null(production_status)) {
    stopifnot(all(production_status %in% c("ABANDONED", "ACTIVE", "CLOSED", "INJECTION", "OTHER", "PRE-PRODUCTION")))
  }
  force(col_select)
  if (!"production_status" %in% col_select && !is.null(production_status)) {
    stop("Filtering on well status requires loading the production_status column")
  }

  .read_one <- function(year, col_select, production_status, filter_fn) {
    path <- file.path(
      DATA_DIR,
      "generated/production/monthly_production",
      paste0("year=", year),
      "file.parquet"
    )
    if (!fs::file_exists(path)) {
      message("Skipping monthly production for ", year)
      return(NULL)
    }

    df <- arrow::read_parquet(path, col_select=col_select)
    if ("month" %in% colnames(df))
      df %<>% dplyr::mutate(prod_date = lubridate::make_date(year=!!year, month=month, day=1L)) %>%
        dplyr::select(-month)
    if ("production_status" %in% col_select) {
      df <- harmonize_production_status(df)
    }
    if (!is.null(production_status)) {
      df %<>% dplyr::filter(.data$production_status %in% !!production_status)
    }
    if (!is.null(filter_fn)) {
      df %<>% filter_fn()
    }
    df
  }
  # Could run in parallel, but takes a lot of memory.
  # furrr::future_map_dfr(years, .read_one, col_select=col_select, production_status=production_status)
  purrr::map_dfr(years, .read_one,
    col_select=col_select, production_status=production_status,
    filter_fn=filter_fn)
}

help_message_header_cols_available <- function() {
  paste(
    "You must specify which columns you want to load. Here are the available columns:",
    '  "api",',  # string
    '  "operator_alias_legacy",',  # string
    '  "operator_company_name",',  # string
    '  "operator_reported",',  # string
    '  "operator_ticker",',  # string
    '  "well_name",',  # string
    '  "well_number",',  # string
    '  "entity_type",',  # string
    '  "county",',  # string
    '  "di_basin",',  # string
    '  "di_play",',  # string
    '  "di_subplay",',  # string
    '  "reservoir",',  # string
    '  "production_type",',  # string
    '  "producing_status",',  # string
    '  "drill_type",',  # string
    '  "first_prod_date",',  # Date
    '  "last_prod_date",',  # Date
    '  "cum_gas",',  # float
    '  "cum_oil",',  # float
    '  "cum_boe",',  # float
    '  "cum_water",',  # float
    '  "cum_mmcfge",',  # float
    '  "cum_bcfge",',  # float
    '  "daily_gas",',  # float
    '  "daily_oil",',  # float
    '  "first_month_oil",',  # float
    '  "first_month_gas",',  # float
    '  "first_month_water",',  # float
    '  "first_6_oil",',  # float
    '  "first_6_gas",',  # float
    '  "first_6_boe",',  # float
    '  "first_6_water",',  # float
    '  "first_12_oil",',  # float
    '  "first_12_gas",',  # float
    '  "first_12_boe",',  # float
    '  "first_12_mmcfge",',  # float
    '  "first_12_water",',  # float
    '  "first_24_oil",',  # float
    '  "first_24_gas",',  # float
    '  "first_24_boe",',  # float
    '  "first_24_mmcfge",',  # float
    '  "first_24_water",',  # float
    '  "first_60_oil",',  # float
    '  "first_60_gas",',  # float
    '  "first_60_boe",',  # float
    '  "first_60_water",',  # float
    '  "first_60_mmcfge",',  # float
    '  "prac_ip_oil_daily",',  # float
    '  "prac_ip_gas_daily",',  # float
    '  "prac_ip_cfged",',  # float
    '  "prac_ip_boe",',  # float
    '  "latest_oil",',  # float
    '  "latest_gas",',  # float
    '  "latest_water",',  # float
    '  "prior_12_oil",',  # float
    '  "prior_12_gas",',  # float
    '  "prior_12_water",',  # float
    '  "last_test_date",',  # Date
    '  "last_flow_pressure",',  # float
    '  "last_whsip",',  # float
    '  "2nd_month_gor",',  # float
    '  "latest_gor",',  # float
    '  "cum_gor",',  # float
    '  "last_12_yield",',  # float
    '  "2nd_month_yield",',  # float
    '  "latest_yield",',  # float
    '  "peak_gas",',  # float
    '  "peak_gas_month_no.",',  # int32
    '  "peak_oil",',  # float
    '  "peak_oil_month_no.",',  # int32
    '  "peak_boe",',  # float
    '  "peak_boe_month_no.",',  # int32
    '  "peak_mmcfge",',  # float
    '  "peak_mmcfge_month_no.",',  # int32
    '  "upper_perforation",',  # float
    '  "lower_perforation",',  # float
    '  "gas_gravity",',  # float
    '  "oil_gravity",',  # float
    '  "completion_date",',  # Date
    '  "well_count",',  # int32
    '  "max_active_wells",',  # int32
    '  "months_produced",',  # int32
    '  "gas_gatherer",',  # string
    '  "oil_gatherer",',  # string
    '  "lease_number",',  # string
    '  "spud_date",',  # Date
    '  "measured_depth_td",',  # float
    '  "true_vertical_depth",',  # float
    '  "gross_perforated_interval",',  # float
    '  "field",',  # string
    '  "state",',  # string
    '  "district",',  # string
    '  "aapg_geologic_province",',  # string
    '  "country",',  # string
    '  "section",',  # string
    '  "township",',  # string
    '  "range",',  # string
    '  "abstract",',  # string
    '  "block",',  # string
    '  "survey",',  # string
    '  "ocs_area",',  # string
    '  "pgc_area",',  # string
    '  "surface_latitude_wgs84",',  # float
    '  "surface_longitude_wgs84",',  # float
    '  "last_12_oil",',  # float
    '  "last_12_gas",',  # float
    '  "last_12_water",',  # float
    '  "entity_id"',  # int32
    sep="\n"
  )
}

load_headers_by_year <- function(years, col_select) {
  paths <- file.path(
    DATA_DIR,
    "generated/production/well_headers",
    paste0("first_prod_year=", years),
    "file.parquet"
  )
  missing_paths <- !fs::file_exists(paths)
  if (any(missing_paths)) {
    message("Skipping headers for years ",
      paste(years[missing_paths], collapse=", "),
      " because the files don't exist.")
    paths <- paths[!missing_paths]
  }
  # col_select doesn't get evaluated here, so it's not counted as missing until
  # load_headers_by_filename
  load_headers_by_filename(paths, col_select)
}

load_headers_by_filename <- function(paths, col_select) {
  if (missing(col_select)) {
    err_msg <- help_message_header_cols_available()
    orig_warn_length <- getOption("warning.length")
    on.exit(options(warning.length = orig_warn_length), add=TRUE)
    options(warning.length = nchar(err_msg) + 30)
    stop(err_msg)
  }
  missing_paths <- !fs::file_exists(paths)
  df <- purrr::map_dfr(paths, arrow::read_parquet, col_select=col_select)
  if ("production_type" %in% colnames(df)) {
    df %<>% harmonize_production_type()
  }
  df
}

elapsed_months <- function(start_date, end_date) {
  stopifnot(length(start_date) == length(end_date) || length(start_date) == 1L || length(end_date) == 1L)
  sd <- as.POSIXlt(start_date)
  ed <- as.POSIXlt(end_date)
  (12L * (ed$year - sd$year)) + (ed$mon - sd$mon)
}

# Define an alternative version of which.max that doesn't return zero-length
# results. Return NA instead.
which_max <- function(x) {
  res <- which.max(x)
  if (length(res) == 0) {
    res <- NA_integer_
  }
  res
}

max_no_warn <- function(x) {
  # max will return -Inf if all arguments are NA, and provide a warning.
  # Instead, I want NA, and no warning.
  suppressWarnings(ans <- max(x, na.rm=TRUE))
  if (is.infinite(ans)) {
    ans <- NA
    # return the right class:
    class(ans) <- class(x)
  }
  ans
}

harmonize_production_type <- function(df) {
  decoder <- tibble::tribble(
    ~raw_value,               ~harmonized_value,
    "(N/A)",                  "OTHER",
    "AIR INJ",                "OTHER",
    "AIR INJECTION",          "OTHER",
    "BIW",                    "OTHER", # ?? only 23 of these
    "BRN",                    "OTHER", # ?? only 1 of these
    "BSW",                    "OTHER", # basic sediment and water
    "CANCEL",                 "OTHER", # only 5
    "CBM",                    "Gas", # coalbed methane
    "CO2",                    "OTHER",
    "CONV",                   "OTHER", # conversion?
    "CONV/PIPE",              "OTHER",
    "CORE",                   "OTHER", # ?? only 3 wells, so not worrying about it
    "DH",                     "OTHER", # ?? only wells, so not worrying about it
    "DISP",                   "INJECTION/DISPOSAL",
    "DISPOSAL",               "INJECTION/DISPOSAL",
    "DRY",                    "OTHER", # ?? only 50 wells, so not worrying about it
    "DS",                     "OTHER", # ?? only 4 wells, so not worrying about it
    "DSP",                    "INJECTION/DISPOSAL", # disposal? (5311 wells)
    "DW",                     "OTHER", # ?? only 13 wells, so not worrying about it
    "ENHANCED RECOVERY INJ",  "INJECTION/DISPOSAL",
    "EOR",                    "INJECTION/DISPOSAL",
    "EXP",                    "OTHER", # ?? only 3 wells, so not worrying about it
    "FS",                     "OTHER", # ?? only 5 wells, so not worrying about it
    "GAS",                    "Gas",
    "GAS CONDENSATE",         "Gas",
    "GAS DISPOSAL",           "INJECTION/DISPOSAL",
    "GAS INJ",                "INJECTION/DISPOSAL",
    "GAS INJECTION",          "INJECTION/DISPOSAL",
    "GAS STORE",              "Gas",
    "GAS/PIPELINE",           "Gas",
    "GD",                     "Gas", # gas development
    "GE",                     "Gas", # gas extraction
    "GS",                     "Gas", # gas? (only 18 wells)
    "HE",                     "OTHER", # ?? < 10 wells, so not worrying about it
    "HELIUM",                 "OTHER",
    "HOUSE GAS",              "Gas",
    "INCOMPLETE",             "OTHER",
    "INJ",                    "INJECTION/DISPOSAL",
    "INJ EOR",                "INJECTION/DISPOSAL",
    "IW",                     "INJECTION/DISPOSAL", # injection well/
    "JUNKED",                 "OTHER",
    "LPG INJECTION",          "INJECTION/DISPOSAL",
    "METH",                   "Gas",
    "MON",                    "OTHER",
    "MONITOR/OBSERVATION",    "OTHER",
    "MS",                     "OTHER", # monitoring storage?
    "NGL",                    "Gas", # Natural gas liquids (7 wells)
    "NITROGEN",               "OTHER",
    "NL",                     "OTHER", # not listed? (N=80)
    "NR",                     "OTHER", # ?? < 10 wells, so not worrying about it
    "O&G",                    "Oil & Gas",
    "OBS",                    "OTHER",
    "OBSERVATION",            "OTHER",
    "OD",                     "Oil", # oil development
    "OE",                     "OTHER", # ?? < 10 wells, so not worrying about it
    "OIL",                    "Oil",
    "OIL (CYCLIC STEAM)",     "Oil",
    "ORPH",                   "OTHER",
    "OTHER",                  "OTHER",
    "PA",                     "OTHER", # plugged and abandoned? (N=66)
    "PRESSURE MAINT.",        "INJECTION/DISPOSAL",
    "SECONDARY RECOVERY INJ", "INJECTION/DISPOSAL",
    "SERVICE",                "OTHER", # ?? < 10 wells, so not worrying about it
    "SKI",                    "OTHER", # ?? < 10 wells, so not worrying about it
    "SOURCE WELL",            "OTHER", # ? 51 wells, so not worrying about it.
    "ST",                     "OTHER", # storage? (N=11)
    "STEAM INJ",              "INJECTION/DISPOSAL",
    "STEAMFLOOD",             "INJECTION/DISPOSAL",
    "STORAGE",                "OTHER",
    "STR",                    "OTHER",
    "STRAT TEST",             "OTHER",
    "SULPHUR",                "OTHER", # ?? < 10 wells, so not worrying about it
    "SVC",                    "OTHER", # ?? < 10 wells, so not worrying about it
    "SW",                     "OTHER", # ?? < 10 wells, so not worrying about it
    "SWD",                    "INJECTION/DISPOSAL", # saltwater disposal
    "SWI",                    "INJECTION/DISPOSAL", # saltwater injection?
    "TA",                     "OTHER", # temporarily abandoned (N=2141)
    "TERTIARY RECOVERY INJ",  "INJECTION/DISPOSAL",
    "TEST WELL",              "OTHER",
    "TH",                     "OTHER", # ?? < 10 wells, so not worrying about it
    "UN",                     "OTHER", # ?? < 10 wells, so not worrying about it
    "UNDERGROUND",            "OTHER", # ?? < 10 wells, so not worrying about it
    "WATER",                  "OTHER",
    "WATER DISPOSAL",         "INJECTION/DISPOSAL",
    "WATER INTAKE",           "OTHER",
    "WATERFLOOD",             "INJECTION/DISPOSAL",
    "WS",                     "OTHER", # water source?
    "WSW",                    "OTHER", # water source well
    "WTR",                    "OTHER",
    "WTR INJ",                "INJECTION/DISPOSAL",
    "WTR SUPPLY",             "OTHER",
    "WTR WELL",               "OTHER",
  )
  harmonize_one_variable(df, decoder, "production_type")
}

load_natural_gas_prices <- function(filename) {
  # NOTE: expected file is "NG_PRI_FUT_S1_D.xls"
  spot_df <- readxl::read_excel(filename, sheet="Data 1", skip=2) %>%
    dplyr::rename(date = Date) %>%
    dplyr::mutate(date = lubridate::date(date)) %>%
    ensure_id_vars(date) %>%
    dplyr::rename(price_spot = `Henry Hub Natural Gas Spot Price (Dollars per Million Btu)`)

  futures_df <- readxl::read_excel(filename, sheet="Data 2", skip=2) %>%
    dplyr::rename(date = Date) %>%
    dplyr::mutate(date = lubridate::date(date)) %>%
    ensure_id_vars(date) %>%
    dplyr::rename(
      price_future_1mo = `Natural Gas Futures Contract 1 (Dollars per Million Btu)`,
      price_future_2mo = `Natural Gas Futures Contract 2 (Dollars per Million Btu)`,
      price_future_3mo = `Natural Gas Futures Contract 3 (Dollars per Million Btu)`,
      price_future_4mo = `Natural Gas Futures Contract 4 (Dollars per Million Btu)`,
    )
  out <- dplyr::full_join(spot_df, futures_df, by="date")
}

load_price_index <- function(filename, base_year = 2018) {
  # From https://fred.stlouisfed.org/series/CPILFENS
  # Consumer Price Index: All Items Less Food and Energy in U.S. City Average, All Urban Consumers
  # Monthly, not seasonally adjusted
  df <- data.table::fread(filename, data.table=FALSE) %>%
    dplyr::rename(date = DATE, price_index = VALUE) %>%
    dplyr::mutate(date = as.Date(date)) # remove data.table's IDate

  # By default, the index is 100 in 1982, which isn't helpful for me.
  # Make it into 2018 dollars.
  base_index <- dplyr::filter(df, lubridate::year(.data$date) == !!base_year) %>%
    purrr::chuck("price_index") %>%
    mean()
  return(dplyr::mutate(df, price_index = price_index / !!base_index))
}
#
# make_prices_real <- function(df) {
#   dplyr::mutate(df,
#     price_spot = price_spot / price_index,
#     price_future_1mo = price_future_1mo / price_index,
#     price_future_2mo = price_future_2mo / price_index,
#     price_future_3mo = price_future_3mo / price_index,
#     price_future_4mo = price_future_4mo / price_index,
#   )
# }

get_complete_wells <- function(header_files) {
  # NOTE: "complete" here means I observe the complete history, NOT that the
  # well has been completed. I should probably change this.
  # Other variables that could be useful:
  # "first_prod_date", "cum_gas", "cum_oil", "di_basin"
  header_cols <- c("entity_id", "first_prod_date", "last_prod_date",
    "completion_date",
    "prac_ip_gas_daily", "county", "production_type", "state")
  headers <- load_headers_by_filename(header_files, col_select=header_cols)
  stopifnot(!anyNA(headers$last_prod_date))
  out <- dplyr::filter(headers, last_prod_date < max(last_prod_date)) %>%
    dplyr::select(-last_prod_date)
  out
}

get_matching_prod <- function(headers, prod_files, prod_cols) {
  # NOTE: the code here doesn't check that prod_files line up with the timeframe
  # implied by the headers data. That's up to you.
  stopifnot("entity_id" %in% colnames(headers))
  prod_cols <- match.arg(prod_cols,
    c("oil", "gas", "daily_avg_oil", "daily_avg_gas"),
    several.ok=TRUE)
  read_cols <- c("entity_id", "month", prod_cols)
  filter_fn <- function(df, ...) {
    dplyr::semi_join(df, ...)
  } %>% purrr::partial(y=headers, by="entity_id")
  out <- load_production_by_filename(prod_files, col_select=read_cols) %>%
    dplyr::semi_join(headers, by="entity_id")
  out
}

check_required_cols <- function(df, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse=", "))
  }
  NULL
}

prep_prices <- function(prices, price_index_file) {
  price_index <- load_price_index(price_index_file)
  out <- dplyr::inner_join(prices, price_index, by="date") %>%
    dplyr::transmute(
      price_spot = price_spot / price_index,
      price_future_1mo = price_future_1mo / price_index,
      price_future_2mo = price_future_2mo / price_index,
      price_future_3mo = price_future_3mo / price_index,
      price_future_4mo = price_future_4mo / price_index,
      year  = lubridate::year(date),
      month = lubridate::month(date)
    ) %>%
    dplyr::group_by(year, month) %>%
    dplyr::summarize_all(mean, na.rm=TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(date = lubridate::make_date(year, month, 1L))
  out
}

# Define convenient versions mean_, max_, min_, and sum_ with na.rm pre-filled
# They also avoid warnings if all NA, and return NA instead of NaN
mean_ <- function(x) {
  suppressWarnings(res <- mean(x, na.rm=TRUE))
  if (is.nan(res)) {
    res <- NA_real_
  }
  res
}

sum_ <- function(x) {
  x <- x[not_na(x)]
  # Special case because sum(c()) == 0, but I want all-NA to give NA
  if (length(x) == 0) {
    return(NA_real_)
  }
  sum(x)
}

max_ <- function(x) {
  suppressWarnings(res <- max(x, na.rm=TRUE))
  if (any(is.infinite(res))) {
    res <- NA_real_
  }
  res
}

min_ <- function(x) {
  suppressWarnings(res <- min(x, na.rm=TRUE))
  if (any(is.infinite(res))) {
    res <- NA_real_
  }
  res
}

sd_ <- function(x) {
  sd(x, na.rm=TRUE)
}

not_na <- Negate(is.na)

noNAs <- function(x) {
  # Return FALSE if x has length == 0 or anyNA(x) is TRUE.
  # This could create issues if you call noNAs on a zero-length vector, but is
  # probably less surprising than the alternatives.
  # (e.g. anyNA(NULL) == TRUE)
  (length(x) > 0) && (!anyNA(x))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

demean_ <- function(x) {
  x - mean_(x)
}

quantile_ <- function(x, probs) {
  stats::quantile(x, probs, na.rm=TRUE, names=FALSE, type=8)
}

geometry_to_lonlat <- function(x) {
  if (any(sf::st_geometry_type(x) != "POINT")) {
    stop("Selecting non-points isn't implemented.")
  }
  coord_df <- sf::st_transform(x, sf::st_crs("+proj=longlat +datum=WGS84")) %>%
    sf::st_coordinates() %>%
    dplyr::as_tibble() %>%
    dplyr::select(X, Y) %>%
    dplyr::rename(longitude = X, latitude = Y)
  out <- sf::st_set_geometry(x, NULL) %>%
    dplyr::bind_cols(coord_df)
  return(out)
}

#' Check that within groups defined by `x`, there's at most one value of `y`
#'
#' @param .data A data.frame
#' @param x Character vector of variables that define the groups
#' @param y Character vector of variables to check for uniqueness with in `x`
#' @return `.data`, invisibly
#'
#' Will raise an error if there's more than one value of `y` within a group
#' defined by a combination of `x`.
#' Does not depend on grouped_df status (data are ungrouped first).
#' NAs are preserved and count as a one value.
check_groups_unique_by <- function(.data, x, y) {
  x <- rlang::syms(x)
  y <- rlang::syms(y)
  if (length(y) > 1) {
    stop("Multiple y variables isn't implemented")
  }
  df <- dplyr::ungroup(.data) %>%
    dplyr::select(!!!x, !!!y) %>%
    dplyr::distinct() %>%
    dplyr::add_count(!!!x, name="n") %>%
    dplyr::filter(n > 1)
  if (nrow(df) > 0) {
    print(head(df, 10))
    stop("Data are *not* unique for ", nrow(df), " groups! See examples above.")
  }
  invisible(.data)
}


write_tex_table <- function(x, filename) {
  stopifnot(inherits(x, "gt_tbl"))
  force(filename)
  gt::as_latex(x) %>% as.character() %>% writeLines(con=filename)
  x
}

#' Provide a function wrapper that loads data if available, otherwise generates and saves it.
#'
#' Note that this is a wrapper (aka adverb or decorator) -- it returns a function.
#' It's a lot like memoise::memoize, but saves a more useful file, and doesn't
#' check that arguments haven't changed.
#' Inspired by econtools in python.
load_or_build <- function(f, filename, envir = environment(f)) {
  stop("This function isn't used yet. Need to think about when we want to provide file name")
  ext <- tolower(tools::file_ext(filename))
  if (ext == "rds") {
    save_fn <- saveRDS
    read_fn <- readRDS
  } else if (ext == "parquet") {
    save_fn <- arrow::write_parquet
    read_fn <- arrow::read_parquet
  } else if (ext == "csv") {
    save_fn <- data.table::fwrite
    read_fn <- function(x) data.table::fread(x, data.table=FALSE)
  } else {
    stop("Sorry, extentions of type '", ext, "' aren't supported yet")
  }
  # The following borrows heavily from memoise
  memo_f <- function(...) {
    mc <- match.call()
    encl <- parent.env(environment())

    if (file.exists(encl$`_filename`)) {
      res <- encl$`_reader`(encl$`_filename`)
    } else {
      mc[[1L]] <- encl$`_f`
      res <- eval(mc, parent.frame())
      encl$`_writer`(res, encl$`_filename`)
    }
    res
  }
  formals(memo_f) <- formals(args(f))
  if (is.null(envir)) {
    envir <- baseenv()
  }
  memo_f_env <- new.env(parent = envir)
  memo_f_env$`_reader` <- read_fn
  memo_f_env$`_writer` <- save_fn
  memo_f_env$`_f` <- f
  memo_f_env$`_filename` <- filename

  environment(memo_f) <- memo_f_env
  # class(memo_f) <- c("memoised", "function")
  memo_f
}

inverse_logit <- function(x) {
  # aka logistic
  # Transform from log odds to probabilities
  # R -> (0, 1)
  1 / (1 + exp(-x))
}
logit <- function(x) {
  # log(x / (1 - x)), but probably more numerically stable
  log(x) - log1p(-x)
}

fill_na <- function(x, val) {
  x[is.na(x)] <- val
  x
}

calc_oil_gas_frac <- function(oil, gas, verbose=FALSE) {
  # Define a continuous measure of how much of this well's production
  # comes from oil. Using DI's simple equivalence 6 mcf = 1 bbl.
  res <- (oil * 6) / (gas + (oil * 6))
  na_frac <- mean(is.na(res))
  if (verbose && na_frac > 0) {
    message("Replacing ", signif(100 * na_frac, 3), "% of oil/gas frac observations with 0")
  }
  ifelse(is.na(res), 0, res)
}

tidy_gsub <- function (x, pattern, replacement, fixed = FALSE) {
  perl <- !fixed
  gsub(pattern, replacement, x, fixed=fixed, perl=perl)
}

make_table_fragment <- function(data, escaped=TRUE, add_comments="") {
  # Identical to gt::to_latex(), except omitting the table start, column
  # headings, and table end. (More flexibility to add those yourself)
  # Built using *internal* functions of gt v0.2.1
  if (digest::sha1(gt::as_latex) != "efc63eab16916eae6fb932b077e839707792ecb7") {
    stop("as_latex function has changed -- the function make_table_fragment may need to be edited")
    # This isn't perfect, of course, because the component functions could also
    # change.
  }
  gt:::stop_if_not_gt(data = data)

  data <- data %>% gt:::build_data(context = "latex")
  table_start <- gt:::create_table_start_l(data = data)
  heading_component <- gt:::create_heading_component(data = data, context = "latex")
  columns_component <- gt:::create_columns_component_l(data = data)
  body_component <- gt:::create_body_component_l(data = data)
  source_notes_component <- gt:::create_source_note_component_l(data = data)
  footnotes_component <- gt:::create_footnotes_component_l(data = data)
  table_end <- gt:::create_table_end_l()

  # Commented out the table start, columns, and end.
  # Another way to do this would be to gsub out the text we don't want.
  text <- paste0(
    # table_start,
    heading_component,
    # columns_component,
    body_component,
    # table_end,
    footnotes_component,
    source_notes_component,
    collapse = ""
  )
  if (any(add_comments != "")) {
    add_comments <- paste("%", paste(add_comments, collapse=", "))
  }

  # gt escapes latex already. If we want non-escaped, we need to reverse.
  # Retain the replacements "~" -> \textasciitilde and "\^" -> \textasciicircum
  # Replace back: \textbackslash -> \, \& -> &, \% -> %, \$ -> $, \# -> #,
  # \_ -> _, \{ -> {, \} -> }
  if (!escaped) {
    text <- text %>%
      tidy_gsub("\\\\([&%$#_{}])", "\\1") %>%
      tidy_gsub("\\textbackslash ", "\\", fixed=TRUE)
  }
  # Remove trailing space (don't need to make it more general whitespace)
  text %<>% tidy_gsub(" \n", "\n", fixed=TRUE) %>%
    paste(
      "%",
      "% DO NOT EDIT! This table is automatically generated. Any changes will be overwritten.",
      "%",
      add_comments,
      ., sep="\n"
    )
  text
}

# Use glue::glue() and rlang::inform() to emit a message. Arguments are the
# inform() default values.
glue_message <- function(
    ...,
    class = NULL,
    .file = NULL,
    .frequency = c("always", "regularly", "once"),
    .frequency_id = NULL
  ) {
  .frequency <- match.arg(.frequency)
  # Need to provde .envir here because otherwise we've changed the frame by
  # calling glue() within another function.
  # It will generate an error to also pass .envir, but that should be rare
  envir = parent.frame()
  rlang::inform(
    glue::glue(..., .envir=envir),
    class=class, .file=.file, .frequency=.frequency, .frequency_id=.frequency_id
  )
}


#' Convert a tidy dataframe of results (e.g. broom output) to a dataframe with
#' formatted estimates and confidence intervals (on different rows)
#'
#'
#' @param df Dataframe with columns named term, estimate, conf_low, and conf_high
#' @param n_sigfig Number of significant figures in output (default 4)
#' @param decimals Number of decimals figures in output (if not NULL, n_sigfig is ignored)
#' @param align String (max 1 character) to make latex alignment easier
#'
#' @examples
#' lm(mpg ~ cyl + am, mtcars) %>%
#'   broom::tidy(conf.int=TRUE) %>%
#'   format_estimate_above_interval()
#' @return A gt object with the estimates and confidence intervals, aligned
format_estimate_above_interval <- function(df, n_sigfig = 4, decimals = NULL, align="") {
  if (!all(c("term", "estimate", "conf_low", "conf_high") %in% names(df))) {
    stop(
      "Data must include variables term, estimate, conf_low, and conf_high.\n",
      "   (These names are for compatibility with broom. Use them even if the ",
      "interval isn't a confidence interval.)"
    )
  }
  N <- nrow(df)
  df %<>% dplyr::ungroup() %>%
    dplyr::select(term, estimate, conf_low, conf_high)
  stopifnot(length(align) == 1, nchar(align) %in% 0:1, N >= 1, !anyNA(df))
  # Do we want to support other fmt_number options?
  if (is.null(decimals)) {
    stopifnot(length(n_sigfig) %in% c(1, N))
    rounding_fn <- purrr::partial(base::signif, digits = n_sigfig)
  } else {
    stopifnot(length(decimals) %in% c(1, N))
    rounding_fn <- purrr::partial(base::round, digits = decimals)
  }
  # Define a little function here to make things easier.
  add_alignment <- function(x, align) {
    x %<>% as.character()
    if (align == "") {
      return(x)
    }
    # add the alignment character after the decimal, if there is one, otherwise
    # at the end
    ifelse(
      grepl(".", x, fixed=TRUE),
      tidy_gsub(x, ".", paste0(".", align), fixed=TRUE),
      paste0(x, align)
    )
  }
  # Round and set sort order
  df %<>% dplyr::mutate_if(is.numeric, rounding_fn) %>%
    dplyr::mutate(term_sort_order = dplyr::row_number())
  ests <- dplyr::mutate(df,
    term_group = term,
    # insert the alignment marker after the decimal point
    estimate = add_alignment(estimate, align),
    row_type = "1 point estimate",
  )
  intervals <- dplyr::mutate(df,
    term_group = term,
    term = "",
    # Make an interval like "[10.00, 11.11]", with the alignment character after
    # the comma.
    estimate = as.character(glue::glue("[{conf_low},{align} {conf_high}]")),
    row_type = "2 conf interval",
  )
  out <- dplyr::bind_rows(ests, intervals) %>%
    dplyr::arrange(term_sort_order, row_type) %>%
    dplyr::select(term_group, row_type, term, estimate)
  out
}

#' @param tidy_results_lst List of results, dataframes, as produced by
#' format_estimate_above_interval()
#' @examples
#' reg1 <- lm(mpg ~ cyl + am, mtcars) %>%
#'   broom::tidy(conf.int=TRUE) %>%
#'   format_estimate_above_interval()
#' reg2 <- lm(mpg ~ cyl + wt + am, mtcars) %>%
#'   broom::tidy(conf.int=TRUE) %>%
#'   format_estimate_above_interval()
#' tab <- list(reg1, reg1, reg2) %>% merge_estimates_df()
merge_estimates_df <- function(tidy_results_lst) {
  # Steps:
  # - check data
  # - figure out the full set of rows
  # - make sure each results df has those rows
  # - sort each df according to the row order in the first df
  # - drop/rename cols to remain unique
  # - merge everything together
  # - run through gt
  stopifnot(length(tidy_results_lst) >= 1)
  is_unique <- purrr::map_lgl(tidy_results_lst, is_id, "term_group", "row_type")
  if (!all(is_unique)) {
    stop(
      "Can't merge non-unique estimates dataframes. (", sum(!is_unique), " of ",
      length(is_unique), " failed)"
    )
  }
  # For each dataframe `i` in the list, rename the "estimate" column estimate_`i`
  # (For clarity later). The !!.y := estimate syntax is necessary because of the
  # way dplyr does nonstandard evaluation.
  est_rename <- paste0("estimate", seq_along(tidy_results_lst))
  tidy_results_lst %<>% purrr::map2(est_rename, ~dplyr::rename(.x, !!.y := estimate))

  full_cross <- purrr::map_dfr(
    tidy_results_lst, ~dplyr::select(., term, term_group, row_type)
    ) %>%
    dplyr::distinct() %>%
    tidyr::expand(term_group, row_type)
  expanded_dfs <- purrr::map(tidy_results_lst, safejoin::safe_full_join,
      y=full_cross, by=c("term_group", "row_type"), check="U V B C L T"
    )
  # Very messy way of filling in the term when the first results df didn't have it
  first_df <- expanded_dfs[[1]] %>%
    dplyr::mutate(term = dplyr::case_when(
      !is.na(term) ~ term,
      row_type == "1 point estimate" ~ term_group,
      TRUE ~ ""))

  # # First, remove the "term" variable from all but the first df
  # # Then join all the tables together (this works as long as tidy_results_lst
  # # contains at least one results dataframe)
  merged_df <- tail(tidy_results_lst, length(tidy_results_lst) - 1) %>%
    purrr::map(~dplyr::select(., -term)) %>%
    purrr::reduce(safejoin::safe_full_join, .init=first_df,
      by=c("term_group", "row_type"), check="U V B C L T")

  # Make a gt
  tab <- gt::gt(merged_df) %>%
    gt::fmt_missing(tidyselect::everything(), missing_text="") %>%
    gt::cols_hide(c("term_group", "row_type"))
  tab
}

clamp <- function(x, lower = -Inf, upper = Inf) {
  stopifnot(lower <= upper)
   x[x < lower] <- lower
   x[x > upper] <- upper
   return(x)
}

#' Extract a piece from each element of a list
#'
#' @param lst Iterable to extract from
#' @param ... Description of pieces to extact. See [purrr::pluck()].
#' @return A list of the extracted elements. Raises an error if the requested
#' sub-element doesn't exist.
#'
#' This behavior is very similar to purrr::map(), but (imo) easier to read and
#' (objectively) stricter.
#'
#' @examples
#' obj1 <- list("a", list(1, elt = "foo"))
#' obj2 <- list("b", list(2, elt = "bar"))
#' x <- list(obj1, obj2)
#' extract_from_each(x, 1) # first element of each
#' extract_from_each(x, 2, "elt") # "elt" element of second element of each
extract_from_each <- function(lst, ...) {
  purrr::map(lst, purrr::chuck, ...)
}

print_pipe <- function(x) {
  # I've defined this enough times for debugging, I guess I should save it.
  print(x)
  x
}

first_of_month <- function(x) {
  lubridate::make_date(lubridate::year(x), lubridate::month(x), 1L)
}

filename_to_model_name <- function(x) {
  stopifnot(length(x) >= 1)
  if (any(basename(x) != "model_fit.rds")) {
    stop(
      "Expected a file named 'model_fit.rds', not\n",
      paste(x, collapse=", ")
   )
  }
  # Convert
  # "data/generated/stan_fits/01_twopart_lognormal_model-bootstrap/model_fit.rds"
  # into list(model="01_twopart_lognormal", prior_only=FALSE, bootstrap=TRUE)
  # Each element of the output has the same length as x.
  d <- basename(dirname(x))
  nm <- stringr::str_match(d, "([0-9]{2}_[^\\-]+)(-prior)?(-bootstrap)?")[, 2, drop=TRUE]
  stopifnot(noNAs(nm))
  bad_model_names <- nm[!nm %in% MODEL_NAMES$all]
  if (length(bad_model_names) > 0) {
    stop("Unrecognized model name parsed from file: ", paste(bad_model_names, collapse=", "))
  }
  nm
}

filename_to_time_period_hr <- function(filename) {
  stopifnot(length(filename) == 1)
  hours_match <- stringr::str_match(filename, "period_(\\d+)_hours")[, 2, drop=TRUE]
  if (is.na(hours_match)) {
    # Do a quick check:
    stopifnot(!stringr::str_detect(filename, "hours"))
    time_H <- 8760 # 1 year
  } else {
    time_H <- as.numeric(hours_match)
    stopifnot(noNAs(time_H))
  }
  time_H
}

approx_equal <- function(x, y) {
  stopifnot(length(x) == length(y) || length(x) == 1 || length(y) == 1)
  (x > y - 1e-10) & (x < y + 1e-10)
}

harmonize_basin_name <- function(df, group_small_CA_basins=FALSE) {
  stopifnot(c("state", "basin") %in% colnames(df))
  out <- dplyr::mutate(df,
    # Insert a space in snake-case basins:
    basin = tidy_gsub(.data$basin, "([a-z])([A-Z])", "\\1 \\2"),
    basin = stringr::str_to_title(basin),
    basin = tidy_gsub(basin, " Basin$", ""),
    # In some places this func, the variable name "state" is the abbreviation and sometimes the full name.
    # This is dirty.
    # San Joaquin is vastly more common than the other options.
    basin = dplyr::if_else(is.na(basin) & .data$state %in% c("CA", "California"), "San Joaquin", basin),
  )

  remaining_na <- dplyr::filter(out, is.na(basin))
  if (nrow(remaining_na) > 0) {
    print(remaining_na)
    stop("There should be no NAs in basin -- see above")
  }
  # Group all other CA basins together because there aren't many sampled wells
  # (by default, don't group)
  if (isTRUE(group_small_CA_basins)) {
    out %<>% dplyr::mutate(
      basin = dplyr::if_else(
        state == "CA" & basin != "San Joaquin",
        "Other California",
        basin
      )
    )
  }
  out
}

read_constants <- function(filename = here::here("code/constants.json")) {
  # See comments in read_constants in outcomes_analysis_common.py for contents.
  const <- jsonlite::fromJSON(filename)
  # Convert these from list to named vector:
  const$TAU_LEVELS <- unlist(const$TAU_LEVELS)
  const$T_LEVELS <- unlist(const$T_LEVELS)
  stopifnot(
    rlang::is_dictionaryish(const$TAU_LEVELS),
    rlang::is_dictionaryish(const$T_LEVELS)
  )
  const
}

escape_regex <- function(x) {
  # There's probably a better way to do this...
  char_to_excape <- "-.,!*+?^$[\\](){}!=|"
  replace_regex <- paste0(
    "(",
    paste0("\\", strsplit(char_to_excape, "")[[1]], collapse="|"),
    ")"
  )
  tidy_gsub(x, replace_regex, "\\\\\\1")
}


check_for_bad_env_vars <- function() {
  # Singularity passes environment variables through to the container. That's a
  # problem for some path variables like RENV_PATHS_ROOT, which won't exist
  # inside the container.
  # in the Snakefile, so here we check if the bad env vars are set and error out
  # with a better error message.
  known_bad <- c(
    "ARROW_THIRDPARTY_DEPENDENCY_DIR",
    "LIBARROW_BINARY",
    "ARROW_HOME",
    "RENV_PATHS_ROOT",
    "R_LIBS",
    "R_LIBS_SITE"
  )
  # R_LIBS_USER is not a problem here because it gets auto-set to something like
  # ~/R/x86_64-conda-linux-gnu-library/4.0, which is fine because '~' becomes
  # the project directory, and then it isn't used.

  in_singularity <- Sys.getenv("SINGULARITY_NAME") != ""
  bad_vars_set <- Sys.getenv(known_bad) != ""
  if (in_singularity && any(bad_vars_set)) {
    err_msg <- paste0(
      "\n\n",
      "Path environment variables are set that should not be. These cause problems ",
      "when running the code in a container.\n",
      "Please re-run Snakemake with\n",
      "snakemake --use-conda --use-singularity --singularity-args='--cleanenv'\n\n",
      "Problem vars: ", paste(known_bad[bad_vars_set], collapse = ", "), "\n\n"
    )
    stop(err_msg, call. = FALSE)
  }
  invisible(NULL)
}

# We want this file to run even if no R packages are installed (because we use
# it in setup_r_library.R), so guard all of these basic setup commands:

if (requireNamespace("magrittr", quietly=TRUE)) {
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
}
if (requireNamespace("rlang", quietly=TRUE)) {
  `%||%` <- rlang::`%||%`
}
if (requireNamespace("here", quietly=TRUE)) {
  source(here::here("code/kdw_package_code.r"))
  DATA_DIR <- here::here("data")
  if (requireNamespace("jsonlite", quietly=TRUE)) {
    MODEL_NAMES <- read_constants()[["MODEL_NAMES"]]
  }
}

elapsed_time_str <- function(start_time) {
  elapsed <- Sys.time() - start_time
  paste(round(elapsed, 2), attr(elapsed, "units"))
}

options(
  stringsAsFactors = FALSE,
  warn = max(1, getOption("warn")),
  scipen = max(5, getOption("scipen"))
)



SnakemakePlaceholder <- setClass("Snakemake",
  slots = c(input = "list", output = "list", params = "list",
    wildcards = "list", threads = "numeric", log = "list", resources = "list",
    config = "list", rule = "character", bench_iteration = "numeric",
    scriptdir = "character", source = "function")
)
