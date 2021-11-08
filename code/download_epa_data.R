`%>%` <- magrittr::`%>%`
options(scipen=99)

#' Download one EPA Envirofacts table
#'
#' @param output The name of the output file to create. Should be
#'   `.../<table_name>.csv.gz` See details.
#' @param verbose Narrate progress? Default false
#' @return The name of the output file, invisibly
#'
#' See table listing here:
#' https://www.epa.gov/enviro/greenhouse-gas-subpart-w-model-all-tables
#' REST API docs: https://www.epa.gov/enviro/web-services
#' Useful tables:
#' All years:
#' - PUB_DIM_FACILITY: "Identifying a collection of data elements containing
#'   detail facility information."
#' - EF_W_EMISSIONS_SOURCE_GHG: "each facility's total reported emissions
#'   foreach industry segment and reporting source reported under Subpart W
#'   across the entire time-series"
#'
#' Reporting year 2015-2018:
#' - EF_W_FACILITY_OVERVIEW
#' - EF_W_ONSHORE_WELLS: Well level info
#' - EF_W_ONS_PROD_MISSING: Seems to be missingness indicators for the above
#'
#' Reporting year 2010-2014:
#' - W_OFFSHORE_SOURCES: "A collection of data elements reported for total
#'   emissions from Offshore Source Emissions."
#' - W_SUB_BASIN: "A collection of data elements reported for each applicable
#'   sub-basin category in the basin."
#' - W_SUBPART_LEVEL_INFORMATION: A collection of data elements containing the
#'   total annual emissions of each greenhouse gas (GHG) listed in Table A-1 of
#'   40 CFR 98 Mandatory Reporting of Greenhouse Gases reported under this
#'   Subpart, expressed in metric tons.
#'
#' Also note, from the site:
#' - All data tables for 2015 and onwards present emissions data in units of
#'   metric tons (mt).
#' - All data tables for 2011-2014 present emissions BOTH in mt (as calculated
#'   by EPA from reported data) AND in units of mt of Carbon Dioxide Equivalent
#'   (CO2e).
#' - For 2011-2012, emissions were reported in units of mt CO2e using GWP's
#'   from the IPCC's Second Assessment Report.
#' - For 2013-2014, emissions were reported in units of mt CO2e using GWP's
#'   from the IPCC's Fourth Assessment Report.
#' - Emissions presented in the table W_SUBPART_LEVEL_INFORMATION span the full
#'   time series (2011 on) and are in units of mt for all years.
download_epa_table <- function(output, verbose=FALSE) {
  output <- output[[1]]
  if (file.exists(output)) {
    if (verbose) {
      message("Output file ", output, " already exists")
    }
    return(output)
  }
  # We'll need these later, so check now.
  # Minor improvement: we're not using processx later, so could change these to
  # system2 calls.
  processx::run("which", "tail")
  processx::run("which", "gzip")
  processx::run("which", "cat")
  table_name <- gsub(".csv.gz", "", basename(output), fixed=TRUE) %>%
    toupper()
  if (verbose) {
    message("Will download ", table_name)
  }
  url_base <- paste0("https://data.epa.gov/efservice/", table_name)

  row_count_response <- xml2::read_xml(paste0(url_base, "/COUNT"))
  row_count <- xml2::as_list(row_count_response)$Envirofacts$RequestRecordCount[[1]] %>%
    as.numeric()

  if (is.null(row_count)) {
    print(row_count_response)
    stop("Failed to parse row count from the xml response above")
  }
  if (verbose) {
    message(table_name, " has ", row_count, " rows. Requesting the full data now")
  }
  # Break into multiple requests for a better experience. Nothing special about
  # 1e5 here, but it's the default number if you don't specify a row count.
  row_starts <- seq(0, row_count, by=1e5)
  row_ends <- c(row_starts[-1], row_count) - 1

  request_urls <- paste0(url_base, "/rows/", row_starts, ":", row_ends, "/CSV")

  if (any(grepl("\\de\\+?\\d", request_urls, perl=TRUE))) {
    print(request_urls)
    stop("Something looks like scientific notation in the request urls")
  }

  download_tempfile <- function(url, verbose) {
    if (verbose) {
      message("Downloading ", url)
    }
    tf <- tempfile()
    curl::curl_download(url, destfile=tf)
    return(tf)
  }
  # Could run this in parallel if necessary:
  tempfiles <- purrr::map_chr(request_urls, download_tempfile, verbose=verbose)

  # Now, reassemble the files, removing the header of all but the first
  for (idx in seq_along(tempfiles)) {
    if (idx > 1) {
      old_tempfile <- tempfiles[idx]
      new_tempfile <- tempfile()
      status_code <- system2("tail", c("--lines=+2", old_tempfile), stdout=new_tempfile)
      if (status_code != 0) {
        stop(
          "Failed to take the tail for file ", idx,
          " which should have been downloaded with this query:\n  ",
          request_urls[idx]
        )
      }
      tempfiles[idx] <- new_tempfile
    }
  }
  if (verbose) {
    message("Combining and zipping files")
  }
  status_code <- system2("cat", c(tempfiles, "| gzip"), stdout=output)
  if (status_code != 0) {
    print(tempfiles)
    stop("Failed to concatenate and zip the files above")
  }
  if (verbose) {
    message(
      "Downloaded data successfully. Saved as:\n  ",
      normalizePath(output, winslash="/", mustWork=FALSE)
    )
  }
  invisible(output)
}

if (!exists("snakemake")) {
  stop("This script is meant to be run by snakemake")
}

download_epa_table(snakemake@output, verbose=TRUE)
