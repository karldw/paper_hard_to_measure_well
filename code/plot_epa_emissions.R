library(ggplot2)
options(warn=2)

suppressMessages(
  here::i_am("code/plot_epa_emissions.R", uuid="5ec8cce3-09c8-4620-b5f9-ce7a70280d02")
)
source(here::here("code/shared_functions.r"))

EPA_CH4_GWP <- 25 # EPA's GWP number for methane
OUR_CH4_GWP <- 29.8 # More comprehensive IPCC number

unzip_files <- function(zf) {
  stopifnot(file.exists(zf))
  td <- tempdir()
  # File for GWP check:
  gwp <- unzip(zf, files="Chapter Text/Executive Summary/Table ES-1.csv", exdir=td)
  # File with emissions:
  emiss <- unzip(zf, files="Chapter Text/Executive Summary/Table ES-2.csv", exdir=td)
  c(gwp = gwp, emiss = emiss)
}

check_epa_gwp <- function(filename) {
  epa_gwp <- data.table::fread(filename, skip=2, header=TRUE) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(Gas == "CH4 a") %>%
    dplyr::transmute(gas = Gas, gwp = as.numeric(tidy_gsub(GWP, ",", "")))
  stopifnot(nrow(epa_gwp) == 1, epa_gwp$gwp == EPA_CH4_GWP)
}

load_emiss <- function(filename) {
  epa_emiss <- data.table::fread(filename, skip=2, header=TRUE)  %>%
    dplyr::as_tibble() %>%
    dplyr::select(-V1) %>%
    dplyr::rename(source = "Gas/Source")

  ch4_row_index <- 36:56
  stopifnot(
    epa_emiss$source[min(ch4_row_index)] == "CH4 c",
    epa_emiss$source[max(ch4_row_index) + 1] == "N2O c"
  )
  epa_ch4 <- dplyr::slice(epa_emiss, ch4_row_index) %>%
    tidyr::pivot_longer(-source, names_to="year", values_to="mmt_co2e") %>%
    dplyr::group_by(source) %>%
    dplyr::filter(any(mmt_co2e != "+")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      source = dplyr::if_else(source == "CH4 c", "total", source),
      year = as.integer(year),
      mmt_co2e = suppressWarnings(as.numeric(mmt_co2e)),
    )


  epa_ch4_total <- dplyr::filter(epa_ch4, source=="total") %>%
    dplyr::select(-source)
  epa_emiss_total <- dplyr::filter(epa_emiss, source == "Total Emissions") %>%
    tidyr::pivot_longer(-source, names_to="year", values_to="mmt_co2e") %>%
    dplyr::mutate(
      year = as.integer(year),
      mmt_co2e = as.numeric(mmt_co2e),
    ) %>%
    dplyr::left_join(dplyr::rename(epa_ch4_total, ch4_mmt_co2e = mmt_co2e), by="year") %>%
    dplyr::mutate(
      mmt_co2e_net_of_ch4 = mmt_co2e - ch4_mmt_co2e
    )

  epa_ch4 %<>% dplyr::mutate(
    mmt_co2e = mmt_co2e * OUR_CH4_GWP / EPA_CH4_GWP,
    source_coarse = dplyr::case_when(
      source == "total" ~ "Total CH4",
      source %in% c("Natural Gas Systems", "Petroleum Systems", "Abandoned Oil and Gas Wells") ~ "Oil and gas",
      source %in% c("Coal Mining", "Abandoned Underground Coal Mines",
        "Petrochemical Production",
        "Stationary Combustion",
        "Mobile Combustion",
        "International Bunker Fuels b") ~ "Other fossil",
      source %in% c("Landfills", "Wastewater Treatment", "Composting") ~ "Waste",
      source %in% c("Enteric Fermentation", "Manure Management",
      "Rice Cultivation", "Field Burning of Agricultural Residues") ~ "Agriculture",
    )
  )
  # Check that we assigned a `source_coarse` to every `source`
  stopifnot(!anyNA(epa_ch4$source_coarse))

  list(epa_ch4 = epa_ch4, epa_emiss_total = epa_emiss_total)
}

make_plots <- function(df_list) {
  col_ordering <- c("Non-CH4 total", "Total CH4",
    "Agriculture", "Waste", "Other fossil", "Oil and gas"
  )

  ch4_2015 <- df_list$epa_ch4 %>%
    ensure_id_vars(source, year) %>%
    dplyr::group_by(source_coarse, year) %>%
    dplyr::summarize(mmt_co2e = sum_(mmt_co2e), .groups="drop") %>%
    dplyr::filter(year == 2015) %>%
    dplyr::select(-year)


  non_ch4_2015 <- df_list$epa_emiss_total %>%
    ensure_id_vars(year) %>%
    dplyr::filter(year == 2015) %>%
    dplyr::transmute(source_coarse = "Non-CH4 total", mmt_co2e = mmt_co2e_net_of_ch4)

  # Using updated values from Alvarez et al. 2018 table s1 total
  alvarez_og_total_ch4 <- 13
  epa_og_total_ch4 <- 8.1
  alvarez_og_prod_co2e <- 7.6 * OUR_CH4_GWP
  additional_emiss_co2e <- (alvarez_og_total_ch4 - epa_og_total_ch4) * OUR_CH4_GWP
  ch4_2015_corrected <- dplyr::mutate(ch4_2015,
      mmt_co2e = dplyr::if_else(
        source_coarse %in% c("Oil and gas", "Total CH4"),
        # Add another source_coarse with updated values from
        # Alvarez et al. 2018 table s1 total
        mmt_co2e + additional_emiss_co2e,
        mmt_co2e
      ),
    x_axis = "Corrected O&G\n(Alvarez et al. 2018)",
  )
  epa_ch4_oil_gas_number <- dplyr::filter(ch4_2015, source_coarse == "Oil and gas")$mmt_co2e
  corrected_increase_pct <- signif(100 * additional_emiss_co2e / epa_ch4_oil_gas_number, 3)
  message("Correction adds ", corrected_increase_pct, "% to EPA estimate")
  corrected_og_total_emiss_pct <- 100 * signif(
    dplyr::filter(ch4_2015_corrected, source_coarse == "Oil and gas")$mmt_co2e /
    dplyr::filter(ch4_2015_corrected, source_coarse == "Total CH4")$mmt_co2e,
    digits=3
  )
  corrected_og_prod_emiss_pct <- 100 * signif(
    alvarez_og_prod_co2e /
    dplyr::filter(ch4_2015_corrected, source_coarse == "Total CH4")$mmt_co2e,
    digits=3
  )
  message("Corrected total: O&G is ", corrected_og_total_emiss_pct, "% of total methane emissions")
  message("Corrected segment: O&G upstream prod is ", corrected_og_prod_emiss_pct, "% of total methane emissions")

  plt_ch4_2015 <- dplyr::mutate(ch4_2015, x_axis = "EPA GHGI for 2015") %>%
    dplyr::bind_rows(ch4_2015_corrected) %>%
    dplyr::filter(source_coarse != "Total CH4") %>%
    dplyr::mutate(
      x_axis = factor(x_axis, levels=c("EPA GHGI for 2015", "Corrected O&G\n(Alvarez et al. 2018)")),
      source_coarse = factor(source_coarse, levels=rev(col_ordering)),
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x=x_axis, y=mmt_co2e, fill=source_coarse)) +
    ggplot2::geom_col(alpha=0.8, color = "black") + # border color
    ggplot2::geom_text(ggplot2::aes(label=source_coarse), position="stack", vjust = 1.6) +
    ggplot2::scale_fill_brewer(palette="Dark2", guide="none") +
    ggplot2::theme_bw() +
    ggplot2::labs(x="", y="MMT CO2e", fill="") +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )

  plt_ghg_2015 <- dplyr::bind_rows(ch4_2015, non_ch4_2015) %>%
    dplyr::mutate(source_coarse = factor(source_coarse, levels=col_ordering)) %>%
    ggplot2::ggplot(ggplot2::aes(x=source_coarse, y=mmt_co2e)) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::labs(x="", y="MMT CO2e")

  tf1 <- tempfile(fileext=".pdf")
  ggplot2::ggsave(
    filename=tf1,
    plot=plt_ch4_2015,
    device = grDevices::cairo_pdf,
    width = 3.15, height = 3.54, units="in"
  )
  reset_datestamp(
    tf1,
    here::here("graphics/epa_emiss_2015_ch4.pdf"),
    category="pdf"
  )
  save_plot(plt_ghg_2015, here::here("graphics/epa_emiss_2015_ghg.pdf"), reproducible=TRUE)
}


unzipped <- unzip_files(here::here("data/epa/EPA_GHG_report_2020_table_data.zip"))
check_epa_gwp(unzipped[["gwp"]])
df_list <- load_emiss(unzipped[["emiss"]])
make_plots(df_list)
