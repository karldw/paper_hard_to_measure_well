
prep_custom_model_data <- function(df, model_name, shift_amount=5) {
  stopifnot(
    inherits(df, "data.frame"),
    shift_amount >= 1,
    shift_amount < 100,
    length(model_name) == 1
  )

  if (model_name %in% MODEL_NAMES$meas_err_models) {
    # See notes above about def of emiss_se_kg_hr_fill_mean
    fm <- brms::brmsformula(
      emiss_kg_hr_fill0 | mi(emiss_se_kg_hr_fill_mean) ~
        asinh(gas_avg_mcfd) + asinh(oil_avg_bbld) + basin +
        prod_oil_frac + asinh(age_yr) + drill_type
    )
  } else {
    # NB: we may take logs of LHS in the stan code (depending on distribution)
    # These asinh transformations are applied for real though.
    fm <- brms::brmsformula(
      emiss_kg_hr_fill0 ~
        asinh(gas_avg_mcfd) + asinh(oil_avg_bbld) + basin +
        prod_oil_frac + asinh(age_yr) + drill_type
    )
  }
  sdata <- brms::make_standata(fm, df)
  # Check that we didn't drop any missings:
  # (in emiss_kg_hr_fill0, non-detections are coded as zero)
  stopifnot(length(sdata$Y) == nrow(df), all(is.finite(sdata$Y)))
  # Define here what counts as a large leak by this shift.
  if (exists("LEAK_SIZE_DEF") && LEAK_SIZE_DEF != shift_amount) {
    stop(
      "LEAK_SIZE_DEF global doesn't match local shift_amount ",
      "(", LEAK_SIZE_DEF, " vs ", shift_amount, ")"
    )
  }
  sdata$shift_amount <- shift_amount

  if (!model_name %in% MODEL_NAMES$meas_err_models) {
    # For models that don't use measurement error, still need to pass "noise" to
    # fit the shared stan data block happy. These values aren't used.
    sdata$noise <- rep_len(0, sdata$N)
  }
  sdata$price <- df$gas_price_per_kg_ch4 %||% stop("price variable is missing")
  stopifnot(noNAs(sdata)) # may be Inf though
  sdata
}


prep_measurement_data <- function(filename) {
  df_list <- readRDS(filename)

  censor_threshold_ground_kg <- 0.04 # from the literature, e.g. Zhang et al (2020)
  censor_threshold_aviris_kg <- 5 # 5-10 in the literature
  methane_mcf_per_kg <- 1 / 18.8916 # Note: this is for pure CH4, not nat gas

  df_list$ground_studies %<>% dplyr::transmute(
    censor_threshold_kg = !!censor_threshold_ground_kg,
    censored = ifelse(emiss_kg_hr <= censor_threshold_kg, "left", "none"),
    emiss_kg_hr_fill_censor = ifelse(censored == "left", censor_threshold_kg, emiss_kg_hr),
    gas_avg_mcfd = gas_production_mcfd,
    basin = basin
  )
  # NB: I'm filling emiss_kg_hr with the censor value, but that's the value brms
  # expects in the modeling.
  # I'm not assuming that these wells actually emit that quantity.
  #
  df_list$jpl_wells_all %<>% dplyr::select(
      #n_overflights, persistence_frac, county,
      emiss_kg_hr, emiss_se_kg_hr,
      basin, drill_type, oil_avg_bbld, gas_avg_mcfd,
      first_60_oil, first_60_gas,
      age_yr, gas_price_per_mcf, well_pad_id, drill_type, gas_frac_methane,
    ) %>%
    dplyr::mutate(
      # n_overflights = fill_na(n_overflights, 0L),
      # persistence_frac = fill_na(persistence_frac, 0),
      censor_threshold_kg = !!censor_threshold_aviris_kg,
      emiss_kg_hr_fill0 = ifelse(is.na(emiss_kg_hr), 0, emiss_kg_hr),
      emiss_kg_hr_fill_censor = ifelse(is.na(emiss_kg_hr), censor_threshold_kg, emiss_kg_hr),
      censored = ifelse(emiss_kg_hr_fill0 <= censor_threshold_kg, "left", "none"),
      # Define a continuous measure of how much of this well's production
      # comes from oil. Using DI's simple equivalence 6 mcf = 1 bbl.
      prod_oil_frac = calc_oil_gas_frac(oil_avg_bbld, gas_avg_mcfd),
    )
  df_list$four_corners_all_wells %<>% dplyr::select(-production_type) %>%
    dplyr::mutate(
      censor_threshold_kg = !!censor_threshold_aviris_kg,
      emiss_kg_hr_fill0 = ifelse(is.na(emiss_kg_hr), 0, emiss_kg_hr),
      emiss_kg_hr_fill_censor = ifelse(is.na(emiss_kg_hr), censor_threshold_kg, emiss_kg_hr),
      censored = ifelse(emiss_kg_hr_fill0 <= censor_threshold_kg, "left", "none"),
      prod_oil_frac = calc_oil_gas_frac(oil_avg_bbld, gas_avg_mcfd),
    )
  df_list$aviris_all <- dplyr::bind_rows(
      df_list$jpl_wells_all,
      df_list$four_corners_all_wells
     ) %>%
     dplyr::mutate(
      # emiss_se_kg_hr isn't reported for the Four Corners data, so fill it in.
      # Assume measurement error is proportional to leak. Could do fancier
      # modeling here, but that seemms a little excessive.
      # Multiply by emiss_kg_hr_fill_censor to make sure it's never NA and
      # avoid stan problems later. The filled values don't matter, but have to
      # be strictly positive
      emiss_se_kg_hr_fill_mean = dplyr::if_else(
        is.na(emiss_se_kg_hr),
        emiss_kg_hr_fill_censor * mean_(emiss_se_kg_hr / emiss_kg_hr),
        emiss_se_kg_hr
      ),
      # Also allow for an all-zero case so we can use the same stan code without
      # allowing measurement error.
      # (Stan is finicky about standard errors that are exactly zero, so make
      # this almost zero.)
      emiss_se_kg_hr_all0 = 0,
      gas_price_per_kg_ch4 = !!methane_mcf_per_kg * gas_price_per_mcf * gas_frac_methane,
     ) %>%
    dplyr::filter(!is.na(oil_avg_bbld), !is.na(gas_avg_mcfd))
  df_list$aviris_trunc <- dplyr::filter(df_list$aviris_all, censored == "none")
  df_list$all_measures <- dplyr::bind_rows(
      df_list$aviris_all, df_list$ground_studies
    ) %>%
    dplyr::select(basin, censored, emiss_kg_hr, gas_avg_mcfd, censor_threshold_kg)
  df_list$all_measured_trunc <- dplyr::filter(df_list$all_measures, censored == "none")

  # check for NAs in aviris_all
  # This isn't an exhaustive list.
  aviris_all_without_na <- dplyr::select(df_list$aviris_all,
    emiss_kg_hr_fill0, emiss_kg_hr_fill_censor, well_pad_id, gas_price_per_mcf,
    basin, drill_type, oil_avg_bbld, gas_avg_mcfd,
    age_yr, well_pad_id, drill_type, gas_frac_methane,
    emiss_se_kg_hr_all0, emiss_se_kg_hr_fill_mean
  )
  cols_with_na <- purrr::map_lgl(aviris_all_without_na, anyNA)
  if (any(cols_with_na)) {
    stop(
      "The following columns should never have NA, but they do  : ",
      paste(colnames(aviris_all_without_na)[cols_with_na], collapse=", ")
    )
  }
  stopifnot(
    anyNA(df_list$aviris_all$emiss_kg_hr),
    !anyNA(df_list$aviris_trunc$emiss_kg_hr)
  )
  df_list
}
