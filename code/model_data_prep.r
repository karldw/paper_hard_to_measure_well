suppressMessages(
  here::i_am("code/model_data_prep.r", uuid="b431dcd9-a312-4e59-a7ab-e68b726534cf")
)
# Put code related to the robustness spec in a different file, so Snakemake
# can track it separately.
source(here::here("code/model_data_prep_robustness.r"))

if (!existsFunction("%>%")) {
  `%>%` <- magrittr::`%>%`
}

prep_custom_model_data <- function(df, model_name, shift_amount=5, robustness_spec_str="main_spec/") {
  reg_spec <- get_robustness_spec(robustness_spec_str)
  rhs_vars <- reg_spec$rhs_vars
  stopifnot(
    length(rhs_vars) == 1,
    is.character(rhs_vars),
    inherits(df, "data.frame"),
    shift_amount >= 1,
    shift_amount < 100,
    length(model_name) == 1
  )


  if (model_name %in% MODEL_NAMES$meas_err_models) {
    # See notes above about def of emiss_se_kg_hr_fill_mean
    fm_string <- paste0(
      "emiss_kg_hr_fill0 | mi(emiss_se_kg_hr_fill_mean) ~ ", rhs_vars
    )
    stop("this model class is no longer supported. If using this type of model, need to re-add brms=2.15.0 as a dependency")
  } else {
    # NB: we may take logs of LHS in the stan code (depending on distribution)
    # These asinh transformations are applied for real though.
    fm_string <- paste0("emiss_kg_hr_fill0 ~ ", rhs_vars)
  }
  .with_dim <- function(x) {
    if (is.null(dim(x))) {
      dim(x) <- length(x)
    }
    x
  }
  # If we were using brms:
  # fm <- brms::brmsformula(fm_string)
  # sdata <- brms::make_standata(fm, df)
  fm <- formula(fm_string)
  model_frame <- model.frame(fm, data = df)
  model_mat <- model.matrix(fm, data = df)
  # substitute for brms::make_standata without all of brms's deps
  colnames(model_mat) <- gsub("[\\(\\) ]", "", colnames(model_mat))
  sdata <- list(
    N = nrow(model_frame),
    Y = .with_dim(model_frame[, 1L, drop=TRUE]),
    K = ncol(model_mat),
    X = model_mat,
    prior_only = 0
  )
  class(sdata) <- c("standata", "list")

  # Check that we didn't drop any missings:
  # (in emiss_kg_hr_fill0, non-detections are coded as zero)
  if (length(sdata$Y) != nrow(df) || !all(is.finite(sdata$Y))) {

    print(robustness_spec_str)
    print(fm_string)
    print(head(model_frame))
    print(sapply(df, anyNA))
    print(length(sdata$Y))
    print(nrow(df))
    print(mean(is.finite(sdata$Y)))
    stop("missing/invalid values for specification printed above")
  }
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
  # For models where gas_avg_mcfd is not in the set of explanatory variables,
  # still need this variable to show up on in sdata because it's used in
  # outcomes_analysis.py
  sdata$gas_avg_mcfd <- df$gas_avg_mcfd %||% stop("gas_avg_mcfd is missing")
  stopifnot(noNAs(sdata)) # may be Inf though
  sdata
}


prep_measurement_data <- function(filename) {
  censor_threshold_aviris_kg <- 5 # 5-10 in the literature
  methane_mcf_per_kg <- 1 / 18.8916 # Note: this is for pure CH4, not nat gas

  # NB: I'm filling emiss_kg_hr with the censor value, but that's the value brms
  # expects in the modeling.
  # I'm not assuming that these wells actually emit that quantity.
  aviris_all <- arrow::read_parquet(filename) %>%
    dplyr::select(
      -production_type, -well_pad_lon, -well_pad_lat, -county, -months_produced,
    ) %>%
    dplyr::mutate(
      censor_threshold_kg = !!censor_threshold_aviris_kg,
      emiss_kg_hr_fill0 = ifelse(is.na(emiss_kg_hr), 0, emiss_kg_hr),
      emiss_kg_hr_fill_censor = ifelse(is.na(emiss_kg_hr), censor_threshold_kg, emiss_kg_hr),
      censored = ifelse(emiss_kg_hr_fill0 <= censor_threshold_kg, "left", "none"),
      # Define a continuous measure of how much of this well's production
      # comes from oil. Using DI's simple equivalence 6 mcf = 1 bbl.
      prod_oil_frac = calc_oil_gas_frac(oil_avg_bbld, gas_avg_mcfd),
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

  # check for NAs in aviris_all
  # This isn't an exhaustive list.
  aviris_all_without_na <- dplyr::select(aviris_all,
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
  stopifnot(anyNA(aviris_all$emiss_kg_hr))
  aviris_all
}


model_prep_paste_with_plus <- function(x) {
  # sort here to ensure consistent ordering when we remove the main spec in model_data_prep_robustness.r
  paste(sort(x), collapse = " + ")
}


#' Return a list specifying what robustness specification we're doing
#'
#' @param robustness_spec_wildcard A string that matches the regex
#' `^(main_spec|robustness[0-9]+)/$`
get_robustness_spec <- function(robustness_spec_wildcard) {

  stopifnot(length(robustness_spec_wildcard) == 1)

  # NOTE: order matters here because it will affect the numbering of the
  # resulting `robustness_idx` (that is, which specification is numbered 17
  # depends on the order of these variables, but all specifications are included
  # eventually.
  main_spec_vars <- c(
    "asinh(gas_avg_mcfd)",
    "asinh(oil_avg_bbld)",
    "basin",
    "prod_oil_frac",
    "asinh(age_yr)",
    "drill_type"
  )

  if (robustness_spec_wildcard == "main_spec/") {
    rhs_vars <- model_prep_paste_with_plus(main_spec_vars)
  } else {
    rhs_vars <- model_prep_get_robustness_rhs(main_spec_vars, robustness_spec_wildcard)
  }
  if (length(rhs_vars) != 1 || anyNA(rhs_vars)) {
    stop(paste0("Bad rhs_vars: '", rhs_vars, "'"))
  }
  # Note: only one element to this list now, but could add other things, like
  # which observations to keep.
  spec <- list(
    rhs_vars = rhs_vars
  )
  return(spec)
}
