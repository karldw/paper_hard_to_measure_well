
model_prep_all_possible_formulas <- function(possible_vars) {
  out <- seq_along(possible_vars) %>%
    purrr::map(
      ~combn(possible_vars, m = ., FUN = model_prep_paste_with_plus)
    ) %>%
    purrr::list_c()
  out
}


model_prep_make_list_robustness_var_opts <- function(robustness_var_set, main_spec_vars) {
  # Take every possible combination of variables, except the main spec
  # and single-variable formulas.
  main_spec_rhs <- model_prep_paste_with_plus(main_spec_vars)

  # These are additional specs that are included individually, but not as every
  # possible combination.
  extra_robustness_specs <- model_prep_robustness_specs_not_full_combo(main_spec_vars)

  # Take every possible combination of variables, except the main spec
  # and single-variable formulas.
  all_combos <- model_prep_all_possible_formulas(robustness_var_set) %>%
    setdiff(main_spec_rhs) %>%
    setdiff(robustness_var_set)

  # Then add on the extra_robustness_specs defined above.
  # (Sort them first for conceptual ease later - they'll be robustness1, ...)
  possible_rhs_vars <- c(extra_robustness_specs, all_combos)
  return(possible_rhs_vars)
}


model_prep_robustness_specs_not_full_combo <- function(main_spec_vars) {
  main_spec_no_basin <- setdiff(main_spec_vars, "basin")
  spec_list <- list(
    c(main_spec_vars, "asinh(tot_count_wells_within_10km)"),
    # When including the sub-basin dummies, have to remove the basin dummies
    # clusters created by dbscan -- see group_into_dbscan_clusters.py
    c(main_spec_no_basin, "cluster"),
    c(main_spec_no_basin, "cluster", "asinh(tot_count_wells_within_10km)")
  )
  return(purrr::map_chr(spec_list, model_prep_paste_with_plus))
}


model_prep_parse_robustness_wildcard <- function(robustness_spec_wildcard) {
  reg_match <- stringr::str_match(
    robustness_spec_wildcard,
    "^(main_spec|robustness[0-9]+)/$"
  )[1, 2]
  robustness_idx <- stringr::str_match(reg_match, "^robustness([0-9]+)$")[1, 2]
  stopifnot(!anyNA(reg_match), !anyNA(robustness_idx))
  return(as.integer(robustness_idx))
}


model_prep_get_robustness_rhs <- function(main_spec_vars, robustness_spec_wildcard) {
  robustness_idx <- model_prep_parse_robustness_wildcard(robustness_spec_wildcard)
  # Have the flexibility to make this a different set, but for now it's not.
  robustness_var_set <- main_spec_vars
  num_robustness_vars <- length(robustness_var_set)
  # Minus 8 cases because we're not including: the cases where the modeling
  # includes no RHS variables (1), only one RHS variable (6), or the main
  # spec (1).
  # This cuts down on computation slightly, while still maintaining the idea
  # of including hypothetical alternative sets of RHS covariates.
  max_idx <- (
    (2 ^ num_robustness_vars)
    - 2
    - num_robustness_vars
    + length(model_prep_robustness_specs_not_full_combo(main_spec_vars))
  )

  if (anyNA(robustness_idx) || robustness_idx < 1 || robustness_idx > max_idx) {
    stop(paste0("Robustness spec must be 'main_spec' or 'robustness{X}', where X is between 1 and ", max_idx))
  }
  all_possible_robustness_rhs <- model_prep_make_list_robustness_var_opts(robustness_var_set, main_spec_vars)
  stopifnot(length(all_possible_robustness_rhs) == max_idx)

  return(all_possible_robustness_rhs[robustness_idx])
}
