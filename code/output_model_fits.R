
source(here::here("code/shared_functions.r"))
source(here::here("code/distribution_model_data_prep.r"))
source(here::here("code/stan_helper_functions.r"))

set.seed(8675309)
# LEAK_SIZE_DEF defines what leak size we're talking about for a "major leak".
# This needs to match the shift_amount in distribution_model_data_prep
LEAK_SIZE_DEF <- 5

if (!exists("snakemake")) {
  message("Using placeholder snakemake")
  TEX_FRAGMENTS <- fs::fs_path(here::here("output/tex_fragments"))
  STAN_FITS <- fs::fs_path(here::here("data/generated/stan_fits"))
  snakemake <- SnakemakePlaceholder(
    input=list(
      # NOTE: the order here (and in the snakemake file) is the order of the
      # columns in the table.
      distribution_fits = c(
        STAN_FITS / "01_twopart_lognormal-bootstrap/model_fit.rds",
        STAN_FITS / "03_twopart_lognormal_meas_err-bootstrap/model_fit.rds",
        # STAN_FITS / "05_twopart_normal_qr/model_fit.rds",
        STAN_FITS / "02_twopart_lognormal_alpha-bootstrap-period_8760_hours/model_fit.rds",
        STAN_FITS / "08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours/model_fit.rds"
      ),
      measurements = here::here("data/generated/methane_measures/matched_wells_all.rds"),
      stan_data_json = c(
        STAN_FITS / "01_twopart_lognormal-bootstrap/stan_data.json",
        STAN_FITS / "03_twopart_lognormal_meas_err-bootstrap/stan_data.json",
        # STAN_FITS / "05_twopart_normal_qr/stan_data.json",
        STAN_FITS / "02_twopart_lognormal_alpha-bootstrap/stan_data.json",
        STAN_FITS / "08_twopart_lognormal_heterog_alpha-bootstrap/stan_data.json"
      )
    ),
    output = list(
      model_prob_leak_plot = "graphics/model_prob_leak_plot.pdf",
      model_cost_vs_q_plot = "graphics/model_cost_vs_q_plot.pdf",
      model_cost_vs_q_dwl_plot = "graphics/model_cost_vs_q_dwl_plot.pdf",
      model_coef_obs_leak  = TEX_FRAGMENTS / "model_parameters_obs_leak.tex",
      model_coef_leak_size = TEX_FRAGMENTS / "model_parameters_leak_size.tex",
      model_coef_footer_obs_leak  = TEX_FRAGMENTS / "model_parameters_footer_obs_leak.tex",
      model_coef_footer_leak_size = TEX_FRAGMENTS / "model_parameters_footer_leak_size.tex",
      # model_cost_alpha_by_leak_size_bin_plot = "graphics/model_cost_alpha_by_leak_size_bin_plot.pdf",
      # model_leak_size_by_leak_size_bin_plot  = "graphics/model_leak_size_by_leak_size_bin_plot.pdf"
      model_prob_size_above_threshold_histogram = "graphics/model_prob_size_above_threshold_histogram.pdf",
      model_cost_alpha_histogram = "graphics/model_cost_alpha_histogram.pdf"
    ),
    threads=1,
    resources=list(mem_mb = 13000)
  )
}

memory_limit(snakemake@resources[["mem_mb"]])
# Set warn = 2 after memory_limit because it may not work
options(scipen = 99, mc.cores=snakemake@threads, warn = 2)


load_samples <- function(parquet_file, n_draws) {
  if (n_draws == "all") {
    return(arrow::read_parquet(parquet_file))
  }
  stopifnot(n_draws >= 1)
  # Load some columns from a parquet file of draws.
  n_col <- arrow::ParquetFileReader$create(parquet_file)$GetSchema()$num_fields
  if (n_col < n_draws) {
    stop("Only ", n_col, " columns are available. (Use n_draws='all' to load all)")
  }
  # Column names are just numbers, "1", "2", ...
  idx <- sample.int(n_col, size=n_draws, replace=FALSE) %>% as.character()
  arrow::read_parquet(parquet_file, col_select=!!idx)
}


# Plot the modeled prob_leak
plot_model_prob_leak <- function(prob_leak_file, outfile) {
  stopifnot(fs::path_file(prob_leak_file) == "prob_leak.parquet")

  # We have a lot of draws that went into the uncertainty. Pick a subset.
  n_draws <- 200
  draws_to_plot <- load_samples(prob_leak_file, n_draws=n_draws)
  # Reset names:
  colnames(draws_to_plot) <- as.character(seq_len(n_draws))
  draws_to_plot %<>% tibble::as_tibble() %>%
    tidyr::pivot_longer(dplyr::everything(), names_to="draw_id", values_to="prob_leak") %>%
    dplyr::mutate(prob_leak_pct = winsorize(prob_leak * 100, trim=c(0, 0.02)))
  # Different lines, but same color for every line.
  colors <- rep_len(RColorBrewer::brewer.pal(n=3, name="Dark2")[1], n_draws)
  plt <- ggplot2::ggplot(draws_to_plot, ggplot2::aes(x=prob_leak_pct, color=draw_id)) +
    ggplot2::stat_density(geom="line", alpha=0.1) +
    ggplot2::scale_color_manual(values=colors, guide="none") +
    ggplot2::theme_bw() +
    ggplot2::labs(x="Probability of leaking (%; winsorized)", y="Density")

  save_plot(plt, outfile, reproducible=TRUE)
}

read_dist_fit <- function(distribution_fit, depvar = c("leak size", "obs leak")) {
  model_name <- filename_to_model_name(distribution_fit)
  stopifnot(model_name %in% MODEL_NAMES$all)
  depvar <- match.arg(depvar)
  varnames <- get_shared_varnames(model_name = model_name)
  if (depvar == "leak size") {
    # Order matters here -- we're going to match up stan coefs and varnames by
    # position!
    stan_coefs <- c("b_y_intercept", "b_y", "sigma_y")
    varnames %<>% c("sigma")
  } else if (depvar == "obs leak") {
    stan_coefs <- c("b_obs_intercept", "b_obs")
    if (model_name == "02_twopart_lognormal_alpha") {
      return(NULL)
    }
  } else {
    stop("programming error")
  }
  if (depvar == "obs leak" && model_name %in% MODEL_NAMES$rhs_ehat_models) {
    varnames %<>% c("e_hat") # e for emissions
  }
  stopifnot(length(distribution_fit) == 1, grepl('intercept', stan_coefs[1]))
  draws <- readRDS(distribution_fit)$draws %>%
    posterior::subset_draws(stan_coefs)
  list(draws = draws, varnames = varnames, model_name = model_name)
}

clean_varnames <- function(x) {
  # This process could obviously be a lot cleaner:
  pretty_names <- c(
    Intercept = "Intercept",
    asinhgas_avg_mcfd = "IHS of gas prod (mcfd)",
    asinhoil_avg_bbld = "IHS of oil prod (bbld)",
    basinSanJoaquin = "Basin: San Joaquin",
    basinSanJuan = "Basin: San Juan",
    basinOtherCalifornia = "Basin: Other California",
    prod_oil_frac = "Oil prod share",
    asinhage_yr = "IHS of age (yr)",
    drill_typeH = "Drill: Horizontal",
    drill_typeU = "Drill: Unknown",
    drill_typeV = "Drill: Vertical",
    drill_typeD = "Drill: Directional",
    sigma = "$\\sigma$",
    e_hat = "$\\hat{e_i}$"
  )
  missing_names <- setdiff(unique(x), names(pretty_names))
  if (length(missing_names) > 0) {
    stop("Missing pretty names for these variables: ", paste(missing_names, collapse=", "))
  }
  pretty_names[x]
}

get_shared_varnames <- function(model_name,
    measurements_file = snakemake@input[["measurements"]]
  ) {
  stopifnot(length(measurements_file) == 1)
  # Load the data the same way it's loaded for the model fitting.
  # These are the names R creates with model.matrix, so they're messy for factors
  # Could memoize this function, but it's already fast.
  coef_names <- prep_measurement_data(measurements_file) %>%
    purrr::chuck("aviris_all") %>%
    prep_custom_model_data(model_name = model_name) %>% #
    purrr::chuck("X") %>%
    colnames()
  coef_names
}

summarize_coefs_once <- function(draws, varnames) {
  stopifnot(posterior::is_draws(draws), is.character(varnames))
  # col_order <- order(names(df))
  # out <- tidyr::pivot_longer(df, dplyr::everything(), names_to="term", values_to="est") %>%
  #   dplyr::group_by(term) %>%
  #   dplyr::summarize(
  #     estimate = signif(mean(est), 3),
  #     conf_low = signif(quantile_(est, 0.025), 2),
  #     conf_high = signif(quantile_(est, 0.975), 2),
  #     .groups="drop"
  #   ) %>%
  #   dplyr::mutate(term = clean_varnames(term))
  est = function(x) signif(mean(x), 3)
  conf95_low  = function(x) signif(quantile(x, 0.025, names=FALSE, type=8), 2)
  conf95_high = function(x) signif(quantile(x, 0.975, names=FALSE, type=8), 2)

  out <- posterior::summarize_draws(draws,
    estimate = est, conf_low = conf95_low, conf_high = conf95_high
  ) %>%
  dplyr::mutate(term = clean_varnames(varnames))
  stopifnot(
    # Check names match expectations. (dplyr already checks lengths)
    grepl("intercept", out$variable[1], ignore.case=TRUE),
    grepl("intercept", varnames[1], ignore.case=TRUE),
    grepl("[1]", out$variable[2], fixed=TRUE),
    grepl("mcfd", varnames[2], fixed=TRUE)
  )
  # Re-sort rows to match the original column order
  # stopifnot(nrow(out) == length(col_order))
  # dplyr::arrange(out, col_order)
  out
}

.bayes_R2 <- function(y, ypred) {
  # Borrowed from brms
  # Subtract y from every column of ypred and multiply by -1
  # (as.array to make sweep a little stricter)
  e <- -1 * sweep(ypred, MARGIN=1, STATS=as.array(y), FUN="-")
  # These are rowVars in the original brms code, but we have transposed the
  # results so each column is a  MCMC draw, and each row as a well.
  # We want to take variance of each column.
  var_ypred <- matrixStats::colVars(ypred)
  var_e <- matrixStats::colVars(e)
  stopifnot(length(var_ypred) == length(var_e))
  # I'm pretty sure this works for binary outcomes too.
  return(var_ypred / (var_ypred + var_e))
}

calc_r2 <- function(model_dir, outcome_name) {
  stopifnot(
    length(model_dir) == 1,
    length(outcome_name) == 1
  )
  stan_data_json <- file.path(model_dir, "stan_data.json")
  y_raw <- jsonlite::read_json(stan_data_json, simplifyVector=TRUE)$Y %||% stop("Missing 'Y' in sdata")
  obs_indicator <- (!is.na(y_raw)) & (y_raw > LEAK_SIZE_DEF)
  generated_dir <- fs::path_dir(stan_data_json) %>% fs::fs_path()
  generated_file <- generated_dir / paste0(outcome_name, ".parquet")

  # pred value matrix, one row per well, one col per draw.
  ypred <- arrow::read_parquet(generated_file) %>% as.matrix()
  if (outcome_name == "prob_leak") {
    y <- as.numeric(obs_indicator)
  } else if (outcome_name == "leak_size_expect") {
    # For leak size, we can only compare the actual observed leak sizes.
    # rows of ypred are wells.
    y <- y_raw[obs_indicator]
    ypred <- ypred[obs_indicator, , drop=FALSE]
  } else {
    stop("Unknown outcome name: ", outcome_name)
  }
  stopifnot(length(y) == nrow(ypred), ncol(ypred) > 1)
  r2_by_draw <- .bayes_R2(y, ypred)
  stopifnot(noNAs(r2_by_draw), r2_by_draw >= 0, r2_by_draw <= 1)
  # could do CI if we wanted.
  mean(r2_by_draw)
}

get_N <- function(model_dir) {
  stopifnot(length(model_dir) == 1)
  stan_data_json <- file.path(model_dir, "stan_data.json")
  N <- jsonlite::read_json(stan_data_json, simplifyVector=FALSE)$N
  N
}

calc_outcome_mean <- function(model_dir, depvar) {
  stopifnot(length(model_dir) == 1, length(depvar) == 1)
  stan_data_json <- file.path(model_dir, "stan_data.json")
  y_raw <- jsonlite::read_json(stan_data_json, simplifyVector=TRUE)$Y %||% stop("Missing 'Y' in sdata")
  obs_indicator <- (!is.na(y_raw)) & (y_raw > LEAK_SIZE_DEF)
  if (depvar == "obs leak") {
    y <- as.numeric(obs_indicator)
  } else if (depvar == "leak size") {
    y <- y_raw[obs_indicator]
  } else {
    stop("Unknown outcome name: ", depvar)
  }
  mean(y)
}

write_coefs_table <- function(snakemake, depvar) {
  # fit_info is a list of lists. Outer list has one element per file in
  # distribution_fits. For each of those, inner list has elements `draws`,
  # `varnames`, and `model_name`
  if (depvar == "leak size") {
    outfile <- snakemake@output$model_coef_leak_size %||% stop("missing outfile")
  } else if (depvar == "obs leak") {
    outfile <- snakemake@output$model_coef_obs_leak  %||% stop("missing outfile")
  } else {
    stop("bad depvar")
  }
  fit_files <- snakemake@input$distribution_fits %||% stop("missing fit files")
  fit_info <- purrr::map(fit_files, read_dist_fit, depvar=depvar) %>%
    purrr::compact()
  model_names <- extract_from_each(fit_info, "model_name") %>% as.character()
  varnames_lst <- extract_from_each(fit_info, "varnames")
  summary_lst <- extract_from_each(fit_info, "draws") %>%
    purrr::map2(varnames_lst, summarize_coefs_once)
  tab <- summary_lst %>%
    purrr::map(format_estimate_above_interval, align="@") %>%
    merge_estimates_df() %>%
    make_table_fragment(escaped=FALSE, add_comments = model_names) %>%
    writeLines(outfile)
  write_coef_footer(snakemake, depvar)
  invisible(summary_lst)
}

write_coef_footer <- function(snakemake, depvar) {
  model_dirs <- dirname(snakemake@input$stan_data_json)
  model_names <- file.path(model_dirs, "model_fit.rds") %>% filename_to_model_name()
  if (depvar == "leak size") {
    outfile <- snakemake@output$model_coef_footer_leak_size
    outcome_name <- "leak_size_expect"
  } else if (depvar == "obs leak") {
    outfile <- snakemake@output$model_coef_footer_obs_leak
    outcome_name <- "prob_leak"
    # No b_obs coef for this model, so no footer for this model:
    model_dirs <- model_dirs[model_names != "02_twopart_lognormal_alpha"]
  } else {
    stop("bad depvar")
  }
  .make_table_line <- function(x) {
    paste(paste(x, collapse = " & "), "\\\\")
  }
  r2 <- purrr::map_dbl(
    model_dirs,
    calc_r2,
    outcome_name=outcome_name
  ) %>% signif(2)
  n <- purrr::map_int(model_dirs, get_N)
  depvar_mean <- purrr::map_dbl(
    model_dirs,
    calc_outcome_mean,
    depvar=depvar
  ) %>% signif(3)
  waic <- rep_len("-", length(model_dirs))
  text_lines <- c(
    .make_table_line(c("$N$", n)),
    .make_table_line(c("$R^2$", r2)),
    # .make_table_line(c("WAIC", waic)),
    .make_table_line(c("Dep. var. mean", depvar_mean))
  )
  writeLines(text_lines, outfile)
}

which_quantile <- function(x, prob) {
  # Get the index of `x` for the value that's closest to the `prob` quantile.
  # Return value has length 1, even if there are ties
  # Analagous to which.min
  stopifnot(length(prob) == 1, length(x) >= 1)
  q <- quantile_(x, prob)
  which.min(abs(x - q))
}

plot_model_cost_param <- function(model_dir, outfile_cost, outfile_dwl) {
  # Create two plots here:
  # 1. Cost param with uncertainty
  # 2. Shaded DWL
  model_dir <- fs::fs_path(model_dir)
  model_name <- filename_to_model_name(model_dir / "model_fit.rds")
  time_period_hr <- filename_to_time_period_hr(model_dir)

  stopifnot(
    model_name %in% MODEL_NAMES$cost_coef_models,
    length(outfile_cost) == 1,
    length(outfile_dwl) == 1
  )
  cost_param_A_file <- model_dir / "cost_param_A.parquet"
  cost_param_alpha_file <- model_dir / "cost_param_alpha.parquet"
  leak_size_expect_file <- model_dir / "leak_size_expect.parquet"
  sdata <- jsonlite::read_json(model_dir / "stan_data.json", simplifyVector=TRUE)

  gas_frac_ch4 <- 0.95
  ch4_kg_per_mcf <- 18.8916
  # price in sdata is $ per kg CH4
  price_per_mcf <- median(sdata$price) * ch4_kg_per_mcf / gas_frac_ch4
  fee_per_ton_co2e <- 5

  ton_co2e_per_ton_ch4 <- 29.8
  ton_per_kg <- 1 / 1000
  fee_per_mcf <- (
    fee_per_ton_co2e
    * ton_co2e_per_ton_ch4
    * ton_per_kg
    * ch4_kg_per_mcf
    * gas_frac_ch4
  )
  scm_per_kg <- 2 # ($58/ton CO2e; kinda low!_
  scm_per_mcf <- scm_per_kg * ch4_kg_per_mcf * gas_frac_ch4



  # For 02_twopart_lognormal_alpha, the coef is constant, so it doesn't matter
  # that we take the median vs whatever. For 08_twopart_lognormal_heterog_alpha,
  # the coef varies by well, so aggregation potentially matters
  # We'll find the well with the median predicted leak size in draw 1 (arbitrary),
  # and then plot that well's parameters.
  # Doing it this way, instead of taking the median in each draw, gives a better
  # sense of the uncertainty, plotting the variation at one well, rather than
  # the quantile.
  leak_size_kg_per_hr <- arrow::read_parquet(leak_size_expect_file) %>%
    as.matrix()
  well_idx <- which_quantile(leak_size_kg_per_hr[, 1], 0.5)
  cost_param_A <- as.matrix(arrow::read_parquet(cost_param_A_file))[well_idx, , drop=FALSE] %>%
    as.vector()
  cost_param_alpha <- as.matrix(arrow::read_parquet(cost_param_alpha_file))[well_idx, , drop=FALSE] %>%
    as.vector()
  leak_size_kg_per_hr <- as.matrix(arrow::read_parquet(leak_size_expect_file))[well_idx, , drop=FALSE] %>%
    as.vector()
  stopifnot(length(cost_param_A) == 4000, length(cost_param_alpha) == 4000)

  leak_size_kg <- leak_size_kg_per_hr * time_period_hr
  leak_size_mcf <- leak_size_kg / ch4_kg_per_mcf / gas_frac_ch4
  line_private <- price_per_mcf
  line_policy <- (price_per_mcf + fee_per_mcf)
  line_optimal <- (price_per_mcf + scm_per_mcf)

  # To find the uncertainty values we want to plot, find the indexes of the 2.5%
  # and 97.5% MC at q=0.99
  # (I think the function is monotonic in this way, so it doesn't matter that
  # we're only doing one point.)
  cost_param <- dplyr::tibble(
    A = cost_param_A,
    alpha = cost_param_alpha,
    mc_point = A * (1 - 0.99) ^ alpha / !!leak_size_mcf,
  )
  draw_percentile = c(2.5, 50, 97.5)
  draw_idx_of_interest <- c(
    low = which_quantile(cost_param$mc_point, 0.025),
    med = which_quantile(cost_param$mc_point, 0.5),
    high = which_quantile(cost_param$mc_point, 0.975)
  )

  # We're going to plot the y-axis in dollars per mcf for two reasons:
  # 1. It makes the policy lines stable, since otherwise they would vary with the e
  # 2. It gives the audience numbers on a scale they might be familar with (eg. $/mcf commodity price)
  # (Could do $/CO2e for similar reasons)

  q_seq <- seq(0.95, 0.99999, length.out=400)
  to_plot <- dplyr::tibble(
      draw_percentile = draw_percentile,
      # Plot quantiles across draws
      alpha = cost_param$alpha[draw_idx_of_interest],
      A     = cost_param$A[draw_idx_of_interest],
      leak_size_mcf = leak_size_mcf[draw_idx_of_interest],
      # Make a list with a copy of q_seq for each value of draw_idx_of_interest
      q = purrr::map(draw_idx_of_interest, ~q_seq),
    ) %>%
    tidyr::unnest(cols="q") %>%
    ensure_id_vars(q, draw_percentile) %>%
    dplyr::mutate(
      marg_cost_per_mcf = A * (1 - q) ^ alpha / leak_size_mcf,
      draw_percentile_fct = factor(draw_percentile),
    ) %>%
    dplyr::filter(marg_cost_per_mcf <= !!line_optimal * 1.2)

  text_labels <- dplyr::tibble(
    text = c("Private", "Policy", "Social"),
    y = c(line_private, line_policy, line_optimal) + c(1.4, 1.7, 1.7),
    x = min(to_plot$q),
  )

  to_plot_median <- dplyr::filter(to_plot, draw_percentile == 50)
  stopifnot(nrow(dplyr::distinct(to_plot_median, A, alpha)) == 1)
  intersection_y <- c(line_private, line_policy, line_optimal) * unique(to_plot_median$leak_size_mcf)
  # unique here because there are duplicates for every q
  intersection_x <- 1 - (intersection_y / unique(to_plot_median$A)) ^ (1 / unique(to_plot_median$alpha))

  dwl_region_orig <- to_plot_median %>%
    dplyr::filter(dplyr::between(.data$q, intersection_x[1], intersection_x[3]))
  dwl_region_with_policy <- to_plot_median %>%
    dplyr::filter(dplyr::between(.data$q, intersection_x[2], intersection_x[3]))
  dwl_region_colors <- c("gray", "gray")

  plt_cost <- ggplot2::ggplot(to_plot) +
    ggplot2::geom_line(ggplot2::aes(q, marg_cost_per_mcf, color=draw_percentile_fct)) +
    ggplot2::geom_hline(yintercept=line_private, linetype="dashed", alpha=0.7) +
    ggplot2::geom_hline(yintercept=line_optimal, linetype="dashed", alpha=0.7) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values=c("gray", "black", "gray"), guide="none") +
    ggplot2::geom_text(ggplot2::aes(x=x, y=y, label=text), data=dplyr::filter(text_labels, text != "Policy"), hjust=0, size=3) +
    ggplot2::labs(
      x="Prob no leak (q)",
      y="",
      subtitle="$ / mcf"
    )
  plt_dwl <- ggplot2::ggplot(to_plot_median, ggplot2::aes(q, marg_cost_per_mcf)) +
    ggplot2::geom_line() +
    ggplot2::geom_line(data=to_plot, alpha=0) + # add to plot to make the axis the same for both graphs.
    ggplot2::geom_ribbon(ggplot2::aes(ymin=marg_cost_per_mcf, ymax=line_optimal), data=dwl_region_orig,        fill=dwl_region_colors[1], alpha=0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=marg_cost_per_mcf, ymax=line_optimal), data=dwl_region_with_policy, fill=dwl_region_colors[2], alpha=0.8) +
    ggplot2::geom_hline(yintercept=line_private, linetype="dashed", alpha=0.7) +
    ggplot2::geom_hline(yintercept=line_policy, linetype="dashed", alpha=0.7) +
    ggplot2::geom_hline(yintercept=line_optimal, linetype="dashed", alpha=0.7) +
    ggplot2::geom_text(ggplot2::aes(x=x, y=y, label=text), data=text_labels, hjust=0, size=3) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x="Prob no leak (q)",
      y="",
      subtitle="$ / mcf"
    )

  save_plot(plt_cost, outfile_cost, reproducible=TRUE)
  save_plot(plt_dwl, outfile_dwl, reproducible=TRUE)
}

load_data_by_leak_size_bin <- function(model_dir) {
  # This function is only used by plot_model_by_leak_size_bin, but makes some
  # big temp objects, so put it in a separate function
  model_dir <- fs::fs_path(model_dir)
  model_name <- filename_to_model_name(model_dir / "model_fit.rds")
  time_period_hr <- filename_to_time_period_hr(model_dir)
  stopifnot(model_name %in% MODEL_NAMES$cost_coef_models)
  cost_param_alpha_file <- model_dir / "cost_param_alpha.parquet"
  leak_size_draw_file <- model_dir / "leak_size_draw.parquet"
  leak_size_kg_per_hr <- arrow::read_parquet(leak_size_draw_file)

  leak_size_means_by_well <- leak_size_kg_per_hr %>% as.matrix() %>% rowMeans()
  size_group_ids <- tibble::tibble(
    well_id = seq_along(leak_size_means_by_well),
    size_group = factor(cut(leak_size_means_by_well, breaks=5, labels=FALSE)),
  )
  leak_size_long <- leak_size_kg_per_hr %>%
    dplyr::mutate(well_id = dplyr::row_number()) %>%
    # Melt is a few seconds faster than pivot_longer here
    data.table::as.data.table() %>%
    data.table::melt(id.vars="well_id", variable.name="draw_id", value.name="leak_size")
  data.table::setkey(leak_size_long, well_id, draw_id)
  cost_param_alpha_long <- cost_param_alpha_file %>%
    arrow::read_parquet() %>%
    dplyr::mutate(well_id = dplyr::row_number()) %>%
    data.table::as.data.table() %>%
    data.table::melt(id.vars="well_id", variable.name="draw_id", value.name="alpha")
  data.table::setkey(cost_param_alpha_long, well_id, draw_id)
  to_plot <- merge(leak_size_long, cost_param_alpha_long, by=c("well_id", "draw_id"), all=TRUE) %>%
    merge(size_group_ids, by="well_id", all=TRUE) %>%
    dplyr::as_tibble() %>%
    dplyr::select(-well_id)
  stopifnot(noNAs(to_plot))

  return(to_plot)
}

plot_model_by_leak_size_bin <- function(model_dir, outfile_alpha, outfile_leak) {
  # Code to generate:
  # "graphics/model_cost_alpha_by_leak_size_bin_plot.pdf",
  # "graphics/model_leak_size_by_leak_size_bin_plot.pdf",
  # Currently not run because it's slow and I went with a simpler graph.
  stopifnot(
    length(model_dir) == 1,
    length(outfile_alpha) == 1,
    length(outfile_leak) == 1
  )
  # Divide wells into quintiles of (mean) leak size
  # Then plot:
  # 1. ridgeplots of leak size (showing variation both across wells and draws within group)
  # 2. ridgeplots of alpha
  to_plot <- load_data_by_leak_size_bin(model_dir)

  plt_alpha <- ggplot2::ggplot(to_plot,
    ggplot2::aes(
      x = -winsorize(alpha, c(0.01, 0)),
      y = size_group,
      height = stat(density)
    )) +
    ggridges::geom_density_ridges(stat = "binline", bins = 100, scale = 0.95, draw_baseline = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Abatement elasticity by leak size quintile",
      x = "Marginal cost elasticity (−α)",
      y = "Leak size quintile"
    )
  plt_leak_size <- ggplot2::ggplot(to_plot, ggplot2::aes(x=leak_size, y=size_group, height = stat(density))) +
    ggridges::geom_density_ridges(stat = "binline", bins = 100, scale = 0.95, draw_baseline = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = "Leak size distribution by quintile of leak size mean",
      x = "Leak size (kg/hr)",
      y = "Leak size quintile"
    )

  save_plot(plt_alpha, outfile_alpha, reproducible=TRUE)
  save_plot(plt_leak_size, outfile_leak, reproducible=TRUE)
}


plot_alpha_histogram <- function(model_dir, outfile_alpha) {
  stopifnot(
    length(model_dir) == 1,
    length(outfile_alpha) == 1
  )
  cost_param_alpha_file <- fs::fs_path(model_dir) / "cost_param_alpha.parquet"
  cost_param_alpha_mean_by_well <- cost_param_alpha_file %>%
    arrow::read_parquet() %>%
    as.matrix() %>%
    rowMeans()
  plt_alpha <- tibble::tibble(alpha = cost_param_alpha_mean_by_well) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = -winsorize(alpha, c(0.01, 0)),
      y = ..density..
    )) +
    ggplot2::geom_histogram(bins=100) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Abatement elasticity",
      subtitle = "Mean across draws for each well",
      x = "Marginal cost elasticity (−α)",
      y = "Density"
    )
  save_plot(plt_alpha, outfile_alpha, reproducible=TRUE)
}

plot_prob_size_above_threshold_histogram <- function(model_dir, outfile) {
  on.exit(options(warn=2))
  options(warn=1)
  stopifnot(
    length(model_dir) == 1,
    length(outfile) == 1
  )
  pq_file <- fs::fs_path(model_dir) / "prob_size_above_threshold.parquet"
  mean_by_well <- pq_file %>%
    arrow::read_parquet() %>%
    as.matrix() %>%
    rowMeans()
  plt <- tibble::tibble(prob_size_above_threshold = mean_by_well) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = prob_size_above_threshold,
      y = ..density..
    )) +
    ggplot2::geom_histogram(bins=100) +
    ggplot2::theme_bw() +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(
      title = "Probability leak size is greater than detection threshold",
      x = "Pr(e > 100 | X)",
      y = "Density"
    )
  save_plot(plt, outfile, reproducible=TRUE)
}


model_to_plot <- "08_twopart_lognormal_heterog_alpha"
model_dir_to_plot <- snakemake@input$distribution_fits[grepl(model_to_plot, snakemake@input$distribution_fits)] %>%
  dirname() %>%
  fs::fs_path()
stopifnot(length(model_dir_to_plot) == 1)

# plot_model_by_leak_size_bin(
#   model_dir_to_plot,
#   snakemake@output$model_cost_alpha_by_leak_size_bin_plot,
#   snakemake@output$model_leak_size_by_leak_size_bin_plot
# )

plot_prob_size_above_threshold_histogram(
  model_dir_to_plot,
  snakemake@output$model_prob_size_above_threshold_histogram
)

plot_alpha_histogram(
  model_dir_to_plot,
  snakemake@output$model_cost_alpha_histogram
)

plot_model_cost_param(
  model_dir_to_plot,
  snakemake@output$model_cost_vs_q_plot,
  snakemake@output$model_cost_vs_q_dwl_plot
)

plot_model_prob_leak(
  model_dir_to_plot / "prob_leak.parquet",
  outfile = snakemake@output$model_prob_leak_plot
)

write_coefs_table(snakemake, "obs leak")
write_coefs_table(snakemake, "leak size")
