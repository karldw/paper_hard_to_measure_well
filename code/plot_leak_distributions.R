

source(here::here("code/shared_functions.r"))

options(scipen=5, SOURCE_DATE_EPOCH=0)


# Note on variable names:
# mcfd is mcf per day (thousands of cubic feet per day), but is reported monthly
# bbld is barrels per day, but is measured monthly.
# There are other variables that report the average on days when production was
# happening, but I don't use those (maybe I should?).



make_pairs_plots <- function(jpl_sites_matched) {
  if (!rlang::is_installed("GGally")) {
    stop("GGally package required for this function")
  }
  pseudo_log <- scales::pseudo_log_trans() # basically asinh

  plt_matched <- jpl_sites_matched %>%
    dplyr::mutate_at(dplyr::vars(oil_avg_bbld, gas_avg_mcfd, emiss_kg_hr, emiss_se_kg_hr), list(log=pseudo_log$transform)) %>%
    sf::st_drop_geometry() %>%
    GGally::ggpairs(
      ggplot2::aes(color=production_type, alpha=0.7),
      columns = c("emiss_kg_hr_log", "oil_avg_bbld_log", "gas_avg_mcfd_log", "persistence_frac"),
      progress=FALSE
    ) +
    ggplot2::theme_bw()
  save_plot(plt_matched, here::here("graphics/pairs_plot_matched_wells_jpl_di.pdf"), scale=3)

}

make_plots <- function(wells_all, jpl_sites_matched, ground_studies) {
  # GAS_SCALE <- ggplot2::scale_x_continuous(trans="pseudo_log", breaks=c(0, 10, 100, 1000, 10000, 100000))
  # EMISS_SCALE <- ggplot2::scale_y_continuous(trans="pseudo_log", breaks=c(0, 10, 30, 100, 1000))

  make_plot_qq <- function(df, add_line=TRUE, ...) {
    df <- dplyr::filter(df, emiss_kg_hr > 0)

    params <- as.list(MASS::fitdistr(df$emiss_kg_hr, "log-normal")$estimate)
    # The `...` here is just so I can pass in a color aesthetic in the combined case
    plt <- ggplot2::ggplot(df, ggplot2::aes(sample=emiss_kg_hr, ...)) +
      ggplot2::geom_qq(distribution=stats::qlnorm, dparams=params) +
      ggplot2::scale_x_continuous(
        trans="pseudo_log",
        limits=c(-0.2, 3000),
        breaks=c(0, 10, 30, 100, 1000)
      ) +
      ggplot2::scale_y_continuous(
        trans="pseudo_log",
        limits=c(-0.2, 3000),
        breaks=c(0, 10, 30, 100, 1000)
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x="Theoretical distribution (kg/hr)",
        y="Observed distribution (kg/hr)"
      )
    if (add_line) {
      plt <- plt + ggplot2::geom_qq_line(distribution=stats::qlnorm, dparams=params)
    }
    plt
  }
  make_plot_emiss_prod_point_density <- function(df, adjust=0.1) {
    min_emiss <- min_(dplyr::filter(df, emiss_kg_hr > 0)$emiss_kg_hr)
    plt <- df %>%
      dplyr::mutate(emiss_kg_hr_filled =
        dplyr::if_else(is.na(emiss_kg_hr), min_emiss - (0.1 * min_emiss), emiss_kg_hr)) %>%
      dplyr::filter(gas_avg_mcfd > 0) %>%
      ggplot2::ggplot(ggplot2::aes(x=gas_avg_mcfd, y=emiss_kg_hr_filled)) +
      ggpointdensity::geom_pointdensity(adjust=adjust, alpha=0.8) +
      # options are "magma", "plasma", "viridis", or "cividis"
      ggplot2::scale_color_viridis_c(option="inferno")+
      ggplot2::geom_hline(yintercept=min_emiss) +
      ggplot2::scale_x_continuous(
        trans="pseudo_log",
        breaks=c(0, 10^(1:6)),
        limits=c(0, 1.5e6)
      ) +
      ggplot2::scale_y_continuous(trans="pseudo_log", breaks=c(0, 10, 30, 100, 1000)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position="none") +
      ggplot2::labs(
        x="Average gas production (mcf/mo)",
        y="Measured emissions (kg/hr)"#,
        # color="Well\ncount"
      )
    plt
  }
  #
  # plt_ecdf_jpl <- ggplot2::ggplot(jpl_sites_matched, ggplot2::aes(x=emiss_kg_hr)) +
  #   ggplot2::stat_ecdf(geom="step") +
  #   ggplot2::scale_x_log10() +
  #   ggplot2::theme_bw()

  plt_obs_count_jpl <- (make_plot_emiss_prod_point_density(wells_all) +
    ggplot2::labs(
      title="Most wells have no JPL-observed emissisons"
    )) %>%
    save_plot(here::here("graphics/observation_count_jpl_flights.pdf"))

  plt_qq_jpl <- (make_plot_qq(jpl_sites_matched) +
    ggplot2::labs(
      title="Observed JPL measures are distributed log-normal"
    ))
    save_plot(plt_qq_jpl, here::here("graphics/jpl_flights_qq_plot.pdf"))

  ground_studies <- dplyr::rename(ground_studies,
    gas_avg_mcfd = gas_production_mcfd
  )
  plt_obs_count_ground <- (make_plot_emiss_prod_point_density(ground_studies, adjust=0.05) +
    ggplot2::labs(
      title="Ground measures have many fewer zeros"
    )) %>%
    save_plot(here::here("graphics/observation_count_ground_studies.pdf"))

  plt_qq_ground <- (make_plot_qq(ground_studies) +
    ggplot2::labs(
      title="Ground-based measurements"
    )) %>%
    save_plot(here::here("graphics/ground_studies_qq_plot.pdf"))

  to_plot_combined_qq <- dplyr::bind_rows(
      dplyr::mutate(ground_studies, src = "Ground studies"),
      sf::st_drop_geometry(dplyr::mutate(jpl_sites_matched, src = "JPL flights")),
    )
  plot_qq_combined <- make_plot_qq(to_plot_combined_qq, add_line = TRUE, color=src) +
    # ggplot2::facet_grid(rows="src") +
    ggplot2::labs(
      title="Comparison QQ plots"
    )
  save_plot(plot_qq_combined, here::here("graphics/combined_ground_jpl_qq_plot.pdf"))

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
    dplyr::transmute(
      emiss_kg_hr = .data$emiss_kg_hr,
      detect_emiss = .data$emiss_kg_hr > 0,
      src = "Ground studies",
    )

  df_list
}

plot_study_comparison <- function(df_list) {
  df_list %<>% standardize_columns_for_plotting()
  # Here's an version with QQ plots. Haven't run it because the function
  # quantile_inflated_censored_lognormal still needs work.
  # plt <- dplyr::filter(combined_df, src %in% c("California", "Four Corners")) %>%
  #   ggplot2::ggplot(ggplot2::aes(sample=emiss_kg_hr, color=src)) +
  #   ggplot2::geom_qq(distribution=quantile_inflated_censored_lognormal, dparams=qq_params) +
  #   ggplot2::geom_qq_line(distribution=quantile_inflated_censored_lognormal, dparams=qq_params) +
  #   ggplot2::theme_bw() +
  #   ggplot2::scale_color_brewer(palette="Dark2")

  ground_studies_fct <- c(
    "Ground studies" = "Ground studies",
    "Ground studies censored at 10 kg/hr" = "Ground cens\nat 10 kg/hr",
    "Ground studies censored at 5 kg/hr" = "Ground cens\nat 5 kg/hr"
  )

  to_plot <- df_list %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      # Define a separate color variable so we can make all the ground study
      # bars and lines share a color later.
      src_color = dplyr::if_else(
        src %in% names(ground_studies_fct),
        names(ground_studies_fct)[1], src
      ),
      src = factor(src,
        levels = c("California", "Four Corners", "Lyon et al.", names(ground_studies_fct)),
        labels = c("California", "Four Corners", "Lyon et al.", unname(ground_studies_fct))
    ))
  to_plot_remote <- df_list[c("jpl_wells_all", "four_corners_all_wells", "lyon")] %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(src = factor(src, levels = c("California", "Four Corners", "Lyon et al.")))


  # I want to use the same colors for the same
  # colors <- RColorBrewer::brewer.pal(n=4, name="Dark2")
  plt_density_no_titles <- dplyr::filter(to_plot, emiss_kg_hr > 0) %>%
    dplyr::filter(!src %in% c(ground_studies_fct[2:3], "Theory")) %>%
    ggplot2::ggplot(ggplot2::aes(x=emiss_kg_hr, color=src_color)) +
    ggplot2::geom_density(adjust=0.7) +
    ggplot2::theme_bw() +
    # ggplot2::scale_color_manual(values=colors[1:3]) +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept=5, alpha=0.3, linetype="dashed") +
    ggplot2::theme(
      strip.text.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(color="white", fill="white"),
      legend.position=c(0.2, 0.85)
    ) +
    ggplot2::labs(
      x="Methane emissions (kg/hr)",
      y="Density",
      color=""
    )

  plt_density_with_titles <- plt_density_no_titles +
    ggplot2::labs(
      title="Ground studies detect smaller leaks than AVIRIS-NG",
      caption=paste(
        "Dashed line denotes detection threshold.",
        sep="\n"
      )
    ) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))

  plt_density_remote_only  <- dplyr::filter(to_plot_remote, emiss_kg_hr > 0) %>%
    ggplot2::ggplot(ggplot2::aes(x=emiss_kg_hr, color=src)) +
    ggplot2::geom_density(adjust=0.7) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept=c(5, 10), alpha=0.1, linetype="dashed") +
    ggplot2::annotate("rect", xmin=5, xmax=10, ymin=-Inf, ymax=Inf, alpha = 0.3, fill="gray") +
    ggplot2::annotate("text", x=7, y=1.2, label="Detect\nlimit", size=2.5) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(color="white", fill="white"),
      legend.position=c(0.82, 0.85)
    ) +
    ggplot2::labs(
      x="Methane emissions (kg/hr)",
      y="Density",
      color=""
    )

  save_plot(plt_density_no_titles,
    here::here("graphics/leak_comparison_sizes.pdf"),
    scale=0.7,
  )
  save_plot(plt_density_with_titles,
    here::here("graphics/leak_comparison_sizes_extra_titles.pdf"),
    scale=1.3,
  )
  save_plot(plt_density_remote_only,
    here::here("graphics/leak_comparison_sizes_remote_only.pdf"),
    scale=0.7,
  )
  pct_detect_ground <- signif(mean(df_list$ground_studies$detect_emiss) * 100, 2)
  caption_ground_studies_detects <- glue::glue(
    "Without censoring, {pct_detect_ground}% of ground ",
    "observations are positive."
  )
  plt_frac_detections <- to_plot %>%
    dplyr::filter(!src %in% c("Ground studies", "Theory")) %>%
    dplyr::group_by(src, src_color) %>%
    dplyr::summarize(pct_detection = 100 * mean(detect_emiss), .groups="drop") %>%
    ggplot2::ggplot(ggplot2::aes(x=src, y=pct_detection, fill=src_color)) +
    ggplot2::geom_bar(stat="identity", alpha=0.8) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_brewer(palette="Dark2") +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::labs(
      y="Pct with detected emissions",
      x=""
    ) +
    ggplot2::theme(legend.position="none")

  plt_frac_detections_remote_only <- to_plot_remote %>%
    dplyr::group_by(src) %>%
    dplyr::summarize(pct_detection = 100 * mean(detect_emiss), .groups="drop") %>%
    ggplot2::ggplot(ggplot2::aes(x=src, y=pct_detection, fill=src)) +
    ggplot2::geom_bar(stat="identity", alpha=1) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_brewer(palette="Dark2") +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::labs(
      y="Pct with detected emissions",
      x=""
    ) +
    ggplot2::theme(legend.position="none")

  save_plot(plt_frac_detections,
    here::here("graphics/leak_comparison_fraction_with_detections.pdf"),
    scale=0.7,
  )
  save_plot(plt_frac_detections_remote_only,
    here::here("graphics/leak_comparison_fraction_with_detections_remote_only.pdf"),
    scale=0.7,
  )
  df_list
}


if (!exists("snakemake")) {
  snakemake <- SnakemakePlaceholder(
    input = list(
      # well_pad_crosswalk = "data/generated/production/well_pad_crosswalk_1970-2018.parquet",
      # headers = glue::glue("data/generated/production/well_headers/first_prod_year={year}/file.parquet", year=1990:2018),
      # prod = glue::glue("data/generated/production/monthly_production/year={year}/file.parquet", year=1990:2018),
      # alvarez_2018 = "data/studies/alvarez_etal_2018/aar7204_Database_S1.xlsx",
      # omara_2018 = "data/studies/omara_etal_2018/Omara_etal_SI_tables.csv",
      # duren_2019_plumes = "data/studies/duren_etal_2019/Plume_list_20191031.csv",
      # duren_2019_sites = "data/studies/duren_etal_2019/41586_2019_1720_MOESM3_ESM.xlsx",
      # lyon_etal_2016_locations = "data/studies/lyon_etal_2016/es6b00705_si_005.xlsx",
      # lyon_etal_2016_measures = "data/studies/lyon_etal_2016/es6b00705_si_004.xlsx",
      # frankenberg_etal_2016_sources = "data/studies/frankenberg_etal_2016/AVNG_sources_all2.xlsx",
      # frankenberg_etal_2016_measures = "data/studies/frankenberg_etal_2016/FourCorners_AV_NG_detections_Werner.xlsx",
      # zavala_araiza_2018 = "data/studies/zavala-araiza_etal_2018/elementa-6-284-s1.xlsx",
      # citygate_price = "data/generated/eia/citygate_price.parquet"

      cleaned_matched_obs = "data/generated/methane_measures/matched_wells_all.rds"
    ),
    output = list(
      plot_obs_count_jpl = "graphics/observation_count_jpl_flights.pdf",
      plot_jpl_flights_qq = "graphics/jpl_flights_qq_plot.pdf",
      lyon_etal_2016 = "data/generated/methane_measures/lyon_etal_2016.parquet",
      aviris_match_fraction_dropped = "output/tex_fragments/aviris_match_fraction_dropped.tex"
    ),
    threads = 4,
    resources = list(mem_mb = 7000),
    rule = ""
  )
}

main <- function(snakemake, make_extra_plots=FALSE) {

  #
  # jpl_wells_all = jpl_wells_all,
  # four_corners_all_wells = four_corners_all_wells,
  # lyon = lyon,
  # ground_studies = ground_studies
  df_list <- readRDS(snakemake@input[["cleaned_matched_obs"]])
  plot_study_comparison(df_list)

  if (make_extra_plots) {
    # Make plots of JPL measures and ground studies.
    make_plots(jpl_wells_all, jpl_sites_matched, ground_studies)
    # Make plots of distributions across studies.
  }


  # Create output data:
  # To add to this list, make sure the list names match names in snakemake@output
  # list(lyon_etal_2016 = read_lyon_etal_2016(snakemake@input)) %>%
    # write_datasets(snakemake@output)
}

arrow::set_cpu_count(snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])

main(snakemake)
