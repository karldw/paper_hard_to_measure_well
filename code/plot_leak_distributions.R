
suppressMessages(
  here::i_am("code/plot_leak_distributions.R", uuid="4b72a9c9-8a15-4a0c-9c9b-92eb6518c6a6")
)

source(here::here("code/shared_functions.r"))
source(here::here("code/model_data_prep.r"))

options(scipen=5)

# Note on variable names:
# mcfd is mcf per day (thousands of cubic feet per day), but is reported monthly
# bbld is barrels per day, but is measured monthly.
# There are other variables that report the average on days when production was
# happening, but I don't use those (maybe I should?).


standardize_columns_for_plotting <- function(df_list, censor_threshold=5) {
  # NOTE: in the output of emiss_kg_hr, zero means emissions could have been
  # detected and quantified but were not. NA means emissions could not have been
  # quantified.
  # Also, this function will drop any dataframes that are not explicitly included.
  out <- list()

  out$jpl_wells_all <- df_list$jpl_wells_all %>%
    dplyr::transmute(
      emiss_kg_hr = dplyr::if_else(is.na(emiss_kg_hr), 0, emiss_kg_hr),
      detect_emiss = emiss_kg_hr > 0,
      src = "California",
    )
  out$four_corners_all_wells <- df_list$four_corners_all_wells %>%
    dplyr::transmute(
      emiss_kg_hr = dplyr::if_else(is.na(emiss_kg_hr), 0, emiss_kg_hr),
      detect_emiss = emiss_kg_hr > 0,
      src = "Four Corners",
    )
  out$lyon <- df_list$lyon %>%
    dplyr::transmute(
      emiss_kg_hr = NA_real_,
      detect_emiss = detect_emiss == 1,
      src = "Lyon et al.",
    )
  out$ground_studies_censored_5kgh <- df_list$ground_studies %>%
    dplyr::transmute(
      emiss_kg_hr = dplyr::if_else(.data$emiss_kg_hr > 5, .data$emiss_kg_hr, 0),
      detect_emiss = .data$emiss_kg_hr > 0,
      src = "Ground studies censored at 5 kg/hr",
    )
  out$ground_studies_censored_10kgh <- df_list$ground_studies %>%
    dplyr::transmute(
      emiss_kg_hr = dplyr::if_else(.data$emiss_kg_hr > 10, .data$emiss_kg_hr, 0),
      detect_emiss = .data$emiss_kg_hr > 0,
      src = "Ground studies censored at 10 kg/hr",
    )

  out$ground_studies <- df_list$ground_studies %>%
    dplyr::transmute(
      emiss_kg_hr = .data$emiss_kg_hr,
      detect_emiss = .data$emiss_kg_hr > 0,
      src = "Ground studies",
    )

  out
}

plot_study_comparison <- function(df_list) {
  df_list %<>% standardize_columns_for_plotting()

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
  plt_density_no_titles <- dplyr::filter(to_plot, emiss_kg_hr > 0) %>%
    dplyr::filter(!src %in% c(ground_studies_fct[2:3], "Theory")) %>%
    ggplot2::ggplot(ggplot2::aes(x=emiss_kg_hr, color=src_color, linetype=src_color)) +
    ggplot2::geom_density(adjust=0.7) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept=5, alpha=0.3, linetype="dashed") +
    ggplot2::theme(
      strip.text.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(color="white", fill="white"),
      legend.position=c(0.2, 0.85)
    ) +
    ggplot2::scale_linetype_manual(values = c(1, 2, 3)) +
    ggplot2::labs(
      x="Methane emissions (kg/hr)",
      y="Density",
      linetype="",
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
    ggplot2::ggplot(ggplot2::aes(x=emiss_kg_hr, color=src, linetype=src)) +
    ggplot2::geom_density(adjust=0.7) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::scale_linetype_manual(values = c(1, 2)) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept=c(5, 10), alpha=0.1, linetype="dashed") +
    ggplot2::annotate("rect", xmin=5, xmax=10, ymin=-Inf, ymax=Inf, alpha = 0.3, fill="gray") +
    ggplot2::annotate("text", x=7, y=1.2, label="Lower\ndetect\nlimit", size=2.5) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(color="white", fill="white"),
      legend.position=c(0.82, 0.85)
    ) +
    ggplot2::labs(
      x="Methane emissions (kg/hr)",
      y="Density",
      linetype="",
      color=""
    )

  save_plot(plt_density_no_titles,
    here::here("graphics/figureA05_leak_comparison_sizes.pdf"),
    scale_mult=0.7,
  )
  # save_plot(plt_density_with_titles,
  #   here::here("graphics/leak_comparison_sizes_extra_titles.pdf"),
  #   scale_mult=1.3,
  # )
  save_plot(plt_density_remote_only,
    here::here("graphics/figure02_leak_comparison_sizes_remote_only.pdf"),
    scale_mult=0.7,
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
    ggplot2::ggplot(ggplot2::aes(x=src, y=pct_detection, color=src_color, linetype=src_color)) +
    ggplot2::geom_bar(stat="identity", alpha=0.8, fill="white") +
    ggplot2::theme_bw() +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::scale_linetype_manual(values = c(1, 2, 3, 4, 4)) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::labs(
      y="Pct with detected emissions",
      x=""
    ) +
    ggplot2::theme(legend.position="none")

  # Make Lyon et al gray here b/c it's not in the accompanying density plot.
  remote_only_colors <- c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "#646464")
  plt_frac_detections_remote_only <- to_plot_remote %>%
    dplyr::group_by(src) %>%
    dplyr::summarize(pct_detection = 100 * mean(detect_emiss), .groups="drop") %>%
    ggplot2::ggplot(ggplot2::aes(x=src, y=pct_detection, color=src, linetype=src)) +
    ggplot2::geom_bar(stat="identity", alpha=0.8, fill="white") +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = remote_only_colors) +
    ggplot2::scale_linetype_manual(values = c(1, 2, 3)) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::labs(
      y="Pct with detected emissions",
      x=""
    ) +
    ggplot2::theme(legend.position="none")

  save_plot(plt_frac_detections,
    here::here("graphics/figureA05_leak_comparison_fraction_with_detect.pdf"),
    scale_mult=0.7,
  )
  save_plot(plt_frac_detections_remote_only,
    here::here("graphics/figure02_leak_comparison_fraction_with_detect_remote_only.pdf"),
    scale_mult=0.7,
  )
  df_list
}

if (!exists("snakemake")) {
  message("This script is meant to be run with snakemake. Using a placeholder.")
  snakemake <- SnakemakePlaceholder(
    input = list(
      cleaned_matched_obs = "data/generated/methane_measures/matched_wells_all.parquet",
      ground_studies = "data/generated/methane_measures/ground_studies.parquet",
      lyon_etal_2016 = "data/generated/methane_measures/lyon_etal_2016.parquet"
    ),
    output = list(
      "graphics/figureA05_leak_comparison_sizes.pdf",
      "graphics/figure02_leak_comparison_sizes_remote_only.pdf",
      # "graphics/leak_comparison_sizes_extra_titles.pdf",
      "graphics/figureA05_leak_comparison_fraction_with_detect.pdf",
      "graphics/figure02_leak_comparison_fraction_with_detect_remote_only.pdf"
    ),
    threads = 4,
    resources = list(mem_mb = 7000),
    rule = ""
  )
}

main <- function(snakemake) {
  df_list <- prep_measurement_data_extra(
    snakemake@input[["cleaned_matched_obs"]],
    snakemake@input[["ground_studies"]]
  )
  df_list$lyon <- arrow::read_parquet(snakemake@input[["lyon_etal_2016"]])
  plot_study_comparison(df_list)
  invisible(NULL)
}

arrow::set_cpu_count(snakemake@threads)
memory_limit(snakemake@resources[["mem_mb"]])

main(snakemake)
