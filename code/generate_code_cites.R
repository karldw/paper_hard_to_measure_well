source(here::here("code/shared_functions.r"))

cite_attached_packages(snakemake@output[[1]], c(
  "base", # R base
  "arrow",
  "brms",
  "broom",
  "cmdstanr",
  "curl",
  "data.table",
  "digest",
  "dplyr",
  "forcats",
  "fs",
  "furrr",
  "ggplot2",
  "future",
  "glue",
  "here",
  "igraph",
  "jsonlite",
  "lubridate",
  "magrittr",
  "matrixStats",
  "nngeo",
  "posterior",
  "processx",
  "purrr",
  "RColorBrewer",
  "readxl",
  "rlang",
  "safejoin",
  "sf",
  "stringr",
  "tibble",
  "tidyr",
  "tidyselect",
  "units",
  "unix"
))