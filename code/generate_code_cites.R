suppressMessages(
  here::i_am("code/generate_code_cites.R", uuid="a15d2768-12b9-457a-ba23-e160281710d5")
)
source(here::here("code/shared_functions.r"))

cite_attached_packages(snakemake@output[[1]], c(
  "base", # R base
  "arrow",
  "cmdstanr",
  "curl",
  "data.table",
  "dplyr",
  "fs",
  "furrr",
  "future",
  "ggplot2",
  "ggpointdensity",
  "glue",
  "gt",
  "here",
  "igraph",
  "jsonlite",
  "lubridate",
  "magrittr",
  "matrixStats",
  "nngeo",
  "posterior",
  "powerjoin",
  "processx",
  "purrr",
  "RColorBrewer",
  "readxl",
  "remotes",
  "rlang",
  "sf",
  "stringr",
  "tibble",
  "tidyr",
  "tidyselect",
  "units",
  "unix",
  "xml2"
))
