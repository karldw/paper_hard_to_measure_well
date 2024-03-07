
# This script sets up the R library. It only needs to be run once.

# Things are a little weird here with the conda environment, so it's worth
# explaining what happens. First, Snakemake asks Conda to create a new
# environment based on code/envs/r_scripts.yml. That environment includes R,
# most of the R packages, and a bunch of system dependencies.
# We'll use remotes to install the remaining packages.

# That means:
# (1) the install script only has to be run once, rather than every time the
#     Snakemake is run, and
# (2) you need to re-run this script if `code/envs/r_scripts.yml` or this script
#     are changed (Snakemke should detect this for you).

source("code/shared_functions.r") # can't use here::here() here

if (isTRUE(getOption("repos") == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Install a couple of packages that weren't available in conda-forge

# Installation can fail with only a warning (not an error), but there are also
# unimportant errors (remotes can incorrectly warn about the lack of compiler).
# So load the packages to make sure they work.
remotes::install_version("nngeo", "0.4.3", dependencies=FALSE, upgrade="never")
loadNamespace("nngeo")
remotes::install_version("unix", "1.5.2", dependencies=FALSE, upgrade="never")
loadNamespace("unix")
remotes::install_version("powerjoin", "0.1.0", dependencies=FALSE, upgrade="never")
loadNamespace("powerjoin")


# Note: the official way to install cmdstanr is to use their custom repo
# https://mc-stan.org/r-packages. However, install_version fails.
# remotes::install_version(
#   "cmdstanr", "0.5.3",
#   dependencies=FALSE,
#   upgrade="never",
#   repos=c("https://mc-stan.org/r-packages", getOption("repos"))
# )
remotes::install_github(
  "stan-dev/cmdstanr", ref="v0.5.3",
  dependencies=FALSE,
  upgrade="never",
)
loadNamespace("cmdstanr")

# Check arrow because it's possible to have a sucessfull install that misses
# most of the features.
if (! (arrow::arrow_available() && arrow::arrow_with_dataset() && arrow::arrow_with_parquet())) {
  stop("Arrow did not install successfully. Please reinstall.")
}
