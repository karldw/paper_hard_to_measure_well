
# This script sets up the R library. It only needs to be run once.

# Things are a little weird here with the container, so it's worth explaining
# what happens. First, Snakemake asks Conda to create a new environment based on
# code/envs/r_scripts.yml. That environment includes R, cmdstan, and a bunch of
# system dependencies. It doesn't include renv or the rest of the R packages.
# When we run `source(".Rprofile")`, renv will bootstrap itself, then
# we'll use renv::restore to install the rest of the packages.

# Installation happens inside the container, but the project directory is
# bind-mounted in the container, so the files are saved in the project folder.
# That means:
# (1) the install script only has to be run once, rather than every time the
#     container is launched,
# (2) the installed packages will be compatible with the container's OS, but
#     might not be compatible with the host system, and
# (3) you need to re-run this script if either `code/envs/r_scripts.yml` or
#     `renv.lock` are changed (Snakemke should detect this for you).

# Snakemake will *not* load renv's .Rprofile by default, because the starting
# directory is .snakemake/scripts.
# https://github.com/snakemake/snakemake/issues/54
# We could load it manually if we wanted witn `source(".Rprofile")`, but doing
# so creates path issues with the renv auto-loader:
# https://rstudio.github.io/renv/articles/docker.html#handling-the-renv-autoloader
# So we're going to use renv, but not source .Rprofile.

source("code/shared_functions.r") # can't use here::here() here

check_for_bad_env_vars()
Sys.setenv(
  # These settings will ask the Arrow installer to download and build Arrow and
  # its third-party components, except S3 support
  LIBARROW_MINIMAL = "false",
  LIBARROW_BINARY = "false",
  ARROW_S3 = "OFF", # don't need it
  ARROW_JEMALLOC = "OFF", # avoid errors, maybe from memory_limit
  ARROW_WITH_BZ2 = "OFF", # avoid https://issues.apache.org/jira/browse/ARROW-14210
  ARROW_R_DEV = "true",
  # As of 0.14, renv doesn't support parallel restore, so the best thing we can
  # do is try to allow each package to compile in parallel.
  # https://github.com/rstudio/renv/issues/459
  MAKEFLAGS = paste0("-j", snakemake@threads)
)

# Installing sf and lwgeom can have issues with linked libraries.
# https://r-spatial.github.io/sf/index.html
# rgdal would require one too, but we're not using it.
conda_prefix <- Sys.getenv("CONDA_PREFIX")
config_args <- paste0("--with-proj-include=", conda_prefix, "/include/ --with-proj-lib=", conda_prefix, "/lib/")

options(
  renv.config.cache.symlinks = FALSE,
  renv.settings.use.cache = FALSE,
  renv.config.install.verbose = FALSE,
  renv.config.mran.enabled = FALSE,
  renv.config.rspm.enabled = FALSE,
  # renv supports configure.args like this
  configure.args.sf = config_args,
  configure.args.lwgeom = config_args
)
if (isTRUE(getOption("repos") == "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

renv::restore(clean=TRUE)

if (! (arrow::arrow_available() && arrow::arrow_with_dataset() && arrow::arrow_with_parquet())) {
  message("Arrow did not install successfully. Purging for easier reinstall.")
  renv::purge("arrow")
  stop("arrow failed to install", call. = FALSE)
}
