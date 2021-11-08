options(warn = 2)
source(here::here("code/shared_functions.r"))
source(here::here("code/stan_helper_functions.r"))
check_cmdstan() # not the same as cmdstanr::check_cmdstan_toolchain

m = cmdstanr::cmdstan_model(
  snakemake@input[[1]],
  include_paths=here::here("code/stan_models"),
  quiet=FALSE,
  cpp_options=list(PRECOMPILED_HEADERS = "false")
)
