options(warn = 2L)
suppressMessages(
  here::i_am("code/setup_cmdstan.R", uuid="17812907-2b57-464f-8ffc-3a7289938d6b")
)

source(here::here("code/stan_helper_functions.r"))
source(here::here("code/version_info_cmdstan.txt"))


check_cmdstan_path()
conda_prefix_env <- Sys.getenv("CONDA_PREFIX") # check_cmdstan_path() checks CONDA_PREFIX

build_flags <- list(
  # Original options:
  # (add these here so we can set append=TRUE and safely re-run this code.)
  "TBB_CXX_TYPE" = "gcc",
  "TBB_INTERFACE_NEW" = "true",
  "TBB_INC" = paste0(conda_prefix_env, "/include/"),
  "TBB_LIB" = paste0(conda_prefix_env, "/lib/"),
  "PRECOMPILED_HEADERS" = "false",
  # end of original options
  # https://blog.mc-stan.org/2022/08/03/options-for-improving-stan-sampling-speed/
  "STAN_CPP_OPTIMS" = "true",
  "CXXFLAGS_OPTIM+= -march=native -mtune=native",
  "CXXFLAGS_OPTIM_TBB+= -march=native -mtune=native",
  "CXXFLAGS_OPTIM_SUNDIALS+= -march=native -mtune=native",
  "CXXFLAGS_WARNINGS+= -w", # suppress all warnings
  "STANC3_VERSION" = CMDSTAN_VERSION
)

cmdstanr::cmdstan_make_local(cpp_options=build_flags, append=FALSE)
# Rebuild with these flags to avoid flag mismatch issues.
# e.g. https://discourse.mc-stan.org/t/recommended-compiler-flags-makes-rstan-model-crash/25689
cmdstanr::rebuild_cmdstan()
# Set flags again in case they were overwritten.
cmdstanr::cmdstan_make_local(cpp_options=build_flags, append=FALSE)
