options(warn = 2)
suppressMessages(
  here::i_am("code/compile_stan.R", uuid="9f61a9f3-519a-4abc-acc5-03c3c98d9d7b")
)
source(here::here("code/shared_functions.r"))
source(here::here("code/stan_helper_functions.r"))
check_cmdstan_path()

print(here::here("code/stan_models"))
# Note that setup_cmdstan.R has set some cpp flags.
m = cmdstanr::cmdstan_model(
  snakemake@input[[1]],
  stanc_options=list("O1"),
  quiet=FALSE,
  force_recompile=TRUE
)
