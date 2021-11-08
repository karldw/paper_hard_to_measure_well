# Run one stan model
# Driven by snakemake. See stan models in code/stan_models/



source(here::here("code/shared_functions.r"))
source(here::here("code/stan_helper_functions.r"))
source(here::here("code/distribution_model_data_prep.r"))
BOOTSTRAP_ITER <- 100 # If we're bootstrapping, how many iter?
# 50 for good accuracy with reasonable re-run times, but ultimately ~100 is preferable

options(scipen=10, warn=1, mc.cores=snakemake@threads)
set.seed(7)
mem_limit <- snakemake@resources[['mem_mb']] %||% NA_integer_
try(memory_limit(mem_limit)) # Sometimes fails on Windows?

# Check paths and cmdstan
fs::dir_create(unique(fs::path_dir(unlist(snakemake@output))))
check_cmdstan()

model_name <- snakemake@wildcards[["model_name"]]
stan_file <- snakemake@input[["stan_file"]]

# Load data
sdata <- snakemake@input[["measurements"]] %>%
  prep_measurement_data() %>%
  purrr::pluck("aviris_all") %>%
  prep_custom_model_data(model_name = model_name)

# Figure out what kind of operation we're doing.
# a. Standard fitting (run stan model, record parameter draws)
# b. Bootstrap fitting (run stan model with a lot of chains, each one draws gets
#    its own bootstrap sample, record parameter draws)
# c. Generate values using results from either (a) or (b), record
bootstrap_str <- snakemake@wildcards[["bootstrap"]] %||% ""
if (endsWith(stan_file, "generate.stan")) {
  # Note that even if we're using bootstrap results, we want "generate" mode
  operate_mode <- "generate"
} else if (bootstrap_str != "") {
  operate_mode <- "bootstrap"
} else {
  operate_mode <- "standard"
}
prior_only_str <- snakemake@wildcards[["prior_only"]] %||% ""
sdata$prior_only <- prior_only_str != ""
sdata$bootstrap <- operate_mode == "bootstrap"

# Parse time interval definition and add to sdata.
# Only useful for cost_coef_models, but code runs for all of them.
time_period_str <- snakemake@wildcards[["time_period"]] %||% ""
if (time_period_str == "") {
  time_period_hr <- 1
} else {
  # time_period_hr will be an integer here, but that doesn't seem like a big issue
  time_period_hr <- tidy_gsub(time_period_str, "-period_(\\d+)_hours", "\\1") %>%
    as.numeric()
  stopifnot(noNAs(time_period_hr))
}
sdata$time_period_hr <- time_period_hr


if (operate_mode == "generate") {
  sdata_to_write <- unclass(sdata)
  # Record the variable names explicitly (they're otherwise dropped when writing)
  # asinh(gas_avg_mcfd) is recorded as asinhgas_avg_mcfd
  sdata_to_write$X_varnames <- colnames(sdata$X) %||% stop("missing X_varnames")
  # write_stan_json is defined in stan_helper_functions, and allows writing
  # otherwise disallowed variables by naming them in skip_conversion
  write_stan_json(
    sdata_to_write,
    snakemake@output$stan_data_json,
    skip_conversion = "X_varnames"
  )
  # We're going to use future-based parallelism to read the generated files
  # (Stan mananges its own parallel execution). Future will choose the number of
  # threads based on mc.cores, which we set above.
  # multicore here is fine, since we're running this in a linux singularity
  # container (multicore doesn't work on windows (excl WSL) or in rstudio)
  future::plan("multicore")

  existing_fit <- readRDS(snakemake@input[["existing_fit"]])[["fit"]]
  csv_files <- existing_fit$output_files()
  if (!all(file.exists(csv_files))) {
    print(csv_files)
    stop(
      "Expected CSV output files listed for ", model_name, " don't exist. ",
      "Have you run the model (on this computer)?"
    )
  }
  arg_list <- list(
    algorithm = "generate",
    fitted_params = existing_fit
  )
  generate_row_count <- sdata$N
  output_files <- snakemake@output[setdiff(names(snakemake@output), "stan_data_json")] %>%
    purrr::compact() %>%
    purrr::as_vector()
  output_dir <- fs::fs_path(unique(dirname(output_files)))
  if (model_name %in% MODEL_NAMES$cost_coef_models) {
    # Only these cost coef models put out these cost_param files, and it's
    # easier to add them here than have two separate snakemake rules.
    output_files <- c(output_files,
      cost_param_A = output_dir / "cost_param_A.parquet",
      cost_param_alpha = output_dir / "cost_param_alpha.parquet"
    )
  }
  output_names <- names(output_files)
  stopifnot(
    length(output_dir) == 1,
    endsWith(output_files, "parquet"),
    length(output_names) == length(output_files),
    length(unique(output_names)) == length(output_names),
    arrow::codec_is_available("zstd"), # Check that write_parquet call will work
    all(output_names != "")
  )
} else if (operate_mode == "bootstrap") {
  # NOTE: we'll end up with (BOOTSTRAP_ITER * iter_sampling) draws, which
  # can be a lot, but it's helpful to have more draws here because it makes
  # the R-hat checks pass and I'm willing to be a little slow to avoid bugs.
  # We'll thin draws later, but the R-hat checks don't pass if we thin here.
  arg_list <- list(
    algorithm = "sampling",
    iter_warmup = 500,
    iter_sampling = 400,
    show_messages = FALSE,
    adapt_delta = 0.9,
    refresh = 0,
    validate_csv = FALSE, # see cmdstanr #280
    chains = BOOTSTRAP_ITER # number of bootstrap iterations
  )
  unclean_division <- ((arg_list$chains * arg_list$iter_sampling) %% 4000 != 0) ||
    (4000 %% arg_list$chains != 0)
  if (unclean_division) {
      stop(
        "Thinning draws will be weirder and less efficient if either:\n",
        "    a. BOOTSTRAP_ITER * iter_sampling %% 4000 != 0 OR\n",
        "    b. 4000 %% BOOTSTRAP_ITER != 0\n",
        "Please fix."
      )
  }
  generate_row_count <- -1 # -1 signals "not generate mode"
  output_files <- snakemake@output[[1]]
} else {
  # Otherwise, we're running HMC to fit the model.
  # Can add args here, if necessary, to pass to the cmdstanr methods
  arg_list <- list(
    # avoid expensive mistakes where we hit 100% max treedepth at max_treedepth = 10
    max_treedepth = 9,
    adapt_delta = 0.9,
    algorithm = "sampling"
  )
  generate_row_count <- -1 # -1 signals "not generate mode"
  output_files <- snakemake@output[[1]]
}


compiled <- cmdstanr::cmdstan_model(stan_file, include_paths=dirname(stan_file))
# Note that fit_stan_model sets a seed by default
fit <- fit_stan_model(model = compiled, sdata = sdata, args = arg_list)
# In generate mode, the program either errors out or gives results.
# In modeling mode, need to actually check the HMC diagnostics.
successful_fit <- (operate_mode == "generate") || (sdata$prior_only) ||
  check_cmdstan_diagnostics(fit, compare_chains = operate_mode != "bootstrap")

if (successful_fit) {
  # Guard against improperly using the results if they're wrong.
  save_stan_results(fit, output_files, generate_row_count = generate_row_count)
}
