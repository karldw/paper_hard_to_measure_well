# Run one stan model
# Driven by snakemake. See stan models in code/stan_models/

suppressMessages(
  here::i_am("code/run_stan.R", uuid="01458a06-b77f-4e52-bc41-3a0b7fe0a5cd")
)


source(here::here("code/shared_functions.r"))
source(here::here("code/stan_helper_functions.r"))
source(here::here("code/model_data_prep.r"))

if (!exists("snakemake")) {
  warning("using placeholder snakemake")

  snakemake <- SnakemakePlaceholder(
    input = list(
      measurements = here::here("data/generated/methane_measures/matched_wells_all.parquet"),
      stan_file = here::here("code/stan_models/09_twopart_lognormal_heterog_alpha_model.stan"),
      data_prep_script = here::here("code/model_data_prep.r"),
      compiled = here::here("code/stan_models/09_twopart_lognormal_heterog_alpha_model"),
      r_lib = here::here("scratch/snakemake_flags/setup_r_library"),
      cmdstan = here::here("scratch/snakemake_flags/setup_cmdstan")
    ),
    output = list(here::here("scratch/stan_fits/main_spec/09_twopart_lognormal_heterog_alpha-period_8760_hours/model_fit.rds")),
    params = list(),
    log = list(""),
    wildcards = list(
      robustness_spec = "main_spec/",
      model_name = "09_twopart_lognormal_heterog_alpha",
      prior_only = "",
      bootstrap = "",
      time_period = "-period_8760_hours"
    ),
    resources = list(tmpdir = "/tmp", attempt_count = 1, mem_mb=4000),
    threads = 1L
  )
}

log_all_output(snakemake@log)


BOOTSTRAP_ITER <- 100 # If we're bootstrapping, how many iter?
# 50 for good accuracy with reasonable re-run times, but ultimately ~100 is preferable

options(mc.cores=snakemake@threads)
set.seed(7)
mem_limit <- snakemake@resources[['mem_mb']] %||% NA_integer_
try(memory_limit(mem_limit)) # Sometimes fails on Windows?

# Check paths and cmdstan
fs::dir_create(unique(fs::path_dir(unlist(snakemake@output))))
check_cmdstan_path()

model_name <- snakemake@wildcards[["model_name"]]
stan_file <- snakemake@input[["stan_file"]]
robustness_spec <- snakemake@wildcards[["robustness_spec"]] %||% stop("need robustness_spec")

# Load data
sdata <- snakemake@input[["measurements"]] %>%
  prep_measurement_data() %>%
  prep_custom_model_data(
    model_name = model_name,
    robustness_spec_str = robustness_spec
  )

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
  # multicore here might run into issues on windows (excl WSL) or in rstudio)
  if (future::supportsMulticore()) {
    future::plan("multicore")
  } else {
    future::plan("sequential")
  }

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
    # Note: these very high counts are slow. They're sufficient (but probably
    # not necessary) to avoid R-hat errors in stan.
    iter_warmup = 1000,
    iter_sampling = 3000,
    show_messages = FALSE,
    adapt_delta = get_adapt_delta(snakemake@resources[["attempt_count"]]),
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
    iter_warmup = 1000,
    iter_sampling = 3000,
    # avoid expensive mistakes where we hit 100% max treedepth at max_treedepth = 10
    max_treedepth = 9,
    adapt_delta = get_adapt_delta(snakemake@resources[["attempt_count"]]),
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

if (robustness_spec == "main_spec/") {
  if (operate_mode == "bootstrap") {
    should_compare_chains <- FALSE
  } else {
    should_compare_chains <- TRUE
  }
  skip_checks <- character(0)
} else {
  should_compare_chains <- FALSE
  skip_checks <- "rhat"
}

successful_fit <- (operate_mode == "generate") || (sdata$prior_only) ||
  check_cmdstan_diagnostics(
    fit,
    compare_chains = should_compare_chains,
    skip_checks = skip_checks
  )

if (successful_fit) {
  # Guard against improperly using the results if they're wrong.
  save_stan_results(fit, output_files, generate_row_count = generate_row_count)

  # Generate depends on the CSVs from running the model, so need to clean both
  # at the end of the generate step.
  if (operate_mode == "generate") {
    clean_stan_csvs(model_fit = existing_fit, generate_fit = fit)
  }
}


# If we're logging, stop
if (sink.number() > 0) {
  sink(file=NULL, type="output")
  sink(file=NULL, type="message")
}
