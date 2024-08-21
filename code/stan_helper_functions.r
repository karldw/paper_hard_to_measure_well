Sys.setenv(CMDSTANR_NO_VER_CHECK = "true")


check_cmdstan_path <- function(verbose=TRUE) {
  source(here::here("code/version_info_cmdstan.txt"))

  # cmdstan is provided in the r_scripts.yml, so this script only has to do a
  # couple of checks and set some flags.

  cmdstan_env <- Sys.getenv("CMDSTAN")
  conda_prefix_env <- Sys.getenv("CONDA_PREFIX")
  if (cmdstan_env == "") {
    stop("CMDSTAN env variable should be set", call.=FALSE)
  }
  if (conda_prefix_env == "") {
    stop("CONDA_PREFIX env variable should be set", call.=FALSE)
  }
  if (!startsWith(cmdstan_env, here::here())) {
    stop(paste0(
      "CMDSTAN env variable should be inside the project, but was not:\n",
      "- CMDSTAN:  '", cmdstan_env, "'\n",
      "- This dir: '", here::here(), "'"
    ), call.=FALSE)
  }
  if (!startsWith(cmdstan_env, conda_prefix_env)) {
    stop(paste0(
      "CMDSTAN should be a subdirectory of CONDA_PREFIX\n",
      "- CMDSTAN:      '", cmdstan_env, "'\n",
      "- CONDA_PREFIX: '", conda_prefix_env, "'"
    ), call.=FALSE)
  }
  installed_version <- cmdstanr::cmdstan_version()
  if (installed_version != CMDSTAN_VERSION) {
    stop(paste0(
      "Installed version (", installed_version, ") ",
      "does not match expected version (", CMDSTAN_VERSION, ")."
    ), call.=FALSE)
  }

  if (verbose) {
    message(cmdstanr::cmdstan_path())
    cmdstanr::check_cmdstan_toolchain()
  }
  invisible(NULL)
}

call_with_args <- function(fn, arg_list, verbose=FALSE) {
  # First slice out unused args.
  allowed_args <- formalArgs(fn) %>% setdiff("...")
  args <- arg_list[allowed_args] %>% purrr::compact()
  stopifnot(rlang::is_dictionaryish(arg_list), rlang::is_dictionaryish(args))
  if (verbose) {
    ignored_arg_names <- setdiff(names(arg_list), allowed_args)
    if (length(ignored_arg_names) > 0) {
      message("Not using arguments: ", paste(ignored_arg_names, collapse=", "))
    }
  }
  # Then call the function with those args
  do.call(fn, args)
}

fit_stan_model <- function(model, sdata, args) {
  # Avoid warnings created by the underlying cmdstanr code.
  orig_warnPartialMatchArgs <- getOption("warnPartialMatchArgs")
  options(warnPartialMatchArgs = FALSE)
  on.exit(options(warnPartialMatchArgs = orig_warnPartialMatchArgs), add=TRUE)

  class(sdata) <- "list"
  args$data <- sdata
  args$output_dir <- here::here("scratch/stan_output")
  if (args$data$prior_only %||% FALSE) {
    args$output_dir <- here::here("scratch/stan_output/prior_only")
  }
  dir.create(args$output_dir, recursive=TRUE, showWarnings=FALSE)

  # Fill in some default args
  # (unset args will use cmdstan defaults)
  defaults <- list(
    chains = 4L,
    seed = 8675309 # different than brm
  )
  # If the seed isn't already set, set it in the defaults.
  # (I want to set a seed for each run_stan call because I want results that
  # are reproducible on a per-call basis.)
  for (d in names(defaults)) {
    # Number of brackets matter here. The code below should work even if
    # defaults[d] is NULL.
    # (Circle 8.1.56 of http://www.burns-stat.com/pages/Tutor/R_inferno.pdf)
    if (is.null(args[[d]])) {
      args[d] <- defaults[d]
    }
  }
  algorithm <- args$algorithm %||% "sampling"
  args$algorithm <- NULL
  if (algorithm == "sampling") {
    if (all(c("warmup", "iter") %in% names(args))) {
      args$iter_sampling <- args$iter - args$warmup
      args$iter_warmup <- args$warmup
      args$warmup <- NULL
      args$iter <- NULL
    } else if ("iter" %in% names(args)) {
      args$iter_warmup <- floor(args$iter / 2)
      args$iter_sampling <- args$iter - args$iter_warmup
      args$iter <- NULL
    } else if ("warmup" %in% names(args)) {
      args$iter_warmup <- args$warmup
      args$iter_sampling <- args$warmup
      args$warmup <- NULL
    }
    fit <- do.call(model$sample, args)
  } else if (algorithm %in% c("fullrank", "meanfield")) {
    # vb does not support parallel execution
    fit <- call_with_args(model$variational, args)
  } else if (algorithm == "optimize") {
    fit <- model$optimize(sdata)
  } else if (algorithm == "generate") {
    fit <- call_with_args(model$generate_quantities, args)
  } else {
    stop("Algorithm '", algorithm, "' is not supported.")
  }
  fit
}

fit_metadata <- function(fit) {
  files <- fit$output_files()
  purrr::map(files, cmdstanr::read_cmdstan_csv, variables="") %>%
    extract_from_each("metadata")
}

#' Extract draws with tidyselect syntax
#' For Stan arrays/vectors/matrices, all elements are selected.
extract_draws <- function(fit, ...) {
  # Avoid warnings created by the underlying posterior code.
  orig_warnPartialMatchArgs <- getOption("warnPartialMatchArgs")
  options(warnPartialMatchArgs = FALSE)
  on.exit(options(warnPartialMatchArgs = orig_warnPartialMatchArgs), add=TRUE)

  # NOTE: this would be simpler with fit$draws(), but because different chains
  # might have different parameters, fit$draws() fails.
  # See https://github.com/stan-dev/cmdstanr/issues/280
  files <- fit$output_files()
  # Read each CSV. From each set of results, pluck out `metadata` element.
  # From each `metadata` element, pluck out `model_params`, then find the set
  # intersection of all param names.
  all_model_param <- fit_metadata(fit) %>%
    extract_from_each("model_params") %>%
    purrr::reduce(intersect)
  stopifnot(
    length(files) >= 1,
    length(files) == length(fit$output_files(include_failed=TRUE)),
    is.character(all_model_param),
    length(all_model_param) >= 1
  )
  # Remove numbering for vector/matrix param
  available <- tidy_gsub(all_model_param, "\\[\\d+\\]$", "") %>% unique()
  # vars_select is "questioning" as of July 2020. Hopefully it's not deprecated,
  # because the replacement eval_select works much worse here.
  # strict=FALSE because we do things like -tidyselect::starts_with("temp"),
  # and we don't want the code to fail if no variables start with "temp"
  selected <- tidyselect::vars_select(available, ..., .strict=FALSE)
  # Recreate the default behavior where fit$draws() returns everything.
  if (length(selected) == 0) {
    selected <- available
  }
  message("  Using files:\n  - ", paste0(files, collapse="\n  - "))
  draws_lst <- purrr::map(files, cmdstanr::read_cmdstan_csv, variables=selected) %>%
    extract_from_each("post_warmup_draws")
  # bind_draws doesn't have a method for a list, but does handle ... args.
  # Use do.call to translate (could also use do_call or purrr::lift_dl). They're
  # all the same idea.
  draws <- do.call(posterior::bind_draws, c(draws_lst, along="chain"))
  draws
}

compare_models <- function(fit_list, criterion=c("loo", "waic", "kfold"), ...) {
  stop("This function is no longer used. If re-enabling, add loo=2.4.1 and brms=2.15.0 to the list of dependencies.")
  criterion <- match.arg(criterion)
  # Minor workaround because brms::loo_compare doesn't accept a list of models
  criterion_fn <- getExportedValue("brms", criterion)
  stopifnot(rlang::is_dictionaryish(fit_list))
  stat_list <- purrr::map(fit_list, criterion_fn, ...)
  loo::loo_compare(stat_list)
}

compile_all_models <- function() {
  d <- here::here("code/stan_models")
  files <- c(
    list.files(d, "*model.stan", full.names=TRUE),
    list.files(d, "*generate.stan", full.names=TRUE)
  )
  for (f in files) {
    message(basename(f), ": ", appendLF=FALSE)
    cmdstanr::cmdstan_model(f, include_paths=d)
  }
  invisible(NULL)
}

cmdstan_diagnose <- function(files, compare_chains=TRUE) {
  stopifnot(length(files) >= 1)
  # This is adapted from the cmdstanr cmdstan_diagnose tool, but allowing for
  # diagnosing files individually or as a group. (It's also quieter)
  target_exe = file.path("bin", cmdstanr:::cmdstan_ext("diagnose"))
  cmdstanr:::check_target_exe(target_exe)
  if (compare_chains || length(files) == 1) {
    all_res <- processx::run(
      command = target_exe,
      args = files,
      wd = cmdstanr::cmdstan_path(),
      echo_cmd = FALSE,
      echo = FALSE,
      error_on_status = TRUE
    )
    # Length-1 character vector
    run_log <- all_res$stdout
  } else {
    # Run this function on each file individually (useful for bootstrapping)
    run_log <- purrr::map_chr(files, cmdstan_diagnose)
  }
  return(run_log)
}

check_cmdstan_diagnostics <- function(fit, quiet=FALSE, compare_chains=TRUE, skip_checks=character(0)) {
  # Searching the output text seems terrible, it's otherwise a hassle to extract
  # all four diagnostics.
  # Sure would be nice if I didn't have to search output strings...

  # Note on treedepth:
  # Used to regex search for "Treedepth satisfactory for all transitions"
  # However, some of the robustness specs have treedepth issues.
  # E.g.
  #   1 of 1000 (0.10%) transitions hit the maximum treedepth limit of 9, or 2^9 leapfrog steps.
  #   Trajectories that are prematurely terminated due to this limit will result in slow exploration.
  #   For optimal performance, increase this limit.
  # And treedepth is an efficiency check, not a correctness check, so I'm no
  # longer checking it.
  fit_method <- fit_metadata(fit) %>% purrr::chuck(1L, "method")
  if (fit_method == "sample") {
    regex_good_results <- c(
      divergence = "No divergent transitions found",
      # treedepth = "Treedepth satisfactory for all transitions",
      energy = "E-BFMI satisfactory",
      rhat = "R-hat values satisfactory all parameters"
    )
  } else if (fit_method == "variational") {
    regex_good_results <- c(
      ess = "Effective sample size satisfactory",
      rhat = "R-hat values satisfactory all parameters",
      no_problems = "Processing complete, no problems detected"
    )
  } else {
    warning("Note: don't know how to check results from method ", fit_method)
    return(TRUE)
  }
  regex_good_results <- regex_good_results[setdiff(names(regex_good_results), skip_checks)]
  diagnose_res <- cmdstan_diagnose(
    fit$output_files(include_failed=FALSE),
    compare_chains = compare_chains
  )
  # cmdstan_diagnose returns a vector (see above). check_once checks one element
  # returning TRUE if all the regex results pass.
  # (By default there is only one element, but will be more if there are
  # multiple chains and compare_chains == FALSE)
  check_once <- function(txt, regex_good_results) {
    all(purrr::map_lgl(regex_good_results, base::grepl, x=txt, fixed=TRUE))
  }
  res <- purrr::map_lgl(diagnose_res, check_once, regex_good_results)
  # Recall that if compare_chains == TRUE, length(res) == 1, but if we're
  # running a lot of short bootstrap chains, we need to accept some error.
  good_enough <- mean(res) > 0.75
  if (!quiet && !good_enough) {
    # cat() out the ones that had any errors:
    purrr::walk(diagnose_res[!res], cat)
    if (compare_chains) {
      msg <- "\nFit performed poorly. See details above."
    } else {
      msg <- glue::glue(
        "\nFit performed poorly in {sum(!res)} of {length(res)} chains. See details above."
      )
    }
    warning(msg)
  }
  good_enough
}

fread_transpose <- function(filename, idx, skip, select, target_draws_per_file) {
  # Pull the draws, which initially have one row per MCMC iteration and one
  # column per observation, then transpose so columns are iterations.
  # data.table already uses only 1 thread for forked processes, but set
  # nThread=1 here in case, e.g., we change to multiprocess parallelism
  mat <- data.table::fread(filename, skip=skip, select=select, nThread=1L) %>%
    base::t()
  col_idx_to_keep <- thin_out_vector(seq_len(ncol(mat)), target_draws_per_file)
  mat <- mat[, col_idx_to_keep, drop=FALSE]
  # Assign sequential numbers (1:1000 when idx == 1, 1001:2000 etc)
  n_draws <- ncol(mat)
  colnames(mat) <- as.character(seq_len(n_draws) + ((idx - 1) * n_draws))
  tibble::as_tibble(mat, .name_repair = "minimal")
}

save_stan_results <- function(fit, output_files, generate_row_count) {
  stopifnot(is.finite(generate_row_count))
  # Are we saving model results (coefficients) or predicted values?
  if (generate_row_count == -1) {
    output_files %<>% unlist()
    if (length(output_files) != 1) {
      print(output_files)
      stop("Expected exactly 1 output file, but got the ones above")
    }
    # Pull out the parameter draws, dropping things that start with temp, or are
    # exactly "Ymi" or "lp__" (it's fine if any of these are not present)
    draws <- extract_draws(fit, -tidyselect::starts_with("temp"), -Ymi, -lp__)
    if (anyNA(draws)) {
      stop("Draws should not be NA!")
    }
    saveRDS(list(fit=fit, draws=draws), output_files)
  } else {
    # Make a small helper function so we can iterate over variables easily.
    # Variables in the generate output map 1:1 to file names.
    helper_fn <- function(output, varname, fit) {
      start_time <- Sys.time()
      files <- fit$output_files(include_failed=FALSE)
      if (length(files) != length(fit$output_files(include_failed=TRUE))) {
        stop("Some stan chains failed.")
      }
      one_file <- files[1]
      all_col_names <- data.table::fread(one_file, skip=varname, nrows=0) %>%
        names()
      col_select <- which(grepl(varname, all_col_names, fixed=TRUE))
      stopifnot(length(col_select) > 100)

      # Run this reading in parallel (this is the slow part of this function -
      # about 1.9 seconds per file). furrr will use whatever the current
      # future::plan is. This is set in run_stan.R.
      # Also note that we could run in parallel at a higher level - over
      # output_files - but doing so would use more memory.
      # skip=varname means skip the stan headers until you find a line that
      # contains the varname.
      df <- furrr::future_map2_dfc(
        files, seq_along(files), fread_transpose,
        skip=varname, select=col_select,
        target_draws_per_file = ceiling(4000 / length(files)) # thin out draws
      )

      stopifnot(names(df) == as.character(seq_len(ncol(df))))
      if (anyNA(df)) {
        stop("Draws should never be NA, but they were for ", varname)
      }
      if (nrow(df) != generate_row_count) {
        stop(glue::glue("Num rows {nrow(df)} doesn't match expected {generate_row_count}"))
      }
      # Use dictionary encoding if there aren't many unique values (this doesn't
      # really matter, but makes for faster writes and smaller files in the
      # other case where values are all unique.)
      use_dictionary <- length(unique(df[["1"]])) < 100
      # NB: iterating columns in R-arrow parquet is really slow.
      # See https://issues.apache.org/jira/browse/ARROW-9557
      # We'll be working in python.
      # Use zstd instead of snappy here because it's about 10% smaller.
      # Compression higher than 5 doesn't give size improvements, but is slower.
      arrow::write_parquet(df, output, version="2.0",
        use_dictionary=use_dictionary, compression="zstd", compression_level=5
      )
      message(glue::glue("  - {output} (in {elapsed_time_str(start_time)})"))
      output
    }
    output_names <- names(output_files)
    stopifnot(
      length(output_names) == length(output_files),
      all(output_names != "")
    )
    purrr::walk2(output_files, output_names, helper_fn, fit=fit)
  }
}

thin_out_vector <- function(x, target) {
  stopifnot(length(target) == 1, target >= 1, length(x) >= 1)
  target <- ceiling(target)
  if (length(x) <= target) {
    return(x)
  }
  # Two passes here:
  # 1. Take every nth element of x
  # 2. Drop any remainder
  take_every_nth <- length(x) %/% target
  keep_idx <- seq.int(from=1, to=length(x), by=take_every_nth)
  if (length(keep_idx) > target) {
    keep_idx <- keep_idx[seq_len(target)]
  }
  x[keep_idx]
}

thin_out_draws <- function(draws, target = 4000) {
  stopifnot(posterior::is_draws_array(draws), target > 0)
  # If we have more draws than we want, first thin, then truncate.
  # See thin_out_vector for a function that does the same for a vector.
  # Result is guaranteed to be at least target (when num draws > desired)
  # and less than or equal to target + nchains(draws)
  num_draws <- posterior::ndraws(draws)
  if (num_draws > target) {
    # thin_draws extracts every nth draw, where n has to be an integer.
    # integer division here will result in an integer that's the lower bound of
    # the number we would want if we could do fractional thinning.
    thin <- num_draws %/% target
    draws %<>% posterior::thin_draws(thin)
    # If thin didn't get us all the way there, truncate the remainder of each
    # chain.
    num_draws <- posterior::ndraws(draws)
    if (num_draws > target) {
      trunc_index <- (target %/% posterior::nchains(draws)) + 1L
      # Check here to avoid indexing outside
      if (trunc_index < posterior::niterations(draws)) {
        draws <- draws[seq_len(trunc_index), , ]
        num_draws <- posterior::ndraws(draws)
      }
    }
    stopifnot(
      num_draws >= target,
      num_draws <= num_draws + posterior::nchains(draws)
    )
  }
  draws
}

# `write_stan_json` is copied directly from cmdstanr::write_stan_json, with the
# addition of the `skip_conversion` argument.
# Stan can't handle the variables that would trigger an error below, but
# sometimes what I want here is the output, as if I had called
# `cmdstanr::write_stan_json()`, but with the option to add additional variable
# names in `skip_conversion`.
write_stan_json <- function (data, file, skip_conversion=charater(0)) {
  if (!is.character(file) || !nzchar(file)) {
    stop("The supplied filename is invalid!", call. = FALSE)
  }
  for (var_name in setdiff(names(data), skip_conversion)) {
    var <- data[[var_name]]
    if (!(
      is.numeric(var) ||
      is.factor(var) ||
      is.logical(var) ||
      is.data.frame(var) ||
      is.list(var))
    ) {
      stop("Variable '", var_name, "' is of invalid type.", call. = FALSE)
    }
    if (is.logical(var)) {
      mode(var) <- "integer"
    }
    else if (is.data.frame(var)) {
      var <- data.matrix(var)
    }
    else if (is.list(var)) {
      var <- list_to_array(var, var_name)
    }
    data[[var_name]] <- var
  }
  jsonlite::write_json(
    data,
    path = file,
    auto_unbox = TRUE,
    factor = "integer",
    digits = NA,
    pretty = TRUE
  )
}

#' Set the adapt_delta Stan parameter based on the attempt count
#'
#' @param attempt_count The number of tries this rule has been tried
#'   (starts from 1).
#' @return A real scalar in interval (0, 1). See Stan docs for
#'   `adapt_delta`. Value weakly increases with `attempt_count`.
get_adapt_delta <- function(attempt_count) {
  stopifnot(
    rlang::is_scalar_integerish(attempt_count),
    noNAs(attempt_count),
    attempt_count >= 1
  )
  if (attempt_count <= 4) {
    adapt_delta <- c(0.9, 0.95, 0.99, 0.995)[attempt_count]
  } else {
    adapt_delta <- 0.99
  }
  stopifnot(!anyNA(adapt_delta))
  if (attempt_count > 1) {
    glue_message("After past failures, using adapt_delta={adapt_delta} for attempt {attempt_count}")
  }
  adapt_delta
}

clean_stan_csvs <- function(model_fit, generate_fit) {
  to_delete <- c(
    model_fit$output_files(include_failed=TRUE),
    generate_fit$output_files(include_failed=TRUE)
  )
  fs::file_delete(to_delete)
}
