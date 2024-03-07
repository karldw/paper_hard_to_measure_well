# Code copied from https://github.com/karldw/kdw.junk/tree/faad2210010660a911f50302b05de6a74e149cf3
# This file is licensed under Mozilla Public License 2.0
# https://github.com/karldw/kdw.junk/blob/faad2210010660a911f50302b05de6a74e149cf3/LICENSE

suppressMessages(
  here::i_am("code/kdw_package_code.r", uuid="b379d3aa-1ee7-4a8a-aec4-ce91856e6c9d")
)
## From winsorize.r:

#' Winsorize one vector
#'
#' @param x Vector to winsorize (numeric)
#' @param trim Fraction of data to trim (if length == 1, same at top and bottom,
#'   if length == 2, trim different fractions)
#' @param na.rm Remove NAs? Fails if NAs are present and this argument is FALSE
#' @return x, but winsorized
#'
#' (this version deals with point masses!)
#' @export
winsorize <- function(x, trim = 0.01, na.rm = FALSE) {
  stopifnot(purrr::is_atomic(x), purrr::is_atomic(trim), length(trim) %in% c(1, 2))
  if (any(trim < 0) || any(trim > 0.5)) {
    stop("trimming must be reasonable")
  }
  if (length(trim) == 1) {
    trim <- c(trim, trim)
  }
  # Things like date can be winsorized, but dplyr::between is unhappy about it.
  # set to a basic class first
  if (typeof(x) %in% c("integer", "double") & !inherits(x, "double") & !inherits(x, "integer")) {
    orig_class <- class(x)
    x <- unclass(x)
    must_reclass <- TRUE
  } else {
    must_reclass <- FALSE
  }
  # quantile approach borrowed from https://www.r-bloggers.com/winsorization/
  # but improved by using type = 1
  lim <- stats::quantile(x, probs = c(trim[1], 1 - trim[2]),
    type = 1, names = FALSE, na.rm = na.rm)
  x[x < lim[1]] <- lim[1]
  x[x > lim[2]] <- lim[2]

  if (must_reclass) {
    class(x) <- orig_class
  }
  x
}


## From reproducible_graphics.r"

#' Save a ggplot plot
#'
#' @param plot The plot created by [ggplot2::ggplot()]
#' @param filename The filename(s) to save (save type depends on extension)
#' @param scale_mult A scale multiplier on the size. Defaults to 1; bigger
#' numbers use a larger canvas. (If length 2, multiply width and height)
#' @param bg The background color, passed to the output device. The default
#' is "transparent". If set to "transparent", the plot will be modified to make
#' the `panel.background`, `plot.background`, `legend.background`, and
#' `legend.box.background` transparent as well. Set it to "white" to retain
#' the normal ggplot behavior.
#' @param device The device to use. Default depends on filename extension. Uses
#' cairo_pdf devices when available. Use "tex" or "tikz" to save with [tikzDevice::tikz()].
#' @param reproducible Logical. Should we try to make the plot reproducible by
#' resetting the embedded timestamp? Defaults to false unless the `SOURCE_DATE_EPOCH`
#' environment variable or `SOURCE_DATE_EPOCH` R option is set. If `reproducible`
#' is `TRUE` and `SOURCE_DATE_EPOCH` isn't set, the timestamp we reset to is
#' 1970-01-01 00:00:00 UTC. Other sources of non-reproducibility aren't handled.
#' Requires system `sed` command.
#'
#' @return The plot (invisibly)
#' @seealso [ggplot2::ggsave()] https://reproducible-builds.org/docs/source-date-epoch/
#'
#' Note that creating reproducible outputs currently depends on the system
#' command `sed`, which isn't installed by default on Windows.
#' This implementation of reproducible outputs is incompatible with multi-file
#' outputs, like tikz.
#'
#' @export
save_plot <- function(plot, filename, scale_mult = 1, bg = "transparent", device=NULL, reproducible=NULL) {
  force(plot)
  if (length(filename) > 1) {
    for (fl in filename) {
      save_plot(
        plot=plot, filename=fl,
        scale_mult=scale_mult, bg=bg, device=device, reproducible=reproducible
      )
    }
    return(invisible(plot))
  }
  stopifnot(dir.exists(dirname(filename)), length(scale_mult) %in% 1:2)
  if (identical(bg, "transparent")) {
    plot <- plot + ggplot2::theme(
      # Borrowed from https://stackoverflow.com/a/41878833
      panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
      plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
      strip.background      = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.key            = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )
  }
  if (is.null(device)) {
    dev_list <- list(
      pdf = function(filename, ...) grDevices::cairo_pdf(filename=filename, ...),
      tex = function(filename, ...) tikzDevice::tikz(filename=filename, ...),
      tikz = function(filename, ...) tikzDevice::tikz(filename=filename, ...)
    )
    device <- dev_list[[tolower(tools::file_ext(filename))]]
  }
  # By default, check if R option or system variable SOURCE_DATE_EPOCH is set
  # and can be parsed as a POSIXct. If so, we'll try to reset the file timestamp
  # to that datetime. If `reproducible` is just `TRUE`, then we'll reset to
  # 1970-01-01 00:00:00 UTC
  if (is.null(reproducible)) {
    reproducible <- !is.null(read_source_date_epoch())
  }
  if (reproducible && (!is_sed_available())) {
    warning("Reproducible graphics (currently) depend on having a system sed command")
    reproducible <- FALSE
  }
  if (reproducible) {
    # Make a tempfile with the same extension as filename so type can be parsed
    ggsave_file <- tempfile(fileext=paste0(".", tools::file_ext(filename)))
    dev_category <- get_dev_category(filename, device)
    if (dev_category == "tikz") {
      rlang::abort(glue::glue("Creating reproducible multi-part files, like {dev_category}, is not supported"))
    }
  } else {
    ggsave_file <- filename
  }
  # Save in the ratio of a beamer slide.
  # This aspect ratio works pretty well for normal latex too
  if (length(scale_mult) == 1) {
    scale_mult <- c(scale_mult, scale_mult)
  }
  ggplot2::ggsave(filename = ggsave_file, plot = plot,
    width = 6.3 * scale_mult[1], height = 3.54 * scale_mult[2], units = "in",
    device = device, bg = bg)
  if (reproducible) {
    reset_datestamp(infile=ggsave_file, outfile=filename, category=dev_category)
  }
  invisible(plot)
}


read_source_date_epoch <- function() {
  source_data_epoch <- getOption("SOURCE_DATE_EPOCH") %||% # First choice
    Sys.getenv("SOURCE_DATE_EPOCH")
  datetime <- as.POSIXct(as.numeric(source_data_epoch), origin="1970-01-01", tz="UTC")
  if (is.na(datetime)) {
    datetime <- NULL
  }
  datetime
}


get_dev_category <- function(filename, device) {
  if (is.function(device)) {
    if (identical(device, grDevices::pdf) || identical(device, grDevices::cairo_pdf)) {
      category <- "pdf"
    } else if (identical(device, grDevices::cairo_ps)) {
      category <- "cairo_ps"
    } else if (identical(device, grDevices::jpeg)) {
      category <- "jpeg"
    } else if (requireNamespace("tikzDevice", quietly=TRUE) && identical(device, tikzDevice::tikz)) {
      category <- "tikz"
    } else {
      category <- "uncorrected"
    }
    if (category != "uncorrected") {
      return(category)
    }
  }
  ext <- tolower(tools::file_ext(filename))
  if (any(c(device, ext) %in% c("pdf", "cairo_pdf"))) {
    category <- "pdf"
  } else if (any(c(device, ext) == "cairo_ps")) {
    category <- "cairo_ps"
  } else if (any(c(device, ext) %in% c("jpg", "jpeg"))) {
    category <- "jpeg"
  } else if (any(c(device, ext) %in% c("tex", "tikz"))) {
    category <- "tikz"
  } else {
    category <- "uncorrected"
  }
  category
}


reset_datestamp <- function(infile, outfile, category) {
  `%||%` <- rlang::`%||%`
  if (category == "uncorrected") {
    file.rename(infile, outfile)
    return(outfile)
  }
  datestring <- switch(category,
    cairo_ps = "%%CreationDate: ",
    jpeg = "%%CreationDate: ",
    pdf = "  /CreationDate ",
    tikz = "% Created by tikzDevice "
  )
  timeformat <- switch(category,
    cairo_ps = "%a %b %d %H:%M:%S %Y",
    jpeg = "%a %b %d %H:%M:%S %Y",
    pdf = "(D:%Y%m%d%H%M%S-00'00)",
    tikz = "%Y-%m-%d %H-%M-%S"
  )
  desired_datetime <- read_source_date_epoch() %||% # First choice
    as.POSIXct(0, origin="1970-01-01", tz="UTC") # third choice
  inregex <- paste0(datestring, ".*")
  outregex <- paste0(
    datestring,
    strftime(desired_datetime, format=timeformat, tz="UTC")
  )
  substitute_text(infile, outfile, inregex, outregex)
}

is_sed_available <- function() {
  suppressWarnings(rc <- system2("sed", "--version", stderr=FALSE, stdout=FALSE))
  can_use_sed <- (rc == 0) && (get_os() != "win")
  can_use_sed
}

substitute_text <- function(infile, outfile, inregex, outregex) {
  # This could definitely be done in R, but it seems like a pain.
  stopifnot(is.character(infile), length(infile) == 1, is.character(outfile),
    length(outfile) == 1, infile != outfile
  )
  if (any(grepl('"', c(inregex, outregex), fixed=TRUE))) {
    stop("This substitution function isn't designed to handle regex that include double quotes. Sorry.")
  }
  if (!is_sed_available()) {
    stop("Text substitution depends on having the system `sed` command available")
  }
  sed_regex <- paste0(
    '"s/',
    gsub("/", "\\/", inregex, fixed=TRUE),
    "/",
    gsub("/", "\\/", outregex, fixed=TRUE),
    '/"'
  )
  # cat(sed_regex)
  # cat("\n")
  # For reasons that aren't clear to me, stdout=outfile doesn't work.
  rc <- system2("sed", c(sed_regex, infile, " > ", outfile))
  if (rc != 0) {
    stop("sed command failed!\n",
      "  Regex:  ", sed_regex, "\n",
      "  Input:  ", infile, "\n",
      "  Output: ", outfile, "\n"
    )
  }
  if ((!file.exists(outfile)) || (file.size(outfile) == 0 && file.size(infile) > 0)) {
    stop("Output file ", outfile, " was not created successfully")
  }
  invisible(outfile)
}








## From convenience.r:

#' Get the current operating system
#'
#' @return "win", "mac", or "linux".
#' Only works for Windows, Mac, and Linux. (Non-Mac unix systems are reported as Linux.)
#' @export
get_os <- function() {
  if (.Platform$OS.type == "windows") {
    "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "mac"
  } else if (.Platform$OS.type == "unix") {
    "linux"  # sorry BSD
  } else {
    stop("Unknown OS")
  }
}

#' Install packages if not already installed
#'
#' @param pkg_list Vector of packages to install
#' @param verbose Narrate progress?
#' @return NULL
#'
#' @export
install_lazy <- function(pkg_list, verbose = TRUE) {
  installed_packages <- utils::installed.packages()[, 1]
  need_to_install <- setdiff(pkg_list, installed_packages)
  already_installed <- pkg_list[pkg_list %in% installed_packages]
  for (pkg in need_to_install) {
    try(utils::install.packages(pkg), silent=TRUE)
  }
  if (verbose) {
    message("Already installed:")
    print(already_installed)
    newly_installed <- need_to_install[need_to_install %in% utils::installed.packages()]
    if (length(newly_installed) > 0) {
      message("Newly installed:")
      print(newly_installed)
    }
  }
  failed_to_install <- setdiff(need_to_install, utils::installed.packages())
  if (length(failed_to_install) > 0) {
    warning("Failed to install these packages:\n  ", paste(failed_to_install))
  }
  invisible()
}


#' Rename columns in a table
#'
#' @param .tbl Table to rename
#' @param .vars A named character vector, where the names are the new column
#'   names and the values are existing column names.
#' @param strict Should the function raise an error if existing column names
#'   can't be found? (Default TRUE)
#' @return The same .tbl, with some renamed columns
#'
#' Note that this function is the same as `colnames()<-` for in-memory
#' data.frames, but also works for remote tbls. It's similar to pandas'
#' `.rename` method.
#'
#' @examples
#' cols_to_rename <- c(cyl2 = "cyl")
#' rename_cols(mtcars, cols_to_rename)
#' @export
rename_cols <- function(.tbl, .vars, strict = TRUE) {
  tbl_names <- colnames(.tbl)
  old_names <- unname(.vars)
  if ((! purrr::is_bare_character(.vars)) ||
    (length(.vars) == 0) ||
    length(names(.vars)) != length(.vars)) {
    stop("Must provide a named character vector of variables to rename. The form should be c(\"new_name\" = \"old_name\")")
  }
  if (anyDuplicated(unname(.vars))) {
    stop("The original names should not be duplicated")
  }
  # Get the index in tbl_names that we're going to rename
  # Will be NA if missing
  rename_idx <- match(old_names, tbl_names)
  if (anyNA(rename_idx)) {
    if (strict) {
      stop("Variables not present to rename:\n  ", paste(old_names[is.na(rename_idx)], collapse = ", "))
    }
    rename_idx <- rename_idx[!is.na(rename_idx)]
  }
  .renaming_fn <- function(x) {
    # new name is stored in the names attribute
    names(.vars[old_names == x])
  }
  out <- dplyr::rename_at(.tbl,
     .vars = dplyr::vars(rename_idx),
     .funs = list(.renaming_fn))
  out
}


#' Get or set memory limits
#'
#' For Linux (or BSD), this function calls [unix::rlimit_as()]
#' For Windows, this function calls [utils::memory.limit()]
#' For Mac OS X, no limiting is available.
#'
#' @param size numeric. Request a new limit, in MiB.
#' @return A vector with the (new) limit, in MiB.
#' @seealso \link[base]{Memory-limits} for other limits.
#' @export
memory_limit <- function(size = NA) {
  os <- get_os()
  if (os == "win") {
    if (is.null(size)) {
      size <- NA
    }
    limit <- utils::memory.limit(size)
  } else if (os == "linux") {
    if (!requireNamespace("unix", quietly=TRUE)) {
      stop("Limiting memory on linux requires the 'unix' package")
    }
    if (is.null(size) || is.na(size)) {
      size <- NULL
    } else {
      size <- size * (1024^2) # rlimit_as expects bytes
    }
    limit <- unix::rlimit_as(size)$cur
  } else {
    warning("Sorry, memory limits on OS X are not supported")
    limit <- NULL
  }
  invisible(limit)
}

#' Make names nicer to work with
#'
#' @param x A character vector of names
#' @return A transformed version of those names
#' @seealso [make.names()] and [tibble::tibble()]'s `.name_repair` argument
#'
#' Resulting names are guaranteed to be unique, and will almost certainly be
#' syntactic.
#'
#' @examples
#' make_better_names(c("Country", "GDP $M", "Coast.Length"))
#' #> [1] "country" "gdp_mn"  "coast_length"
#' # Note that the guarantee to make the names unique can lead to some surprises
#' # for example, "a_and_b" becomes "a_and_b_3" in this case:
#' make_better_names(c("a and b", "a-and-b", "a.and.b", "a_and_b"))
#' #> c("a_and_b", "a_and_b_1", "a_and_b_2", "a_and_b_3")
#' # Here's a way to have a bad time:
#' make_better_names(c("", "x", "X", "x_1"))
#' @export
make_better_names <- function(x) {
  `.` <- NULL # make R CMD CHECK happy
  better_names <- gsub("%", "pct", x, fixed=TRUE) %>%
    gsub("$M", "mn", ., fixed=TRUE) %>%
    gsub("$B", "bn", ., fixed=TRUE) %>%
    gsub("$T", "tn", ., fixed=TRUE) %>%
    make.names() %>%
    tolower() %>%
    gsub(".", "_", ., fixed=TRUE) %>%
    gsub("_+", "_", ., perl=TRUE) %>%
    gsub("^_|_$", "", ., perl=TRUE)
  loop_count <- 0L
  while (anyDuplicated(better_names) != 0) {
    loop_count <- loop_count + 1L
    if (loop_count > 100L) {
      stop("Failed to make names unique!")
    }
    better_names <- gsub(".", "_", make.names(better_names, unique=TRUE), fixed=TRUE)
  }
  better_names
}


## From uniqueness.r:

#' Do the claimed variables identify rows?
#'
#' Just like Stata's isid. For normal tables, this runs faster if data.table is
#' installed.
#'
#' @param df A dataframe to test
#' @param ... Variable names, following [dplyr::select] rules
#' @param notifier A function to report conditions you wouldn't want in an ID variable
#'  (Defaults to [base::warning]. Other reasonable options might be [base::stop] to
#'   escalate issues or [base::force] to not report them.)
#' @return TRUE/FALSE for ID-ness
#' @examples
#' is_id(mtcars, cyl)  # FALSE
#' is_id(Loblolly, Seed) # FALSE
#' is_id(Loblolly, Seed, age) # TRUE
#' vars <- c("Seed", "age")
#' is_id(Loblolly, vars) # TRUE
#' @export
is_id <- function(df, ..., notifier = base::warning) {
  UseMethod("is_id")
}

#' @export
is_id.sf <- function(df, ..., notifier = base::warning) {
  df <- sf::st_drop_geometry(df)
  NextMethod()
}

#' @export
is_id.data.frame <- function(df, ..., notifier = base::warning) {
  df_names <- colnames(df)
  # eval_select checks if columns are missing
  claimed_id_vars <- df_names[tidyselect::eval_select(rlang::expr(c(...)), df)]

  stopifnot(is.character(claimed_id_vars), length(claimed_id_vars) > 0,
            is.function(notifier))
  df_id_cols_only <- dplyr::ungroup(dplyr::select(df, tidyselect::all_of(claimed_id_vars)))
  id_cols_with_na <- purrr::map_lgl(df_id_cols_only, anyNA)
  if (any(id_cols_with_na)) {
    err_msg <- paste("ID variables cannot be NA. Problem variables:",
      paste(colnames(id_cols_with_na)[id_cols_with_na], collapse = ", "), sep = "\n")
    notifier(err_msg)
    return(FALSE)
  }
  total_row_count <- nrow(df_id_cols_only)
  if (total_row_count == 0) {
    notifier("No rows!")
    return(FALSE)
  }

  # Timing considerations:
  # - anyDuplicated from data.table is fastest by far
  # - dplyr::distinct() is good, and faster than `count` when there's one column
  # - dplyr::count() is good, and can be better than `distinct` when there are
  #   multiple columns and the data are unique
  # - base::anyDuplicated is slow
  # Tests on this data:
  # df <- purrr::map_dfr(1:10, ~dplyr::mutate(nycflights13::flights, rep_group = .))
  if (requireNamespace("data.table", quietly=TRUE)) {
    ids_are_unique <- anyDuplicated(data.table::as.data.table(df_id_cols_only)) == 0
  } else {
    ids_are_unique <- nrow(dplyr::distinct(df_id_cols_only)) == total_row_count
  }
  return(ids_are_unique)
}

#' @export
is_id.tbl_lazy <- function(df, ..., notifier = base::warning) {
  `.` <- NULL # make R CMD CHECK happy.
  df <- dplyr::collect(df, 0L)
  df_names <- colnames(df)

  # eval_select checks if columns are missing
  claimed_id_vars <- df_names[tidyselect::eval_select(rlang::expr(c(...)), df)]

  stopifnot(is.character(claimed_id_vars), length(claimed_id_vars) > 0,
            is.function(notifier))

  any_vars <- dplyr::any_vars  # does nothing except satisfy R CMD CHECK
  df_id_cols_only <- dplyr::ungroup(dplyr::select(df, tidyselect::all_of(claimed_id_vars)))
  df_nas <- dplyr::filter_all(df_id_cols_only, any_vars(is.na(.)))
  df_nas_nrow <- nrow(utils::head(df_nas, 1), force = TRUE)
  # If the df_nas table has any rows, at least one ID variable contains NAs
  if (df_nas_nrow > 0) {
    notifier("ID variables cannot be NA.")
    return(FALSE)
  }

  total_row_count <- nrow(df_id_cols_only, force = TRUE)
  if (total_row_count == 0) {
    notifier("No rows!")
    return(FALSE)
  }
  nrow_distinct <- nrow(dplyr::distinct(df_id_cols_only), force = TRUE)

  nrow_distinct == total_row_count
}


#' Raise an error if the claimed variables don't uniquely identify rows.
#'
#' Calls [is_id] (with warnings as errors), then returns the original data if [is_id]
#' returns `TRUE`.
#'
#' @param df A dataframe to test
#' @param ... Passed to [is_id]
#' @return Original `df`
#' @export
ensure_id_vars <- function(df, ...) {
  if (! isTRUE(is_id(df, ..., notifier = base::stop))) {
    claimed_id_vars <- tidyselect::vars_select(colnames(df), ...)
    stop("Variables don't uniquely identify rows: ",
         paste(claimed_id_vars, collapse = ", "))
  }
  df
}



## From cite_packages.r



#' Cite currently attached pacakges
#'
#' @param outfile Write results to a file. The default, NULL, doesn't write a file.
#' @param include Extra packages to include, beyond those that have been
#'  attached (ones that you've called `library()` for.)
#'  The default, "base", includes a citation for R itself.
#' @return A list of citations, invisibly
#'
#' Note that some packages provide multiple citations. This function will pick only the
#' first one they provided.
#'
#' @examples
#' cite_attached_packages()
#' @export
cite_attached_packages <- function(outfile = NULL, include = "base") {
  prev_enc <- getOption("encoding")
  on.exit(options(encoding = prev_enc), add = TRUE)
  options(encoding = "utf8")
  cite_one <- function(pkg, ...) {
    cite <- utils::toBibtex(utils::citation(pkg, ...))
    if (identical(pkg, "base")) {
      pkg <- "r_core"
    }

    # Handle cases where authors list multiple citations (only pick the first)
    cite_lines <- paste(cite)
    empty_lines <- which(cite_lines == "")
    bracket_lines <- which(cite_lines == "}")  # lone bracket on a line
    if (length(empty_lines) + 1 != length(bracket_lines)) {
      stop("Multi-citation detection failed for ", pkg)
    }
    if (length(empty_lines) > 0) {
      if (! all((empty_lines - 1) %in% bracket_lines)) {
        stop("Multi-citation detection failed for ", pkg)
      }
      # Take only the first cite by taking all the lines up until the first
      # empty line (which we just confirmed is after a closing bracket)
      cite <- cite[seq_len(bracket_lines[1])]
      class(cite) <- "Bibtex"
    }
    cite[1] <- gsub("(@[a-zA-Z]+{),", paste0("\\1", pkg, ","), cite[1], perl = TRUE)

    # Also adjust to the biblatex "@software" and set a version
    # cite[1] <- gsub("@Manual{", "@software{", cite[1], fixed = TRUE)
    version_row <- which(names(cite) == "note")
    version_text <- gsub("R package version ", "",
                         cite[version_row], fixed = TRUE)
    cite[version_row] <- version_text
    # names(cite)[version_row] <- "version"

    return(cite)
  }

  # Don't auto-include RevoUtilsMath (can still be manually included with the include param)
  # If no packages have been attached, otherPkgs will be NULL, which is fine.
  pkg_names <- setdiff(names(utils::sessionInfo()$otherPkgs), "RevoUtilsMath")

  include <- include[! is.na(include)]
  pkg_names <- c(pkg_names, include)

  cites <- lapply(pkg_names, cite_one)

  if (!is.null(outfile)) {
    cites_text <- paste(lapply(cites, paste, collapse = "\n"), collapse = "\n\n")
    writeLines(cites_text, outfile, useBytes = TRUE)
  }
  invisible(cites)
}



## From spatial.r


#' Union intersecting geometries in x
#
#' @param x An simple features object (`sf`, `sfg`, or `sfg`)
#' @return The union of intersecting geometries
#' Specifically, union all of the intersecting groups in x
#' https://gis.stackexchange.com/a/323067
#' @examples
#' library(sf)
#' sq = function(pt, sz = 1) {
#'  st_polygon(list(rbind(c(pt - sz), c(pt[1] + sz, pt[2] - sz), c(pt + sz),
#'    c(pt[1] - sz, pt[2] + sz), c(pt - sz))))
#' }
#' x = st_sf(box = 1:6, st_sfc(sq(c(4.2,4.2)), sq(c(0,0)), sq(c(1, -0.8)),
#'   sq(c(0.5, 1.7)), sq(c(3,3)), sq(c(-3, -3))))
#' st_union_intersection(x)
#' @seealso [sf::st_union()]
#' @export
st_union_intersection <- function(x) {
  if (!requireNamespace("igraph", quietly=TRUE)) {
    stop("st_union_intersection requires the igraph package")
  }
  if (!requireNamespace("sf", quietly=TRUE)) {
    stop("st_union_intersection requires the sf package")
  }
  if (length(x) == 0) {
    return(x)
  }
  # Helper function:
  union_by_index <- function(idx, geom) {
    sf::st_union(geom[idx])
  }
  # doesn't matter if x is already a geometry
  x <- sf::st_geometry(x)
  unioned_list <- sf::st_intersects(x) %>%
    igraph::graph_from_adj_list() %>%
    igraph::components() %>%
    igraph::groups() %>%
    lapply(union_by_index, geom=x)
  do.call(c, unioned_list)
}
