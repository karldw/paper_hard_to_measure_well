suppressMessages(
  here::i_am("code/generate_code_cites.R", uuid="a15d2768-12b9-457a-ba23-e160281710d5")
)
source(here::here("code/shared_functions.r"))




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


extract_entry_names_from_bib <- function(filename) {
  lines <- readLines(filename)
  # Find things that look like @article{___, or @manual{___,
  matches <- stringr::str_match(lines, r"[^@[A-Za-z]+\{(.+),]")[, 2]
  matches[!is.na(matches)]
}

make_code_cite_bib <- function(outfile) {
  # This should be done before copy_all() so the bib_to_include file gets copied.
  code_bib_files <- here::here(c("output/software_cites_r.bib", "output/software_cites_python.bib"))

  bib_names <- lapply(code_bib_files, extract_entry_names_from_bib) |>
    unlist()
  stopifnot(anyDuplicated(bib_names) == 0)
  cite_lines <- paste0("\\nocite{", bib_names, "}")
  writeLines(cite_lines, outfile)
  invisible(outfile)
}


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
  "ggpattern",
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


make_code_cite_bib(snakemake@output[[2]])
