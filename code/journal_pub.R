suppressMessages(
here::i_am("code/journal_pub.R", uuid="45c011d0-4728-4906-a78e-f5fdc65dc2f9")
)

OUTDIR <- here::here("output/single_dir_submission_files")
TEX_FRAGMENTS_DIR <- here::here("output/tex_fragments")

make_empty_outdir <- function() {
  if (fs::dir_exists(OUTDIR)) {
    fs::dir_delete(OUTDIR)
  }
  fs::dir_create(OUTDIR, recurse=FALSE)
}

check_length <- function(x) {
  too_long <- x[nchar(x) > 64]
  if (length(too_long) > 0) {
    print(too_long)
    stop("Paths printed above are too long.")
  }
  invisible(x)
}


copy_all <- function() {
  # note that paper.tex and the files in individual_figures.tex are handled separately.
  # software_cites_generated.tex"
  copy_from_output <- c(
    "data_cites.bib",
    "define_acronyms.tex",
    "methane_measurement_refs.bib",
    "paper_appendix_shared_preamble.tex",
    "refs.bib"
  )

  copy_from_graphics <- c(
    "ORCIDiD_iconvector.pdf",
    "figure01_epa_emiss_2015_ch4.pdf",
    "figure02_leak_comparison_fraction_with_detect_remote_only.pdf",
    "figure02_leak_comparison_sizes_remote_only.pdf",
    "figure03_model_cost_vs_q_dwl_plot.pdf",
    "figure04_outcomes_fee_frac=1pct.pdf",
    "figure05_outcomes_dwl_emis_frac=1pct.pdf"
  )

  check_length(c(copy_from_output, copy_from_graphics))

  fs::file_copy(
    here::here("output", copy_from_output),
    OUTDIR
  )

  fs::file_copy(
    here::here("graphics", copy_from_graphics),
    OUTDIR
  )

}


TEX_FRAGMENTS_REGEX <- paste0(
  # Match either \import{tex_fragments/}{ or \input{../tex_fragments/
  r"[(\\import\{tex_fragments/\}\{|\\input\{\.\./tex_fragments/)]",
  # Match the tex fragment filename, then a closing `}`
  r"[(.+\.tex)\}]"
)


replace_one_tex_fragment <- function(x) {
  message(x)
  match <- stringr::str_match(x, TEX_FRAGMENTS_REGEX)
  val <- readLines(fs::path(TEX_FRAGMENTS_DIR, match[1, 3]), warn=FALSE) |>
    trimws(whitespace = "[ \t\r\n%]", which="right") |>
    paste0(collapse="\n")
  val
}


replace_each_tex_fragment <- function(input_file, output_file) {
  lines <- readLines(input_file)
  # First, run replace_one_tex_fragment to replace each matching line
  # (non-matches are unchanged), then also replace
  #\import{individual_figures/}{...} with \input{...} because things are now
  # all in the same directory.
  outlines <- stringr::str_replace_all(lines, TEX_FRAGMENTS_REGEX, replace_one_tex_fragment) |>
    stringr::str_replace_all(
      r"[\\import\{individual_figures/\}\{(.+\.tex)\}]",
      r"[\\input\{\1\}]"
    )
  # note that this will write with "\n" line endings, regardless of input
  writeLines(outlines, output_file)
}


replace_all_tex_fragments <- function() {
  # Note: I know this is like a dramatically worse version of knitr or sweave. Sorry. This was the easiest way to meet requirements I didn't initially design for.
  individual_figures <- c(
    "figure04_expected_fee_1pct.tex",
    "figure05_policy_outcomes_frac=1pct.tex",
    "table01_well_summary_stats.tex"
  )
  check_length(individual_figures)
  message("Replacements:")
  replace_each_tex_fragment(here::here("output", "paper.tex"), fs::path(OUTDIR, "paper.tex"))
  purrr::walk2(
    here::here("output/individual_figures", individual_figures),
    fs::path(OUTDIR, individual_figures),
    replace_each_tex_fragment
  )
}


main <- function() {
  make_empty_outdir()
  copy_all()
  replace_all_tex_fragments()
  invisible(NULL)
}


main()