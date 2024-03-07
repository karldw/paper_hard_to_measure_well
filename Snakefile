import os
from pathlib import Path

# Set some variables that define workflow execution in code rather than the
# command line: (note that the syntax has changed over time - need a leading
# underscore now; see https://github.com/snakemake/snakemake/issues/2434)
workflow._use_conda = True
workflow._conda_frontend = "mamba"
# Note: use_singularity = False is the default, but making that explicit here
# because using Singularity leads to weird segfaults with CmdStan, depending on
# the compliation details.
workflow._use_singularity = False
# --cleanenv because singularity passes existing environment variables through
# to the container. That's a problem for some path variables and some build args
# workflow.singularity_args = "--cleanenv"

# The snakemake rules below already include the code as input
workflow._rerun_triggers = frozenset({"mtime", "input", "param"})

# Note on threads:
# The number listed is the max snakemake will schedule. Snakemake will scale
# things down if fewer are available. Snakemake will also not actually enforce
# the thread count; that's up the the called program.

# There's probably a clever way to write this snakemake rule, but this works.
LATEX_CMD = "latexmk -norc -pdflua -cd -bibtex-cond1 -silent {input.tex}"
TEX_FRAGMENTS = Path("output/tex_fragments")
INDIVIDUAL_FIGS = Path("output/individual_figures")
SNAKEMAKE_FLAGS = Path("scratch/snakemake_flags")
SNAKEMAKE_FLAGS.mkdir(parents=True, exist_ok=True)
# Stan fits are things like saved coefficients and big parquet files of MCMC draws
# (this is different than the raw CSVs Stan writes to scratch/stan_output)
STAN_FITS = Path("scratch/stan_fits/")
STAN_FITS.mkdir(parents=True, exist_ok=True)
POLICY_OUTCOMES = Path("data/generated/policy_outcomes/")
POLICY_OUTCOMES.mkdir(parents=True, exist_ok=True)


# Retry failed Stan jobs up to twice with higher adapt_delta. This is slower,
# but avoids transient and/or some fit errors, at the cost of trying things twice if there's an irresolvable issue.
# This is like snakemake's --retries, but is only applied to the Stan rules,
# since errors in other code won't be fixed by retrying.
STAN_RETRY_ATTEMPTS = 3

if not Path("data/drillinginfo").is_dir():
    raise FileNotFoundError("Missing contents in data/")

# numbers (basic, and decimals require 0)
number_regex = r"\d+\.?\d*"
RE_AUDIT_tauT = f'({"|".join([tau + "-" + T for tau in ("low", "med", "high") for T in ("1week", "3month")])}|{number_regex})'

# If using a container, a reasonable one might be:
# container: "docker://condaforge/mambaforge:22.9.0-3"

# Suggestion for running (though for some reason, with script this doesn't always use the conda environments? YMMV)
# script --quiet --return --log-out "run_$(date --iso-8601).log" --append --command "ionice -c3 nice -n19 snakemake --cores 8"
# less --RAW-CONTROL-CHARS --chop-long-lines "$(ls -1 --sort=time *.log | head --lines 1)"

# Remove R paths that would interfere.
import os
import shutil
os.environ.pop("RENV_PATHS_ROOT", None)
os.environ.pop("R_LIBS", None)
os.environ.pop("R_LIBS_USER", None)
os.environ.pop("R_LIBS_SITE", None)

def toolbox_available_but_unused():
    """Check if we're in a setup where the computer is set up to use toolbox,
    but we're not using it.
    https://github.com/containers/toolbox
    There are other containers that this doesn't test for, but this was the test
    that was relevant during development.
    """
    which_toolbox = shutil.which("toolbox")
    try:
        toolbox_path = os.environ["TOOLBOX_PATH"]
    except KeyError:
        toolbox_path = ""
    res = which_toolbox is not None and toolbox_path == ""
    return res

if toolbox_available_but_unused():
    raise OSError("\nDid you forget to run Snakemake inside toolbox?\n")


rule all:
    input:
        "output/paper.pdf",
        # "output/all_figures.pdf",


rule generate_di_parquet_files:
    input:
        # Note: We have data on other states, but no need to process it.
        monthly = expand(
            "data/drillinginfo/entity_production_monthly/{area}_Producing_Entity_Monthly_Production.CSV.gz",
            area = ("CA_Kern", "CA_notKern", "CO", "NM")
        ),
        headers = expand(
            "data/drillinginfo/entity_production_headers/{area}_Production_Headers.CSV.gz",
            area = ("CA_Kern", "CA_notKern", "CO", "NM")
        ),
    output:
        monthly = directory("data/generated/production/monthly_production"),
        headers = directory("data/generated/production/well_headers"),
    conda: "code/envs/py_scripts.yml"
    threads: 6
    script:
        "code/convert_well_info_to_parquet.py"


rule group_into_well_pads:
    # The years used are parsed from the output filename. Different years in
    # year_range will result in different well pad definitions.
    input:
        headers = rules.generate_di_parquet_files.output.headers,
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/group_into_well_pads.R",
    output:
        "data/generated/production/well_pad_crosswalk_{year_range}.parquet"
    threads: 4
    resources:
        mem_mb = 10_000 # 10 GB
    conda: "code/envs/r_scripts.yml"
    script:
        "code/group_into_well_pads.R"


rule match_jpl_measurements:
    input:
        # NOTE: you can use different well pad defs by changing the years here.
        # See the group_into_well_pads rule.
        well_pad_crosswalk = "data/generated/production/well_pad_crosswalk_1970-2018.parquet",
        prod = rules.generate_di_parquet_files.output.monthly,
        headers = rules.generate_di_parquet_files.output.headers,
        alvarez_2018 = "data/studies/alvarez_etal_2018/aar7204_Database_S1.xlsx",
        omara_2018 = "data/studies/omara_etal_2018/Omara_etal_SI_tables.csv",
        duren_2019_plumes = "data/studies/duren_etal_2019/Plume_list_20191031.csv",
        duren_2019_sites = "data/studies/duren_etal_2019/41586_2019_1720_MOESM3_ESM.xlsx",
        lyon_etal_2016_measures = "data/studies/lyon_etal_2016/es6b00705_si_004.xlsx",
        lyon_etal_2016_locations = "data/studies/lyon_etal_2016/es6b00705_si_005.xlsx",
        frankenberg_etal_2016_sources = "data/studies/frankenberg_etal_2016/AVNG_sources_all2.xlsx",
        frankenberg_etal_2016_measures = "data/studies/frankenberg_etal_2016/FourCorners_AV_NG_detections_Werner.xlsx",
        zavala_araiza_2018 = "data/studies/zavala-araiza_etal_2018/elementa-6-284-s1.xlsx",
        nat_gas_prices = "data/generated/nat_gas_prices_by_basin.parquet",
        # Include the script so edits are picked up.
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/match_jpl_measurements.R",
    output:
        # plot_obs_count_jpl = "graphics/observation_count_jpl_flights.pdf",
        # plot_jpl_flights_qq = "graphics/jpl_flights_qq_plot.pdf",
        lyon_etal_2016 = "data/generated/methane_measures/lyon_etal_2016.parquet",
        ground_studies = "data/generated/methane_measures/ground_studies.parquet",
        cleaned_matched_obs = "data/generated/methane_measures/matched_wells_all.parquet",
        aviris_match_fraction_dropped = TEX_FRAGMENTS / "aviris_match_fraction_dropped.tex",
    threads: 1
    resources:
        mem_mb = 4000 # 4 GB
    conda: "code/envs/r_scripts.yml"
    script:
        "code/match_jpl_measurements.R"


rule group_into_dbscan_clusters:
    """
    Assign cluster IDs using dbscan. This is only used for robustness, so do it
    outside the main path of match_jpl_measurements
    """
    input:
        cleaned_matched_obs = "data/generated/methane_measures/matched_wells_all.parquet",
        script = "code/group_into_dbscan_clusters.py",
    output:
        matched_wells_with_dbscan = "data/generated/methane_measures/matched_wells_with_dbscan.parquet",
        map_dbscan_clusters = "graphics/map_dbscan_clusters.png",
    threads: 1
    conda: "code/envs/py_scripts.yml"
    script:
        "code/group_into_dbscan_clusters.py"


rule plot_leak_distributions:
    input:
        cleaned_matched_obs = "data/generated/methane_measures/matched_wells_all.parquet",
        ground_studies = "data/generated/methane_measures/ground_studies.parquet",
        lyon_etal_2016 = "data/generated/methane_measures/lyon_etal_2016.parquet",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/plot_leak_distributions.R",
    output:
        "graphics/leak_comparison_sizes.pdf",
        "graphics/leak_comparison_sizes_remote_only.pdf",
        "graphics/leak_comparison_sizes_extra_titles.pdf",
        "graphics/leak_comparison_fraction_with_detections.pdf",
        "graphics/leak_comparison_fraction_with_detections_remote_only.pdf",
    threads: 1
    resources:
        mem_mb = 7000 # 7 GB
    conda: "code/envs/r_scripts.yml"
    script:
        "code/plot_leak_distributions.R"


rule convert_inkscape:
    input:
        "graphics/drawing_{draw_name}.svg"
    output:
        pdf = "graphics/drawing_{draw_name}.pdf",
        pdf_tex = "graphics/drawing_{draw_name}.pdf_tex"
    shell:
        # TODO: these inkscape options are deprecated. Fix.
        "inkscape {input} --export-filename={output.pdf} --export-latex --without-gui --export-dpi=300 --export-area-drawing --export-pdf-version='1.5' --export-background-opacity=0.0"


rule generate_code_cites_r:
    # Python code cites are manual.
    # (There are also some manual tweaks to output/software_cites_r.bib to fix
    # typos in the package citations.)
    input:
        "code/generate_code_cites.R",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
    output:
        "output/software_cites_r.bib",
    conda: "code/envs/r_scripts.yml"
    script:
        "code/generate_code_cites.R"


rule predict_lyon:
    # Note: this rule supports the in-text claim that leaks are hard to predict,
    # but we don't use the numberical output, and this rule isn't part of a
    # regular run.
    input:
        "data/generated/methane_measures/lyon_etal_2016.parquet",
    conda:  "code/envs/predict_lyon_lightgbm.yml"
    threads: 6
    script:
        # "code/predict_lyon_with_tpot.py"
        "code/predict_lyon_with_lightgbm.py"


rule summary_stats_tables:
    input:
        headers = rules.generate_di_parquet_files.output.headers,
        prod = rules.generate_di_parquet_files.output.monthly,
        well_pad_crosswalk = "data/generated/production/well_pad_crosswalk_1970-2018.parquet",
        matched_leakage = "data/generated/methane_measures/matched_wells_all.parquet",
        ground_studies = "data/generated/methane_measures/ground_studies.parquet",
        lyon_etal_2016 = "data/generated/methane_measures/lyon_etal_2016.parquet",
        data_prep_script = "code/model_data_prep.r",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/summary_stats.R",
    output:
        # Order matters for both of these (we'll index by position later)
        well_summary_stats = [
            TEX_FRAGMENTS / "well_summary_stats_aviris_matched.tex",
            TEX_FRAGMENTS / "well_summary_stats_lyon.tex",
            TEX_FRAGMENTS / "well_summary_stats_all_2018.tex",
        ],
        covariate_balance_by_leak = TEX_FRAGMENTS / "well_covariate_balance.tex",
        covariate_balance_by_flyover = TEX_FRAGMENTS / "well_covariate_balance_by_flyover.tex",
        # TODO: num_rows_all_2018.tex is number of well pads in the well pad
        # crosswalk, which we've cut down to only the states we're using.
        # Either expand that or rename the file
        obs_counts = [
            TEX_FRAGMENTS / "num_rows_aviris_matched.tex",
            TEX_FRAGMENTS / "num_rows_aviris_with_leak.tex",
            TEX_FRAGMENTS / "num_rows_lyon.tex",
            TEX_FRAGMENTS / "num_rows_all_2018.tex",
        ],
        percent_wells_with_gas = TEX_FRAGMENTS / "summ_percent_wells_nonzero_gas_in_us.tex",
    conda: "code/envs/r_scripts.yml"
    script:
        "code/summary_stats.R"


rule setup_r_library:
    # Read description in setup_r_library.R
    input:
        "code/envs/r_scripts.yml",
        "code/setup_r_library.R",
    output:
        touch(SNAKEMAKE_FLAGS / "setup_r_library"),
    conda: "code/envs/r_scripts.yml"
    # Sets MAKEFLAGS for parallel compilation
    threads: 4
    script:
        "code/setup_r_library.R"


rule setup_cmdstan:
    input:
        "code/version_info_cmdstan.txt",
        "code/setup_cmdstan.R",
        SNAKEMAKE_FLAGS / "setup_r_library",
    output:
        touch(SNAKEMAKE_FLAGS / "setup_cmdstan"),
    conda: "code/envs/r_scripts.yml"
    threads: 4
    script:
        "code/setup_cmdstan.R"


rule compile_stan:
    input:
        stan_file = "code/stan_models/{model_name}_{model_or_generate}.stan",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        cmdstan = SNAKEMAKE_FLAGS / "setup_cmdstan",
        script = "code/compile_stan.R",
    output:
        # this would be {filename}.exe on windows
        "code/stan_models/{model_name}_{model_or_generate}"
    conda: "code/envs/r_scripts.yml"
    resources:
        mem_mb = 4000
    wildcard_constraints:
        model_or_generate = "(model|generate)"
    threads: 1
    script:
        "code/compile_stan.R"


rule compile_models:
    input:
        expand(
            "code/stan_models/{model_name}_{model_or_generate}",
            model_name = [
                "08_twopart_lognormal_heterog_alpha",
                "09_twopart_lognormal_heterog_alpha",
            ],
            model_or_generate = ["model", "generate"],
        )


def stan_max_threads(wildcards):
    """How many threads can stan possibly use?"""
    if wildcards.bootstrap == "-bootstrap":
        threads = 100
    else:
        threads = 4
    return threads


rule stan_fit_model:
    input:
        measurements = "data/generated/methane_measures/matched_wells_with_dbscan.parquet",
        stan_file = "code/stan_models/{model_name}_model.stan",
        data_prep_script = "code/model_data_prep.r",
        compiled =  "code/stan_models/{model_name}_model",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        cmdstan = SNAKEMAKE_FLAGS / "setup_cmdstan",
    output:
        # Note: the run_stan.R code pulls this file by *index*
        STAN_FITS / "{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}/model_fit.rds",
    # Note: in non-bootstrap stan, we'll only use 4 threads (not doing within-chain threading)
    # As always, snakemake will scale this down if fewer cores are available.
    threads: stan_max_threads
    wildcard_constraints:
        # robustness_spec is either "main_spec/" or something like "robustness1/"
        robustness_spec = r"(main_spec|robustness[0-9]+)/",
        # This is a little messy, but here we say anything other than a "-",
        # "/" or "." is part of the model name, which prevents the model_name
        # regex from being greedy and taking over the "-prior" text.
        # To add other param, just put them after the prior,
        model_name = r"[^\-/\.]+",
        prior_only = "(-prior)?", # only accept "-prior" or ""
        bootstrap = "(-bootstrap)?",
        time_period = r"(-period_\d+_hours)?",
    retries: STAN_RETRY_ATTEMPTS
    resources:
        attempt_count = lambda wildcards, attempt: attempt
    log: "scratch/logs/cmdstan/stan_fit_model_{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}.log",
    conda: "code/envs/r_scripts.yml"
    script:
        "code/run_stan.R"


rule stan_generate:
    input:
        measurements = "data/generated/methane_measures/matched_wells_with_dbscan.parquet",
        stan_file = "code/stan_models/{model_name}_generate.stan",
        data_prep_script = "code/model_data_prep.r",
        existing_fit = STAN_FITS / "{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}/model_fit.rds",
        compiled = "code/stan_models/{model_name}_generate",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        cmdstan = SNAKEMAKE_FLAGS / "setup_cmdstan",
    output:
        # Note: The code will use these names (except stan_data_json) to look
        # for variables in the stan output. Therefore, all must be named.
        # The cost_coef models will additionally write cost_param_alpha.parquet
        # and cost_param_A.parquet (special cases in run_stan.R)
        leak_size_draw   = STAN_FITS / "{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}/leak_size_draw.parquet",
        leak_size_expect = STAN_FITS / "{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}/leak_size_expect.parquet",
        prob_leak        = STAN_FITS / "{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}/prob_leak.parquet",
        stan_data_json   = STAN_FITS / "{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}/stan_data.json",
        prob_size_above_threshold = STAN_FITS / "{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}/prob_size_above_threshold.parquet",
    threads: stan_max_threads
    wildcard_constraints:
        robustness_spec = r"(main_spec|robustness[0-9]+)/",
        # Match a string of at least one character, matching any characters except
        # "-", "/", and "."
        model_name = r"[^\-/\.]+",
        prior_only = "(-prior)?", # only accept "-prior" or ""
        bootstrap = "(-bootstrap)?",
        time_period = r"(-period_\d+_hours)?",
    retries: STAN_RETRY_ATTEMPTS
    log: "scratch/logs/cmdstan/stan_generate_{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}.log",
    conda: "code/envs/r_scripts.yml"
    script:
        "code/run_stan.R"


rule plot_epa_emissions:
    """Plot EPA's emissions numbers"""
    input:
        epa_data = "data/epa/EPA_GHG_report_2020_table_data.zip",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/plot_epa_emissions.R"
    output:
        ch4_2015 = "graphics/epa_emiss_2015_ch4.pdf",
        ghg_2015 = "graphics/epa_emiss_2015_ghg.pdf",
    conda: "code/envs/r_scripts.yml"
    script:
        "code/plot_epa_emissions.R"


rule tex_to_docx:
    input:
        tex = "output/{filename}.tex",
        pdf = "output/{filename}.pdf" # force the dependencies
    output:
        "output/{filename}.docx"
    shell:
        "pandoc {input.tex} -o {output}"


rule pdf_to_grayscale:
    input:
        pdf = "{filename}.pdf",
    output:
        "{filename}_grayscale.pdf",
    shell:
        """
        gs -o "{filename}_grayscale.pdf -sDEVICE=pdfwrite -sColorConversionStrategy=Gray -dProcessColorModel=/DeviceGray -dCompatibilityLevel=1.4 -dNOPAUSE -dBATCH {filename}"
        """


def outcomes_analysis_inputs(wildcards):
    """
    Some of the models, those that directly estimate cost parameters, will vary
    with the `time_period` wildcard. The other models don't depend on `time_period`
    at all. I'd like to avoid running the models when unnecessary, so this
    function tells Snakemake that for only models  in `cost_coef_models`, the
    fit input files depend on the value of `time_period`.

    In all cases, the *output* of the `outcomes_analysis` rule depends on the
    time_period.

    Snakemake docs:
    https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-input-functions
    """
    cost_coef_models = {
        "02_twopart_lognormal_alpha",
        "04_twopart_lognormal_meas_err_alpha",
        "06_twopart_normal_qr_alpha",
        "08_twopart_lognormal_heterog_alpha",
        "09_twopart_lognormal_heterog_alpha",
    }
    all_models = {
        "01_twopart_lognormal",
        "02_twopart_lognormal_alpha",
        "03_twopart_lognormal_meas_err",
        "04_twopart_lognormal_meas_err_alpha",
        "05_twopart_normal_qr",
        "06_twopart_normal_qr_alpha",
        "07_twopart_normal_qr_meas_err",
        "09_twopart_lognormal_heterog_alpha",
    }
    model_name = wildcards.model_name
    if model_name not in all_models:
        raise ValueError(f"Unknown model name: {model_name}")

    if model_name in cost_coef_models:
        time_period = wildcards.time_period
    else:
        time_period = ""
    fit_dir = STAN_FITS / f"{wildcards.robustness_spec}{model_name}{wildcards.prior_only}{wildcards.bootstrap}{time_period}"

    input_files = {
        "leak_size_draw":   fit_dir / "leak_size_draw.parquet",
        "leak_size_expect": fit_dir / "leak_size_expect.parquet",
        "prob_leak":        fit_dir / "prob_leak.parquet",
        "prob_size_above_threshold": fit_dir / "prob_size_above_threshold.parquet",
        "stan_data_json":   fit_dir / "stan_data.json",
        "script":           "code/outcomes_analysis.py",
        "helper":           "code/outcomes_analysis_helpers.py",
        "constants":        "code/constants.json",
    }
    # Repack the dictionary because snakemake requires strings, not `Path`s
    # input_files = {k: str(v) for k, v in input_files.items()}
    return input_files


def outcomes_analysis_max_threads(wildcards):
    """How many threads can outcomes_analysis.py reasonably use?"""
    if wildcards.audit_rule in {"target_x", "target_e_low", "target_e_high"}:
        threads = 20
    else:
        threads = 2
    return threads


rule outcomes_analysis:
    # See comments in outcomes_analysis_inputs
    input: unpack(outcomes_analysis_inputs)
    output:
        results_summary = POLICY_OUTCOMES /
            "{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}" /
            "audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.parquet",
    wildcard_constraints:
        robustness_spec = r"(main_spec|robustness[0-9]+)/",
        model_name = r"[^\-/\.]+",
        bootstrap = "(-bootstrap)?",
        prior_only = "(-prior)?", # only accept "-prior" or ""
        time_period = r"(-period_\d+_hours)?",
        audit_rule= "(none|uniform|remote_low|remote_high|target_x|target_e_low|target_e_high)",
        # Note: could also make these numbers, if we wanted more flexibility.
        # The downside is harder to read / less consistent filenames.
        audit_amount="(0pct|1pct|10pct|optimal-100usd|optimal-600usd)",
        # For audit_tauT, accept anything like "low-1week", "med-3month", or "3.2"
        audit_tauT=RE_AUDIT_tauT
    threads: outcomes_analysis_max_threads
    resources:
        # note: this memory number is a rough guess, and not enforced in each thread.
        # (each thread takes something like 100 MB)
        mem_mb = 3000 # 3 GB
    conda: "code/envs/py_scripts.yml"
    log: "scratch/logs/outcomes_analysis/{robustness_spec}{model_name}{prior_only}{bootstrap}{time_period}_audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.log",
    script:
        "code/outcomes_analysis.py"

## Begin a bunch of rules that take output from outcomes_analysis and make outputs

# rule policy_shadow_price_plot:
#     input:
#         outcome_summaries = expand(
#             POLICY_OUTCOMES / "main_spec" /
#             "08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
#             "audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT=med-3month.parquet",
#             audit_rule=["uniform", "target_x", "target_e_high"],
#             audit_amount=["1pct", "10pct"],
#         ),
#         r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
#         script = "code/policy_shadow_price_plot.R",
#         constants = "code/constants.json",
#         policy_output_helper_functions = "code/policy_output_helper_functions.r",
#     output:
#         shadow_price_plot = "graphics/audit_shadow_price_plot.pdf",
#     threads: 1
#     resources:
#         mem_mb = 4000,
#     conda: "code/envs/r_scripts.yml"
#     script:
#         "code/policy_shadow_price_plot.R"


rule policy_robustness_plots:
    input:
        # Note: hard-coding the 1-59 range (59 = 2 ^ num_robustness_vars - 2 - num_robustness_vars + 3
        # where we're omitting specs that are equal to the main spec, or have
        # less than 2 rhs variables, and adding 3 non-combinatorial specs.
        # (and remember range() doesn't include the right endpoint)
        # not using bootstrap version
        outcome_summaries = expand(
            POLICY_OUTCOMES / "robustness{idx}" /
            "09_twopart_lognormal_heterog_alpha-period_8760_hours" /
            "audit_outcome_summary_rule={{audit_rule}}_frac={{audit_amount}}_tauT={{audit_tauT}}.parquet",
            idx = range(1, 60)
        ) + [
            POLICY_OUTCOMES / "main_spec" /
            "09_twopart_lognormal_heterog_alpha-period_8760_hours" /
            "audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.parquet",
        ],
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/policy_robustness_plots.R",
        # constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r",
        robustness_spec_script = "code/model_data_prep_robustness.r",
    output:
        plot_file_dwl = "graphics/robustness_audit_result_dwl_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.pdf",
        plot_file_emiss = "graphics/robustness_audit_result_emission_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.pdf",
    threads: 1
    resources:
        mem_mb = 4000,
    wildcard_constraints:
        audit_rule = "(none|uniform|remote_low|remote_high|target_x|target_e_low|target_e_high)",
        audit_amount = "(0pct|1pct|10pct|optimal-100usd|optimal-600usd)",
        audit_tauT = RE_AUDIT_tauT,
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_robustness_plots.R"


def linspace_int(start, stop, num=50):
    # np.linspace and np.round without numpy (assuming endpoint=True)
    # Rounding the results because the numbers get stored as filenames, and
    # it's messy to have filenames like
    # audit_outcome_summary_rule=remote_low_frac=0pct_tauT=24687.272727272728.parquet"
    return [round(start + x * (stop - start) / (num - 1)) for x in range(num)]


rule policy_aggregate_abatement_plot:
    input:
        outcome_summaries = expand(
            POLICY_OUTCOMES / "main_spec" /
            "09_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT={audit_tauT}.parquet",
            # Here we're using SOCIAL_COST_METHANE_PER_KG=2 and time_H=8760
            audit_tauT = linspace_int(start=0, stop=2 * 8760 * 1.5, num=100)
        ),
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/policy_aggregate_abatement_plot.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r",
    output:
        aggregate_abatement_plot = "graphics/aggregate_abatement_curve.pdf",
    threads: 1
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_aggregate_abatement_plot.R"


def policy_tables_generator_inputs(wildcards):
    """Figure out what input files we need for the tables

    Expected fee:
    - uniform
    - target_x
    - target_e_low
    - target_e_high
    For all of these, do audit_amount and tau_T from the wildcards
    (tau of high, med, and low are requested separately)

    Outcome emissions / DWL:
    - uniform
    - target_x
    - target_e_low
    - target_e_high
    - remote_low (set audit_amount = 0pct)
    - remote_high (set audit_amount = 0pct)
    For all of these, do T from the wildcards, and allow tau to be high, med, or low
    For the first 4, pull audit_amount from the wildcards
    """

    # Which policies go into the table?
    # This isn't set in stone; could add or remove if you want.
    if wildcards.outcome_type == "expected_fee":
        audit_policies = ["uniform", "target_x", "target_e_low", "target_e_high", "remote_low", "remote_high"]

    elif wildcards.outcome_type in ("outcome_dwl_emis"):
        audit_policies = ["uniform", "target_x", "target_e_low", "target_e_high", "remote_low", "remote_high"]
    else:
        raise ValueError("Unknown outcome_type: ", wildcards.outcome_type)

    if wildcards.audit_tau == "all-":
        tau_list = ["low-", "med-", "high-"]
    else:
        # Here audit_tau includes the trailing "-"
        assert wildcards.audit_tau.endswith("-")
        tau_list = [wildcards.audit_tau]
    input_file_part = []
    for pol in audit_policies:
        if pol in ["none", "remote_low", "remote_high"]:
            amt = "0pct"
        else:
            amt = wildcards.audit_amount
        for tau in tau_list:
            tauT = f"{tau}{wildcards.audit_T}"
            input_filename = f"audit_outcome_summary_rule={pol}_frac={amt}_tauT={tauT}.parquet"
            # Avoid bug that was happening with f-strings adding spaces (this
            # was solved by running Snakemake with Python 3.11 instead of 3.12,
            # but unclear why.)
            assert " " not in input_filename
            input_file_part.append(input_filename)

    model_dir = POLICY_OUTCOMES / "main_spec" / "09_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours"
    input_files = {
        "outcome_summaries": [model_dir / f for f in input_file_part],
    }
    return input_files


rule policy_tables_generator:
    input:
        unpack(policy_tables_generator_inputs),
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/policy_tables_generator.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r",
    output:
        tables = [
            TEX_FRAGMENTS / "{outcome_type}_frac={audit_amount}_tauT={audit_tau}{audit_T}.tex",
        ]
    wildcard_constraints:
        outcome_type = "(expected_fee|outcome_dwl_emis)",
        audit_tau = r"(low|med|high|all)\-",
        audit_T = r"\d+(week|month)",
    threads: 1
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_tables_generator.R"

# These are organized by table output. For instance, one table in the text
# is from individual_figures/expected_fee_1pct_1week.tex, which in turn
# includes the three sub-tables for different tau levels.
# We could make this shorter with expand(), but it's this way for clarity.
rule policy_table_expected_fee_1pct:
    # Order doesn't matter
    input:
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=high-1week.tex",
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=med-1week.tex",
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=low-1week.tex",
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=high-3month.tex",
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=med-3month.tex",


rule policy_table_outcomes_dwl_emiss_1pct:
    # This isn't a real rule that gets run. We just use it as a cleaner input list.
    # Order doesn't matter
    input:
        TEX_FRAGMENTS / "outcome_dwl_emis_frac=1pct_tauT=low-1week.tex",
        TEX_FRAGMENTS / "outcome_dwl_emis_frac=1pct_tauT=med-1week.tex",
        TEX_FRAGMENTS / "outcome_dwl_emis_frac=1pct_tauT=high-1week.tex",
        TEX_FRAGMENTS / "outcome_dwl_emis_frac=1pct_tauT=med-3month.tex",
        TEX_FRAGMENTS / "outcome_dwl_emis_frac=1pct_tauT=high-3month.tex",



rule policy_snippets_generator:
    """
    Note that this goes one at a time, which is pretty inefficient, since each
    one gets its own R session, but it works okay.
    """
    input:
        outcome_summaries = (POLICY_OUTCOMES / "main_spec" /
            "09_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
            "audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.parquet",
        ),
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/policy_snippets_generator.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r",
    output:
        # We end up running this multiple times per policy, since we want
        # different snippets. It's slightly slower, but means we don't end up
        # with a ton of extra files.
        snippets = TEX_FRAGMENTS / "OUTCOME={outcome}_RULE={audit_rule}_FRAC={audit_amount}_tauT={audit_tauT}{CI}.tex",
    # log: "scratch/logs/policy_snippets/OUTCOME={outcome}_RULE={audit_rule}_FRAC={audit_amount}_tauT={audit_tauT}.log",
    wildcard_constraints:
        audit_rule = "(none|remote_low|remote_high|uniform|target_x|target_e_low|target_e_high)",
        # Note: currently the optimal-* form isn't supported.
        audit_amount="(0pct|1pct|10pct|optimal-100usd|optimal-600usd)",
        audit_tauT = RE_AUDIT_tauT,
        CI = "(_CI)?",
    threads: 1
    resources:
        mem_mb = 2500, # not enforced; see note in script
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_snippets_generator.R"


# Same as tables, but make figures instead
rule policy_graphs_generator:
    input:
        outcome_summaries = expand(
            POLICY_OUTCOMES / "main_spec" /
            "09_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
            "audit_outcome_summary_rule={audit_rule}_frac={{audit_amount}}_tauT={audit_tauT}.parquet",
            audit_rule=["uniform", "target_x", "target_e_high", "target_e_low"],
            audit_tauT=["low-1week", "low-3month", "med-1week", "med-3month",  "high-1week", "high-3month", "2500", "7500"],
        ) +
        expand(
            POLICY_OUTCOMES / "main_spec" /
            "09_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
            "audit_outcome_summary_rule={audit_rule}_frac=0pct_tauT={audit_tauT}.parquet",
            audit_rule=["remote_low", "remote_high"],
            audit_tauT=["low-1week", "low-3month", "med-1week", "med-3month",  "high-1week", "high-3month", "2500", "7500"],
        ),
        r_lib = "scratch/snakemake_flags/setup_r_library",
        script = "code/policy_graphs_generator.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r"
    output:
        plot_fee = "graphics/outcomes_fee_frac={audit_amount}.pdf",
        plot_dwl_emis_ordinal = "graphics/outcomes_dwl_emis_frac={audit_amount}.pdf",
        plot_dwl_emis_cardinal = "graphics/outcomes_dwl_emis_frac={audit_amount}_cardinal.pdf",
    wildcard_constraints:
        audit_amount="(0pct|1pct|10pct|optimal-100usd|optimal-600usd)",
    threads: 1
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_graphs_generator.R"


rule policy_snippets_generator_robustness:
    input:
        # Note: using the non-bootstrap version here.
        outcome_summaries = (
            POLICY_OUTCOMES /
            "{robustness_spec}09_twopart_lognormal_heterog_alpha-period_8760_hours" /
            "audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.parquet",
        ),
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/policy_snippets_generator.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r",
        robustness_spec_script = "code/model_data_prep_robustness.r",
    output:
        # We end up running this multiple times per policy, since we want
        # different snippets. It's slightly slower, but means we don't end up
        # with a ton of extra files.
        snippets = TEX_FRAGMENTS / "{robustness_spec}OUTCOME={outcome}_RULE={audit_rule}_FRAC={audit_amount}_tauT={audit_tauT}{CI}.tex",
    # log: "scratch/logs/policy_snippets/OUTCOME={outcome}_RULE={audit_rule}_FRAC={audit_amount}_tauT={audit_tauT}.log",
    wildcard_constraints:
        # Note that here robustness_spec does not allow for main_spec.
        # This rule is just for robustness generation.
        robustness_spec = r"(robustness[0-9]+)/",
        audit_rule = "(none|remote_low|remote_high|uniform|target_x|target_e_low|target_e_high)",
        # Note: currently the optimal-* form isn't supported.
        audit_amount="(0pct|1pct|10pct|optimal-100usd|optimal-600usd)",
        audit_tauT = RE_AUDIT_tauT,
        CI = "(_CI)?",
    threads: 1
    resources:
        mem_mb = 2500,
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_snippets_generator.R"


rule policy_required_snippets:
    # This isn't a real rule that gets run. We just use it as a cleaner input list.
    input:
        [TEX_FRAGMENTS / f for f in [

        "OUTCOME=fee_p75_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=fee_p99_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",

        "OUTCOME=fee_mean_RULE=uniform_FRAC=1pct_tauT=med-1week.tex",
        "OUTCOME=fee_mean_RULE=target_x_FRAC=1pct_tauT=med-1week.tex",
        "OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex",
        "OUTCOME=emission_tonCO2e_RULE=none_FRAC=0pct_tauT=high-3month.tex",

        "OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=2298.tex",
        "OUTCOME=emission_reduce_pct_RULE=target_e_high_FRAC=1pct_tauT=2298.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=2298.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_x_FRAC=1pct_tauT=med-1week.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=uniform_FRAC=1pct_tauT=med-1week.tex",
        "OUTCOME=net_private_cost_per_mcf_pct_price_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex",

        "OUTCOME=welfare_gain_pct_RULE=uniform_FRAC=1pct_tauT=med-1week.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=med-1week.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=low-1week.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=low-1week.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=high-3month.tex",
        # Note: these vs_uniform ones also read the "uniform" rule output.
        # We could have Snakemake track and require this, but haven't.
        "OUTCOME=welfare_gain_pct_vs_uniform_RULE=target_x_FRAC=1pct_tauT=low-1week.tex",
        "OUTCOME=welfare_gain_pct_vs_uniform_RULE=target_x_FRAC=1pct_tauT=high-3month.tex",
        "OUTCOME=welfare_pct_of_target_e_low_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=welfare_pct_of_target_e_low_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex",
        ]]


rule policy_required_snippets_robustness:
    # This isn't a real rule that gets run. We just use it as a cleaner input list.
    input:
        [TEX_FRAGMENTS / f"robustness{idx}" / f
        for idx in range(1, 53)
        for f in [
        # "OUTCOME=fee_p10_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=fee_p10_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        # "OUTCOME=fee_p90_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=fee_p90_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        # "OUTCOME=fee_mean_RULE=uniform_FRAC=1pct_tauT=med-3month.tex",
        # "OUTCOME=fee_mean_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        # "OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=2500.tex",
        # "OUTCOME=emission_tonCO2e_RULE=none_FRAC=0pct_tauT=high-3month.tex",
        # "OUTCOME=emission_reduce_pct_RULE=target_e_high_FRAC=1pct_tauT=2500.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        # "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex",
        # "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=uniform_FRAC=1pct_tauT=med-3month.tex",
        # "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=2500.tex",
        # "OUTCOME=net_private_cost_per_mcf_pct_price_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex",
        # "OUTCOME=welfare_gain_pct_RULE=uniform_FRAC=1pct_tauT=med-3month.tex",
        # "OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        # "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex",
        # "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        ]
        ]


rule policy_audit_gains_rel_plot:
    input:
        outcome_summaries = expand(
            POLICY_OUTCOMES / "main_spec" /
            "09_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
            "audit_outcome_summary_rule={audit_rule}_frac={{audit_amount}}_tauT={audit_tauT}.parquet",
            audit_rule=["uniform", "target_x", "target_e_high"],
            audit_tauT=["low-1week", "low-3month", "med-1week", "med-3month",  "high-1week", "high-3month"],
        ),
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/policy_audit_gains_rel_plot.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r",
    output:
        "graphics/audit_gains_rel_plot_{out_measure}_frac={audit_amount}{simple}.pdf",

    wildcard_constraints:
        out_measure = "(dwl|emis|fee_per_kg_mean|fee_per_kg_med|fee_per_kg_p[0-9]{2})",
        audit_amount="(0pct|1pct|10pct|optimal-100usd|optimal-600usd)",
        simple = "(_simple)?"
    threads: 1
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_audit_gains_rel_plot.R"

## End of policy outputs


rule output_model_fits:
    input:
        distribution_fits = [
            STAN_FITS / "main_spec" /
            "09_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours/model_fit.rds",
        ],
        measurements = "data/generated/methane_measures/matched_wells_with_dbscan.parquet",
        stan_data_json = [
            STAN_FITS / "main_spec" /
            "09_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours/stan_data.json",
        ],
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/output_model_fits.R",
    output:
        model_prob_leak_plot = "graphics/model_prob_leak_plot.pdf",
        model_cost_vs_q_plot = "graphics/model_cost_vs_q_plot.pdf",
        model_cost_vs_q_dwl_plot = "graphics/model_cost_vs_q_dwl_plot.pdf",
        model_coef = TEX_FRAGMENTS / "model_parameters.tex",
        model_coef_footer  = TEX_FRAGMENTS / "model_parameters_footer.tex",
        model_cost_alpha_histogram = "graphics/model_cost_alpha_histogram.pdf",
        model_prob_size_above_threshold_histogram = "graphics/model_prob_size_above_threshold_histogram.pdf",
    threads: 1
    resources:
        mem_mb = 13000 # 13 GB
    conda: "code/envs/r_scripts.yml"
    script:
        "code/output_model_fits.R"


rule prep_natural_gas_prices:
    input:
        snl_price_file = "data/snl/snl_gas_price_series.xls",
        price_index_file = "data/fred/CPILFENS.csv",
        basin_mapping = "code/snl_nat_gas_basin_match.csv",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/prep_nat_gas_price_data.R",
    output:
        nat_gas_prices = "data/generated/nat_gas_prices_by_basin.parquet",
    threads: 1
    resources:
        mem_mb = 2000,
    conda: "code/envs/r_scripts.yml"
    script:
        "code/prep_nat_gas_price_data.R"


rule plot_natural_gas_prices:
    # This is a separate script because prep_natural_gas_prices is upstream of
    # some slow code, and I don't want to re-run that every time I change the plot
    input:
        nat_gas_prices = "data/generated/nat_gas_prices_by_basin.parquet",
        matched_leakage = "data/generated/methane_measures/matched_wells_all.parquet",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/plot_natural_gas_prices.R",
    output:
        nat_gas_price_timeseries = "graphics/nat_gas_price_timeseries.pdf",
        nat_gas_price_histogram  = "graphics/nat_gas_price_histogram.pdf",
    threads: 1
    resources:
        mem_mb = 2000,
    conda: "code/envs/r_scripts.yml"
    script:
        "code/plot_natural_gas_prices.R"


rule bea_industry_facts:
    input:
        gdp_file = "data/bea/gdp_by_industry.xls",
        deprec_file = "data/bea/current_cost_depreciation_of_private_fixed_assets_by_industry.xls",
        price_index_file = "data/fred/CPILFENS.csv",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/bea_industry_facts.R",
    output:
        net_value = TEX_FRAGMENTS / "oil_gas_industry_net_value.tex",
    threads: 1
    resources:
        mem_mb = 1800
    conda: "code/envs/r_scripts.yml"
    script:
        "code/bea_industry_facts.R"


rule plot_well_operator_stats:
    input:
        well_pad_crosswalk = "data/generated/production/well_pad_crosswalk_1970-2018.parquet",
        matched_leakage = "data/generated/methane_measures/matched_wells_all.parquet",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/plot_well_operator_stats.R",
    output:
        well_pads_per_operator_cdf = "graphics/well_pads_per_operator_cdf.pdf",
        wells_per_pad_histogram = "graphics/wells_per_pad_histogram.pdf"
    threads: 1
    resources:
        mem_mb = 4000
    conda: "code/envs/r_scripts.yml"
    script:
        "code/plot_well_operator_stats.R"


rule zip_for_replication:
    input:
        "code/zip_for_replication.py",
    output:
        public_zip = "data/replication/replication_public.zip",
        drillinginfo_zip = "data/replication/replication_drillinginfo.zip",
    threads: 1
    conda: None
    script:
        "code/zip_for_replication.py"


# Outputs, using the above as inputs.
rule figures:
    input:
        tex = "output/all_figures.tex",
        includes = [
            "graphics/leak_comparison_sizes.pdf",
            "graphics/game_tree_audit_covariates.tikz",
            "graphics/game_tree_audit_leak.tikz",
            "graphics/game_tree_audit_leak_censor.tikz",
            "graphics/nat_gas_price_timeseries.pdf",
            rules.summary_stats_tables.output,
            # rules.outcomes_analysis_outputs.output,
            # "graphics/audit_shadow_price_plot.pdf",
            "graphics/aggregate_abatement_curve.pdf",
            "graphics/audit_gains_rel_plot_dwl_frac=1pct.pdf",
            "graphics/audit_gains_rel_plot_emis_frac=1pct.pdf",
            rules.policy_table_expected_fee_1pct.input,
            rules.policy_table_outcomes_dwl_emiss_1pct.input,
            # rules.policy_table_policy_outcomes_frac_1pct_T_3month.output,
            rules.output_model_fits.output,
            rules.plot_epa_emissions.output,
            rules.plot_well_operator_stats.output,
            "graphics/robustness_audit_result_dwl_rule=target_e_high_frac=1pct_tauT=med-3month.pdf",
            "graphics/robustness_audit_result_emission_rule=target_e_high_frac=1pct_tauT=med-3month.pdf",
            "graphics/outcomes_dwl_emis_frac=1pct.pdf",
            "graphics/outcomes_fee_frac=1pct.pdf",
        ],
    output:
        "output/all_figures.pdf",
    shell:
        LATEX_CMD


rule slides:
    input:
        tex = "output/2020-12-07_EI_RIP_slides.tex",
        includes = rules.figures.input.includes,
        drawings = [
            "graphics/drawing_cost_curve_dwl_no_policy.pdf",
            "graphics/drawing_cost_curve_dwl_with_policy.pdf",
        ],
    output:
        "output/2020-12-07_EI_RIP_slides.pdf",
    shell:
        LATEX_CMD


rule paper:
    input:
        tex = "output/paper.tex",
        includes = rules.figures.input.includes,
        snippets = rules.policy_required_snippets.input,
        # snippets_robust = rules.policy_required_snippets_robustness.input,
        snippets2 = TEX_FRAGMENTS / "oil_gas_industry_net_value.tex",
        bibs = [
            "output/refs.bib",
            "output/methane_measurement_refs.bib",
            "output/software_cites_r.bib",
            "output/software_cites_python.bib",
        ],
        misc = [
            "output/define_acronyms.tex",
            "graphics/marks_2018_fig4a.png",
            "graphics/natural_gas_leakage_percentages_marks_fig1.png",
            "graphics/frankenberg_etal_2016_fig4-4.jpg",
            "graphics/frankenberg_etal_2016_fig1.jpg",
            "graphics/ORCIDiD_iconvector.pdf",
        ],
    output:
        "output/paper.pdf",
    shell:
        LATEX_CMD
