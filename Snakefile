import os
from pathlib import Path
# Can set these in code rather than the command line:
workflow.use_conda = True
workflow.use_singularity = True
# --cleanenv because singularity passes existing environment variables through
# to the container. That's a problem for some path variables and some build args
workflow.singularity_args = "--cleanenv"
# Retry failed jobs once. This can avoid erroring out for transient problems,
# at the cost of trying things twice if there's a genuine issue.
workflow.restart_times = 1

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
STAN_FITS = Path("data/generated/stan_fits/")

if not Path("data/drillinginfo").is_dir():
    raise FileNotFoundError("Missing contents in data/")

# numbers (basic, and decimals require 0)
number_regex = r"\d+\.?\d*"
RE_AUDIT_tauT = f'({"|".join([tau + "-" + T for tau in ("low", "med", "high") for T in ("1week", "3month")])}|{number_regex})'

# Use this container image for all computations -- even on other operating
# systems, computations will run in a linux virtual machine.
container: "docker://continuumio/miniconda3:4.10.3"
# https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#ad-hoc-combination-of-conda-package-management-with-containers
# Everything except latex and placeholder (`touch`) rules get run in the container.


rule all:
    input:
        "output/paper.pdf",
        # "output/all_figures.pdf",


rule all_models:
    input:
        # Force all currently working models:
        STAN_FITS / "01_twopart_lognormal/model_fit.rds",
        STAN_FITS / "01_twopart_lognormal-bootstrap/model_fit.rds",
        STAN_FITS / "02_twopart_lognormal_alpha/model_fit.rds",
        STAN_FITS / "02_twopart_lognormal_alpha-bootstrap/model_fit.rds",
        STAN_FITS / "03_twopart_lognormal_meas_err/model_fit.rds",
        STAN_FITS / "03_twopart_lognormal_meas_err-bootstrap/model_fit.rds",
        STAN_FITS / "04_twopart_lognormal_meas_err_alpha/model_fit.rds",
        STAN_FITS / "05_twopart_normal_qr/model_fit.rds",
        STAN_FITS / "07_twopart_normal_qr_meas_err/model_fit.rds",
        # STAN_FITS / "06_twopart_normal_qr_alpha_model.rds",  100% max treedepth
        STAN_FITS / "08_twopart_lognormal_heterog_alpha/model_fit.rds",

        # And their generated quantities:
        STAN_FITS / "01_twopart_lognormal/leak_size_draw.parquet",
        STAN_FITS / "01_twopart_lognormal-bootstrap/leak_size_draw.parquet",
        STAN_FITS / "02_twopart_lognormal_alpha/leak_size_draw.parquet",
        STAN_FITS / "02_twopart_lognormal_alpha-bootstrap/leak_size_draw.parquet",
        STAN_FITS / "03_twopart_lognormal_meas_err/leak_size_draw.parquet",
        # STAN_FITS / "03_twopart_lognormal_meas_err-bootstrap/leak_size_draw.parquet",  # Failing
        STAN_FITS / "04_twopart_lognormal_meas_err_alpha_generate/leak_size_draw.parquet",
        # STAN_FITS / "04_twopart_lognormal_meas_err_alpha_generate-bootstrap/leak_size_draw.parquet", # Failing
        STAN_FITS / "08_twopart_lognormal_heterog_alpha/leak_size_draw.parquet",
        STAN_FITS / "08_twopart_lognormal_heterog_alpha-bootstrap/leak_size_draw.parquet",


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


rule well_pad_mapping:
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
        # See the well_pad_mapping rule.
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
        # lyon_etal_2016 = "data/generated/methane_measures/lyon_etal_2016.parquet",
        cleaned_matched_obs = "data/generated/methane_measures/matched_wells_all.rds",
        aviris_match_fraction_dropped = TEX_FRAGMENTS / "aviris_match_fraction_dropped.tex",
    threads: 1
    resources:
        mem_mb = 4000 # 4 GB
    conda: "code/envs/r_scripts.yml"
    script:
        "code/match_jpl_measurements.R"


rule plot_leak_distributions:
    input:
        cleaned_matched_obs = "data/generated/methane_measures/matched_wells_all.rds",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/plot_leak_distributions.R",
    output:
        "graphics/leak_comparison_sizes.pdf",
        "graphics/leak_comparison_sizes_remote_only.pdf",
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
    container: None
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
        matched_leakage = "data/generated/methane_measures/matched_wells_all.rds",
        data_prep_script = "code/distribution_model_data_prep.r",
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
    input:
        "renv.lock",
        "code/envs/r_scripts.yml",
        "code/setup_r_library.R",
    output:
        touch(SNAKEMAKE_FLAGS / "setup_r_library"),
    conda: "code/envs/r_scripts.yml"
    # Sets MAKEFLAGS for parallel compilation
    threads: 4
    script:
        "code/setup_r_library.R"


rule compile_stan:
    input:
        stan_file = "code/stan_models/{filename}.stan",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
    output:
        # this would be {filename}.exe on windows
        "code/stan_models/{filename}"
    conda: "code/envs/r_scripts.yml"
    resources:
        mem_mb = 2000
    threads: 1
    script:
        "code/compile_stan.R"


rule stan_fit_model:
    input:
        measurements = "data/generated/methane_measures/matched_wells_all.rds",
        stan_file = "code/stan_models/{model_name}_model.stan",
        data_prep_script = "code/distribution_model_data_prep.r",
        stan_version = "code/version_info_cmdstan.txt",
        compiled =  "code/stan_models/{model_name}_model",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
    output:
        # Note: the run_stan.R code pulls this file by *index*
        STAN_FITS / "{model_name}{prior_only}{bootstrap}{time_period}/model_fit.rds",
    # Note: in non-bootstrap stan, we'll only use 4 threads (not doing within-chain threading)
    # As always, snakemake will scale this down if fewer cores are available.
    threads: 100
    wildcard_constraints:
        # This is a little messy, but here we say anything other than a "-" is
        # part of the model name, which prevents the model_name regex from being
        # greedy and taking over the "-prior" text.
        # To add other param, just put them after the prior,
        model_name = r"[^\-]+",
        prior_only = "(-prior)?", # only accept "-prior" or ""
        bootstrap = "(-bootstrap)?",
        time_period = r"(-period_\d+_hours)?",
    conda: "code/envs/r_scripts.yml"
    script:
        "code/run_stan.R"


rule stan_generate:
    input:
        measurements = "data/generated/methane_measures/matched_wells_all.rds",
        stan_file = "code/stan_models/{model_name}_generate.stan",
        data_prep_script = "code/distribution_model_data_prep.r",
        stan_version = "code/version_info_cmdstan.txt",
        existing_fit = STAN_FITS / "{model_name}{prior_only}{bootstrap}{time_period}/model_fit.rds",
        compiled = "code/stan_models/{model_name}_generate",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
    output:
        # Note: The code will use these names (except stan_data_json) to look
        # for variables in the stan output. Therefore, all must be named.
        # The cost_coef models will additionally write cost_param_alpha.parquet
        # and cost_param_A.parquet (special cases in run_stan.R)
        leak_size_draw   = STAN_FITS / "{model_name}{prior_only}{bootstrap}{time_period}/leak_size_draw.parquet",
        leak_size_expect = STAN_FITS / "{model_name}{prior_only}{bootstrap}{time_period}/leak_size_expect.parquet",
        prob_leak        = STAN_FITS / "{model_name}{prior_only}{bootstrap}{time_period}/prob_leak.parquet",
        stan_data_json   = STAN_FITS / "{model_name}{prior_only}{bootstrap}{time_period}/stan_data.json",
        prob_size_above_threshold = STAN_FITS / "{model_name}{prior_only}{bootstrap}{time_period}/prob_size_above_threshold.parquet",
    threads: 100
    wildcard_constraints:
        # Match a string of at least one character, matching any character except
        # "-"
        model_name = r"[^\-]+",
        prior_only = "(-prior)?", # only accept "-prior" or ""
        bootstrap = "(-bootstrap)?",
        time_period = r"(-period_\d+_hours)?",
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


# rule render_dot:
#     """Use graphviz's `dot` command to render a file
#     Depends on having `dot` on the PATH. Works for any file type dot supports,
#     including pdf, png, and svg.
#     """
#     input:
#         "graphics/{basename}.dot"
#     output:
#         "graphics/{basename}.{ext}"
#     wildcard_constraints:
#         ext = "(pdf|svg|svgz|png|ps)"
#     container: None
#     shell:
#         "dot {input} -T {wildcards.ext} -o {output}"


rule tex_to_docx:
    input:
        tex = "output/{filename}.tex",
        pdf = "output/{filename}.pdf" # force the dependencies
    output:
        "output/{filename}.docx"
    shell:
        "pandoc {input.tex} -o {output}"


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
    }
    all_models = {
        "01_twopart_lognormal",
        "02_twopart_lognormal_alpha",
        "03_twopart_lognormal_meas_err",
        "04_twopart_lognormal_meas_err_alpha",
        "05_twopart_normal_qr",
        "06_twopart_normal_qr_alpha",
        "07_twopart_normal_qr_meas_err",
        "08_twopart_lognormal_heterog_alpha",
    }
    model_name = wildcards.model_name
    if model_name not in all_models:
        raise ValueError(f"Unknown model name: {model_name}")

    if model_name in cost_coef_models:
        time_period = wildcards.time_period
    else:
        time_period = ""
    fit_dir = STAN_FITS / f"{model_name}{wildcards.prior_only}{wildcards.bootstrap}{time_period}"

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


rule outcomes_analysis:
    # See comments in outcomes_analysis_inputs
    input: unpack(outcomes_analysis_inputs)
    output:
        results_summary = STAN_FITS /
            "{model_name}{prior_only}{bootstrap}{time_period}" /
            "audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.parquet",
    wildcard_constraints:
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
    threads: 1
    resources:
        mem_mb = 3000 # 3 GB
    conda: "code/envs/py_scripts.yml"
    log: "scratch/logs/outcomes_analysis/{model_name}{prior_only}{bootstrap}{time_period}_audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT={audit_tauT}.log",
    script:
        "code/outcomes_analysis.py"

## Begin a bunch of functions that take output from outcomes_analysis and make outputs

rule policy_shadow_price_plot:
    input:
        outcome_summaries = expand(
            STAN_FITS / "08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
            "audit_outcome_summary_rule={audit_rule}_frac={audit_amount}_tauT=med-3month.parquet",
            audit_rule=["uniform", "target_x", "target_e_high"],
            audit_amount=["1pct", "10pct"],
        ),
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/policy_shadow_price_plot.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r",
    output:
        shadow_price_plot = "graphics/audit_shadow_price_plot.pdf",
    threads: 1
    resources:
        mem_mb = 2000,
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_shadow_price_plot.R"


def linspace(start, stop, num=50):
    # np.linspace without numpy (assuming endpoint=True)
    return [start + x * (stop - start) / (num - 1) for x in range(num)]

rule policy_aggregate_abatement_plot:
    input:
        outcome_summaries = expand(
            STAN_FITS / "08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
            "audit_outcome_summary_rule=remote_low_frac=0pct_tauT={audit_tauT}.parquet",
            # Here we're using SOCIAL_COST_METHANE_PER_KG=2 and time_H=8760
            audit_tauT = linspace(start=0, stop=2 * 8760 * 1.5, num=100)
        ),
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/policy_aggregate_abatement_plot.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r",
    output:
        aggregate_abatement_plot = "graphics/aggregate_abatement_curve.pdf",
    threads: 1
    resources:
        mem_mb = 2000,
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
    if wildcards.table_type == "expected_fee":
        audit_policies = ["uniform", "target_x", "target_e_low", "target_e_high"]

    elif wildcards.table_type in ("outcome_dwl", "outcome_emis"):
        audit_policies = ["uniform", "target_x", "target_e_low", "target_e_high", "remote_low", "remote_high"]
    else:
        raise ValueError("Unknown table_type: ", wildcards.table_type)

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
            input_file_part.append(
                f"audit_outcome_summary_rule={pol}_frac={amt}_tauT={tauT}.parquet"
            )
    model_dir = STAN_FITS / "08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours"
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
            TEX_FRAGMENTS / "{table_type}_frac={audit_amount}_tauT={audit_tau}{audit_T}.tex",
        ]
    wildcard_constraints:
        table_type = "(expected_fee|outcome_dwl|outcome_emis)",
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

rule policy_table_expected_fee_1pct_1week:
    input:
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=high-1week.tex",
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=med-1week.tex",
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=low-1week.tex",
    # Note: the individual_figures file of the same name is hand-written, but
    # it's nice to give snakemake a way to decide whether to demand these files.
    container: None
    output:
        touch(SNAKEMAKE_FLAGS / "expected_fee_1pct_1week")

rule policy_table_expected_fee_1pct_3month:
    input:
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=high-3month.tex",
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=med-3month.tex",
        TEX_FRAGMENTS / "expected_fee_frac=1pct_tauT=low-3month.tex",
    # Note: the individual_figures file of the same name is hand-written, but
    # it's nice to give snakemake a way to decide whether to demand these files.
    container: None
    output:
        touch(SNAKEMAKE_FLAGS / "expected_fee_1pct_3month")

# Not currently used:
rule policy_table_expected_fee_frac_optimal_100usd_T_1week:
    input:
        TEX_FRAGMENTS / "expected_fee_frac=optimal-100usd_tauT=med-1week.tex",
    container: None
    output:
        touch(SNAKEMAKE_FLAGS / "expected_fee_optim_100")
rule policy_table_expected_fee_frac_optimal_600usd_T_1week:
    input:
        TEX_FRAGMENTS / "expected_fee_frac=optimal-600usd_tauT=med-1week.tex",
    container: None
    output:
        touch(SNAKEMAKE_FLAGS / "expected_fee_optim_600")

rule policy_table_policy_outcomes_frac_1pct_T_1week:
    # This isn't a real rule that gets run. We just use it as a cleaner input list.
    input:
        TEX_FRAGMENTS / "outcome_emis_frac=1pct_tauT=all-1week.tex",
        TEX_FRAGMENTS / "outcome_dwl_frac=1pct_tauT=all-1week.tex",
    container: None
    output:
        touch(SNAKEMAKE_FLAGS / "policy_outcomes_frac=1pct_T=1week")

rule policy_table_policy_outcomes_frac_1pct_T_3month:
    input:
        TEX_FRAGMENTS / "outcome_emis_frac=1pct_tauT=all-3month.tex",
        TEX_FRAGMENTS / "outcome_dwl_frac=1pct_tauT=all-3month.tex",
    container: None
    output:
        touch(SNAKEMAKE_FLAGS / "policy_outcomes_frac=1pct_T=3month")

rule policy_snippets_generator:
    input:
        outcome_summaries = (STAN_FITS /
            "08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
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
        mem_mb = 2500,
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_snippets_generator.R"

rule policy_required_snippets:
    # This isn't a real rule that gets run. We just use it as a cleaner input list.
    input:
        [TEX_FRAGMENTS / f for f in [

        "OUTCOME=fee_p10_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=fee_p10_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",

        "OUTCOME=fee_p90_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=fee_p90_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",

        "OUTCOME=fee_mean_RULE=uniform_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=fee_mean_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",

        "OUTCOME=emission_tonCO2e_RULE=none_FRAC=0pct_tauT=high-3month.tex",

        "OUTCOME=emission_reduce_pct_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",

        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=emission_reduce_tonCO2e_RULE=uniform_FRAC=1pct_tauT=med-3month.tex",

        "OUTCOME=net_private_cost_per_mcf_pct_price_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex",

        "OUTCOME=welfare_gain_pct_RULE=uniform_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=med-3month.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex",
        "OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex",
        ]]
    container: None
    output:
        touch(SNAKEMAKE_FLAGS / "policy_required_snippets")

rule policy_audit_gains_rel_plot:
    input:
        outcome_summaries = expand(
            STAN_FITS / "08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours" /
            "audit_outcome_summary_rule={audit_rule}_frac={{audit_amount}}_tauT={audit_tauT}.parquet",
            audit_rule=["uniform", "target_x", "target_e_high"],
            audit_tauT=["low-1week", "med-1week", "med-3month",  "high-1week", "high-3month"],
        ),
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/policy_audit_gains_rel_plot.R",
        constants = "code/constants.json",
        policy_output_helper_functions = "code/policy_output_helper_functions.r",
    output:
        "graphics/audit_gains_rel_plot_{out_measure}_frac={audit_amount}{simple}.pdf",

    wildcard_constraints:
        out_measure = "(dwl|emis|fee_per_kg_p90|fee_per_kg_mean|fee_per_kg_med|fee_per_kg_p10)",
        audit_amount="(0pct|1pct|10pct|optimal-100usd|optimal-600usd)",
        simple = "(_simple)?"
    threads: 1
    resources:
        mem_mb = 3000,
    conda: "code/envs/r_scripts.yml"
    script:
        "code/policy_audit_gains_rel_plot.R"

## End of policy outputs


rule output_model_fits:
    input:
        distribution_fits = [
            STAN_FITS / "01_twopart_lognormal-bootstrap/model_fit.rds",
            STAN_FITS / "03_twopart_lognormal_meas_err/model_fit.rds", # bootstrapped measurement error isn't working yet
            # STAN_FITS / "05_twopart_normal_qr/model_fit.rds",
            STAN_FITS / "02_twopart_lognormal_alpha-bootstrap-period_8760_hours/model_fit.rds",
            STAN_FITS / "08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours/model_fit.rds",
        ],
        measurements = "data/generated/methane_measures/matched_wells_all.rds",
        stan_data_json = [
            STAN_FITS / "01_twopart_lognormal-bootstrap/stan_data.json",
            STAN_FITS / "03_twopart_lognormal_meas_err/stan_data.json",
            # STAN_FITS / "05_twopart_normal_qr/stan_data.json",
            STAN_FITS / "02_twopart_lognormal_alpha-bootstrap-period_8760_hours/stan_data.json",
            STAN_FITS / "08_twopart_lognormal_heterog_alpha-bootstrap-period_8760_hours/stan_data.json",
        ],
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/output_model_fits.R",
    output:
        model_prob_leak_plot = "graphics/model_prob_leak_plot.pdf",
        model_cost_vs_q_plot = "graphics/model_cost_vs_q_plot.pdf",
        model_cost_vs_q_dwl_plot = "graphics/model_cost_vs_q_dwl_plot.pdf",
        model_coef_obs_leak = TEX_FRAGMENTS / "model_parameters_obs_leak.tex",
        model_coef_leak_size = TEX_FRAGMENTS / "model_parameters_leak_size.tex",
        model_coef_footer_obs_leak  = TEX_FRAGMENTS / "model_parameters_footer_obs_leak.tex",
        model_coef_footer_leak_size = TEX_FRAGMENTS / "model_parameters_footer_leak_size.tex",
        # model_cost_alpha_by_leak_size_bin_plot = "graphics/model_cost_alpha_by_leak_size_bin_plot.pdf",
        # model_leak_size_by_leak_size_bin_plot  = "graphics/model_leak_size_by_leak_size_bin_plot.pdf",
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
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/plot_natural_gas_prices.R",
    output:
        nat_gas_price_timeseries = "graphics/nat_gas_price_timeseries.pdf",
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
        matched_leakage = "data/generated/methane_measures/matched_wells_all.rds",
        r_lib = SNAKEMAKE_FLAGS / "setup_r_library",
        script = "code/plot_well_operator_stats.R",
    output:
        well_pads_per_operator_cdf = "graphics/well_pads_per_operator_cdf.pdf",
        wells_per_pad_histogram = "graphics/wells_per_pad_histogram.pdf"
    threads: 1
    resources:
        mem_mb = 2000
    conda: "code/envs/r_scripts.yml"
    script:
        "code/plot_well_operator_stats.R"


rule zip_for_replication:
    input:
        "code/zip_for_replication.py",
    output:
        public_zip = "data/replication/replication_public.zip",
        drillinginfo_zip = "data/replication/replication_drillinginfo.zip",
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
            "graphics/audit_shadow_price_plot.pdf",
            "graphics/aggregate_abatement_curve.pdf",
            "graphics/audit_gains_rel_plot_dwl_frac=1pct.pdf",
            "graphics/audit_gains_rel_plot_emis_frac=1pct.pdf",
            "graphics/audit_gains_rel_plot_emis_frac=1pct_simple.pdf",
            rules.policy_table_expected_fee_1pct_1week.output,
            rules.policy_table_expected_fee_1pct_3month.output,
            rules.policy_table_policy_outcomes_frac_1pct_T_1week.output,
            rules.policy_table_policy_outcomes_frac_1pct_T_3month.output,

            rules.output_model_fits.output,
            rules.plot_epa_emissions.output,
            rules.plot_well_operator_stats.output,
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
    container: None
    shell:
        LATEX_CMD


rule paper:
    input:
        tex = "output/paper.tex",
        includes = rules.figures.input.includes,
        snippets = rules.policy_required_snippets.output,
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
    container: None
    output:
        "output/paper.pdf",
    shell:
        LATEX_CMD
