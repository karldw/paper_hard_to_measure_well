
# Information Matters: Feasible Policies for Reducing Methane Emissions

Authors: Karl Dunkle Werner [![ORCID logo](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-0523-7309) and [Wenfeng Qiu](https://wenfengqiu.com/)


## Overview

The code in this replication package builds an analysis dataset using observed methane leaks and a database of wells. It then estimates a structural model of well leaks and performs counterfactual policy simulations, plotting the results.

The code is orchestrated by [Snakemake](https://snakemake.github.io/), and takes the authors about 18 hours to run.


## Abbreviated steps to replicate

1. Read this README
1. Make sure you have an appropriate computing resources (see below)
1. Unzip the replication files.
1. If the data is saved somewhere outside the project folder, mount or link a copy inside the project folder. (Useful for development only)
1. Install Mamba and Snakemake (see below)
1. Check results into git for later comparison (optional)
1. Run Snakemake
1. Check results (optional)

Putting it all the together:

```sh
# 3. Unzip
mkdir methane_replication # or whatever you want
cd methane_replication
unzip path/to/replication_public.zip -d .
unzip path/to/replication_drillinginfo.zip -d .

# 4. OPTIONAL
# If the data is saved somewhere outside the project folder, mount
# a copy inside the project folder.
# This is only useful if the data are stored somewhere *outside*
# the project folder. You may need to change these paths to fit
# your situation
data_drive="$HOME/Dropbox/data/methane_abatement"
project_dir="$(pwd)"
cd "$project_dir"
mkdir -p "scratch"  # optional - feel free to have this as a link to somewhere else
ln -s "$data_drive" "data"


# 5. Install Mamba and Snakemake
# If mamba is not already installed, follow instructions here:
# https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
mamba env create --name snakemake --file code/envs/install_snakemake.yml
mamba activate snakemake
snakemake --version
# Should show version, not an error

# 6. Check in files to Git for later comparison (optional)
git init
git add .
git commit -m "Initial replication files"
snakemake --delete-all-output # can add --dry-run if you want to see what will happen first
# Note: --delete-all-output does not clean the conda environments, nor the
# additional packages that get added to them.

# 7. Run Snakemake to create all outputs
# (this takes about a day with 12-core CPU)
# snakemake --dry-run to see what will be run
snakemake --cores 12
# Alternate version
# script --quiet --return --log-out "run_$(date --iso-8601).log" --append --command "ionice -c3 nice -n19 snakemake --cores 8"

# 8. Check results (optional)
# Check everything into git, rerun snakemake, and verify results are the same.
git add .
git commit -m "Replication run"
# Results may not be binary-identical, but should not be meaningfully different.
# Stan lists reasons results might differ:
# https://mc-stan.org/docs/reference-manual/reproducibility.html
# If re-run on exactly the same system, results will *mostly* be reproduced, but
# a few of the PDF graphs are not identical run-to-run.
git diff
```


## Setup

### Hardware requirements

Parts of this program require significant amounts of memory and disk space.
Most parts also benefit from having multiple processors available.
(The slow parts parallelize well, so speed should increase almost linearly with processor count.)

The tasks that require significant memory are noted in the Snakemake file (see the `mem_mb` fields).
The highest requirement for any task is 13 GB, though most are far lower.
(These could be overstatements; we haven't tried that hard to find the minimum memory requirements for each operation.)

The programs also use about 120 GB of on-disk storage in `scratch/` in addition to the ~5 GB of input data and ~1 GB output data.

Running the whole thing takes roughly 15-20 hours on a 2022 computer with 12 CPU cores\.

R and Python packages are downloaded from the internet - a reasonably quick internet connection will be helpful. After installation, these packages take about 4 GB of space.


### Operating system

The code was developed on Linux, and should mostly work on other operating systems. However, because most parts of the code have not been tested on Windows or macOS, the code would likely require minor tweaks. To replication on Windows, we recommend running the code inside the [Windows Subsystem for Linux](https://ubuntu.com/desktop/wsl) (WSL).

The code is known to work with Ubuntu Linux 24.04.

Stan compiles some custom programs; the system must be allowed to run these compiled programs. By default, this is not a problem, but some computers with additional security may not allow it.


### Software

This project uses [Snakemake](https://snakemake.readthedocs.io) (v7.32.4), [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (v1.5.8), and Conda (24.7.1) to manage dependencies.
- To get started, first install Mamba (mini or full-sized). Mamba will bring Conda along. Most any version of Mamba should be fine.
    - The code can also use conda, without using mamba. If doing that, edit the `_conda_frontend` line in the `Snakefile` file from `"mamba"` to `"conda"`.
- Then use Mamba to install Snakemake from the file `code/envs/install_snakemake.yml`. This will install the correct version.

`mamba env create --name snakemake --file code/envs/install_snakemake.yml`

**Run all other commands in that activated `snakemake` environment**.

If you close the terminal window, you need to re-run `mamba activate snakemake` before running the rest of the commands. (You can name the environment something other than `snakemake` if you want.)

Note that other versions of Snakemake might work, but Snakemake has a lot of breaking changes from version to version. Versions of Snakemake >=7.8.0 and < 8.0.0 should probably be fine. This code will _not_ work with Snakemake 8 or newer.

Snakemake will automatically install the R and Python packages needed, along with cmdstan. These are specified in `code/envs/r_scripts.yml` and `code/envs/py_scripts.yml`

The code relies on other tools, which are assumed to be available on the system.


- LuaLaTeX, `latexmk`, and `bibtex` as distributed in TeX Live 2022 or 2023 (newer versions should be fine) to compile the `paper.tex` and `online_appendix.tex` PDFs.
- The tex files use some fonts (Libertinus) that are widely distributed, but may not be installed by default. If you get an error about missing the font, [download it](https://github.com/alerque/libertinus/releases) and unzip it to somewhere your system can find fonts (`~/.fonts/Libertinus` works on Linux).
- Ghostscript (10.02.1 used; versions from 9 onward should be fine).
    - If Ghostscript is not available (`gs` on the command line), the code will still run, but the colors in the PDF will be slightly different.
- `sha256sum` (version 9.4 used; any version should be fine). If `sha256sum` is not available, the `sha256_sum_drillinginfo` rule will fail and Snakemake will be unhappy, but only for the `zip_for_replication` rule.
- `curl` was used (as documented in `code/download_prices.sh`) to download a few files from EIA and FRED, but these downloads are included in the replication package.



#### What does Snakemake do?

Snakemake uses rules to generate outputs and manages the code environments to make it all work.
Snakemake uses mamba to install packages. All of this is handled transparently as the rules are run.

It can be useful to run `snakemake --dry-run` to see the planned jobs.

Snakemake keeps track of what needs to run and what doesn't. If something goes wrong midway through, snakemake will see that some outputs are up-to-date and others aren't, and won't re-run the things that don't need it.

Snakemake tries to manage computing resources (RAM and CPU) based on the claims made in the Snakemake rules. Snakemake does not enforce these limits; some of the code does.

In Snakemake, we have code that unsets these environment variables, with the aim of avoiding path confusion:
`RENV_PATHS_ROOT`, `R_LIBS`, `R_LIBS_USER`, `R_LIBS_SITE`, `CMDSTAN`.



### Files and data

The code uses relative paths and R's [here](https://cran.r-project.org/web/packages/here/index.html) package. It should not be necessary to edit the paths in any code.


#### Accessing data

We need make sure the code can access the right files.
There are two ways this can be done, the straightforward way and the way Karl does it.

##### Recommended file access

Straightforward approach:
Unzip the replication files, either interactively or with the commands below.

```sh
mkdir methane_replication # or whatever you want
cd methane_replication
unzip path/to/replication_public.zip -d .
unzip path/to/replication_drillinginfo.zip -d .
```

##### Alternative file access

Less straightforward, arguably better for development

- Store the `data` and/or `scratch` folders somewhere else (e.g. `data` in Dropbox).
- Create symbolic links (on Windows, 'junctions') to point to the `data` and `scratch` folders.
- If you don't like symlinks, [bind mounts](https://unix.stackexchange.com/questions/198590/what-is-a-bind-mount) are fine too.

Example:
```sh
ln -s path/to/dropbox/data/this_project data
ln -s path/to/scratch_drive scratch
```


#### File structure

All files in `output/tex_fragments` (except the README), `data/generated`, `data/replication`, and `scratch/` are auto-generated and can safely be deleted. All other files in `data/` should not be deleted.
`data/` and `scratch/` are ignored by `.gitignore`.


<details>

<summary><b>
Annotated file listing
</b></summary>

```
.
├── code
│   ├── bea_industry_facts.R: Read and summarize the BEA table. Create outputs/tex_fragments/oil_gas_industry_net_value.tex
│   ├── compile_stan.R: Compile one Stan model file (see models in code/stan_models/)
│   ├── constants.json: Constants (some unused) that the analysis assumes, such as methane GWP and the detection threshold.
│   ├── convert_well_info_to_parquet.py: Read the Drillinginfo csv files (data/drillinginfo/entity_production_*) into parquet files (data/generated/production/monthly_production and well_headers/)
│   ├── download_prices.sh: Download price series from FRED, save data/fred/CPILFENS.csv
│   ├── envs
│   │   ├── install_snakemake.yml: Environment to install Snakemake (see instructions above)
│   │   ├── predict_lyon_lightgbm.yml: Environment used for predicting leak in the Lyon et al. dataset. Used in unreported analysis of data/studies/lyon_etal_2016/.
│   │   ├── py_scripts.yml: Environment used by Snakemake to run Python scripts.
│   │   └── r_scripts.yml: Environment used by Snakemake to run R scripts.
│   ├── generate_code_cites.R: Create latex citations for R code. The script creates output/software_cites_r_raw.bib, which is edited by hand into output/software_cites_r.bib and output/software_cites_generated.tex, which is just a bunch of \nocite{} commands to include.
│   ├── group_into_dbscan_clusters.py: Provide cluster groupings using the DBSCAN algorithm. For a robustness check. Creates data/generated/matched_wells_with_dbscan from parquet data/generated/matched_wells_all.parquet
│   ├── group_into_well_pads.R: Group individual wells into well pads by proximity. Creates data/generated/production/well_pad_crosswalk_1970-2018.parquet.
│   ├── journal_pub.R: Copy and modify the various latex files to be journal-compliant. Saves files in output/single_dir_submission_files/, and inlines a bunch of the files in output/tex_fragments/ into the paper.tex itself.
│   ├── match_jpl_measurements.R: Match well pads with measured methane leaks. Creates data/generated/matched_wells_all.parquet, data/generated/methane_measures/lyon_etal_2016.parquet, data/generated/methane_measures/ground_studies.parquet, and output/tex_fragments/aviris_match_fraction_dropped.tex.
│   ├── model_data_prep.r: Specifications for the main models that will be run in Stan (used by run_stan.R).
│   ├── model_data_prep_robustness.r: Specifications for the robustness check models that will be run in Stan (used by run_stan.R).
│   ├── outcomes_analysis_helpers.py: Provides functions used in the policy optimization of outcomes_analysis.py.
│   ├── outcomes_analysis.py: Run policy simulations. Relies on the estimated cost and leak parameters from Stan (in scratch/stan_fits/), and creates all the files in data/generated/policy_outcomes.
│   ├── outcomes_analysis_tests.py: Tests that the outcomes_analysis.py functions work as expected.
│   ├── output_model_fits.R: Provide summaries of the fitted Stan models. Creates graphics/model_cost_vs_q_dwl_plot.pdf (figure 3), graphics/model_cost_alpha_histogram.pdf (figure A9), outputs/tex_fragments/model_parameters.tex (table A3), and outputs/tex_fragments/model_parameters_footer.tex (table A3).
│   ├── plot_epa_emissions.R: Plot EPA emissions inventory and Alvarez et al comparison in graphics/epa_emiss_2015_ch4.pdf (figure 1).
│   ├── plot_leak_distributions.R: Plot distributions of leak sizes and presence measured in various scientific studies (see data/studies/). Creates graphics/leak_comparison_sizes.pdf (figure A5 left), graphics/leak_comparison_sizes_remote_only.pdf (figure 2 left), graphics/leak_comparison_fraction_with_detections.pdf (figure A5 right), and graphics/leak_comparison_fraction_with_detections_remote_only.pdf (figure 2 right).
│   ├── plot_natural_gas_prices.R: Plot natural gas prices from data/snl/. Creates graphics/nat_gas_price_timeseries.pdf (figure A7) and graphics/nat_gas_price_histogram.pdf (figure A8).
│   ├── plot_well_operator_stats.R: Uses well pad mapping and well data and would create graphics/well_pads_per_operator_cdf.pdf and graphics/wells_per_pad_histogram.pdf. These are not currently used, but code may be of interest.
│   ├── policy_aggregate_abatement_plot.R: Plot aggregate abatement at different fee levels for the set of wells in this sample. Not currently used, but retained because it may be of interest. Would create graphics/aggregate_abatement_curve.pdf.
│   ├── policy_audit_gains_rel_plot.R: Not currently used, but retained because it may be of interest for relative effects. (See policy_graphs_generator.R instead).
│   ├── policy_graphs_generator.R: Plot simulated policy outcomes, using files saved by code/outcomes_analysis.py. Creates graphics/outcomes_fee_frac=1pct.pdf (figure 4), graphics/outcomes_dwl_emis_frac=1pct.pdf (figure 5), and graphics/outcomes_dwl_emis_frac={audit_amount}_cardinal.pdf (figure A10).
│   ├── policy_output_helper_functions.r: Shared functions to deal with the policy outputs. Used by code/policy_tables_generator.R, code/policy_snippets_generator.R, code/policy_aggregate_abatement_plot.R, code/policy_audit_gains_rel_plot.R, code/policy_graphs_generator.R, and code/policy_robustness_plots.R.
│   ├── policy_robustness_plots.R: Plot policy outcomes for the robustness specifications. Creates graphics/robustness_audit_result_dwl_rule=target_e_high_frac=1pct_tauT=med-3month.pdf (figure A11 top) and graphics/robustness_audit_result_emission_rule=target_e_high_frac=1pct_tauT=med-3month.pdf (figure A11 bottom).
│   ├── policy_snippets_generator.R: Pull policy outcomes into various tex snippets that are included in the paper. List:  output/tex_fragments/OUTCOME=fee_p75_RULE=target_x_FRAC=1pct_tauT=med-3month.tex, output/tex_fragments/OUTCOME=fee_p99_RULE=target_x_FRAC=1pct_tauT=med-3month.tex, output/tex_fragments/OUTCOME=fee_mean_RULE=uniform_FRAC=1pct_tauT=med-1week.tex, output/tex_fragments/OUTCOME=fee_mean_RULE=target_x_FRAC=1pct_tauT=med-1week.tex, output/tex_fragments/OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex, output/tex_fragments/OUTCOME=emission_tonCO2e_RULE=none_FRAC=0pct_tauT=high-3month.tex, output/tex_fragments/OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=2298.tex, output/tex_fragments/OUTCOME=emission_reduce_pct_RULE=target_e_high_FRAC=1pct_tauT=2298.tex, output/tex_fragments/OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=2298.tex, output/tex_fragments/OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex, output/tex_fragments/OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex, output/tex_fragments/OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex, output/tex_fragments/OUTCOME=emission_reduce_tonCO2e_RULE=target_x_FRAC=1pct_tauT=med-1week.tex, output/tex_fragments/OUTCOME=emission_reduce_tonCO2e_RULE=uniform_FRAC=1pct_tauT=med-1week.tex, output/tex_fragments/OUTCOME=net_private_cost_per_mcf_pct_price_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_RULE=uniform_FRAC=1pct_tauT=med-1week.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=med-1week.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=low-1week.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=low-1week.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=high-3month.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_vs_uniform_RULE=target_x_FRAC=1pct_tauT=low-1week.tex, output/tex_fragments/OUTCOME=welfare_gain_pct_vs_uniform_RULE=target_x_FRAC=1pct_tauT=high-3month.tex, output/tex_fragments/OUTCOME=welfare_pct_of_target_e_low_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex, output/tex_fragments/OUTCOME=welfare_pct_of_target_e_low_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex
│   ├── policy_tables_generator.R: Summarize the simulated policy outcomes in text form. Creates the pieces that make up tables A4 and A5. List: output/tex_fragments/expected_fee_frac=1pct_tauT=high-1week.tex, output/tex_fragments/expected_fee_frac=1pct_tauT=high-3month.tex, output/tex_fragments/expected_fee_frac=1pct_tauT=low-1week.tex, output/tex_fragments/expected_fee_frac=1pct_tauT=low-3month.tex, output/tex_fragments/expected_fee_frac=1pct_tauT=med-1week.tex, output/tex_fragments/expected_fee_frac=1pct_tauT=med-3month.tex, output/tex_fragments/outcome_dwl_emis_frac=1pct_tauT=high-1week.tex, output/tex_fragments/outcome_dwl_emis_frac=1pct_tauT=high-3month.tex, output/tex_fragments/outcome_dwl_emis_frac=1pct_tauT=low-1week.tex, output/tex_fragments/outcome_dwl_emis_frac=1pct_tauT=low-3month.tex, output/tex_fragments/outcome_dwl_emis_frac=1pct_tauT=med-1week.tex, output/tex_fragments/outcome_dwl_emis_frac=1pct_tauT=med-3month.tex
│   ├── prep_nat_gas_price_data.R: Pre-process the SNL natural gas price data (see data/snl/). Creates data/generated/nat_gas_prices_by_basin.parquet.
│   ├── run_stan.R: Run a cmdstan model in either MCMC fitting mode or generate mode. Uses models saved in code/stan_models. Saves output in scratch/stan_output/ (the raw CSVs) and scratch/stan_fits/ (model files model_fit.rds and stan_data.json, generated files cost_param_alpha.parquet, cost_param_A.parquet, leak_size_draw.parquet, leak_size_expect.parquet,  prob_leak.parquet, prob_size_above_threshold.parquet).
│   ├── setup_cmdstan.R: After the R library is set up, re-compile stan. (See notes below.) The Stan files are saved in the conda folder - remove and reinitialize the conda environment if you want a completely fresh install.
│   ├── setup_r_library.R: Install additional R packages into the R conda environment.
│   ├── shared_functions.r: Many shared functions used by other R scripts.
│   ├── snl_nat_gas_basin_match.csv: Mapping between the SNL basin names and my well regions. Used by prep_nat_gas_price_data.R. (Despite being a csv, this was created by hand, and so is tracked as code.)
│   ├── stan_helper_functions.r: Functions used by run_stan.R
│   ├── stan_models
│   │   ├── 09_twopart_lognormal_heterog_alpha_generate.stan: Create the generated outputs. See run_stan.R
│   │   └── 09_twopart_lognormal_heterog_alpha_model.stan: Fit the model. See run_stan.R
│   ├── summary_stats.R: Create summary stats for wells (tables 1 and A2) and some in-text tex_fragments. Uses processed data from drillinginfo and from various studies. Output files: output/tex_fragments/table01_well_summary_stats_aviris_matched.tex, output/tex_fragments/table01_well_summary_stats_lyon.tex, output/tex_fragments/table01_well_summary_stats_all_2018.tex, output/tex_fragments/tableA2_well_covariate_balance.tex, output/tex_fragments/table01_num_rows_aviris_matched.tex, output/tex_fragments/table01_num_rows_aviris_with_leak.tex, output/tex_fragments/table01_num_rows_lyon.tex, output/tex_fragments/table01_num_rows_all_2018.tex, output/tex_fragments/intext_summ_percent_wells_nonzero_gas_in_us.tex.
│   ├── tests.R: Tests for a couple of functions in shared_functions.r.
│   ├── version_info_cmdstan.txt: Version of cmdstan used. This is a txt file because it's checked in both R and Python code. There will be an error if it doesn't match the version specified in envs/r_scripts.yml.
│   └── zip_for_replication.py: Create the replication files.
├── .condarc: Configuration file for conda (sets strict channel priority to avoid version match issues)
├── data: See notes below, in the data availability section.
│   ├── bea
│   │   ├── current_cost_depreciation_of_private_fixed_assets_by_industry.xls
│   │   ├── gdp_by_industry.xls
│   │   └── source_notes.yml
│   ├── drillinginfo
│   │   ├── docs
│   │   │   ├── data_dict_production.pdf
│   │   │   ├── data_dict_wells.pdf
│   │   │   ├── Drillinginfo Well dictionary.pdf
│   │   │   ├── reference_eur_calculation.pdf
│   │   │   └── reference_recompletion.pdf
│   │   ├── entity_production_headers
│   │   │   ├── CA_Kern_Production_Headers.CSV.gz
│   │   │   ├── CA_notKern_Production_Headers.CSV.gz
│   │   │   ├── CO_Production_Headers.CSV.gz
│   │   │   └── NM_Production_Headers.CSV.gz
│   │   ├── entity_production_monthly
│   │   │   ├── CA_Kern_Producing_Entity_Monthly_Production.CSV.gz
│   │   │   ├── CA_notKern_Producing_Entity_Monthly_Production.CSV.gz
│   │   │   ├── CO_Producing_Entity_Monthly_Production.CSV.gz
│   │   │   └── NM_Producing_Entity_Monthly_Production.CSV.gz
│   │   └── SHA256SUM
│   ├── eia
│   │   └── emission_annual.xls
│   ├── epa
│   │   ├── egrid2018_summary_tables.pdf
│   │   └── EPA_GHG_report_2020_table_data.zip
│   ├── fred
│   │   └── CPILFENS.csv
│   ├── snl
│   │   ├── gas_prices_download_notes.txt
│   │   ├── na_gas_methodology.pdf
│   │   ├── snl_gas_price_series.xls
│   │   └── snl_gas_volume_series.xls
│   └── studies
│       ├── alvarez_etal_2018
│       │   ├── aar7204_Database_S1.xlsx
│       │   ├── aar7204_Database_S2.zip
│       │   ├── Alvarez_et_al_2018_Science_methane_supply_chain.pdf
│       │   └── data_source_notes.txt
│       ├── duren_etal_2019
│       │   ├── 41586_2019_1720_MOESM3_ESM.pdf
│       │   ├── 41586_2019_1720_MOESM3_ESM.xlsx
│       │   ├── ANG_L1B_L2_Data_Product_Readme_v02.txt
│       │   ├── AVIRIS-NG Flight Lines - AVIRIS-NG Flight Lines.csv
│       │   ├── AVIRIS-NG Flight Lines - AVIRIS-NG Flight Lines.csv.gz
│       │   ├── CH4_Plume_AVIRIS-NG.pdf
│       │   ├── Duren_etal_2019_Nature_California_superemitters.pdf
│       │   ├── Plume_list_20191031.csv
│       │   └── source_notes.txt
│       ├── frankenberg_etal_2016
│       │   ├── AVNG_sources_all2.xlsx
│       │   ├── FourCorners_AV_NG_detections_Werner.xlsx
│       │   └── source_notes.md
│       ├── lyon_etal_2016
│       │   ├── es6b00705_si_003.xlsx
│       │   ├── es6b00705_si_004.xlsx
│       │   ├── es6b00705_si_005.xlsx
│       │   ├── es6b04303_si_002.xlsx
│       │   ├── Lyon_etal_2016_EST_aerial_survey_of_production_sites.pdf
│       │   └── Lyon_etal_2016_EST_aerial_survey_of_production_sites_supplement.pdf
│       ├── omara_etal_2018
│       │   ├── Omara_etal_SI_tables.csv
│       │   ├── Readme_OmaraEtAl2018_GriddedInventoryKmzFile.docx
│       │   └── US2015_35kmBy35km_GriddedCH4EmissFor2015NGProdSites.kmz
│       └── zavala-araiza_etal_2018
│           └── elementa-6-284-s1.xlsx
├── .gitignore
├── graphics
│   ├── epa_emiss_2015_ch4.pdf: See plot_epa_emissions.R above (figure 1)
│   ├── frankenberg_etal_2016_fig1.jpg: Pulled from paper (figure A6)
│   ├── frankenberg_etal_2016_fig4-4.jpg: Pulled from paper (figure A6)
│   ├── game_tree_audit_covariates.tikz: Created by hand (figure A1)
│   ├── game_tree_audit_leak_censor.tikz: Created by hand (figure A3)
│   ├── game_tree_audit_leak.tikz: Created by hand (figure A2)
│   ├── leak_comparison_fraction_with_detections.pdf: See plot_leak_distributions.R above (figure A5)
│   ├── leak_comparison_fraction_with_detections_remote_only.pdf: See plot_leak_distributions.R above (figure 2)
│   ├── leak_comparison_sizes.pdf: See plot_leak_distributions.R above (figure A5)
│   ├── leak_comparison_sizes_remote_only.pdf: See plot_leak_distributions.R above (figure 2)
│   ├── model_cost_alpha_histogram.pdf: See output_model_fits.R above (figure A9)
│   ├── model_cost_vs_q_dwl_plot.pdf: See output_model_fits.R above (figure 3)
│   ├── nat_gas_price_histogram.pdf: See plot_natural_gas_prices.R above (figure A8)
│   ├── nat_gas_price_timeseries.pdf: See plot_natural_gas_prices.R above (figure A7)
│   ├── natural_gas_leakage_percentages_marks_fig1.png: Pulled from paper (figure A4)
│   ├── ORCIDiD_iconvector.pdf: ORCID logo (converted from svg)
│   ├── outcomes_dwl_emis_frac=1pct_cardinal.pdf: See policy_graphs_generator.R above (figure A10)
│   ├── outcomes_dwl_emis_frac=1pct.pdf: See policy_graphs_generator.R above (figure 4)
│   ├── outcomes_fee_frac=1pct.pdf: See policy_graphs_generator.R above (figure 5)
│   ├── robustness_audit_result_dwl_rule=target_e_high_frac=1pct_tauT=med-3month.pdf: See policy_robustness_plots.R above (figure A11)
│   └── robustness_audit_result_emission_rule=target_e_high_frac=1pct_tauT=med-3month.pdf: See policy_robustness_plots.R above (figure A11)
├── output
│   ├── data_cites.bib: Citations for data.
│   ├── define_acronyms.tex: Acronyms, using latex glossaries package and accsup to make small caps come out as uppercase.
│   ├── individual_figures: These are tex files that wrap the graphics and tables.
│   │   ├── tableA02_covariate_balance_comparison_detect_leak.tex: Table A2
│   │   ├── tableA01_Cusworth_etal_2019_table2.tex: Table A1, using a simplified version of table 2 in Cusworth et al. (2019).
│   │   ├── expected_fee_1pct_app_table.tex: Table A4
│   │   ├── expected_fee_1pct.tex: Figure 4
│   │   ├── frankenberg_measurement_figs.tex: Figure A6
│   │   ├── game_tree_target_e_censor.tex: Figure A3
│   │   ├── game_tree_target_e_no_censor.tex: Figure A2
│   │   ├── game_tree_target_x.tex: Figure A1
│   │   ├── model_parameters.tex: Table A3
│   │   ├── nat_gas_price_timeseries.tex: Figure A7
│   │   ├── policy_outcomes_dwl_emis_app_table_frac=1pct.tex: Table A5
│   │   ├── policy_outcomes_frac=1pct.tex: Figure 5
│   │   ├── README.md
│   │   └── table01_well_summary_stats.tex: Table 1
│   ├── latexmkrc: latexmk configuration file that allows online_appendix.tex to reference the figure numbers in paper.tex
│   ├── methane_measurement_refs.bib: Citations for methane measurements
│   ├── online_appendix.tex: Appendix file
│   ├── paper_appendix_shared_preamble.tex: Latex preamble, shared between paper and appendix
│   ├── paper.tex: Main paper
│   ├── refs.bib: Main citations
│   ├── software_cites_generated.tex: Software citation \nocite commands for R
│   ├── software_cites_python.bib: Python citations
│   ├── software_cites_r.bib: R citations (edited version of software_cites_r_raw.bib)
│   ├── software_cites_r_raw.bib: Auto-generated version of R citations (see code/generate_code_cites.R)
│   └── tex_fragments: These are small tex files that are auto-included in the paper or appendix. The sources are listed above in the code folder.
│       ├── aviris_match_fraction_dropped.tex
│       ├── expected_fee_frac=1pct_tauT=high-1week.tex
│       ├── expected_fee_frac=1pct_tauT=high-3month.tex
│       ├── expected_fee_frac=1pct_tauT=low-1week.tex
│       ├── expected_fee_frac=1pct_tauT=low-3month.tex
│       ├── expected_fee_frac=1pct_tauT=med-1week.tex
│       ├── expected_fee_frac=1pct_tauT=med-3month.tex
│       ├── model_parameters_footer.tex
│       ├── model_parameters.tex
│       ├── num_rows_all_2018.tex
│       ├── num_rows_aviris_matched.tex
│       ├── num_rows_aviris_with_leak.tex
│       ├── num_rows_lyon.tex
│       ├── oil_gas_industry_net_value.tex
│       ├── outcome_dwl_emis_frac=1pct_tauT=high-1week.tex
│       ├── outcome_dwl_emis_frac=1pct_tauT=high-3month.tex
│       ├── outcome_dwl_emis_frac=1pct_tauT=low-1week.tex
│       ├── outcome_dwl_emis_frac=1pct_tauT=low-3month.tex
│       ├── outcome_dwl_emis_frac=1pct_tauT=med-1week.tex
│       ├── outcome_dwl_emis_frac=1pct_tauT=med-3month.tex
│       ├── OUTCOME=emission_reduce_pct_RULE=target_e_high_FRAC=1pct_tauT=2298.tex
│       ├── OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=2298.tex
│       ├── OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex
│       ├── OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex
│       ├── OUTCOME=emission_reduce_tonCO2e_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=emission_reduce_tonCO2e_RULE=target_x_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=emission_reduce_tonCO2e_RULE=uniform_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=emission_tonCO2e_RULE=none_FRAC=0pct_tauT=high-3month.tex
│       ├── OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=2298.tex
│       ├── OUTCOME=fee_mean_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=fee_mean_RULE=target_x_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=fee_mean_RULE=uniform_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=fee_p75_RULE=target_x_FRAC=1pct_tauT=med-3month.tex
│       ├── OUTCOME=fee_p99_RULE=target_x_FRAC=1pct_tauT=med-3month.tex
│       ├── OUTCOME=net_private_cost_per_mcf_pct_price_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex
│       ├── OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-1week.tex
│       ├── OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=high-3month.tex
│       ├── OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=low-1week.tex
│       ├── OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=welfare_gain_pct_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex
│       ├── OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=high-3month.tex
│       ├── OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=low-1week.tex
│       ├── OUTCOME=welfare_gain_pct_RULE=target_x_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=welfare_gain_pct_RULE=uniform_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=welfare_gain_pct_vs_uniform_RULE=target_x_FRAC=1pct_tauT=high-3month.tex
│       ├── OUTCOME=welfare_gain_pct_vs_uniform_RULE=target_x_FRAC=1pct_tauT=low-1week.tex
│       ├── OUTCOME=welfare_pct_of_target_e_low_RULE=target_e_high_FRAC=1pct_tauT=med-1week.tex
│       ├── OUTCOME=welfare_pct_of_target_e_low_RULE=target_e_high_FRAC=1pct_tauT=med-3month.tex
│       ├── README.md
│       ├── summ_percent_wells_nonzero_gas_in_us.tex
│       ├── well_covariate_balance.tex
│       ├── well_summary_stats_all_2018.tex
│       ├── well_summary_stats_aviris_matched.tex
│       └── well_summary_stats_lyon.tex
├── README.md
└── Snakefile: The file that organizes the analysis and tells Snakemake what to run. Start here to understand the flow.

```

</details>


## Data sources and availability

The data in this study come from a variety of sources, with the sources in bold providing the central contribution.

- Scientific studies (except as noted, all from the published papers and supplementary information on the journal website). Files are saved in `data/studies/`
    - [Alverez et al. (2018)](https://dx.doi.org/10.1126/science.aar7204)
        - data/studies/alvarez_etal_2018/aar7204_Database_S1.xlsx is used in the "ground studies" figures (e.g. A5). This paper's supplementary info results are plotted in figure 1. The other files are not used directly.
    - **[Duren et al. (2019)](https://dx.doi.org/10.1038/s41586-019-1720-3)**
        - This study provides the emissions from California wells.
        - data/studies/duren_etal_2019/Plume_list_20191031.csv and data/studies/duren_etal_2019/41586_2019_1720_MOESM3_ESM.xlsx provides data on methane plumes and their associated matched sites respectively. Both are available on the journal website.
        - data/studies/duren_etal_2019/AVIRIS-NG Flight Lines - AVIRIS-NG Flight Lines.csv provides data on which areas were flown over (this same file is used for the Frankenberg el al. flight coverage). AVIRIS-NG flight lines downloaded 2020-02-25 from https://docs.google.com/spreadsheets/d/1qbCEFO90rfEUMGLENY-qr6xmRjPk2JoGxcFFd6zO7kc/edit#gid=2050285191. Data portal: https://avirisng.jpl.nasa.gov/dataportal/. Flight lines also available here: https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1727
    - **[Frankenberg et al. (2016)](https://dx.doi.org/10.1073/pnas.1605617113)**
        - This study provides the emissions from Four Corners (New Mexico and Colorado) wells.
        - data/studies/frankenberg_etal_2016/AVNG_sources_all2.xlsx and data/studies/frankenberg_etal_2016/FourCorners_AV_NG_detections_Werner.xlsx provide data on sources and detected plumes. The former is available on the journal website. The latter was received by email from the authors when we asked for more detail on the detected methane plumes.
    - [Lyon et al. (2016)](https://dx.doi.org/10.1021/acs.est.6b00705)
        - This study provides data on detections for a wide variety of sites. It informs the Lyon et al. bars in figures 2 and A5, as well as panel B of table 1.
        - data/studies/lyon_etal_2016/es6b00705_si_004.xlsx provides data on their measures (note that leak sizes are not quantified). data/studies/lyon_etal_2016/es6b00705_si_005.xlsx provides data on sources. Both are available on the journal website.
        - The other files are informative, but not used directly.
    - [Omara et al. (2018)](https://dx.doi.org/10.1021/acs.est.8b03535)
        - data/studies/omara_etal_2018/Omara_etal_SI_tables.csv is used in the "ground studies" figures (e.g. A5). It is available from the journal website.
        - The other files are informative, but not used directly.
    - [Zavala-Araiza et al. (2018)](https://dx.doi.org/10.1525/elementa.284)
        - data/studies/zavala-araiza_etal_2018/elementa-6-284-s1.xlsx is used in the "ground studies" figures (e.g. A5). It is available from the journal website.
- US Agencies
    - [BEA](https://apps.bea.gov/iTable) (`data/bea/`)
        - Used for value added statistics
        - Files: current_cost_depreciation_of_private_fixed_assets_by_industry.xls and gdp_by_industry.xls
    - [EIA](https://www.eia.gov/electricity/data/state/emission_annual.xls): (`data/eia/emission_annual.xls`)
        - Used for emissions comparison statistic
        - File: emission_annual.xls
    - [EPA](https://www.epa.gov/ghgreporting/data-sets) (`data/epa/`)
        - Used for emissions inventory statistics
        - Files: egrid2018_summary_tables.pdf and EPA_GHG_report_2020_table_data.zip
    - [St. Louis Federal Reserve](https://fred.stlouisfed.org/series/CPILFENS) (`data/fred/CPILFENS.csv`)
        - Used for inflation adjustment
        - File: CPILFENS.csv
- Data providers
    - SNL: prices at trading hubs (`data/snl/`)
        - Data from the SNL "Commodity Charting" tool on 2020-10-26. See details in `data/snl/gas_prices_download_notes.txt`
        - Data file snl_gas_price_series.xls. Files na_gas_methodology.pdf and snl_gas_volume_series.xls may be useful references
    - **Enverus** (formerly Drillinginfo): Well production and characteristics (`data/drillinginfo/`)
        - Files in data/drillinginfo/docs/ may provide a useful reference.
        - Files in data/drillinginfo/entity_production_headers/ provide per-well data. In Drillinginfo terminology, "entity" can be a well or a lease, depending on the state and time. For this sample, they're always wells. The data are segmented by region to handle data download limits (CA_Kern_Production_Headers.CSV.gz, CA_notKern_Production_Headers.CSV.gz, CO_Production_Headers.CSV.gz, and NM_Production_Headers.CSV.gz)
        - Files in data/drillinginfo/entity_production_monthly provide well-by-month data. This data is used to provide information on the well operations (mainly production) _when the well was observed._ These follow the same geographic segmentation. (CA_Kern_Producing_Entity_Monthly_Production.CSV.gz, CA_notKern_Producing_Entity_Monthly_Production.CSV.gz, CO_Producing_Entity_Monthly_Production.CSV.gz, and NM_Producing_Entity_Monthly_Production.CSV.gz)
        - SHA256SUM provides hashes of these *.CSV.gz files.

All datasets are included in the `replication_public.zip` file, except the Enverus data, which is included in a separate `replication_drillinginfo.zip` file.
The Enverus data was extracted using their map interface; it is highly unlikely that another researcher would be able to download exactly the same files. Please contact Karl Dunkle Werner for access.
We believe our Enverus data subset can be shared with people who have access to the Enverus entity production headers and entity production monthly datasets.





## Development docs


### Old, Windows-specific notes

In the past, we worked to run the code on Windows without usign WSL. This subsection has notes on some Windows-specific issues we encountered. These are out of date, and are provided as a starting point to help debugging, but **we strongly recommend using WSL if running on Windows.**


Snakemake was unable to activate a conda environment at runtime. Therefore, we turned `use_conda` off, manually installed the appropriate R packages, and manually enabled the conda environment specified below.

Additionally, the ipopt package has more installation challenges on Windows. Instructions to install it manually are below.


- After creating the environment specified below, install cyipopt [manually](https://github.com/matthias-k/cyipopt#from-source-on-windows).
    - Download and extract ipopt from the link
    - Download and extract cyipopt
    - Copy the ipopt folders into cyipopt
    - Activate the environment
    - `python ../cyipopt/setup.py install` (assuming you saved the cyipopt directory next to the `methane_abatement` directory)




<details>
<summary>
Old Snakemake environment used on Windows
</summary>

```yml
# Snakemake is refusing to create environments on Windows, so use this env instead
# DO NOT set use_conda
# Usage:
# conda activate snakemake
# snakemake --cores 6 --keep-going > logfile.log 2>&1
name: snakemake
channels:
  - conda-forge
  - defaults
  - bioconda
dependencies:
  - cython=0.29 # dep for cyipopt manual install
  - future=0.18 # dep for cyipopt manual install
  - numexpr=2.7
  - numpy=1.21
  - pandas=1.3
  - pyarrow=5.0
  - python=3.9
  - pyprojroot=0.2
  - setuptools=50.3 # dep for cyipopt manual install
  - scipy=1.7
  - six=1.15 # dep for cyipopt manual install
  - snakemake-minimal=5.27.4
  - tqdm=4.62 # parallel progress bars

```

</details>




### Running R

The code is executed in a conda R environment managed by Snakemake. Sometimes you just want to run the R code manually. To activate the environment for debugging purposes, you can use mamba to figure out the name and activate it.

```sh
mamba env list
# pick the relevant entry, e.g. .snakemake/mamba/abcdef012345
mamba activate <long path listed above>
unset R_LIBS_USER
# Check that this result points to the environment specific to this project
Rscript -e '.libPaths()'
```

### Cmdstan install notes

There are several ways to install cmdstan. The code will manage it for you (via conda-forge and then re-compiled locally). See `code/setup_cmdstan.R`. Here are some notes about what issues I ran into with different approaches.

Cmdstan is available in mamba (via conda-forge). Pinning a particular version can be finicky if pinning other packages (e.g. R). The combination in `code/envs/r_scripts.yml` works, but some other combinations don't.
Installing via mamba brings pre-compiled binaries, which won't work if some other build flags are used (e.g. `-march=native`). Not all models are affected -- it seems like only ones that use SIMD instructions.
See https://discourse.mc-stan.org/t/recommended-compiler-flags-makes-rstan-model-crash/25689

Note that a mamba install brings the appropriate build tools, which is nice.


### Cleaning up

Many of the files generated during the build can be removed with `snakemake --dryrun --delete-all-output` (remove the `--dryrun` if you're sure).
Some of the files are not tracked by Snakemake and need to be cleaned up in other ways:

```sh
snakemake --dryrun --delete-all-output
snakemake --list-conda-envs
# See the two environment paths
mamba env remove -p <first path>
mamba env remove -p <second path>

# Code tries to clean up, but will leave CSVs if there are errors:
rm scratch/stan_output/*.csv
```
