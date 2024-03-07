
# Information Matters: Feasible Policies for Reducing Methane Emissions

Authors: Karl Dunkle Werner [![ORCID logo](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-0523-7309) and [Wenfeng Qiu](https://wenfengqiu.com/)


## Steps to replicate (TL;DR)

1. Read this README
1. Make sure you have an appropriate OS (Linux or WSL2) and the necessary computing resources (see below)
1. Unzip the replication files.
1. If the data is saved somewhere outside the project folder, mount or link a copy inside the project folder. (Useful for development only)
1. Install Mamba and Snakemake (see below)
1. Run Snakemake
1. Check results

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

# 6. Run Snakemake to create all outputs
# (this takes about a day with 12-core CPU)
/usr/bin/time -v snakemake
# snakemake --dry-run to see what will be run


# 7. Check results (optional and slow)
# Check everything into git, rerun snakemake, and verify results are the same.
git init
git add .
git commit -m "Replication run 1"
snakemake --delete-all-output
rm -r scratch/*
rm -r .snakemake/conda
snakemake --use-conda
# Results should be binary-identical if everything worked correctly
# (except software_cites_r.bib, which has some manual edits, and a few of the PDF graphs)
git diff
```


## Setup

### Operating system

The code was developed on Linux, and should mostly work on other operating systems. Requirements include:
(1) a working compiler (provided by mamba), and
(2) some of the commands use unix tools like `sed`, `tail`, and `cat`.
The compiler is a firm requirement; the unix tools could probably be avoided if you really wanted to.

This code used to use Singularity as a container mechanism, but since Singularity can only be installed on Linux, this required either (a) directly using Linux, (b) using WSL2, or (3) running a virtual machine. None of those are great for people who just want to run code.

Additionally, we got to a point where there were weird segfault in CmdStan when compiled with certain (common) build flags like `-march=native`, possibly related to [this issue](https://github.com/apptainer/singularity/issues/845). So.. that was a pain.


### Software

This project uses [Snakemake](https://snakemake.readthedocs.io) (v7.32.4) and [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (v4.10.3) to manage dependencies.
- To get started, first install Mamba (mini or full-sized).
- Then use Mamba to install Snakemake  from the file `code/envs/install_snakemake.yml` (in the replication zipfile).

See code above.

Run all other commands in that activated environment.
If you close the terminal window, you need to re-run `mamba activate snakemake` before running the rest of the commands.
These downloads can be large.

Note that other versions of Snakemake might work, but Snakemake has a lot of breaking changes from version to version. Versions of Snakemake >=7.8.0 and < 8.0.0 should probably be fine.
Mamba has fewer breaking changes; the newest should be fine.


#### What does Snakemake do?

Snakemake uses rules to generate outputs and manages the code environments to make it all work.
Snakemake uses mamba to install packages. All of this is handled transparently as the rules are run.

It can be useful to run `snakemake --dry-run` to see the planned jobs.

Snakemake keeps track of what needs to run and what doesn't. If something goes wrong midway through, snakemake will see that some outputs are up-to-date and others aren't, and won't re-run the things that don't need it.

Snakemake tries to manage computing resources (RAM and CPU) based on the claims made in the Snakemake rules. Snakemake does not enforce these limits; some of the code does.

In Snakemake, we have code that unsets these environment variables, with the aim of avoiding path confusion:
`RENV_PATHS_ROOT`, `R_LIBS`, `R_LIBS_USER`, `R_LIBS_SITE`, `CMDSTAN`.



### Files and data

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

- Store the `data` and `scratch` folders somewhere else (e.g. `data` in Dropbox).
- Create symbolic links (on Windows, 'junctions') to point to the `data` and `scratch` folders.
- If you don't like symlinks, [bind mounts](https://unix.stackexchange.com/questions/198590/what-is-a-bind-mount) are fine too.

Example:
```sh
ln -s path/to/dropbox/data/this_project data
ln -s path/to/scratch_drive scratch
```


#### File structure

All files in `output/tex_fragments`, `data/generated`, and `scratch/` are auto-generated and can safely be deleted. All other files in `data/` should not be deleted.
*Some* files in `graphics/` are auto-generated, but the ones that are in the replication zipfile are not.
`data/` and `scratch/` are ignored by `.gitignore`.

### Other

- The PDF outputs are built with Latexmk and LuaLaTeX.
    - For size reasons, LuaLaTeX is not included in the set of software managed by mamba. The `paper` job, which runs `latexmk` might fail if it's not installed on your computer. All the outputs up to that point will be present.
    - The tex files use some fonts that are widely distributed, but may not be installed by default.


## Computing resources for a full run

In addition to the programs above, parts of this program require significant amounts of memory and disk space.
Most parts also benefit from having multiple processors available.
(The slow parts parallelize well, so speed should increase almost linearly with processor count.)

The tasks that require significant memory are noted in the Snakemake file (see the `mem_mb` fields).
The highest requirement for any task is 10 GB, though most are far lower.
(These could be overstatements; we haven't tried that hard to find the minimum memory requirements for each operation.)
The programs also use about 80 GB of storage in `scratch/` in addition to the ~10 GB of input data and ~8 GB output data.

Running the whole thing takes 23 hours on a computer with 4 CPUs.
According to `/usr/bin/time -v`, it uses 22:45 of wall time and 81:45 of user time.
Maximum resident set size is (allegedly) 3.62 GiB (this seems low).


## Data sources and availability

The data in this study come from a variety of sources, with the sources in bold providing the central contribution.

- Scientific studies (except as noted, all from the published papers and supplementary information)
    - Alverez et al. (2018)
    - **Duren et al. (2019)**
    - **Frankenberg et al. (2016)**
        - Includes data received by email from the authors
    - Lyon et al. (2016)
    - Omara et al. (2018)
    - Zavala-Araiza et al. (2018)
- US Agencies
    - BEA
    - EIA
    - EPA
    - St. Louis Federal Reserve
- Data providers
    - SNL: prices at trading hubs
    - **Enverus** (formerly Drillinginfo): Well production and characteristics

All datasets are included in the replication_public.zip file, except the Enverus data.
I believe my Enverus data subset can be shared with people who have access to the Enverus entity production headers and entity production monthly datasets.




## Development docs

#### Windows-specific instructions

** Note: these instructions worked for an older version of the code, but may no longer work on Windows. Try running in WSL2 instead. **

- After installing mamba, use `code/snakemake_environment_windows.yml` to create `snakemake` environment (will error if you already have one named `snakemake`)
    - `mamba env create -f code/snakemake_environment_windows.yml`
- Install cyipopt [manually](https://github.com/matthias-k/cyipopt#from-source-on-windows).
    - Download and extract ipopt
    - Download and extract cyipopt
    - Copy the ipopt folders into cyipopt
    - Activate `snakemake` environment
    - `python ../cyipopt/setup.py install` (assuming you saved the cyipopt directory next to the `methane_abatement` directory)


##### Running on Windows
- Activate the newly created `snakemake` environment, and do not set `--use_conda` when running Snakemake.
- There will be some warnings.

#### Overleaf

Connect to Overleaf by Git. See [details here](https://www.overleaf.com/learn/how-to/How_do_I_push_a_new_project_to_Overleaf_via_git%3F).
```sh
git remote add overleaf https://git.overleaf.com/5d6e86df6820580001f6bdfa
git checkout master
git pull overleaf master --allow-unrelated-histories
```


### Running R

The code is executed in an R environment managed by Snakemake. To activate the environment for debugging purposes, you can use mamba to figure out the name and activate it.

```sh
mamba env list
# pick the relevant entry, e.g. .snakemake/mamba/abcdef012345
mamba activate <long path listed above>
unset R_LIBS_USER
Rscript -e '.libPaths()' # should point to the env specific to this project
```

### Cmdstan install notes

There are several ways to install cmdstan. The code will manage it for you (via conda-forge and then re-compiled locally). See `code/setup_cmdstan.R`. Here are some notes about what issues I ran into with different approaches.

#### Mamba

cmdstan is available in mamba (via conda-forge). Pinning a particular version can be finicky if pinning other packages (e.g. R)
Installing via mamba brings pre-compiled binaries, which won't work if some other build flags are used (e.g. -march=native). Not all models are affected - it seems like only ones that use SIMD instructions.
See https://discourse.mc-stan.org/t/recommended-compiler-flags-makes-rstan-model-crash/25689

Note that a mamba install brings the appropriate build tools, which is nice.

#### Tarball

I was running into issues where Eigen was unhappy with bounds checks for some reason. This seemed to depend on the build flags used, so might be avoidable.


### Cleaning up

Many of the files generated during the build can be removed with `snakemake --dryrun --delete-all-output` (remove the `--dryrun` if you're sure).
Some of the files are not tracked by Snakemake and need to be cleaned up in other ways:

```sh
snakemake --list-conda-envs
# See the two environment paths
mamba env remove -p <first path>
mamba env remove -p <second path>

# Code tries to clean up, but will leave CSVs if there are errors:
rm scratch/stan_output/*.csv
```
