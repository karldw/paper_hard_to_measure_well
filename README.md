
# Hard to Measure Well: Can Feasible Policies Reduce Methane Emissions?

Authors: Karl Dunkle Werner [![ORCID logo](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-0523-7309) and [Wenfeng Qiu](https://wenfengqiu.com/)


## Steps to replicate (TL;DR)

1. Read this README
1. Make sure you have an appropriate OS (Linux or WSL2) and the necessary computing resources (see below)
1. Unzip the replication files.
1. If the data is saved somewhere outside the project folder, mount a copy inside the project folder. (Useful for development only)
1. Install Conda and Snakemake (see below)
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
# This is only necessary if the data are stored somewhere *outside*
# the project folder. You may need to change these paths to fit
# your situation
data_drive="$HOME/Dropbox/data/methane_abatement"
scratch_drive="$HOME/scratch/methane_abatement"
project_dir="$(pwd)"
mkdir -p "$scratch_drive" "$project_dir/data" "$project_dir/scratch"
sudo mount --bind "$data_drive"    "$project_dir/data"
sudo mount --bind "$scratch_drive" "$project_dir/scratch"


# 6. Install Conda and Snakemake
# If conda is not already installed, follow instructions here:
# https://docs.conda.io/en/latest/miniconda.html
conda env create --name snakemake --file code/envs/install_snakemake.yml
conda activate snakemake
snakemake --version
singularity --version
# Should show versions, not an error


# 7. Run Snakemake to create all outputs
# (this takes about a day with 4 CPU)
/usr/bin/time -v snakemake
# snakemake --dry-run to see what will be run


# 8. Check results (optional and slow)
# Check everything into git, rerun snakemake, and verify results are the same.
git init
git add .
git commit -m "Replication run 1"
snakemake --delete-all-output
rm -r scratch/*
rm -r .snakemake/conda
snakemake --use-conda --use-singularity --singularity-args='--cleanenv'
# Results should be binary-identical if everything worked correctly
# (except software_cites_r.bib, which has some manual edits)
git diff
```


## Setup

### Operating system

This code uses Singularity. You don't have to install it yourself, but you do have to be on an operating system where it can be installed. Good options are any recent version of Linux or Windows WSL2 (but not WSL1).

On macOS, or on Windows outside WSL2, things are more difficult. [One approach](https://sylabs.io/guides/3.5/admin-guide/installation.html#installation-on-windows-or-mac) is to install Vagrant, use Vagrant to create a virtual machine, and run everything inside that virtual machine. Good luck.

For more detail, see [Singularity's installation docs](https://singularity.hpcng.org/admin-docs/3.7/installation.html) (only the pre-install requirements; conda will [install Singularity](https://github.com/conda-forge/singularity-feedstock) for you)


### Software

This project uses [Snakemake](https://snakemake.readthedocs.io) (v6.8.0) and [Conda](https://docs.conda.io/en/latest/miniconda.html) (v4.10.3) to manage dependencies.
- To get started, first install Conda (mini or full-sized).
- Then use Conda to install Snakemake and Singularity from the file `install_snakemake.yml` (in the replication zipfile).

In a terminal:

```sh
conda env create --name snakemake --file code/envs/install_snakemake.yml
```

Run all other commands in that activated environment.
If you close the terminal window, you need to re-run `conda activate snakemake` before running the rest of the commands.
These downloads can be large.


#### What does Snakemake do?

Snakemake uses rules to generate outputs and manages the code environment to make it all work.

In particular, we're following a pattern Snakemake calls an [Ad-hoc combination of Conda package management with containers](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#ad-hoc-combination-of-conda-package-management-with-containers).

Snakemake uses Singularity (an alternative to Docker) to run code in a virtual environment, and uses conda to install packages.
All of this is handled transparently as the rules are run.

It can be useful to run `snakemake --dry-run` to see the planned jobs.

Snakemake keeps track of what needs to run and what doesn't. If something goes wrong midway through, snakemake will see that some outputs are up-to-date and others aren't, and won't re-run the things that don't need it.

The Snakefile is set up to retry failing jobs once, to avoid issues where temporary issues cause the build to fail (e.g. "Error creating thread: Resource temporarily unavailable").
If you would rather not restart failed jobs, remove the line `workflow.restart_times = 1` from `Snakefile`.
Note that Snakemake will still stop after failing twice (it will not run other jobs).


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
- Create your own [bind mounts](https://unix.stackexchange.com/questions/198590/what-is-a-bind-mount) to point to the `data` and `scratch` folders.
(See an example in `code/bind_mount_folders.sh`)

For people familiar with Singularity:
Note that `$SINGULARITY_BIND` doesn't work, because it's not used until the Singularity container is running, so Snakemake thinks files are missing.

For people familiar with symlinks:
Using symlinks (in place of bind mounts) do not work here, because Singularity will not follow them.



#### File structure

All files in `output/tex_fragments`, `data/generated`, and `scratch/` are auto-generated and can safely be deleted. All other files in `data/` should not be deleted.
*Some* files in `graphics/` are auto-generated, but the ones that are in the replication zipfile are not.
`data/` and `scratch/` are ignored by `.gitignore`.

### Other

- The PDF outputs are built with Latexmk and LuaLaTeX.
    - For size reasons, LuaLaTeX is not included in the set of software managed by conda. The `paper` job, which runs `latexmk` might fail if it's not installed on your computer. All the outputs up to that point will be present.
    - The tex files use some fonts that are widely distributed, but may not be installed by default.
- Note that the code depends on [`moodymudskipper/safejoin`](https://github.com/moodymudskipper/safejoin) which is *a different package* than `safejoin` on CRAN.
`moodymudskipper/safejoin` will be [renamed](https://github.com/moodymudskipper/safejoin/issues/44).
    - In case the original author deletes the repository, a copy is [here](https://github.com/karldw/safejoin).


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

These notes are modestly outdated, and aren't useful for replication.

### Other installation instructions

#### Installing Stan

1. Download and extract CmdStan
1. Add these lines to a file named `local` in the CmdStan `make/` directory. (Create `local` if it doesn't already exist)
```
O_STANC=3
STANCFLAGS+= --O --warn-pedantic
STAN_THREADS=true
STAN_CPP_OPTIMS=true
STANC3_VERSION=2.27.0  # change this version to match the downloaded cmdstan version
```
1. Edit user environment variable `CMDSTAN` to the folder (e.g. `~/.R/cmdstan-2.27.0`)
1. Windows only: `mingw32-make install-tbb` (even if `make` is installed)
1. Follow prompts, including adding the TBB dir to your path (Windows only)
1. Run `make build` and `make examples/bernoulli/bernoulli` (see install instructions)
1. Installation tries to download `stanc` (because compilation is a hassle), but sometimes I've had to download manually from https://github.com/stan-dev/stanc3/releases

#### Windows-specific instructions

- After installing conda, use `code/snakemake_environment_windows.yml` to create `snakemake` environment (will error if you already have one named `snakemake`)
    - `conda env create -f code/snakemake_environment_windows.yml`
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
