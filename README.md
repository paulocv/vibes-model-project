# VIBES: A Multi-Scale Modeling Approach Integrating Within-Host and Between-Hosts Dynamics in Epidemics

![Banner: VIBES - Viral dynamics Individual-Based Epidemic Simulator](docs/assets/VIBES_GH_banner.png)

A multi-scale modeling framework for respiratory infectious diseases
that explicitly integrates viral dynamics based on patient-level data with population-level transmission on a data-driven network of social contacts.

This repository holds the code used for the following publications:

> [Ventura, P. C., Dam Jeong, Y., Litvinova, M., Kummer, A. G., Iwami, S., Yu, H., ... & Ajelli, M. (2026). VIBES: A multiscale modeling approach integrating within-host and between-hosts dynamics in epidemics. Proceedings of the National Academy of Sciences, 123(13), e2523055123.](https://doi.org/10.1073/pnas.2523055123)

> [Ventura, P. C., Kummer, A. G., Jeong, Y. D., Mhade, S., Iwami, S., Merler, S., ... & Ajelli, M. (2026). The impact of the distributions of the latent and infectious periods on epidemic dynamics. Research Square preprint.](https://doi.org/10.21203/rs.3.rs-11154983/v1)

This repository is licensed under the permissive GNU General Public License, but please cite the original publication if you use the code in any scientific project.

The core Viral Dynamics Individual-Based Epidemic Simulator (VIBES) code lives on a [separate repository](https://github.com/paulocv/vibes-core) and used here as a submodule. For this reason, you need to git clone this repo with the `--recurse-submodules` option.

Jump to:

- [Configuration guide](#configuration-guide)
- [General usage](#general-usage)
- [Steps for reproducing the results](#steps-for-reproducing-the-results)
- [Documentation of output files](#documentation-of-output-files)

# Configuration guide

## Requisites

You must have the following software installed in order to run the VIBES core code:

- git
- A C++ compiler
- The [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)

Additionally, to reproduce all results from the papers, you may need:

- [uv](https://docs.astral.sh/uv/getting-started/installation/), which manages the python interpreter and libraries. The project is pinned to python 3.9 and to the library versions used in the development (see [pyproject.toml](pyproject.toml)); uv downloads the right python version automatically, so you do not need to install it yourself.

## Downloading the necessary files

Clone this repository and its submodules using git:

```bash
git clone --recurse-submodules https://github.com/paulocv/vibes-model-project.git
```

## Compiling the VIBES engine

For more details, please visit the VIBES core engine repository: [https://github.com/paulocv/vibes-core](https://github.com/paulocv/vibes-core)

On a Unix system (Linux, MacOS), you can use the `compile_main.sh` script that comes with VIBES to quickly compile the engine, assuming you have the g++ compiler and GSL installed. From the root of this repository ("`cd vibes-model-project`" if needed), run:

```bash
./vibes-core/compile_main.sh
```

If compilation is successful, the command should return with no messages after a few seconds.

To check that configuration and compilation were successful, run a simple VIBES simulation using this command from the root of this repository:

```bash
./main inputs/quick_test.in outputs/quick_test
```

You should see some feedback messages from the engine on your terminal if successful, or an error message if something went wrong (e.g. missing a data file, compilation unsuccessful, failure to link with GSL). 

## Configuring the python environment 

The core VIBES engine is written in C++, but for this project we use python scripts for running multiple simulations, postprocessing outputs and visualizing results. 

The python environment is managed with [uv](https://docs.astral.sh/uv/). The exact interpreter (python 3.9) and package versions used in the development are declared in [pyproject.toml](pyproject.toml) and locked in [uv.lock](uv.lock). From the root of this repository, run:

```shell
uv sync
```

This creates a `.venv` directory with all required packages (it may take a few minutes). You can then either activate the environment and use `python` normally:

```shell
source .venv/bin/activate
```

or prefix any python command with `uv run`, without activating anything (e.g. `uv run python run_beta_grid.py ...`). The python commands in the rest of this document are written assuming the environment is activated.

To check if the environment creation was successful, you can run VIBES for multiple infectiousness values using:

```shell
python run_beta_grid.py inputs/quick_test.in outputs/tests/beta_grid_test
```

> Note: All python scripts must be called from the root of this repository, since they use paths relative to it (e.g. `./main` and `inputs/`).


# General usage

## Running VIBES once

After the VIBES engine is compiled, you can run it to produce a sequence of independent executions with the same parameters by calling the `./main` executable with a path to the input file and the output directory:

```bash
./main [path_to_input_file] [path_to_output_dir] [extra_args]
```

Full documentation on running the engine can be found in the VIBES core repository: [github.com/paulocv/vibes-core](https://github.com/paulocv/vibes-core)

## Running sets of simulations

The `run_beta_grid.py` script is used to call the VIBES engine for a predefined grid of relative infection risk values (the master `beta` parameter). For this script to work, you first need to compile the core VIBES engine from its C++ code into the `./main` executable.

Having the python environment configured (see [Configuring the python environment](#configuring-the-python-environment)), from the root directory of this repository, run:

```
python run_beta_grid.py [path_to_input_file] [path_to_output_dir] [extra-args]
```

Where [path_to_input_file] should be replaced by the path to an parameters file in the `inputs` directory, [path_to_output_dir] should be the path to the desired folder where results should be stored (preferably inside `outputs`), and [extra-args] can be replaced by command line arguments that are accepted either by the `run_beta_grid.py` script or the C++ VIBES engine.


## Postprocessing simulations

### Time series from a single simulation
Time-dependent metrics can be calculated from the results of each call to the VIBES core engine (the `./main` executable) using the `analysis/postprocess_time_series.py` script. Assuming your terminal is in the root of the repository, call:
```shell
python analysis/postprocess_time_series.py [path_to_output_dir]
```

Where [path_to_output_dir] should be replaced by the path informed to the `./main` script, or to one of the simulations created by `run_beta_grid.py`. This will aggregate multiple stochastic runs with the same set of parameters.

### Scalar metrics from a set of simulations

You can postprocess simulation sets created with `run_beta_grid.py` using the `analysis/postprocess_beta_grid.py`. From the root of this repository, run:

```shell
python analysis/postprocess_beta_grid.py [path_to_output_dir]
```

Where in this case [path_to_output_dir] is the output path informed to `run_beta_grid.py`, which should contain multiple simulation output subdirectories prefixed with `sim_`. The prefix can be modified with the `-s` or `--sim-prefix` command-line argument. For a list of all arguments, you can run:

```shell
python analysis/postprocess_beta_grid.py --help
```
Or check the source code.  


[//]: # (## Plotting results)


# Steps for reproducing the results

A full guide for reproducing published results from this code is available in the [`docs/reproducibility_guide.md](docs/reproducibility_guide.md) file.

Specific guide for each paper:

- ["VIBES: A multiscale modeling approach integrating within-host and between-hosts dynamics in epidemics"](docs/reproducibility_guide.md#first-paper)
- ["The impact of the distributions of the latent and infectious periods on epidemic dynamics"](docs/reproducibility_guide.md#second-paper)


# Documentation of output files

- Outputs from the core VIBES engine (through calling `./main`) are documented here: [https://github.com/paulocv/vibes-core/blob/main/docs/output_files.md](https://github.com/paulocv/vibes-core/blob/main/docs/output_files.md)

- Outputs from the python scripts from this repository are documented here: [docs/outputs.md](docs/outputs.md)
