# VIBES: A Multi-Scale Modeling Approach Integrating Within-Host and Between-Hosts Dynamics in Epidemics

![Banner: VIBES - Viral dynamics Individual-Based Epidemic Simulator](docs/assets/VIBES_GH_banner.png)

_General details of the project and the model_

Jump to:

- [Configuration guide](#configuration-guide)
- [General usage](#general-usage)
- [Steps for reproducing the results](#steps-for-reproducing-the-results)

# Configuration guide

## Requisites

You must have the following software installed in order to run the VIBES core code:

- git
- A C++ compiler
- The [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)

Additionally, to reproduce all results from the papers, you may need:

- python (version 3.9.7 was used in the development)
- Required python libraries (see [requirements.in](/requirements.in))

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

To reliably run the python scripts with the same package versions as used in development, we recommend creating an environment and installing packages using the [requirements.txt](requirements.txt)) file. Having `python` and `pip` installed, run:

```shell
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The third command may take a few minutes to run. To check if the environment creation was successful, you can run VIBES for multiple infectiousness values using:

```shell
python run_beta_grid.py inputs/quick_test.in outputs/tests/beta_grid_test
```

# General usage

## Running sets of simulations

The `run_beta_grid.py` script is used to call the VIBES engine for a predefined grid of relative infection risk values (the master `beta` parameter). For this script to work, you first need to compile the core VIBES engine from its C++ code into the `./main` executable.

Having a python interpreter configured with the required libraries, from the root directory of this repository, run:

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

This guide shows how to re-run and postprocess simulations published in [Ventura et al. 2025, arXiv](https://arxiv.org/abs/2508.13354).

## Running and processing main analysis simulations

In the paper, we ran a beta grid for the VIBES model without interventions and with 100% isolation of symptomatic infections, both with 50 independent executions for each value of beta. To reproduce these simulations, run the following commands on a terminal open in the root of this repository:
```bash
python run_beta_grid.py inputs/model_paper_in/vl_model_baseline.in outputs/model_paper_out/vl_model_baseline -n 50
python run_beta_grid.py inputs/model_paper_in/vl_model_symp-isol.in outputs/model_paper_out/vl_model_symp-isol -n 50
```

This should create the `outputs/model_paper_out` directory and populate it with results from the simulations. To postprocess these simulations and obtain the data files that were used for the main text Figures 2 and 3, run:


# Documentation of output files

- Single execution outputs: https://github.com/paulocv/vibes-core/blob/main/docs/output_files.md
- 
