# VIBES: A Multi-Scale Modeling Approach Integrating Within-Host and Between-Hosts Dynamics in Epidemics

_General details of the project and the model_

# Configuration guide

## Requisites

You must have the following software installed in order to run the VIBES core code:

- git
- A C++ compiler
- The [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)

Additionally, to reproduce all results from the papers, you may need:

- python (version 3.9.7 was used in the development)
- Required python libraries (see [requirements.in](/requirements.in)

## Downloading necessary files

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

## Running a test simulation

To check that the configuration was successful, run a simple VIBES simulation using this command from the root of this repository:

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

# Steps for reproducing the results

## Running sets of simulations

The `run_beta_grid.py` script is used to call the VIBES engine for a predefined grid of relative infection risk values (the master `beta` parameter). For this script to work, you first need to compile the core VIBES engine from its C++ code into the `./main` executable.

Having a python interpreter configured with the required libraries, from the root directory of this repository, run:

```
python run_beta_grid.py [path_to_input_file] [path_to_output_dir] [extra-args]
```

Where [path_to_input_file]

## Postprocessing simulations

## Plotting results