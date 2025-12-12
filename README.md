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
- Required python libraries (see [requirements.txt](/requirements.txt))

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

To check thtat the configuration was successful, run a simple VIBES simulation using this command from the root of this repository:

```bash
./main inputs/quick_test.in outputs/quick_test
```

You should see some feedback messages from the engine on your terminal if successful, or an error message if something went wrong (e.g. missing a data file, compilation unsuccessful, failure to link with GSL). 


# Steps for reproducing the results

## Running sets of simulations

The `run_beta_grid.py` script is used to call the VIBES engine for a predefined grid of relative infection risk values (the master `beta` parameter). For this script to work, you need to first compile the C++ core VIBES engine.

Having a python interpreter configured with the required libraries, from the root directory of this repository, run:

```
python run_beta_grid.py inputs/[path_to_input_file] outputs/[path_to_output_dir]
```

Where

## Postprocessing simulations

## Plotting results