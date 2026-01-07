## Outputs from `run_beta_grid.py` 

Running the core engine through the `run_beta_grid.py` python wrapper will produce the following outputs:

### `index.csv`

A comma-separated file, containing each value of the master `beta` coefficient (`beta_coef`) for which the simulator was called, along with the relative path to the output folder for that value of beta (`sim_dir`). Usually the simulation output folder is named `sim_0000`, with numbers increased sequentially.

### `sim_0000`, `sim_0001`, ...

Each directory contains the outputs from the `./main` (VIBES engine) with a different beta value, whose contents are documented [here](https://github.com/paulocv/vibes-core/blob/main/docs/output_files.md).


## Outputs from `postprocess_beta_grid.py`

### `seqsets_analysis.csv`

Results of posptprocessing each set of simulations with the same parameters. The file contains a markup header that can be parsed with `toolbox.file_tools.read_config_file()`, and a csv content separated from the header by a `-----` line. 

#### Header variables

- `main_dir` = Path to the output directory as informed to `postprocess_beta_grid.py`.
- `expw_cases_start` = Threshold number of cumulative incident infections that mark the start of the exponential window. 
- `expw_cases_end` = Threshold number of cumulative incident infections that mark the end of the exponential window.
- `expw_min_steps` = Mininum number of time steps required within the exponential window.
- `exp_fit_min_rsquare` = Lowest r-square obtained from exponentially fitting the incidence within the exponential window. Used to detect possible bad fit.
- `quantile_q` = Probability used to determine the lower quantile of each metric, essentially for calculating confidence intervals. The upper quantile is given by `1 - quantile_q`.
- `steps_gp` = Number of time steps that are grouped to compose one day. Should be 1. / `dt`, where `dt` is the length of a time step in days.


#### CSV content variables


## Outputs from `postprocess_time_series.py`

TODO
