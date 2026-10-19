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

This data frame contains the estimates of metrics for each set of independent VIBES executions with the same parameters. Most metrics have a mean, median, standard deviation, lower quantile and upper quantile estimates with the respective suffixes: `_mean`, `_median`, `_std`, `_quant_low`, and `_quant_high`. 

Many variables are exported within this process. We document here the most relevant ones that were used in our scientific publications, along with their preffixes:

- `avg_gentime_asymp` = Average generation time of asymptomatic infectors
- `avg_gentime_h` = Average generation time of infections in households
- `avg_gentime` = Overall average generation time
- `avg_gentime_sympt` = Average generation time of symptomatic infectors
- `avg_effective_infectious_period` = Average effective (realized) infectious period, defined as the time between start of infectious period and each infection.
- `avg_latent_period_infected` = Average latent period over the ensemble of infected individuals (each infectee counts as one point).
- `avg_latent_period_infector` = Average latent period over the ensemble of _infectors_. The infector of each infector-infectee pair counts as one point, meaning that the latent period of each individual is accounted _r_ times (_r_ = number of secondary infections).
- `avg_serialint` = Overall average serial interval
- `empir_repnum` = Reproduction number calculated from infection chains whose primary infection occurred within the exponential window
- `exp_growth` = Ggrowth rate, calculated from the exponential fit of the incidence in the exponential phase.
- `exp_repnum` = Reproduction number calculated from the exponential growth rate
- `final_frac_r` = Fraction of removed individuals by the end of the simulation
- `frac_infecs_c` = Fraction of all infections that occurred in the community
- `frac_infecs_h` = Fraction of all infections that occurred in households
- `frac_infecs_s` = Fraction of all infections that occurred in schools
- `frac_infecs_w` = Fraction of all infections that occurred in workplaces
- `presymp_transm` = Fraction of presymptomatic transmission
- `window_gentime` = Average generation time for infections that occurrred within the exponential window
- `window_serialint` = Serial interval for infections that occurrred within the exponential window

### `seqsets_hist_[**].csv`

CSV files with the `seqsets_hist_` preffix contain histograms of some metrics, aggregated for each set of independent executions:

- `seqsets_hist_r.csv` = Number of secondary infections
- `seqsets_hist_tg.csv` = Generation time
- `seqsets_hist_serialint.csv` = Serial interval
- `seqsets_hist_tass.csv` = Time of secondary infection relative to time of infector symptom onset. Negative values represent pre-symptomatic transmission.
- `seqsets_hist_avg_r_vs_eff_infec_period` = Average reproduction number grouped by the length of the effective infectious period. This is not a histogram in the strict sense, but a grouped average.
- `seqsets_hist_avg_r_vs_latent_period` = Average reproduction number grouped by the length of the latent period. This is not a histogram in the strict sense, but a grouped average.

Generation time histograms are further disaggregated between `sympt` (from symptomatic infectors), `asympt` (from asymptomatic infectors), `window` (within the exponential window) and `total` (overall). Secondary cases and TASS are also disaggregated between `total` and `window`.


### `exec_samples.csv`

This file contains selected metrics from individual VIBES executions, excluding those that did not produce a valid exponential window. Each row contains metrics from a single execution.

Variables: 

- `exs_label` = Name of the set of independent executions with the same parameters.
- `i_exec` = ID of the execution within the same `exs_label`.
- `infec_peak` = Peak amplitude in the fraction of infectious individuals.
- `infec_peak_t` = Time of the peak (in steps) in the fraction of infectious individuals.
- `final_frac_r` = Fraction of removed individuals by the end of the simulation.
- `incidence_peak_t` = Time of the peak (in steps) of the incidence of new infections per time step.
- `incidence_peak_shifted_t` = Time of the peak (in steps) of the incidence of new infections per time step, shifted by the exponential window start.
- `incidence_daily_peak_shifted_t` = Time of the peak (in steps) of the incidence of new infections per day, shifted by the exponential window start.

## Outputs from `postprocess_time_series.py`

This script produces processed time series from the raw VIBES engine outputs.

### `a_time_series.csv` 

General time series statistics from the set of independent VIBES executions with the same parameters.
Most metrics have a mean, median, standard deviation, lower quantile and upper quantile estimates with the respective suffixes: `_mean`, `_median`, `_std`, `_quant_low`, and `_quant_high`.

Most time series are indexed by the id of the time step, `i_t`, which is also shifted such that `i_t = 0` represents the time step in which the exponential window started. 

Among the many variables exported by this script, we document here the most relevant ones that were used in our scientific publications, along with their preffixes:

- `new_E` = Number of new infections per time step.
- `newE_days` = Number of new infections per day, where `i_t` is in days instead of time step.
- `R` = Reproduction number calculated from transmission chains whose primary infection occurred at the given time step.
- `R_h` = Reproduction number of secondary infections in households, calculated from transmission chains whose primary infection occurred at the given time step.
- `R_s` = Same as `R_h` for schools.
- `R_w` = Same as `R_w` for workplaces.
- `R_c` = Same as `R_c` for the community.
- `T_g` = Generation time calculated from infections that occurred at the given time step.

Metrics with the `smooth_` prefix represent the 7-day moving average of the corresponding metric.
