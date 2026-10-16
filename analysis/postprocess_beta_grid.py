"""
Calculate scalar metrics from a set of VIBES simulations.
Each metric is a statistical property of each simulation (a set of
independent executions).
This procedure also calculate histograms for the generation time, serial
interval, and symptom-to-infection time (TASS).

Multiprocessing is implemented for the executions of each given set.
Each set is iterated sequentially.
Run "calc_exec_seqsets.py --help" for program usage.

To implement new metrics, please mind the following procedures:
* Add a boolean to determine its calculation (or use an existing one)
* Update the dependency graph in function 'create_dependent_bools'
* Implement its calculation for a single execution in function 'calc_single_execution'
    * If the metric depends on the validity of the exponential window, manually set it to None if expw_valid is False.
* Implement its aggregation over executions in the `calc_execset_statistics` function.

Usage
-----
The program has only one mandatory argument, which is the path to the
beta grid output directory. Example:
```
python analysis/postprocess_beta_grid.py outputs/tests/beta_grid_test
```

For a complete list of arguments, run:
```
python analysis/postprocess_beta_grid.py --help
```
Check also the variables `params` and `bool_switches` in the source code
for further customization.
"""

import argparse
import glob
from datetime import datetime

# import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys
import time
import warnings

from collections import defaultdict

import scipy.stats
from pathos.multiprocessing import ProcessPool

sys.path.append("./")
from analysis_tools import ExecSet, get_time_window_by_thresholds, calc_expon_repnum_from_tghist, LAYER_ID, \
    FILE_SECTION_BREAK
from analysis_tools import LAYER_CHAR, INVALID_SERIAL_TIME, NUM_LAYERS
from toolbox.file_tools import read_config_file, write_config_string, seconds_to_hhmmss

# Global parameters
DEFAULT_SET_PREFIX = "sim"
DEFAULT_EXEC_PREFIX = "run"
DEFAULT_OUTPUT_PREFIX = "seqsets_analysis"
DEFAULT_HIST_PREFIX = "seqsets_hist"
DEFAULT_CURVE_LENGTH = 81  # Max length of transmissibility curves, for runs that do not show this on header.
DEFAULT_NUM_THREADS = 4   # Number of concurrent processes
PATH_SEP = os.path.sep

# Options
pd.set_option('display.max_rows', 8300)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 150)

# Suppresses pandas performance warning for fragmented dataframe. There's no performance loss in this use case.
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def main():

    # Define which metrics should be or not calculated
    # ------------------------------------------------
    bool_switches = defaultdict(
        lambda: False,
        # --- Reproduction numbers
        calc_empir_repnum=True,
        calc_exp_repnum=True,

        # --- Generation time
        calc_total_gentime=True,
        calc_window_gentime=True,
        calc_sympt_gentime=True,

        # --- Serial interval
        calc_total_serialint=True,
        calc_frac_neg_serialint=True,
        calc_window_serialint=True,

        # --- Histograms
        calc_tg_hist=True,
        calc_R_hist=True,
        calc_serialint_hist=True,
        calc_window_tg_hist=True,  # Overridden if calc_exp_repnum is True
        calc_window_R_hist=True,
        calc_window_serialint_hist=True,
        calc_sympt_tg_hist=True,  # Tg by source symptom status
        calc_tass_hist=True,
        calc_presymp_transm_samples=False,  # Individual pre-symptomatic transmission fraction

        # --- Time series properties
        calc_infec_peak=True,
        calc_layer_absence=True,  # num_nodes - presence in a given layer.

        # --- Misc
        calc_attack_rates=False,  # Not implemented
        calc_presymp_transm=True,
        calc_test_behav_rates=True,
        calc_latent_periods=True,

        export=True,
        plot=True,  # Not used yet
    )

    # Numeric post-processing parameters
    # ----------------------------------
    params = dict(
        # Exponential phase detection (also used to determine whether there was a significant outbreak)
        expw_cases_start=1000,
        expw_cases_end=7000,
        expw_min_steps=round(0.0 / 0.25),
        exp_fit_min_rsquare=0.98,  # Discards exponential window measures whose fit R^2 s under this value

        # Statistics over executions
        quantile_q=0.025,
        steps_gp=4,  # Number of steps to group to compose a "day". This is set independently of dt.
    )

    _set_cap = None  # Cuts the number of execution sets. None to run all

    # PREAMBLE
    # --------
    args: argparse.Namespace = parse_arguments()
    num_threads = int(args.num_threads)
    create_dependent_bools(bool_switches)
    execset_list = create_list_of_execsets(args.main_dir, set_prefix=args.set_prefix, exec_prefix=args.exec_prefix)
    num_sets = len(execset_list)

    # MAIN CALCULATION LOOP
    # ---------------------
    pool = ProcessPool(ncpus=num_threads)
    final_df = pd.DataFrame()
    final_samples_df_list = list()
    hists = {key: list() for key in [
        "tg", "R", "serialint", "tass", "avg_r_vs_latent_period", "avg_r_vs_eff_infec_period"
    ]}  # Stores lists of histograms

    print("Importing and calculating individual stats...")

    # Warn about some relevant modifications
    if not bool_switches["export"]:
        print("******************\nWARNING: EXPORT = FALSE\n******************\n")
    if _set_cap is not None:
        print("******************\nREMINDER: _set_cap is active\n******************\n")

    _set_cap = len(execset_list) if _set_cap is None else _set_cap
    for i_set, exset in enumerate(execset_list[:_set_cap]):
        print(f"Exec set = \"{exset.label}\" ({i_set + 1} of {num_sets})")

        # --- Calculation of metrics from single executions
        def task(input_tup):
            i_exec, exec_dir = input_tup
            data = import_single_execution(exec_dir, bool_switches)
            return calc_single_execution_metrics(data, bool_switches, params, i_exec)

        xt0 = time.time()
        exec_result_list: list[ExecResults] = pool.map(task, enumerate(exset.exec_dirs))
        xtf = time.time()

        print("\tExec metrics time  ........................ = {:0.3f}s".format(xtf - xt0))
        metrics_cpu_time = sum(res.cpu_time for res in exec_result_list)
        print("\tCPU time calculating metrics .............. = {:0.3f}s ({:0.3f}s per thread)"
              "".format(metrics_cpu_time, metrics_cpu_time / num_threads))

        # --- Aggregation: statistics over the set of executions
        xt0 = time.time()
        stat_df, hist_dict, out_samples_df = calc_execset_statistics(exec_result_list, bool_switches, params, exset)
        xtf = time.time()
        print("\tExec time calculating statistics over runs  = {:0.3f}s".format(xtf - xt0))
        # print(stat_df.T)

        # --- Storage of each result into containers for all simulations
        final_df = final_df.append(stat_df, ignore_index=True)
        final_samples_df_list.append(out_samples_df)
        for key, l in hists.items():
            # Stores the calculated hists. Use None for non-existent ones, as some simulations may lack some hists.
            if key in hist_dict:
                l.append(hist_dict[key])
            # l.append(hist_dict[key] if key in hist_dict else None)


    # Set index to the exec set labels
    final_df.index = [exs.label for exs in execset_list[:_set_cap]]
    final_df.sort_index(axis=1, inplace=True)  # Sort columns alphabetically
    print(final_df.T)

    # Same with the samples df, a multi-indexed dataframe
    final_samples_df = pd.concat(
        final_samples_df_list,
        keys=[exs.label for exs in execset_list[:_set_cap]],
        names=["exs_label", "i_exec"],
    )

    # EXPORT OUTPUTS
    # --------------
    if bool_switches["export"]:
        main_fname = os.path.join(args.main_dir, args.output_prefix + ".csv")
        if not os.path.exists(os.path.dirname(main_fname)):
            os.mkdir(os.path.dirname(main_fname))

        # Main df with scalar metrics and statistics
        export_main_dict_with_header(main_fname, final_df, args, params)

        # Df with full sets of execution metric samples
        _fname = os.path.join(args.main_dir, "exec_samples.csv")
        final_samples_df.to_csv(_fname)

        # Histograms
        for file_suffix, hist_key in [
            ("tg", "tg"),
            ("r", "R"),
            ("serialint", "serialint"),
            ("tass", "tass"),
            ("avg_r_vs_latent_period", "avg_r_vs_latent_period"),
            ("avg_r_vs_eff_infec_period", "avg_r_vs_eff_infec_period"),
        ]:
            if len(hists[hist_key]) > 0:
                export_hist_df(make_hist_fname(file_suffix, args), hists[hist_key], execset_list)


# ----------------------------------------------------------------------------------------------------------------------
# PROGRAM OPERATIONAL FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------
def parse_arguments():
    """Define and interpret command line arguments."""
    parser = argparse.ArgumentParser(description="Process a sequence of execution sets.")

    parser.add_argument("main_dir")
    parser.add_argument("-s", "--set-prefix", "--sim-prefix", default=DEFAULT_SET_PREFIX,
                        help="Expected prefix of the simulations (sets of executions).")
    parser.add_argument("-e", "--exec-prefix", default=DEFAULT_EXEC_PREFIX,
                        help="Expected prefix of the single execution (run) directories.")
    parser.add_argument("-o", "--output-prefix", default=DEFAULT_OUTPUT_PREFIX,
                        help="Prefix of the output dataframe file with calculated metrics.")
    parser.add_argument("--hist-prefix", default=DEFAULT_HIST_PREFIX,
                        help="Prefix of the output histogram files.")
    parser.add_argument("-t", "--num-threads", "--ncpus", default=DEFAULT_NUM_THREADS,
                        help="Number of threads (workers) to use for the import and pre-analysis.")

    return parser.parse_args()


def create_dependent_bools(bool_switches: defaultdict):
    """Based on provided boolean switches, create other dependent ones
    based on a predefined dependency graph.
    Changes are made in plcae to bool_switches.
    """
    b = bool_switches  # Just a name alias to shorten lines

    # Primary (input) bool switches that influence each other
    # -------------------------------------------------------

    if b["calc_exp_repnum"]:
        b["calc_window_tg_hist"] = True

    # Dependent bool switches
    # -----------------------

    # --- Defines whether to filter node_event_df by valid serial intervals in the exp window.
    if b["calc_window_serialint"] or b["calc_window_serialint_hist"]:
        b["calc_window_serialint_event_df"] = True

    # --- Defines whether to calculate the exponential window
    if any((b[k] for k in ["calc_empir_repnum", "calc_window_gentime", "calc_window_serialint_event_df",
                           "calc_tg_hist", "calc_R_hist", "calc_serialint_hist",
                           "calc_window_tg_hist", "calc_window_R_hist", "calc_window_serialint_hist"])):
        b["calc_window_event_df"] = True

    if b["calc_exp_growth"] or b["calc_exp_repnum"]:
        b["calc_exp_fit"] = True
        b["calc_window_case_series"] = True

    # --- Data sliced in the exponential window
    if b["calc_window_event_df"] or b["calc_window_case_series"] or b["calc_window_tg_hist"]:
        b["calc_exp_window"] = True
        b["calc_case_series"] = True

    # --- Serial interval and symptom-related metrics
    if b["calc_window_serialint"] or b["calc_window_serialint_hist"]:
        b["calc_window_serialint_event_df"] = True

    if b["calc_total_serialint"] or b["calc_frac_neg_serialint"] or b["calc_window_serialint_event_df"]:
        b["filter_valid_serialint"] = True

    if any(b[k] for k in ["calc_presymp_transm", "calc_tass_hist", "calc_sympt_tg_hist", "calc_presymp_transm_samples"]):
        b["calc_symptsrc_event_df"] = True

    # --- Defines which files to import
    # Unused switches
    b["import_r_tg_series"] = False
    b["import_sympt_series"] = False

    # Precalculated Tg histogram (overall)
    if b["calc_tg_hist"]:
        b["import_tg_hist"] = True

    # Node event data
    if any((b[k] for k in ["filter_valid_serialint", "calc_window_event_df", "calc_symptsrc_event_df",
                           "calc_sympt_gentime", "calc_test_behav_rates"])):
        b["import_node_event_data"] = True

    # status_series file
    if b["calc_case_series"] or b["calc_infec_peak"]:
        b["import_status_series"] = True

    if b["calc_attack_rates"]:
        b["import_attack_rates"] = True

    if b["calc_layer_absence"]:
        b["import_presence_series"] = True


def create_list_of_execsets(main_dir: str, set_ids=None, exec_ids=None,
                            set_prefix=DEFAULT_SET_PREFIX, exec_prefix=DEFAULT_EXEC_PREFIX,
                            set_fmt="{:04d}", exec_fmt="{:05d}"):
    """Find paths of the execution sets in the given folder, and then
    the path of each execution from each set.
    """

    # Define the paths of the execution sets.
    # ---------------------------------------
    full_set_prefix = os.path.join(main_dir, set_prefix)

    if set_ids is None:
        # Finds all the subdirectories matching the set pattern
        execset_paths: list = glob.glob(full_set_prefix + "_*" + PATH_SEP)
        execset_paths.sort()
    else:
        # Creates the list of subdirectory paths based on the given ids
        execset_paths = [full_set_prefix + "_" + set_fmt.format(i) + PATH_SEP for i in set_ids]

    # Define the paths of the executions in each set
    # ----------------------------------------------
    return [ExecSet(path, label=os.path.basename(path[:-1]),
                    exec_ids=exec_ids, exec_prefix=exec_prefix, exec_fmt=exec_fmt)
            for path in execset_paths]


def write_output_header(args, params):
    """Prepare header of the main scalar outputs file, containing
    the post-processing parameters.
    """
    header = str()

    header += str(datetime.now()) + '\n'
    header += write_config_string({
        "main_dir": args.main_dir,
    })
    header += write_config_string(params)
    header += FILE_SECTION_BREAK

    return header


def export_main_dict_with_header(fname, final_df, args, params):
    """Write all contents of the main scalar outputs file."""
    header = write_output_header(args, params)
    content = final_df.to_csv()

    # --- Exports
    with open(fname, "w") as fp:
        fp.write(header)
        fp.write(content)


def make_hist_fname(suffix, args):
    """Prepare the file name of a histogram output file."""
    return os.path.join(args.main_dir, args.hist_prefix + f"_{suffix}.csv")


def export_hist_df(fname, df, execset_list):
    """Save the contents of a histogram data frame into a file."""
    final_tghist_df = pd.concat(df, axis=1, keys=[exs.label for exs in execset_list])
    final_tghist_df.to_csv(fname)


# ----------------------------------------------------------------------------------------------------------------------
# DATA MANIPULATION FUNCTIONS AND CLASSES
# ----------------------------------------------------------------------------------------------------------------------

class ExecData:
    """
    Class to hold data from a single simulation execution, for specific use in calc_exec_seqsets.py.
    There may be classes with the same name in other modules of this project which are _not_ interchangeable.
    """

    # --- Data containers
    # Imported
    summary_dict: dict = None
    r_tg_df: pd.DataFrame = None
    status_df: pd.DataFrame = None
    sympt_df: pd.DataFrame = None
    presence_df: pd.DataFrame = None
    tg_hist_df: pd.DataFrame = None
    node_event_df: pd.DataFrame = None
    attack_rates_df: pd.Series = None

    # Calculated, views
    new_e_series: pd.Series = None  # Number of new infections (S -> E) in each time step (time series)
    case_series: pd.Series = None  # Accumulated number of cases (since E status)
    window_case_series: pd.Series = None   # Slice of the case_series in the exponential time frame.
    window_node_event_df: pd.DataFrame = None  # Slice of the node_event_df in the exponential time frame.
    node_event_serialint_df: pd.DataFrame = None  # Only valid serial interval entries from node_event_df.
    node_event_tg_df: pd.DataFrame = None  # Only valid own_generation_time values.
    window_node_event_serialint_df: pd.DataFrame = None  # Valid serial interval inside the exponential time frame.
    symptsrc_event_df: pd.DataFrame = None  # Filtered by nodes infected by symptomatic sources.
    asympsrc_event_df: pd.DataFrame = None  # Nodes infected by asymptomatic sources

    # --- Execution parameters and metrics
    # Imported
    dt: float = None
    num_nodes: int = None

    # Calculated
    exp_window: tuple = None  # (i_t_start, i_t_end)


class ExecResults:
    """Class to hold all results and some other parameters from a
    single execution.
    """
    __slots__ = ["metrics", "cpu_time", "dt", "num_nodes", "TMP",
                 "tg_hist", "r_hist", "serint_hist", "w_tg_hist", "w_r_hist",
                 "w_serint_hist", "sympt_tg_hist", "asymp_tg_hist",
                 "tass_hist", "tass_lhists", "max_curve_length",
                 "avg_r_vs_latent_period_df", "avg_r_vs_eff_infec_period_df",
                 "presymp_transm_samples",
                 ]

    def __init__(self):
        self.metrics = dict()

        # Histograms (non-windowed and windowed)
        self.tg_hist: pd.Series = None
        self.w_tg_hist: pd.Series = None  # Windowed generation time histogram
        self.w_r_hist: pd.Series = None  # Windowed secondary cases histogram
        self.w_serint_hist: pd.Series = None  # Windowed serial interval histogram
        self.sympt_tg_hist: pd.Series = None
        self.asymp_tg_hist: pd.Series = None
        self.tass_hist: pd.Series = None  # Time-After-Source-Symptoms histogram
        self.tass_lhists: list[pd.Series] = None  # TASS histograms per layer
        self.avg_r_vs_latent_period_df: pd.DataFrame = None
        self.avg_r_vs_eff_infec_period_df: pd.DataFrame = None
        self.presymp_transm_samples: pd.Series = None  # ....

        # Scalar metrics and parameters
        self.dt = None
        self.num_nodes = None
        self.max_curve_length = None
        self.cpu_time = 0.0
        self.TMP = None
        # self.metrics = defaultdict(lambda: None)


def import_single_execution(exec_dir: str, bool_switches: defaultdict):
    """Load all requested data from a single VIBES execution, given its
    directory path.
    """
    data = ExecData()

    data.summary_dict = read_config_file(os.path.join(exec_dir, "summary.out"))
    data.dt = float(data.summary_dict["step_dt"])
    data.num_nodes = int(data.summary_dict["num_nodes"])

    if bool_switches["import_r_tg_series"]:
        data.r_tg_df = pd.read_csv(os.path.join(exec_dir, "r_and_tg.csv"), skiprows=1, header=0).rename(
            str.strip, axis=1)

    if bool_switches["import_status_series"]:
        data.status_df = pd.read_csv(os.path.join(exec_dir, "status_series.csv"), skiprows=1, header=0).rename(
            str.strip, axis=1)
        data.new_e_series = pd.Series(data.status_df["num_S"].values[:-1] - data.status_df["num_S"].values[1:])

    if bool_switches["import_sympt_series"]:
        data.sympt_df = pd.read_csv(os.path.join(exec_dir, "symptom_series.csv"), skiprows=1, header=0).rename(
            str.strip, axis=1)

    if bool_switches["import_presence_series"]:
        data.presence_df = pd.read_csv(os.path.join(exec_dir, "presence_series.csv"), skiprows=1, header=0).rename(
            str.strip, axis=1)

    if bool_switches["import_tg_hist"]:
        try:
            data.tg_hist_df = pd.read_csv(os.path.join(exec_dir, "tg_hist.csv"), skiprows=1, header=0).rename(
                str.strip, axis=1)
        except FileNotFoundError:
            data.tg_hist_df = None

    if bool_switches["import_node_event_data"]:  # -- Time consuming
        fname = os.path.join(exec_dir, "node_event_data.csv")
        with open(fname, "r") as fp:
            nrows = int(fp.readline().split(",")[1]) + 1
        data.node_event_df = pd.read_csv(fname, skiprows=1, header=0, nrows=nrows, sep=',', engine="c").rename(
            str.strip, axis=1)

    if bool_switches["import_attack_rates"]:
        data.attack_rates_df = pd.read_csv(os.path.join(exec_dir, "attack_rates.csv"), skiprows=1, header=0, sep=',',
                                           engine="c").rename(str.strip, axis=1)

    return data


def calc_single_execution_metrics(data: ExecData, bool_switches: defaultdict, params: dict, i_exec: int):
    """Calculate the metrics of a single VIBES execution.
    This is expected to be the most time-consuming part in the main script.
    """

    res = ExecResults()

    xt0 = time.time()

    # Preamble
    # --------
    # Calculate the time series of cumulative cases
    if bool_switches["calc_case_series"]:
        data.case_series = data.num_nodes - data.status_df["num_S"]

    # Create a filtered node event df with valid generation times
    if bool_switches["import_node_event_data"]:
        data.node_event_tg_df = data.node_event_df[data.node_event_df["own_generation_time"] > 0]

    # Create a filtered node_event_df with only valid serial intervals
    if bool_switches["filter_valid_serialint"]:
        data.node_event_serialint_df = data.node_event_df[data.node_event_df["own_serial_time"] > INVALID_SERIAL_TIME]

    # Create filtered event df for nodes whose infection source was symptomatic.
    if bool_switches["calc_symptsrc_event_df"]:
        if "own_tass" not in data.node_event_tg_df:
            warnings.warn("\nWARNING: node_event_data does not have 'own_tass' column. "
                          "Cannot calculate related metrics.")
            bool_switches["calc_presymp_transm"] = bool_switches["calc_symptsrc_event_df"] = False
            # UPDATE HERE IF MORE METRICS DEPEND ON THIS
        else:
            mask = data.node_event_tg_df["own_tass"] > INVALID_SERIAL_TIME
            data.symptsrc_event_df = data.node_event_tg_df[mask]
            data.asympsrc_event_df = data.node_event_tg_df[~mask]

    # Constant metrics (expected to be the same for all executions)
    res.dt = data.dt
    res.num_nodes = data.num_nodes
    res.max_curve_length = int(data.summary_dict.get("max_curve_length", DEFAULT_CURVE_LENGTH))

    # Ready-to-use metrics (these require no boolean switches)
    # --------------------------------------------------------
    res.metrics["final_num_r"] = int(data.summary_dict["final_num_removed"])
    # Calc. fraction of infections in each layer
    num_infections = int(data.summary_dict["num_infections"])
    for lc in LAYER_CHAR:
        res.metrics["frac_infecs_" + lc] = int(data.summary_dict["num_infections_" + lc]) / num_infections

    res.metrics["num_steps"] = int(data.summary_dict["num_steps"])

    # Exponential window and metrics in it
    # ------------------------------------
    # IMPORTANT: please register new metrics to None if they are stored inside the "if expw_valid" clause.

    if bool_switches["calc_exp_window"]:
        # Calculate window using threshold in the series of cumulative cases
        (i_t_start, i_t_end) = res.metrics["expw_it_start"], res.metrics["expw_it_end"] = \
            get_time_window_by_thresholds(
            data.case_series, params["expw_cases_start"], params["expw_cases_end"], params["expw_min_steps"])

        # Determine if the window is valid
        expw_valid = res.metrics["expw_valid"] = None not in (i_t_start, i_t_end) and i_t_start != i_t_end

        # --- Slices the data by the given window
        # Slice the node event dataframe
        if bool_switches["calc_window_event_df"] and expw_valid:
            data.window_node_event_df = data.node_event_df.loc[  # Tested: faster than with groupby/concat
                (data.node_event_df["infec_step"] >= i_t_start) &
                (data.node_event_df["infec_step"] < i_t_end)
            ]

        # Slice the node event for valid-only serial intervals
        if bool_switches["calc_window_serialint_event_df"]:
            df = data.node_event_serialint_df
            data.window_node_event_serialint_df = df.loc[
                (df["infec_step"] >= i_t_start) & (df["infec_step"] < i_t_end)
            ]

            # Slice the time series of new cases
        if bool_switches["calc_window_case_series"]:
            data.window_case_series = data.case_series[i_t_start:i_t_end] if expw_valid else None

        # --- Windowed metrics
        if bool_switches["calc_empir_repnum"]:
            if expw_valid:
                # Calc overall reproduction number
                res.metrics["empir_repnum"] = data.window_node_event_df["num_second_infecs"].mean()

                # Calc reproduction number of symptomatic/asymptomatic individuals
                mask = data.window_node_event_df["so_step"] >= 0
                sympt_df = data.window_node_event_df.loc[mask]
                asymp_df = data.window_node_event_df.loc[~mask]

                res.metrics["empir_repnum_sympt"] = sympt_df["num_second_infecs"].mean()
                res.metrics["empir_repnum_asymp"] = asymp_df["num_second_infecs"].mean()

            else:
                res.metrics["empir_repnum"] = res.metrics["empir_rempnum_sympt"] = res.metrics["empir_rempnum_asymp"] \
                    = None

        if bool_switches["calc_window_gentime"]:
            df = data.window_node_event_df
            res.metrics["window_gentime"] = \
                data.dt * df["own_generation_time"].mean() if expw_valid else None
            res.metrics["window_gentime_h"] = \
                data.dt * df["own_generation_time"].loc[df["infec_layer"] == 0].mean() if expw_valid else None

        if bool_switches["calc_window_serialint"]:
            if expw_valid:
                df = data.window_node_event_serialint_df
                res.metrics["window_serialint"] = data.dt * df["own_serial_time"].mean()
                res.metrics["window_serialint_h"] = data.dt * df["own_serial_time"].loc[df["infec_layer"] == 0].mean()

            else:
                res.metrics["window_serialint"] = res.metrics["window_serialint_h"] = None

        # --- "Shifted" metrics - Times calculated from the start of the exp. window
        if bool_switches["calc_infec_peak"]:
            if expw_valid:
                # # --- Infectious prevalence peak  (unused for now)
                # infec_tseries = data.status_df["num_Ia"] + data.status_df["num_Is"]
                # res.metrics["infec_peak_shifted_t"] = infec_tseries.idxmax() - i_t_start

                # --- Infection incidence peak
                daily_incidence = data.new_e_series.groupby(lambda i_t: i_t // params["steps_gp"]).sum()
                res.metrics["incidence_peak_shifted_t"] = data.new_e_series.idxmax() - i_t_start
                res.metrics["incidence_daily_peak_shifted_t"] = daily_incidence.idxmax() - (i_t_start // params["steps_gp"])  # This one is in daily units already
                del daily_incidence
            else:
                # res.metrics["infec_peak_shifted_t"] = None  # Unused for now
                res.metrics["incidence_peak_shifted_t"] = None
                res.metrics["incidence_daily_peak_shifted_t"] = None
                pass

        # --- Windowed histograms
        if bool_switches["calc_window_tg_hist"] and expw_valid:  # w_tg_hist already None, so expw_valid goes here
            res.w_tg_hist = calc_tg_hist(data.window_node_event_df["own_generation_time"], False, res.max_curve_length)

        if bool_switches["calc_window_R_hist"] and expw_valid:
            df: pd.DataFrame = data.window_node_event_df
            # Use pandas method for simplicity. R is not bounded, so with pandas is easier to concatenate later.
            res.w_r_hist = df["num_second_infecs"].value_counts(sort=False)

        if bool_switches["calc_window_serialint_hist"] and expw_valid:
            df: pd.DataFrame = data.window_node_event_serialint_df
            # Use pandas method for simplicity to deal with negative serial intervals.
            res.w_serint_hist = df["own_serial_time"].value_counts(sort=False)

        # --- Windowed metrics that require an exponential fit
        # Calculate the exponential fit itself
        if bool_switches["calc_exp_fit"] and expw_valid:
            t_slice = np.arange(data.dt * i_t_start, data.dt * i_t_end, data.dt)
            res.metrics["exp_growth"], res.metrics["exp_nt0"], res.metrics["exp_fit_rsquare"] = \
                calc_exp_fit(t_slice, data.case_series[i_t_start:i_t_end])

        elif bool_switches["calc_exp_fit"]:  # This means that expw_valid = False
            res.metrics["exp_growth"], res.metrics["exp_nt0"], res.metrics["exp_fit_rsquare"] = None, None, None

        # Calculate exponential reproduction number
        if bool_switches["calc_exp_repnum"] and expw_valid:
            # tg_hist = calc_tg_hist(data.window_node_event_df["own_generation_time"])
            # res.metrics["exp_repnum"] = calc_expon_repnum_from_tghist(tg_hist, res.metrics["exp_growth"], data.dt)
            res.metrics["exp_repnum"] = calc_expon_repnum_from_tghist(res.w_tg_hist / res.w_tg_hist.sum(),
                                                                      res.metrics["exp_growth"], data.dt,
                                                                      len_hist=res.max_curve_length)

        elif bool_switches["calc_exp_repnum"]:
            res.metrics["exp_repnum"] = None

    # Time series properties
    # ----------------------
    if bool_switches["calc_infec_peak"]:
        # Time step index w/out displacing and w/out multiplying by dt
        infec_tseries = data.status_df["num_Ia"] + data.status_df["num_Is"]
        res.metrics["infec_peak_t"] = infec_tseries.idxmax()
        res.metrics["infec_peak"] = infec_tseries.loc[res.metrics["infec_peak_t"]]

        daily_incidence = data.new_e_series.groupby(lambda i_t: i_t // params["steps_gp"]).sum()
        res.metrics["incidence_peak"] = data.new_e_series.max()
        res.metrics["incidence_peak_t"] = data.new_e_series.idxmax()
        res.metrics["incidence_daily_peak"] = daily_incidence.max()
        del daily_incidence

    if bool_switches["calc_layer_absence"]:
        # Maximum absence in community
        num_absence_series = data.num_nodes - data.presence_df["presence_c"]
        # res.metrics["absence_c"] = 1.0 - data.presence_df["presence_c"].min() / data.num_nodes
        res.metrics["absence_c"] = num_absence_series.max() / data.num_nodes
        res.metrics["absence_c_sum"] = num_absence_series.sum() / data.num_nodes * data.dt

    # Overall (non-windowed) histograms
    # ---------------------------------
    if bool_switches["calc_tg_hist"]:
        if data.tg_hist_df is not None:  # Tg histogram was read from file.
            # All layers
            res.tg_hist = data.tg_hist_df["freq_tg"]

        else:  # File wasn't opened. Calculate from event data.
            # All layers
            res.tg_hist = pd.Series(calc_tg_hist(data.node_event_tg_df["own_generation_time"], normalize=False,
                                                 minlength=res.max_curve_length))

    if bool_switches["calc_R_hist"]:
        # Pandas method, for simplicity with handling unbounded values.
        res.r_hist = data.node_event_df["num_second_infecs"].value_counts(sort=False)

    if bool_switches["calc_serialint_hist"]:
        res.serint_hist = data.node_event_serialint_df["own_serial_time"].value_counts(sort=False)


    # --- Generation time per source symptoms
    if bool_switches["calc_sympt_tg_hist"]:

        res.sympt_tg_hist = pd.Series(calc_tg_hist(data.symptsrc_event_df["own_generation_time"], normalize=False,
                                                   minlength=res.max_curve_length))
        res.asymp_tg_hist = pd.Series(calc_tg_hist(data.asympsrc_event_df["own_generation_time"], normalize=False,
                                                   minlength=res.max_curve_length))

    if bool_switches["calc_tass_hist"]:
        # Per-window hist
        # # Not Implemented

        # Total (all layers)
        res.tass_hist = data.symptsrc_event_df["own_tass"].value_counts(sort=False)

    if bool_switches["calc_presymp_transm_samples"]:
        num_presymp_transm_sr = (
            (data.symptsrc_event_df["own_tass"] < 0)
            .groupby(data.symptsrc_event_df["infec_src"])
            .sum()
        )
        num_secondary_infec_sr = data.symptsrc_event_df.groupby("infec_src").size()
        res.presymp_transm_samples = pd.DataFrame({
            "num_cases": num_secondary_infec_sr,
            "num_presymp_cases": num_presymp_transm_sr,
        })

    # Other metrics
    # -------------
    # --- Total generation time (during the whole outbreak), overall and per layer
    if bool_switches["calc_total_gentime"]:
        try:
            res.metrics["avg_gentime"] = float(data.summary_dict["avg_gentime"])
            res.metrics["avg_gentime_h"] = float(data.summary_dict["avg_gentime_h"])
        except KeyError:
            # To avoid this error, tg could be calculated from `node_event_data`
            sys.stderr.write(f"HEY, execution {i_exec} does not have \"avg_gentime\" pre-calculated on summary.csv.")
            raise

    # --- Generation time split by symptomatic and asymptomatic sources
    if bool_switches["calc_sympt_gentime"]:

        # # -()- Calculate from own symptom status and sum of secondary generation times (does not require TASS column)
        # df = data.node_event_df
        # # Filter for valid tg and splits between symptomatic and asymptomatic
        # mask = df["so_step"] >= 0  # Mask for symptomatic targets
        # df_s = df[mask]
        # df_a = df[~mask]
        #
        # # Calculate Tg based on target's symptom status
        # res.metrics["avg_gentime_sympt"] = df_s["sum_generation_time"].sum() / df_s["num_second_infecs"].sum()
        # res.metrics["avg_gentime_asymp"] = df_a["sum_generation_time"].sum() / df_a["num_second_infecs"].sum()
        # # -\\-

        # # -()- Calculate from own generation time and symptom status of source infector (REQUIRES TASS)
        res.metrics["avg_gentime_sympt"] = data.symptsrc_event_df["own_generation_time"].mean()
        res.metrics["avg_gentime_asymp"] = data.asympsrc_event_df["own_generation_time"].mean()
        # -\\-

        # --- Two-step (TS) Generation time of secondary infections, based on symptoms of source.
        # res.metrics[""]
        for df, metrics_name in [
            (data.symptsrc_event_df, "avg_gentime_ts_sympt"),
            (data.asympsrc_event_df, "avg_gentime_ts_asymp"),
        ]:
            s = df["num_second_infecs"].sum()
            res.metrics[metrics_name] = df["sum_generation_time"].sum() / s if s else None

    # --- Total serial interval (during the whole outbreak)
    if bool_switches["calc_total_serialint"]:
        # --- OVERALL SERIALINTERVAL (ALL LAYERS)
        res.metrics["avg_serialint"] = data.dt * data.node_event_serialint_df["own_serial_time"].mean()

        # --- PER-LAYER SERIALINTERVAL

        # # -()- Groupby method to calculate in all layers
        # gp = data.node_event_serialint_df.groupby("infec_layer")
        # gp_avg = gp["own_serial_time"].mean() * data.dt
        # for li, lc in enumerate(LAYER_CHAR):
        #     res.metrics[f"avg_serialint_{lc}"] = gp_avg[li] if li in gp_avg else None
        # -//-

        # -()- Simple loc method to calculate only in households
        lc, li = "h", 0  # Household layer char and id
        df = data.node_event_serialint_df[data.node_event_serialint_df["infec_layer"] == li]
        res.metrics[f"avg_serialint_{lc}"] = df["own_serial_time"].mean() * data.dt
        # -//-

    if bool_switches["calc_presymp_transm"]:
        df = data.symptsrc_event_df
        res.metrics["presymp_transm"] = df[df["own_tass"] < 0].shape[0] / df.shape[0] \
            if df.shape[0] > 0 else None

    if bool_switches["calc_frac_neg_serialint"]:
        df = data.node_event_serialint_df
        res.metrics["frac_neg_serialint"] = (df["own_serial_time"] < 0).sum() / df.shape[0] \
            if df.shape[0] > 0 else None

    if bool_switches["calc_test_behav_rates"]:
        df = data.node_event_df
        res.metrics["ptest_rate"] = (df["ptest_step"] >= 0).sum() / num_infections
        res.metrics["behav_rate"] = (df["behav_step"] >= 0).sum() / num_infections
        res.metrics["behav_frac"] = res.metrics["behav_rate"] * num_infections / res.num_nodes  # Same but by total pop.

    # --- Household attack rates
    if bool_switches["calc_attack_rates"]:
        warnings.warn("\nATTACK RATE not implemented. Skipped...")

    # --- Latent period of infected and infectors
    if bool_switches["calc_latent_periods"]:
        df: pd.DataFrame = data.node_event_df
        _has_infector_mask = df["infec_src"] >= 0  # Exclude index cases (-1)

        # --- Calculate latent period
        df["latent_period_days"] = (df["transm_start_step"] - df["infec_step"]) * data.dt
        infector_df = df.set_index("i_node").reindex(df.loc[_has_infector_mask, "infec_src"])
        pd.Series().mean()

        res.metrics["avg_latent_period_infected"] = df["latent_period_days"].mean(skipna=True)
        res.metrics["avg_latent_period_infector"] = infector_df["latent_period_days"].mean(skipna=True)

        # --- Calculate the effective infectious period, which depends on latent period
        # and generation time.

        # # -()- Visit each node and check their secondary infections
        # _num = df["sum_generation_time"] * data.dt - df["num_second_infecs"] * df["latent_period_days"]
        # _den = df["num_second_infecs"]
        # res.metrics["avg_effective_infectious_period"] = (_num.sum() / _den.sum())

        # -()- Visit each node and check their infecion source (same result, but also provides the proper ensemble)
        # OBS: use `.values` to avoid alignment with index (not the same)
        own_tg_sr = df.loc[_has_infector_mask, "own_generation_time"] * data.dt
        src_lp_sr = infector_df["latent_period_days"]
        res.metrics["avg_effective_infectious_period"] = (
                own_tg_sr.values - src_lp_sr.values
        ).mean()

        # --- Histogram-ish: Nr. of secondary infections vs infector latent period (binned for better stats)
        # _bins = np.arange(0, df["latent_period_days"].max() + data.dt, 1.0)  # 1-day bins
        res.avg_r_vs_latent_period_df = (
            df["num_second_infecs"]
            # .groupby(pd.cut(df["latent_period_days"], bins=_bins))
            .groupby(np.ceil(df["latent_period_days"]))  # Round up to the next full integer day
            .agg(["mean", "sum", "size"])
            .sort_index()
        )

        # --- Histogram-ish: Nr. of secondary infections vs effective infectious period (binned for better stats)
        eip_sr = data.dt * df["sum_generation_time"] / df["num_second_infecs"] - df["latent_period_days"]  #  nan if no second. infections
        #  ^ ^ This is the avg. effective infectious period of each infector: Avg_tg - LP. Unit: Days.
        # _bins = np.arange(0, _sr.max() + data.dt, 1.0)  # 1-day bins
        res.avg_r_vs_eff_infec_period_df = (
            df["num_second_infecs"]
            # .groupby(pd.cut(df["latent_period_days"], bins=_bins))
            .groupby(np.ceil(eip_sr))  # Round up to the next full integer day
            .agg(["mean", "sum", "size"])
            .sort_index()
        )

    # Finalization
    # ------------
    xtf = time.time()
    res.cpu_time = xtf - xt0

    return res


def calc_execset_statistics(exec_result_list: list[ExecResults], bool_switches, params, exset: ExecSet):

    # ------------------------------------------------------------------------------------------------------------------
    # PREAMBLE
    # ------------------------------------------------------------------------------------------------------------------
    probe_res = exec_result_list[0]  # Take first result to extract parameters that should be constant
    # keys = list(probe_res.metrics.keys())  # Assumes from the first execution
    dt = probe_res.dt
    num_nodes = probe_res.num_nodes
    num_execs = len(exec_result_list)
    max_curve_length = probe_res.max_curve_length
    b = bool_switches  # Alias

    # Create a single df with all result metrics together.
    samp_df = pd.DataFrame([res.metrics for res in exec_result_list])  # Create from a list of dicts
    # samp_df = pd.DataFrame({key: [res.metrics[key] for res in exec_result_list] for key in keys})  # Equivalent

    stat_df = pd.DataFrame()  # Stores the final metrics
    out_samples_df = pd.DataFrame()  # Stores sample sets that are to be exported as is

    # ------------------------------------------------------------------------------------------------------------------
    # METRIC CALCULATION
    # ------------------------------------------------------------------------------------------------------------------

    # Basic values
    # ------------
    stat_df["num_execs"] = [num_execs]
    stat_df["dt"] = [dt]
    stat_df["num_nodes"] = [num_nodes]

    # Metrics that do not require filters
    # -----------------------------------
    # NONE

    # Metrics filtered by valid outbreaks only
    # ----------------------------------------
    outb_mask = samp_df["expw_valid"]  # Having a window by number of cases is the criterion for valid outbreak.
    stat_df["frac_valid_outb"] = [outb_mask.mean()]
    filt_samp_df = samp_df.loc[outb_mask]

    # --- Layer proportions: just the mean is enough
    for lc in LAYER_CHAR:
        # stat_df[f"frac_infecs_{lc}_mean"] = filt_samp_df[f"frac_infecs_{lc}"]
        metric = f"frac_infecs_{lc}"
        calc_metric_stats(filt_samp_df[metric], metric, params, stat_df)

    # --- Metrics divided by N *after* average
    for metric, metric_pos in [
        ("final_num_r", "final_frac_r"), 2*("infec_peak",), 2*("incidence_peak",),
        2*("incidence_daily_peak", )
    ]:
        calc_metric_stats(filt_samp_df[metric], metric_pos, params, stat_df, 1. / num_nodes)
    # calc_metric_stats(filt_samp_df["final_num_r"], "final_frac_r", params, stat_df, 1. / num_nodes)

    # --- Metrics multiplied by dt *after* average
    for metric, metric_pos in [
        ("num_steps", "out_duration"), 2*("avg_gentime_sympt",), 2*("avg_gentime_asymp",),
        2*("avg_gentime_ts_sympt",), 2*("avg_gentime_ts_asymp",), 2*("infec_peak_t",),
        2*("incidence_peak_t",), 2*("incidence_peak_shifted_t",)
    ]:
        if metric in samp_df:
            calc_metric_stats(filt_samp_df[metric], metric_pos, params, stat_df, dt)

    # --- Metrics with no further normalization
    for metric in ["avg_gentime", "avg_serialint", "frac_neg_serialint", "presymp_transm"] + \
                  [f"avg_gentime_{lc}" for lc in LAYER_CHAR] + \
                  [f"avg_serialint_{lc}" for lc in LAYER_CHAR] + ["ptest_rate", "behav_rate", "behav_frac"] + \
                  ["absence_c", "absence_c_sum"] + \
                  ["incidence_daily_peak_shifted_t"] + \
                  ["avg_latent_period_infected", "avg_latent_period_infector", "avg_effective_infectious_period"] \
            :
        if metric in samp_df:
            series: pd.Series = filt_samp_df[metric]
            calc_metric_stats(series, metric, params, stat_df)

    # === Extra (2025-11-12)
    # --- Keep the ensemble of some metrics
    # Metrics filtered by valid outbreak
    for metric, metric_pos, coef in [
        ("infec_peak", "infec_peak", 1. / num_nodes),
        ("infec_peak_t", "infec_peak_t", dt),
        ("final_num_r", "final_frac_r", 1. / num_nodes),
        ("incidence_peak_t", "incidence_peak_t", dt),
        ("incidence_peak_shifted_t", "incidence_peak_shifted_t", dt),
        ("incidence_daily_peak_shifted_t", "incidence_daily_peak_shifted_t", 1.),
        # ("final_frac_r", 1. / num_nodes)
    ]:
        if metric in filt_samp_df:
            out_samples_df[metric_pos] = filt_samp_df[metric] * coef

    # =========

    # Metrics filtered by valid outbreaks AND quality of exponential fit
    # ------------------------------------------------------------------
    exp_mask = (samp_df["exp_fit_rsquare"] > params["exp_fit_min_rsquare"]) & outb_mask
    stat_df["frac_valid_expfit"] = [exp_mask.mean()]

    # --- Metrics with no further normalization
    for metric in ["empir_repnum", "empir_repnum_sympt", "empir_repnum_asymp", "exp_repnum", "window_gentime",
                   "window_gentime_h", "window_serialint", "window_serialint_h", "exp_growth"]:
        if metric in samp_df:
            calc_metric_stats(samp_df[metric][exp_mask], metric, params, stat_df)

    # --- Metrics multiplied by dt after stats
    # metric_pre = name in samp_df. |  metric_pos = name that will be in stat_df
    for metric_pre, metric_pos in \
            zip(["expw_it_start", "expw_it_end"], ["expw_t_start", "expw_t_end"]):
        if metric_pre in samp_df:
            calc_metric_stats(samp_df[metric_pre][exp_mask], metric_pos, params, stat_df, dt)

    # ------------------------------------------------------------------------------------------------------------------
    # HISTOGRAMS
    # ------------------------------------------------------------------------------------------------------------------
    hist_dict: dict[str, pd.DataFrame] = dict()

    # Generation time
    # ---------------
    if b["calc_window_tg_hist"] or b["calc_tg_hist"]:
        # Initialize
        df = hist_dict["tg"] = pd.DataFrame()

        # --- Windowed
        if b["calc_window_tg_hist"]:
            df["window"] = np.zeros(max_curve_length, dtype=np.int64)
            aggregate_exec_hists_fixed(df, "window", (res.w_tg_hist for res in exec_result_list))

        # --- Overall (non-windowed)
        if b["calc_tg_hist"]:
            df["total"] = np.zeros(max_curve_length, dtype=np.int64)
            aggregate_exec_hists_fixed(df, "total", (res.tg_hist for res in exec_result_list))

        # --- For symptomatic and asymptomatic nodes
        if b["calc_sympt_tg_hist"]:
            df["sympt"] = np.zeros(max_curve_length, dtype=np.int64)
            df["asymp"] = np.zeros(max_curve_length, dtype=np.int64)

            aggregate_exec_hists_fixed(df, "sympt", (res.sympt_tg_hist for res in exec_result_list))
            aggregate_exec_hists_fixed(df, "asymp", (res.asymp_tg_hist for res in exec_result_list))

    # Number of secondary cases
    # -------------------------
    if b["calc_window_R_hist"] or b["calc_R_hist"]:
        # Initialize
        # ----------
        df = hist_dict["R"] = pd.DataFrame()

        # --- Windowed
        if b["calc_window_R_hist"]:
            df["window"] = aggregate_exec_hists_variable((res.w_r_hist for res in exec_result_list))

        # --- Overall (non-windowed)
        if b["calc_R_hist"]:
            df["total"] = aggregate_exec_hists_variable((res.r_hist for res in exec_result_list))

        df.sort_index(inplace=True)

    # Serial interval
    # ---------------
    if b["calc_window_serialint_hist"] or b["calc_serialint_hist"]:
        df = hist_dict["serialint"] = pd.DataFrame()

        # --- Windowed
        if b["calc_window_serialint_hist"]:
            df["window"] = aggregate_exec_hists_variable((res.w_serint_hist for res in exec_result_list))

        # --- Overall (non-windowed)
        if b["calc_serialint_hist"]:
            df["total"] = aggregate_exec_hists_variable((res.serint_hist for res in exec_result_list))

        # Final processing
        df.sort_index(inplace=True)

    # Time after source symptoms (TASS)
    # ---------------------------------
    if b["calc_tass_hist"]:
        df = hist_dict["tass"] = pd.DataFrame()

        # --- All layers
        df["total"] = aggregate_exec_hists_variable((res.tass_hist for res in exec_result_list))

        df.sort_index(inplace=True)

    # Number of secondary cases vs relevant periods (latent, eff. infectious)
    # ------------------------------------------------------
    if b["calc_latent_periods"]:
        # --- Latent period
        # AGGREGATION METHOD: Add all sums and sizes, then divide.
        # This means that all individuals accross all executions have
        # the same weight.
        df_list = [res.avg_r_vs_latent_period_df for res in exec_result_list]
        concat_df = pd.concat(df_list, axis=0)
        df = concat_df[["sum", "size"]].groupby(concat_df.index).sum()
        df["mean"] = df["sum"] / df["size"]

        hist_dict["avg_r_vs_latent_period"] = df.copy()  # Keep only the mean column, which is the actual metric

        # --- Secondary infections by effective generation time
        df_list = [res.avg_r_vs_eff_infec_period_df for res in exec_result_list]
        concat_df = pd.concat(df_list, axis=0)
        df = concat_df[["sum", "size"]].groupby(concat_df.index).sum()
        df["mean"] = df["sum"] / df["size"]

        hist_dict["avg_r_vs_eff_infec_period"] = df.copy()  # Keep only the mean column, which is the actual metric


    # Individual pre-symptomatic infection samples
    # --------------------------------------------
    if b["calc_presymp_transm_samples"] and probe_res.presymp_transm_samples is not None:
        df: pd.DataFrame = pd.concat(
            [r.presymp_transm_samples for r in exec_result_list],
            keys=range(len(exec_result_list)),
            names=["i_exec"]
        ).reset_index().drop(columns=["infec_src"])  # To save space, ditch the infector ID

        # --- Save simulation-wide samples in-loco for simplicity
        fname = os.path.join(exset.dir, "presymp_transm_samples.parquet")

        if b["export"]:
            df.to_parquet(fname, index=False)
            print("Exported presymptomatic samples to '{}'".format(fname))

    return stat_df, hist_dict, out_samples_df


# ----------------------------------------------------------------------------------------------------------------------
# AUX FUNCTIONS THAT COULD BE MOVED TO ANOTHER MODULE (bcareful not to conflict with similar functions in analysis_tools
# ----------------------------------------------------------------------------------------------------------------------


def calc_metric_stats(series: pd.Series, metric: str, params, stat_df: pd.DataFrame, coeff=1.0):
    """Calculate the common statistics for a given metric, storing them into stat_df."""
    stat_df[metric + "_mean"] = [series.mean() * coeff]
    stat_df[metric + "_std"] = [series.std() * coeff]
    stat_df[metric + "_quant_low"] = [series.quantile(q=params["quantile_q"]) * coeff]
    stat_df[metric + "_quant_high"] = [series.quantile(q=1. - params["quantile_q"]) * coeff]


def calc_exp_fit(t_array, y_array):
    t0 = t_array[0]

    # Call the fit routine (linear fit of the logarithm)
    p_coefs, *args = np.polyfit(t_array, np.log(y_array), deg=1, full=True)

    # Parameter intepretation
    growth = p_coefs[0]  # r (growth coefficient)
    y_t0 = np.exp(p_coefs[1] + growth * t0)  # fitted y value at t0, start of the slice

    # Goodness of fit: r^2 = Pearson's r
    fit_r, fit_pval = scipy.stats.pearsonr(y_array, y_t0 * np.exp(growth * (t_array - t0)))
    fit_rsquare = fit_r ** 2

    return growth, y_t0, fit_rsquare


def calc_tg_hist(tg_sample: pd.Series, normalize=True, minlength=None):
    """
    Calculate a histogram of integer non-negative values from a pandas series.
    This function was optimized. The conversion to numpy array early on provides much better performance than
    using pandas.
    """
    # -()- Using Numpy method
    tg_array = tg_sample.to_numpy()
    tg_array = tg_array[tg_array > 0]  # Filter non-positive (invalid) values
    hist = np.bincount(tg_array, minlength=minlength)

    return hist / tg_array.shape[0] if normalize else hist


def aggregate_exec_hists_fixed(df, key, hist_seq, normalize=True):
    """
    Sums the frequencies of a set of histograms. Uses simple '+' operation, but only works when the values are
    the same in all histograms.

    OBS: for some unknown reason, the df and key must be passed separately. Passing df[key] as a single argument
    fails in some cases and I don't know why.
    """
    for hist in hist_seq:
        if hist is not None:
            df[key] += hist

    if normalize:
        norm = df[key].sum()
        df[key] /= norm if norm else 1


def aggregate_exec_hists_variable(hist_seq, normalize=True):
    """
    Sums the frequencies of a set of histograms. Uses dataframe .add operation, which generates a new copy every time,
    but is able to include previously non-existent values at each new sum.
    Fixes a previous version in which only the index values of the first-ever histogram would be considered.

    Returns the aggregated histogram as a series, which can be then incorporated into a collection dataframe.
    """

    series = pd.Series(dtype=np.int64)

    for hist in hist_seq:
        if hist is not None:
            series = series.add(hist, fill_value=0)

    if normalize:
        norm = series.sum()
        series /= norm if norm else 1

    return series


if __name__ == "__main__":
    mt0 = time.time()
    main()
    mtf = time.time()
    xt = mtf - mt0
    print("---------------------------")
    print("Total execution time (main) = {} ({:0.2f}s)".format(seconds_to_hhmmss(xt), xt))
