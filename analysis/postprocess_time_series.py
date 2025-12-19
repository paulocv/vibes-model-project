"""
Calculate time-dependent metrics from one VIBES simulation directory -
this is, a set of independent executions from the same parameters.

Results are exported to a single data frame as csv,
and an InputMarkup summary of scalar metrics.

Steps in the procedures:
- Import and individual execution preprocess (parallelized with Pathos)
- Rearrange all executions into a single hierarchical data frame (each
set of columns is a metric, each column is an
    execution).
- Calculate statistics over executions.
- Export results

Usage
-----
The program only receives one argument, which is the path to the
simulation directory. Example:
```
python analysis/postprocess_time_series.py outputs/tests/quick_test
```

"""

import time
import warnings
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from pathos.pools import ProcessPool
import sys

sys.path.append("./")
from analysis.analysis_tools import filter_invalid_r_and_tg_series, ExecSet, LAYER_CHAR
from toolbox.file_tools import read_config_file, write_config_string


NUM_THREADS = 4
PERIOD_TO_PRINT = 10  # Only print to terminal every this number of executions completed

pd.set_option('display.max_rows', 80)
pd.set_option('display.max_columns', 20)


# ----------------------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
# ----------------------------------------------------------------------------------------------------------------------

def main():

    # PARAMETERS
    # ----------------------------

    bool_switches = {
        "reindex_by_peak": False,  # Has priority over the switches below.
        "reindex_by_cases": True,  # Reindex by first time when cases >= `cases_thres`
        "reindex_by_r1": False,
        "calc_r_tseries": True,
        "calc_tg_tseries": True,
        "calc_newcases_tseries": True,
        "calc_presence_tseries": True,
        "calc_smooth": True,  # Calculate rolling averages in addition to the regular tseries.
        "do_export": True,
        "do_plot": False,
    }

    params = {
        "quantile": 0.025,
        "cases_thres": 500,
        "smooth_steps": round(7 / 0.25),  # For dt = 0.25 days
        "steps_gp": 4,  # Number of time steps that compose one day. Used (independently of dt) to group new cases
    }

    _EXEC_CAP = None  # Set none to calculate all

    df_output_fname = "a_time_series.csv"
    summary_fname = "a_summary_tseries.out"

    # IMPORT AND CALC
    # -------------------------------------------------------------

    try:
        input_path = sys.argv[1]
    except IndexError:
        raise IndexError("Hey, program needs one argument: simulation directory.")

    # Some useful warns
    if not bool_switches["do_export"]:
        print("################################ WARNING: do_export = False")

    exs = ExecSet(input_path)
    pool = ProcessPool(ncpus=NUM_THREADS)

    if _EXEC_CAP:
        print(f"############################### WARNING: _EXEC_CAP is active ({_EXEC_CAP} execs)")

    def task(arg_tuple):
        i_dir, exec_dir = arg_tuple

        # Import command
        data = import_execution_data(exec_dir)

        # Preprocess
        preprocess_execution_data(data, **bool_switches, **params)

        if i_dir % PERIOD_TO_PRINT == 0:
            print(f"{i_dir} ", end="")
        sys.stdout.flush()

        return data

    content = enumerate(exs.exec_dirs[:_EXEC_CAP or len(exs.exec_dirs)])
    print(f"Importing and individual preprocessing... (total = {exs.num})")
    xt0 = time.time()
    data_list: list[ExecData] = pool.map(task, content)
    xtf = time.time()
    print("\n----------------- Time = {:0.2e} s".format(xtf - xt0))

    # --- Gather parameters from the first execution
    data = data_list[0]
    dt = data.dt
    num_nodes = int(data.summary_dict["num_nodes"])
    num_layers = int(data.summary_dict["num_layers"])

    # AGGREGATE AND STATISTICS
    # -------------------------------------------------------------

    print("Arranging (aggregating) executions...")
    full_df = aggregate_executions(data_list, exs.num, **params, **bool_switches)

    print("Calculating statistics...")
    xt0 = time.time()
    stats_df, stats_metrics = calc_stats(full_df, **params, **bool_switches)
    xtf = time.time()
    print("\n----------------- Time = {:0.2e} s".format(xtf - xt0))

    # agg_r_df.hist()  # This creates one for each column. YOu can transpose, but plotting is terribly slow.

    # EXPORT
    # -------------------------------------------------------------
    if bool_switches["do_export"]:
        # --- Main data frame to a csv file
        stats_df.to_csv(os.path.join(input_path, df_output_fname))

        # --- Metrics and processing parameters
        summary_dict = OrderedDict()

        summary_dict["input_path"] = input_path
        summary_dict.update(params)
        summary_dict.update(bool_switches)
        summary_dict["dt"] = dt
        summary_dict["num_nodes"] = num_nodes
        summary_dict["num_layers"] = num_layers

        with open(os.path.join(input_path, summary_fname), "w") as fp:
            fp.write(write_config_string({key: str(val) for key, val in summary_dict.items()}))

    # PLots tests
    # ------------------------

    # ------#
    #
    if bool_switches["do_plot"]:
        # Plot all at once
        metric_prefix_list = ["newE", "R", "Tg"]
        num_axes = len(metric_prefix_list)
        fig, axes = plt.subplots(nrows=num_axes, figsize=(10, 3*num_axes))

        for i_ax, metric_prefix in enumerate(metric_prefix_list):
            ax = axes[i_ax]

            ax.plot(dt * stats_df.index, stats_df[f"{metric_prefix}_mean"])
            ax.fill_between(dt * stats_df.index, stats_df[f"{metric_prefix}_quant_low"],
                            stats_df[f"{metric_prefix}_quant_high"], alpha=0.3)

            if metric_prefix == "R":
                ax.plot([0, 400], [1, 1], "k--")

            ax.set_xlabel("t")
            ax.set_ylabel(f"{metric_prefix}(t)")

            # ax.set_xlim(-280, +200)

        fig.tight_layout()
        plt.show()


# ----------------------------------------------------------------------------------------------------------------------
# CLASSES, BUNCHES
# ----------------------------------------------------------------------------------------------------------------------


class ExecData:
    """
    Class to hold data from a single simulation execution, for specific use in this file.
    There may be classes with the same name in other modules of this project which are _not_ interchangeable.
    """

    # Data containers
    summary_dict: dict = None
    r_tg_df: pd.DataFrame = None
    status_df: pd.DataFrame = None
    sympt_df: pd.DataFrame = None
    # new_e_array: np.ndarray = None
    new_e_series: pd.Series = None
    presence_series: pd.DataFrame = None

    # Execution parameters and metrics
    dt: float = None
    num_nodes: int = None

    i_t_center: int = None


# ----------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------


def import_execution_data(exec_dir):
    """Load data of one VIBES execution from files."""
    data = ExecData()

    # Import data containers
    data.summary_dict = read_config_file(os.path.join(exec_dir, "summary.out"))
    data.r_tg_df = pd.read_csv(os.path.join(exec_dir, "r_and_tg.csv"), skiprows=1, header=0).rename(str.strip, axis=1)
    data.status_df = pd.read_csv(os.path.join(exec_dir, "status_series.csv"), skiprows=1, header=0).rename(
        str.strip, axis=1)
    filter_invalid_r_and_tg_series(data.r_tg_df)
    # data.sympt_df = pd.read_csv(os.path.join(exec_dir, "event_sympt_onset.csv"), skiprows=1, header=0)

    # Optional imports
    try:
        data.presence_series = pd.read_csv(os.path.join(exec_dir, "presence_series.csv"), skiprows=1,
                                           header=0).rename(str.strip, axis=1)
    except FileNotFoundError:
        warnings.warn("PRESENCE SERIES FILE NOT FOUND: '{}'".format(os.path.join(exec_dir, "presence_series.csv")) +
                      "To prevent errors, disable bool_switches['calc_presence_series']")

    # --- Extract important parameters and metrics from summary dictionary
    data.dt = float(data.summary_dict["step_dt"])
    data.num_nodes = int(data.summary_dict["num_nodes"])

    # --- Number of new infections (series)

    # -()- Values method (instead of .to_numpy())
    new_e_array = data.status_df["num_S"].values[:-1] - data.status_df["num_S"].values[1:]
    new_e_array = np.insert(new_e_array, 0, 0)
    # -//-

    data.new_e_series = pd.Series(new_e_array)

    return data


def preprocess_execution_data(
        data: ExecData, reindex_by_peak=False, reindex_by_cases=False,
        reindex_by_r1=False, **kwargs
):
    """Prepare data from one execution for the main calculations.

    This function calculates an "alignment point" for this execution,
    a time point in which all simulations will be aligned (as the "zero"
    reference).
    """

    # Reindexing (time displacement) based on some criterion.
    # -------------------------------------------------------
    will_reindex = True
    if reindex_by_peak:
        # Set the incidence peak as the alignment point with other executions
        # # -()- Simple detection by max()
        # # data.i_t_center = np.argmax(data.new_e_array)
        # data.i_t_center = data.new_e_series.idxmax()

        # -()- Detection using a rolling window (smoothed series)
        smooth_new_e_series = data.new_e_series.rolling(7, center=True).mean()
        data.i_t_center = smooth_new_e_series.idxmax()

    elif reindex_by_cases:
        # Set an incidence threshold crossing as the alignment point with other executions
        cases_series = data.status_df["num_Ia"] + data.status_df["num_Is"] + data.status_df["num_R"]
        data.i_t_center = (cases_series >= kwargs["cases_thres"]).argmax()

    else:
        will_reindex = False

    if will_reindex:
        # --- Data time displacement
        data.r_tg_df.index -= data.i_t_center
        data.status_df.index -= data.i_t_center
        # data.sympt_df.index -= data.i_t_center
        data.new_e_series.index -= data.i_t_center

        if data.presence_series is not None:
            data.presence_series.index -= data.i_t_center


def _aggregate_series_to_df(var_name, data_list: list[ExecData], df_list: list[pd.DataFrame], colname_list):
    """Auxiliary function: append a data frame and its variable to a registry for concatenation."""
    df_list.append(pd.concat([data.r_tg_df[var_name] for data in data_list], axis=1, ignore_index=True))
    colname_list.append(var_name)


def aggregate_executions(data_list: list[ExecData], num_execs=None, calc_r_tseries=False, calc_tg_tseries=False,
                         calc_newcases_tseries=False, calc_presence_tseries=False, steps_gp=4, **kwargs):
    """Join data from all executions into a single dataframe with hierarchical columns."""
    if num_execs is None:
        num_execs = len(data_list)

    df_list = []
    colname_list = []

    # Signature of each df: rows = i_t, cols = i_exec
    if calc_r_tseries:
        _aggregate_series_to_df("R", data_list, df_list, colname_list)
        _aggregate_series_to_df("R_h", data_list, df_list, colname_list)
        _aggregate_series_to_df("R_s", data_list, df_list, colname_list)
        _aggregate_series_to_df("R_w", data_list, df_list, colname_list)
        _aggregate_series_to_df("R_c", data_list, df_list, colname_list)

    if calc_tg_tseries:
        _aggregate_series_to_df("Tg", data_list, df_list, colname_list)
        _aggregate_series_to_df("Tg_h", data_list, df_list, colname_list)
        _aggregate_series_to_df("Tg_s", data_list, df_list, colname_list)
        _aggregate_series_to_df("Tg_w", data_list, df_list, colname_list)
        _aggregate_series_to_df("Tg_c", data_list, df_list, colname_list)

    if calc_newcases_tseries:
        # New cases per time step
        df_list.append(pd.concat([data.new_e_series for data in data_list], axis=1, ignore_index=True))
        colname_list.append("newE")

        # New cases per "day" (a set of time steps)
        df_list.append(pd.concat([
            data.new_e_series.groupby(data.new_e_series.index // steps_gp).sum()  # Series of grouped steps
            for data in data_list], axis=1, ignore_index=True))
        colname_list.append("newE_days")

    if calc_presence_tseries:
        df_list.append(pd.concat([1. - data.presence_series["presence_c"] / data.num_nodes for data in data_list],
                                 axis=1, ignore_index=True))
        colname_list.append("absence_c")

    # --- Further aggregate the data frames of all metrics into one big hierarchical data frame
    return pd.concat(df_list, keys=colname_list, axis=1)


def calc_stats(full_df: pd.DataFrame, calc_r_tseries=False, calc_tg_tseries=False, calc_newcases_tseries=False,
               calc_smooth=False, calc_presence_tseries=False, **kwargs):
    """Perform main calculations over the ensemble of simulation time series.
    Calculations include time series of the reproduction number,
    incidence, generation time, etc., as well as moving (rolling) averages
    of some metrics.
    """
    metrics = dict()
    stats_df = pd.DataFrame(index=full_df.index)

    # REGULAR AVERAGE TIME SERIES
    # ---------------------------

    # --- New cases series
    if calc_newcases_tseries:
        print("New cases...")
        df = full_df["newE"]
        stats_df["newE_mean"] = df.mean(axis=1)
        stats_df["newE_std"] = df.std(axis=1)
        stats_df["newE_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
        stats_df["newE_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

        print("New cases in days...")
        df = full_df["newE_days"]
        stats_df["newE_days_mean"] = df.mean(axis=1)
        stats_df["newE_days_std"] = df.std(axis=1)
        stats_df["newE_days_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
        stats_df["newE_days_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

    # --- Overall reproduction number
    if calc_r_tseries:
        print("R(t)...")
        df = full_df["R"]
        stats_df["R_mean"] = df.mean(axis=1)
        stats_df["R_std"] = df.std(axis=1)
        stats_df["R_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
        stats_df["R_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

    # --- Per-layer reproduction number
    if calc_r_tseries:
        print("R(t) per layer...")
        for lc in LAYER_CHAR:
            df = full_df[f"R_{lc}"]
            stats_df[f"R_{lc}_mean"] = df.mean(axis=1)
            stats_df[f"R_{lc}_std"] = df.std(axis=1)
            stats_df[f"R_{lc}_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
            stats_df[f"R_{lc}_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

    # --- Overall generation time
    if calc_tg_tseries:
        print("Tg(t)...")
        df = full_df["Tg"]
        stats_df["Tg_mean"] = df.mean(axis=1)
        stats_df["Tg_std"] = df.std(axis=1)
        stats_df["Tg_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
        stats_df["Tg_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

        # --- Per-layer genration time
        print("Tg(t) per layer...")
        for lc in LAYER_CHAR:
            df = full_df[f"Tg_{lc}"]
            stats_df[f"Tg_{lc}_mean"] = df.mean(axis=1)
            stats_df[f"Tg_{lc}_std"] = df.std(axis=1)
            stats_df[f"Tg_{lc}_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
            stats_df[f"Tg_{lc}_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

    # --- Community absence (isolation or isolated fraction)
    if calc_presence_tseries:
        print("absence_c(t)...")
        df = full_df["absence_c"]
        stats_df["absence_c_mean"] = df.mean(axis=1)
        stats_df["absence_c_std"] = df.std(axis=1)
        stats_df["absence_c_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
        stats_df["absence_c_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

    # TIME SERIES OF ROLLING AVERAGES
    # -------------------------------
    if calc_smooth:
        rolling = full_df.rolling(kwargs["smooth_steps"], center=False)  # Point i_t is the right edge of window.

        # --- New cases series
        if calc_newcases_tseries:
            print("New cases smooth...")
            df = rolling["newE"].mean()
            stats_df["smooth_newE_mean"] = df.mean(axis=1)
            stats_df["smooth_newE_std"] = df.std(axis=1)
            stats_df["smooth_newE_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
            stats_df["smooth_newE_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

            print("New cases smooth in days...")
            df = rolling["newE_days"].mean()
            stats_df["smooth_newE_days_mean"] = df.mean(axis=1)
            stats_df["smooth_newE_days_std"] = df.std(axis=1)
            stats_df["smooth_newE_days_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
            stats_df["smooth_newE_days_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

        if calc_r_tseries:
            # --- Overall reproduction number
            print("R(t) smooth...")
            df = rolling["R"].mean()
            stats_df["smooth_R_mean"] = df.mean(axis=1)
            stats_df["smooth_R_std"] = df.std(axis=1)
            stats_df["smooth_R_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
            stats_df["smooth_R_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

            # --- Per-layer reproduction number
            print("R(t) per layer smooth...")
            for lc in LAYER_CHAR:
                df = rolling[f"R_{lc}"].mean()
                stats_df[f"smooth_R_{lc}_mean"] = df.mean(axis=1)
                stats_df[f"smooth_R_{lc}_std"] = df.std(axis=1)
                stats_df[f"smooth_R_{lc}_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
                stats_df[f"smooth_R_{lc}_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

        if calc_tg_tseries:
            # --- Overall generation time
            print("Tg(t) smooth...")
            df = rolling["Tg"].mean()
            stats_df["smooth_Tg_mean"] = df.mean(axis=1)
            stats_df["smooth_Tg_std"] = df.std(axis=1)
            stats_df["smooth_Tg_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
            stats_df["smooth_Tg_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

            print("Tg(t) per layer smooth...")
            for lc in LAYER_CHAR:
                df = rolling[f"Tg_{lc}"].mean()
                stats_df[f"smooth_Tg_{lc}_mean"] = df.mean(axis=1)
                stats_df[f"smooth_Tg_{lc}_std"] = df.std(axis=1)
                stats_df[f"smooth_Tg_{lc}_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
                stats_df[f"smooth_Tg_{lc}_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

            # # --- Household generation time
            # print("Tg(t) per layer smooth...")
            # df = rolling["Tg_h"].mean()
            # stats_df["smooth_Tg_h_mean"] = df.mean(axis=1)
            # stats_df["smooth_Tg_h_std"] = df.std(axis=1)
            # stats_df["smooth_Tg_h_quant_low"] = df.quantile(kwargs["quantile"], axis=1)
            # stats_df["smooth_Tg_h_quant_high"] = df.quantile(1. - kwargs["quantile"], axis=1)

    # Etc
    stats_df.index.name = "i_t"
    print(stats_df)

    return stats_df, metrics


if __name__ == "__main__":
    main()
