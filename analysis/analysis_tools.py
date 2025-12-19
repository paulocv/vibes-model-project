"""Utilities for post-processing VIBES simulations.
Toolbox module, not to be called directly (use with import statements).
"""

# import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
# from pathos.pools import ProcessPool
import scipy.stats
import statsmodels.api as sm
# import sys
import glob

from toolbox.file_tools import read_config_file, get_bool_from_dict, HEADER_END
from toolbox.plot_tools import make_axes_seq

MIN_EXP_RSQUARE = 0.99  # Exponential fit r-square (warning raised below this value)
INVALID_SERIAL_TIME = -999  # Nodes with this value for the serial time are invalid
ZERO_EPSILON = 1.E-10  # Threshold for considering zero in some float comparisons.
FILE_SECTION_BREAK = HEADER_END

# Conversion between layer character and layer index
LAYER_CHAR = list("hswc")
LAYER_ID = {lc: i for i, lc in enumerate(LAYER_CHAR)}
NUM_LAYERS = len(LAYER_CHAR)
LAYER_LABELS = ["Household", "School", "Workplace", "Community"]
LAYER_LABELS_PLUR = ["Households", "Schools", "Workplaces", "Community"]


# **********************************************************************************************************************
# ** AUX FUNCTIONS
# **********************************************************************************************************************

def cpp_get_bool_from_dict(d, key, raise_error=False, std=False):
    """Returns a boolean read from a string at an input dictionary.
    True if the string is either "True", "true" or "1", false otherwise.
    """
    return get_bool_from_dict(d, key, ["True", "true", "1"], raise_error, std)


# **********************************************************************************************************************
# ** DATA INPUT
# **********************************************************************************************************************

class InputBunch:
    """Class that contains both input parameters and output structures
    from a single VIBES model run (execution).
    """

    # Imports a single dataset
    def __init__(self, exec_dir: str):
        self.exec_dir = exec_dir
        # --- Summary and parameters
        summary_dict = self.summary_dict = read_config_file(os.path.join(exec_dir, "summary.out"))

        # --- Unpacks essential parameters
        self.dt = float(summary_dict["step_dt"])
        self.num_nodes = int(summary_dict["num_nodes"])
        self.num_steps = int(summary_dict["num_steps"])
        self.final_num_removed = int(summary_dict["final_num_removed"])
        self.num_infections_layer = {lc: int(summary_dict["num_infections_" + lc]) for lc in "hswc"}
        self.transm_gamma_shape = float(summary_dict["transm_gamma_shape"])
        self.transm_gamma_scale = float(summary_dict["transm_gamma_scale"])

        self.avg_gentime = float(summary_dict.get("avg_gentime", -1.0))
        self.avg_gentime_h = float(summary_dict.get("avg_gentime_h", -1.0))
        self.avg_gentime_s = float(summary_dict.get("avg_gentime_s", -1.0))
        self.avg_gentime_w = float(summary_dict.get("avg_gentime_w", -1.0))
        self.avg_gentime_c = float(summary_dict.get("avg_gentime_c", -1.0))
        self.num_infections_h = float(summary_dict["num_infections_h"])
        self.num_infections_s = float(summary_dict["num_infections_s"])
        self.num_infections_w = float(summary_dict["num_infections_w"])
        self.num_infections_c = float(summary_dict["num_infections_c"])

        # --- Boolean switches of output files
        self.export_status_series = cpp_get_bool_from_dict(summary_dict, "export_status_series", std=True)
        self.export_event_logs = cpp_get_bool_from_dict(summary_dict, "export_event_logs", std=False)
        self.export_r_and_tg_series = cpp_get_bool_from_dict(summary_dict, "export_r_and_tg_series", std=True)
        self.export_tg_histogram = cpp_get_bool_from_dict(summary_dict, "export_tg_histogram", std=False)
        self.export_node_event_data = cpp_get_bool_from_dict(summary_dict, "export_node_event_data", std=False)
        self.export_symptom_series = cpp_get_bool_from_dict(summary_dict, "export_symptom_series", std=False)
        self.export_sec_attack_rate = cpp_get_bool_from_dict(summary_dict, "export_sec_attack_rate", std=False)

        # --- Import dataframes from other files
        if self.export_r_and_tg_series:
            self.r_tg_df = pd.read_csv(os.path.join(exec_dir, "r_and_tg.csv"), skiprows=1, header=0
                                       ).rename(str.strip, axis=1)
            filter_invalid_r_and_tg_series(self.r_tg_df)
        else:
            self.r_tg_df = None

        self.status_df = pd.read_csv(os.path.join(exec_dir, "status_series.csv"), skiprows=1, header=0
                                     ).rename(str.strip, axis=1) if self.export_status_series else None

        self.tg_hist_df = pd.read_csv(os.path.join(exec_dir, "tg_hist.csv"), skiprows=1, header=0
                                      ).rename(str.strip, axis=1) if self.export_tg_histogram else None

        self.node_event_df = pd.read_csv(os.path.join(exec_dir, "node_event_data.csv"), skiprows=1, header=0
                                         ).rename(str.strip, axis=1) if self.export_node_event_data else None

        if self.export_event_logs:
            self.event_infec_df = pd.read_csv(os.path.join(exec_dir, "event_infec.csv"), skiprows=1,
                                              header=0).rename(str.strip, axis=1)
            self.event_sympt_df = pd.read_csv(os.path.join(exec_dir, "event_sympt_onset.csv"), skiprows=1,
                                              header=0).rename(str.strip, axis=1)
            self.event_test_df = pd.read_csv(os.path.join(exec_dir, "event_test.csv"), skiprows=1,
                                              header=0).rename(str.strip, axis=1)
        else:
            self.event_infec_df = self.event_sympt_df = None

        self.attack_rate_df = pd.read_csv(os.path.join(exec_dir, "attack_rates.csv"), skiprows=1, header=0,
                                          index_col=0).rename(str.strip, axis=1) if self.export_sec_attack_rate else None

        # --- Other containers and parameters
        self.t_array = np.arange(0.0, self.dt * self.num_steps, self.dt)
        self.s_series = self.status_df["num_S"] if self.export_status_series else None
        self.e_series = self.status_df["num_Ea"] + self.status_df["num_Es"] if self.export_status_series else None
        self.i_series = self.status_df["num_Ia"] + self.status_df["num_Is"] if self.export_status_series else None
        self.case_series = self.e_series + self.i_series + self.status_df["num_R"] \
            if self.export_status_series else None
        self.new_e_array = self.s_series.iloc[:-1].to_numpy() - self.s_series.iloc[1:].to_numpy() \
            if self.export_status_series else None


class ExecSet:
    """
    Class that contains data from a VIBES simulation, which is a set of executions.
    Differ from an ExecSeq as the executions are assumed to belong to a single parameter set.
    """
    def __init__(self, sim_dir, label=None, exec_ids=None, exec_prefix="run", exec_fmt="{:05d}"):
        """

        :param sim_dir: string . Directory of the simulation, where execution subdirectories are expected to be.
        :param label: Any. Label of the execution set.
        :param exec_ids: List of integers. If informed, restricts the executions to the given indexes. Otherwise,
            all execution directories matching the exec_prefix are taken.
        :param exec_prefix: string. Prefix of the single execution directory, without trailing underscore.
        :param exec_fmt: Integer numeric format string for the execution index, as used in the execution directories.
        """
        self.dir = sim_dir

        # Define the list of execution directories.
        if exec_ids is None:
            generic_path = os.path.join(sim_dir, exec_prefix + "_*/")
            self.exec_dirs = glob.glob(generic_path)
        else:
            self.exec_dirs = [os.path.join(sim_dir, exec_prefix + '_' + exec_fmt.format(i)) for i in exec_ids]

        # Check paths
        if len(self.exec_dirs) == 0:
            raise ValueError("Hey, no executions matching the given simulation folder \"{}\" and exec. ids were found."
                             "".format(sim_dir))
        for path in self.exec_dirs:
            if not os.path.exists(path):
                raise FileNotFoundError("Hey, execution directory not found: \"{}\"".format(path))

        # Set other parameters
        self.label = label
        self.num = len(self.exec_dirs)


def filter_invalid_r_and_tg_series(r_tg_df: pd.DataFrame):
    """
    Filters the time series of R (reproduction number) and Tg (generation time) for time steps in which there were
    no primary cases. These are marked as -1 in the layer-specific Tg values. Time steps in which all layers have
    invalid values.

    Changes are made in place.
    """
    global_mask = pd.Series(np.repeat(True, r_tg_df.shape[0]))

    for c_layer in LAYER_CHAR:
        tg_topic = f"Tg_{c_layer}"
        r_topic = f"R_{c_layer}"

        # Determine invalid values
        mask = r_tg_df[tg_topic] < ZERO_EPSILON  # True for _invalid_  values
        global_mask = global_mask & mask

        # Replace invalid values
        r_tg_df[tg_topic].loc[mask] = np.nan
        r_tg_df[r_topic].loc[mask] = np.nan

    # Finally: filter overall R(t) series using validity among all layers
    r_tg_df["R"].loc[global_mask] = np.nan
    r_tg_df["Tg"].loc[global_mask] = np.nan


def import_population_as_df(fname, add_zero=True):
    """Load data from a synthetic population file, as used by VIBES."""
    df = pd.read_table(fname, sep=" ", index_col=0,
                       names=["age", "fam_id", "school.ID", "teacher", "university.ID", "work.ID"])
    if add_zero:  # Include extra node at the beginning of the list, as in simulation code
        zero_df = pd.DataFrame({"age": 0, "fam_id": 0, "school.ID": -1,
                                "teacher": -1, "university.ID": -1, "work.ID": -1}, index=[0])
        df = pd.concat([zero_df, df])
    return df


def import_age_data(fname):
    """Load age-dependent data as used by VIBES."""
    return pd.read_csv(fname, index_col=0, header=0).rename(str.strip, axis=1)


class TransmCurveCollection:
    """A class to hold data from a collection of infectiousness curves.
    """
    df: pd.DataFrame  # a.loc[i_curve][i_step] | Infectiousness over time
    sot_sr: pd.Series  # Series of the SOT (Symptom onsed times)
    num_steps: int
    num_curves: int
    dt: float
    is_sympt: bool


def import_transm_curves(fname, is_sympt, normalize_curves=True) \
        -> TransmCurveCollection:
    """Load transmissibility curves from a file, as used by VIBES."""
    obj = TransmCurveCollection()
    obj.is_sympt = is_sympt

    # --- Read the file header
    with open(fname, "r") as fp:
        items = fp.readline().split(",")
        obj.dt = float(items[0])
        obj.num_steps = int(items[1])
        obj.num_curves = int(items[2])

    # --- Import file bulk content
    obj.df = pd.read_csv(fname, skiprows=1)

    # --- Processes the rows
    t_vals = obj.df.pop("Unnamed: 0")[1:]  # Remove row 0 (time stamps)
    obj.df = obj.df.T   # Rows = i_curve + 1; columns = i_t + 1

    obj.sot_sr = obj.df.pop(0)  # Remove column 0 (SOT) and put into new data frame

    # --- Normalize each curve by its sum.
    if normalize_curves:
        obj.df = obj.df.div(obj.df.sum(axis=1), axis=0)

    # --- Rename stuff to simplify
    obj.df.index.name = "curve"
    obj.df.columns.name = "steps"

    return obj


# **********************************************************************************************************************
# ** ANALYSIS METHODS
# **********************************************************************************************************************


def calc_cumulative_average(a):
    """Returns an array whose element i is the average of a from elements 0 to i."""
    # len_a = len(a)
    return a.cumsum() / np.arange(1, a.shape[0] + 1)


def fit_negative_binomial(data):
    """
    Parameter 'data' must be a sample of integer values, expected to be drawn from a negative binomial distribution.
    """
    exog = np.ones_like(data)
    nb_fit_res = sm.NegativeBinomial(data, exog, method="nb2").fit(start_params=[1.0, 1.0],
                                                                   # method="basinhopping"
                                                                   )

    #     print(nb_fit_res.summary())
    #     print(nb_fit_res.params)

    alpha_fit = nb_fit_res.params[1]  # Overdispersion
    mean_fit = np.exp(nb_fit_res.params[0])  # Mean
    p_fit = 1 / (1 + mean_fit * alpha_fit)  # Probability of success
    n_fit = 1 / alpha_fit  # Desired number of successes to stop

    return n_fit, p_fit, mean_fit, alpha_fit


def calc_avg_gentime_from_hist(xd, hist_key="freq_tg"):
    """Calculate the average generation time starting from a histogram of
    individual-level generation times."""
    if xd.tg_hist_df is None:  # No Tg histogram
        raise ValueError("Hey, attempted to calculate average generation time from histogram, but the execution has"
                         " no Tg histogram dataframe.")

    return xd.dt * np.average(xd.tg_hist_df.index.to_numpy(), weights=xd.tg_hist_df[hist_key].to_numpy())


class ExpFitResults:
    """Results from a fitting attempt of an exponential growth period."""

    def __init__(self, growth, n_t0, fit_rsquare, fit_pval, i_t_start, i_t_end, dt):
        self.growth = growth
        self.n_t0 = n_t0
        self.fit_rsquare = fit_rsquare
        self.fit_pval = fit_pval
        self.i_t_start = i_t_start
        self.i_t_end = i_t_end
        self.dt = dt
        self.func = lambda t: n_t0 * np.exp(growth*(t - dt * i_t_start))


def get_time_window_by_thresholds(series, start_val, end_val, min_steps=0):
    """
    Thresholds are expected to be crossed down-to-up.
    """
    # --- Extraction of the points
    i_t_start = np.argmax(series >= start_val)
    i_t_end = max(np.argmax(series > end_val), i_t_start + min_steps)

    # --- Verification of the results
    len_series = len(series)
    if i_t_end > len_series - 1:  # If min_steps, it is possible that i_t_end goes beyond the time series
        i_t_end = len_series - 1  # Simply clamp to the end of the series
        # raise ValueError("Hey, using min_steps = {} for this time series causes i_t_end to be beyond the end of the "
        #                  "series.".format(min_steps))

    if i_t_start == 0:  # Either not found or found at first position
        if series[0] < start_val:  # Series does not cross the threshold.
            i_t_start = None
    if i_t_end == 0:  # Either not found or found at first position
        if series[0] <= end_val:
            i_t_end = None

    return i_t_start, i_t_end


def calc_repnum_windowed(xd: InputBunch, i_t_start, i_t_end):
    """Calculate the reproduction number for secondary infections within
    a time window of a VIBES simulation.
    """

    if xd.node_event_df is None:
        raise ValueError("Hey, node_event_df is required to calculate windowed reproduction number.It is likely that "
                         "export_status_series == false for this execution.")

    # Select and aggregate the infection events by their time step
    gp = xd.node_event_df.groupby("infec_step")
    window_df = pd.concat([gp.get_group(i_t) for i_t in range(i_t_start, i_t_end) if i_t in gp.groups.keys()])

    r_mean = window_df["num_second_infecs"].mean()

    return r_mean


def calc_exp_growth(xd: InputBunch, i_t_start, i_t_end):
    """Calculate the exponential growth rate, assuming that `i_t_start`
     and `i_t_end` delimit the epidemic exponential growth phase.
     """

    if xd.case_series is None:
        raise ValueError("Hey, case_series is required to calculate windowed exponential growth. It is likely that "
                         "export_status_series == false for this execution.")

    # Do the data slicing
    slice_t = xd.t_array[i_t_start:i_t_end]
    slice_n = xd.case_series.to_numpy()[i_t_start:i_t_end]

    t0, tf = xd.dt * i_t_start, xd.dt * i_t_end  # Slice border in days

    # --- Calculation and verification of the linearized fit
    p_coefs, *args = np.polyfit(slice_t, np.log(slice_n), deg=1, full=True)

    # Parameter intepretation
    growth = p_coefs[0]  # r
    n_t0 = np.exp(p_coefs[1] + growth * t0)  # fitted nr. of cases at t0, start of the slice

    def fit(t):
        return n_t0 * np.exp(growth * (t - t0))

    # Goodness of fit: r^2 = Pearson's r
    fit_rsquare, fit_pval = scipy.stats.pearsonr(slice_n, fit(slice_t))
    fit_rsquare = fit_rsquare ** 2

    return ExpFitResults(growth, n_t0, fit_rsquare, fit_pval, i_t_start, i_t_end, xd.dt)


def calc_windowed_tg_histogram(xd: InputBunch, i_t_start, i_t_end):
    """Calculate the generation time histogram within a time window of
    a VIBES simulation.
    """

    if xd.node_event_df is None:
        raise ValueError("Hey, node_event_df is required to calculate windowed generation time hist. It is likely that "
                         "export_node_event_data == false for this execution.")

    # Select and aggregate the infection events by their time step
    gp = xd.node_event_df.groupby("infec_step")
    window_df = pd.concat([gp.get_group(i_t) for i_t in range(i_t_start, i_t_end) if i_t in gp.groups.keys()])

    if "own_generation_time" in window_df:
        # Extracts Tg (in steps) directly from the data
        tg_steps_array = window_df["own_generation_time"].to_numpy()
        tg_steps_array = tg_steps_array[tg_steps_array > 0]
    else:
        # Use some Pandas magic to match sources and calculate the generation time (in steps). Time-consuming.
        # IMPORTANT : This does NOT exclude index infections (infec_step = -1), neither from source, nor from target.
        tg_steps_array = window_df["infec_step"].to_numpy() - np.fromiter(
            (xd.node_event_df["infec_step"].loc[xd.node_event_df["i_node"] == src].to_numpy() for src in
             window_df["infec_src"]), dtype=int)

    # Construct the histogram (relative frequencies)
    num_events = tg_steps_array.shape[0]

    return np.bincount(tg_steps_array) / num_events if num_events > 0 else None


def calc_expon_repnum_from_tghist(tg_hist, growth, dt, len_hist=None):
    """Calculate the exponential reproduction number, i.e., using
    the exponential growth rate (`growth`) and the generation time
    distribution (`tg_hist`).
    As equation (3.6) of Wallinga, Lipsitch 2007 (https://doi.org/10.1098/rspb.2006.3754).
    Assumes a uniform grid for Tg values: arange(0, (dt*len(tg_hist) - 1).
    Assumes the histogram is not empty (no test!).
    """
    rdt = growth * dt
    if len_hist is None:
        try:
            len_hist = tg_hist.shape[0]
        except AttributeError:
            len_hist = len(tg_hist)

    denominator = (tg_hist[1:] * np.exp(-rdt * np.arange(len_hist - 1))).sum()

    return rdt / (1 - np.exp(-rdt)) / denominator


def calc_index_reparray_from_simuset(sim_dir: str):
    """Calculate the number of secondary infections from one index case
    for multiple executions, taking the simulation directory as an input.

    This function locates all execution (run) directories in that simulation
    and calls `calc_index_reparray_from_execdirs` to calcualte R_index.
    """
    # Get all execution directories
    wildcard = "run_*/"
    exec_dirs = glob.glob(os.path.join(sim_dir, wildcard))
    if len(exec_dirs) == 0:
        raise ValueError("Hey, no execution directories matching '{}' wildcard were found.".format(wildcard))

    return calc_index_reparray_from_execdirs(exec_dirs)


def calc_index_reparray_from_execdirs(exec_dirs):
    """
    Returns an array with the numbers of secondary cases from the index case of a set of simulations.
    Parameter exec_dirs is expected to contain valid paths for execution directories.
    """
    num_execs = len(exec_dirs)
    result = np.empty(num_execs, dtype=int)

    for i_exec, exec_dir in enumerate(exec_dirs):
        df = pd.read_csv(os.path.join(exec_dir, "node_event_data.csv"), skiprows=1, header=0).rename(
            str.strip, axis=1)

        # Assumes the first node with infec_steps < 0 as the only index case.
        result[i_exec] = int(df.loc[df["infec_step"] < 0]["num_second_infecs"])

    return result


def calc_quantiles_from_hist(
        hist_df: pd.DataFrame, quantiles: np.ndarray=None):
    """Calculate quantiles of a pmf distribution given as a dataframe.
    Useful for Generation time, Incubation period and TASS histograms.

    The result is a new data frame having:
    - The columns of `hist_df` as indexes;
    - The required quantile q values as columns;
    - The calculated quantiles as entries.
    """
    # --- Assert `quantiles` argument
    if quantiles is None:
        quantiles = np.array([0.025, 0.25, 0.75, 0.975])

    if not isinstance(quantiles, np.ndarray):
        quantiles = np.array(quantiles)

    if (quantiles < 0).any() or (quantiles > 1).any():
        raise ValueError("Hey, quantiles must be between 0 and 1.")

    # ---
    csum_df = hist_df.cumsum(axis=0)

    sr_list = list()
    for i_q, q in enumerate(quantiles):
        sr_list.append((csum_df > q).idxmax())

    res_df = pd.concat(sr_list, keys=quantiles, axis=1)

    return res_df


def filter_serial_intervals(node_event_df: pd.DataFrame):
    """Returns a view of the pandas data frame containing only nodes with valid serial intervals."""
    return node_event_df.loc[
        (node_event_df["own_serial_time"] != INVALID_SERIAL_TIME) & (node_event_df["infec_layer"] >= 0)
    ]


def calc_serial_interval_tseries(xd: InputBunch, node_event_df: pd.DataFrame = None, need_to_filter=True):
    """
    Constructs a data frame with serial interval time series, for each layer and for the total population.
    Besides the input bunch ```xd```, you can inform the node_event_df, which can be previously filtered by the nodes
    that have valid serial intervals.
    """
    if node_event_df is None:
        node_event_df = xd.node_event_df

    if need_to_filter:
        node_event_df = filter_serial_intervals(node_event_df)
        # node_event_df = node_event_df.loc[
        #     (node_event_df["own_serial_time"] != INVALID_SERIAL_TIME) & (node_event_df["infec_layer"] >= 0)
        # ]

    # Separate by layer and time step
    layer_grouped = node_event_df.groupby(["infec_layer", "infec_step"])
    total_grouped = node_event_df.groupby("infec_step")

    layer_ts_series = layer_grouped["own_serial_time"].mean()

    result_df = pd.DataFrame({"Ts_" + lc: xd.dt * layer_ts_series.loc[i] for i, lc in enumerate(LAYER_CHAR)})
    result_df["Ts"] = xd.dt * total_grouped["own_serial_time"].mean()

    return result_df


def calc_serial_interval_stats(
        xd: InputBunch, node_event_df: pd.DataFrame = None,
        need_to_filter=True
):
    """Calculate mean, standard deviation and the fraction of negative
    values of the serial interval, given individual infection pairs from
    `node_event_df`.
    """
    # --- Data checking
    if node_event_df is None:
        node_event_df = xd.node_event_df

    if need_to_filter:
        node_event_df = filter_serial_intervals(node_event_df)

    # --- Calculations
    if node_event_df.shape[0] == 0:  # Empty dataframe
        print("Empty data frame")
        mean = std = negative_frac = None
    else:
        mean = xd.dt * node_event_df["own_serial_time"].mean()
        std = xd.dt * node_event_df["own_serial_time"].std()
        negative_frac = (node_event_df["own_serial_time"] < 0).sum() / node_event_df["own_serial_time"].shape[0]

    return mean, std, negative_frac


# **********************************************************************************************************************
# ** PLOTTING METHODS
# **********************************************************************************************************************


def make_plot_tables(nplots, nrows, ncols, figsize=(16, 9)):
    """Create a sequence of matplotlib figures with subplots, such that
    a total of `nplots` axes fit into pages with `nrows` and `ncols`
    rows and columns.
    """

    # Calculate table parameters
    plots_per_page = nrows * ncols
    num_pages = (nplots - 1) // plots_per_page + 1
    plots_last_page = nplots - (num_pages - 1) * plots_per_page

    nplots_list = (num_pages - 1) * [plots_per_page] + [plots_last_page]

    figs = []
    all_axes = []
    for nplots in nplots_list:
        fig, axes = make_axes_seq(nplots, ncols, figsize[0], figsize[1] / nrows)
        figs.append(fig)
        all_axes += axes

    return figs, all_axes


def plot_exp_growth_fit_from_xdata(ax, xd: InputBunch, res: ExpFitResults):
    """Visualize results from an exponential growth fit, showing the
    incidence (in log scale), the exponenitial fit function, the
    growth ratio (numeric) and the fit's r square statistics.
    """
    return plot_exp_growth_fit(ax, xd.case_series, res, xd.t_array)


def plot_exp_growth_fit(ax, case_series, res: ExpFitResults, t_array=None, put_axis_labels=True):
    """Visualize results from an exponential growth fit, showing the
    incidence (in log scale), the exponenitial fit function, the
    growth ratio (numeric) and the fit's r square statistics.
    """
    # Get preliminary data
    dt = res.dt
    t_array = np.arange(0.0, dt * len(case_series), dt) if t_array is None else t_array

    # Plot curve data
    ax.plot(t_array, case_series, label="Data")  # Number of cases series

    # Plot exponential fit (and extrapolations)
    broader_slice_t = t_array[max(0, res.i_t_start - 45):res.i_t_end + 45]
    ax.plot(broader_slice_t, res.func(broader_slice_t), lw=2, ls="--", label="Exp. fit")

    # Edges of the data slice
    ax.plot(2 * [res.i_t_start * dt], ax.get_ylim(), "k--")  # Slice edges
    ax.plot(2 * [res.i_t_end * dt], ax.get_ylim(), "k--", label="Slice borders")  # Slice edges

    # Additional info
    ax.text(0.6, 0.25, "growth = {:0.4f}\n$R^2$ = {:0.5f}".format(res.growth, res.fit_rsquare),
            transform=ax.transAxes,
            color="black" if res.fit_rsquare > MIN_EXP_RSQUARE else "red"
            )

    # --- Plot setup
    ax.set_yscale("log")
    ax.set_xlim(dt * (res.i_t_start - 100), dt * (res.i_t_end + 300))
    if put_axis_labels:
        ax.set_xlabel("$t$ (days)")
        ax.set_ylabel("Number of cases")
