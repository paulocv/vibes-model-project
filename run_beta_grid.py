"""
Calls the VIBES simulator with multiple values of the
relative infection risk `beta`.

First argument [mandatory]: input file
Second argument [mandatory]: simulation output directory
`-t` or `--num-threads` [optional]: number of parallel processes to use (default: 4)
Any remaining arguments are forwarded to the VIBES engine program and must be recognized there.
"""

import argparse
from datetime import datetime
import os
import random
import sys
import time

import numpy as np
import pandas as pd
from pathos.pools import ProcessPool

from toolbox.file_tools import read_config_file, write_config_string, seconds_to_hhmmss


class Params:
    std_num_threads = 4
    """Standard number of processes to use for parallel execution. The loop of simulations
    is parallelized over this number of processes."""
    
    i0_sim = 0
    """Integer index of the first simulation."""

    master_beta_list = np.array(
        [0.35, 0.50, 0.67, 0.85, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
    )
    """Default grid of beta values to simulate."""

    vibes_executable = "./main"
    # seed_leap = 1 # If not defined, it will be set according to number of executions.
    do_remove_tmp = True
    """If True, temporary input files created for each simulation are removed after execution."""


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process a sequence of execution sets.")

    parser.add_argument("main_dir")
    parser.add_argument(
        "-t", "--num-threads", default=Params.std_num_threads, type=int,
        help=f"Number of threads to use for this run (default = {Params.std_num_threads})."
    )

    return parser.parse_known_args()


def run_simulation_set(args, other_args):
    """"""
    # --- Handle some parameters
    num_sims = len(Params.master_beta_list)
    if len(sys.argv) < 3:
        raise ValueError(
            "Hey, this program needs to be called with 2 arguments "
            "(input file and output directory)"
        )
    main_input_fname = sys.argv[1]
    main_output_dir = sys.argv[2]

    # --- Read Input file
    main_input_dict = read_config_file(main_input_fname)

    # --- Determine seed leap Find the "-n" argument in program call (number of executions), if present.
    i_numexec = None
    try:
        i_numexec = sys.argv.index("-n")  # Argument overrides input metrics
    except ValueError:
        pass

    if not hasattr(Params, "seed_leap"):
        if i_numexec is None:
            Params.seed_leap = int(main_input_dict["num_execs"])  # Num exec from input file
        else:
            Params.seed_leap = int(sys.argv[i_numexec + 1])  # Num exec from program call

    # Preparation
    # -----------
    input_dict = main_input_dict.copy()
    tmp_input_files = []
    sim_dir_names = []
    for i_sim, beta_coef in enumerate(Params.master_beta_list):
        # Modify all layer betas
        for name in ["beta_h", "beta_s", "beta_w", "beta_c"]:
            input_dict[name] = str(float(main_input_dict[name]) * beta_coef)
        input_dict["original_input_file"] = main_input_fname

        # Modify seed with a seed leap between simulations.
        input_dict["rnd_seed"] = str(int(input_dict["rnd_seed"]) + Params.seed_leap)

        # Make simulation folder and tmp input paths
        i_sim_str = "{:04d}".format(Params.i0_sim + i_sim)
        tmp_input_fname = os.path.join("tmp_inputs", str(random.randint(0, 10000)) + "tmp_sim" + i_sim_str + ".in")

        sim_dir_name = "sim_" + i_sim_str
        sim_dir_names.append(sim_dir_name)
        sim_output_dir = os.path.join(main_output_dir, sim_dir_name)

        # Export new dictionary to tmp input
        os.makedirs(os.path.dirname(tmp_input_fname), exist_ok=True)
        with open(tmp_input_fname, "w") as fp:
            fp.write(write_config_string(input_dict))
        tmp_input_files.append(tmp_input_fname)

    # --- Data frame with information about each simulation
    index_df = pd.DataFrame({
        "beta_coef": Params.master_beta_list,
        "sim_dir": sim_dir_names,
    })

    # Execution
    # ---------
    def task(input_tup):
        # Argument unpacking
        input_fname, sim_dir_name, i_sim = input_tup
        output_path = os.path.join(main_output_dir, sim_dir_name)  # Constructs the full path
        cmd = "{} {} {}".format(Params.vibes_executable, input_fname, output_path) + " " + " ".join(sys.argv[3:])

        print("---------STARTING SIMULATION {} of {} -----".format(i_sim, num_sims))
        print(cmd)
        os.system(cmd)
        print("---------FINISHED SIMULATION {} of {} -----".format(i_sim, num_sims))

    content = list(zip(
        tmp_input_files,
        sim_dir_names,
        range(num_sims),
    ))

    pool = ProcessPool(ncpus=int(args.num_threads))
    result_list = pool.map(task, content)

    # Post-execution / exports
    # ------------------------
    index_df.to_csv(os.path.join(main_output_dir, "index.csv"))

    if Params.do_remove_tmp:
        for fname in tmp_input_files:
            os.remove(fname)


def main():
    args, other_args = parse_arguments()
    run_simulation_set(args, other_args)

    print(args)
    print(other_args)


if __name__ == "__main__":
    xt0 = time.time()
    main()
    xtf = time.time()
    print("--------------------")
    print("Finished at {}".format(datetime.now()))
    print("TOTAL EXECUTION TIME (python script): {} ({}s).".format(seconds_to_hhmmss(xtf - xt0), xtf - xt0))
