"""
Generate the transmissibility curves of the Geometric Parametrization (the traditional SEIR-like model).

Each curve is a rectangular pulse: transmissibility is 1 from the end of the latent period (E -> I) to the end of the
infectious period (I -> R), and 0 otherwise. Both periods are geometrically distributed (Markovian).
Symptom onset is placed a fixed number of steps after the E -> I transition.

Baseline settings (paper):
    dt  = 0.25 d
    eta = dt / 3.75   (E -> I probability per step; mean latent period of 3.75 d)
    mu  = dt / 3.50   (I -> R probability per step; the "gamma" of the paper)
    seed = 12

Usage:
    python transmissibility/generate_geometric_curves.py [OUT_DIR]
OUT_DIR defaults to transmissibility/geometric (relative to the working directory, i.e. run from the repo root).
Writes OUT_DIR/v01_sympt.csv and OUT_DIR/v01_asymp.csv.
"""

import argparse
import os
import time

import numpy as np

MAX_ATTEMPTS = 1000

DT = 0.25  # Time step (days)
NUM_CURVES = 10000
LATENT_DAYS = 3.75  # Mean latent period (days), eta = dt / LATENT_DAYS
INFECTIOUS_DAYS = 3.5  # Mean infectious period (days), mu = dt / INFECTIOUS_DAYS
SOT_BEFORE_ONSET_DAYS = 2.25  # Time (days) between the start of the infectious period and symptom onset
MAX_DURATION_DAYS = 40  # Maximum duration of transmissibility curves
DEFAULT_SEED = 12
DEFAULT_OUT_DIR = "transmissibility/geometric"  # Relative to the working directory (run from the repo root)


class Params:
    def __init__(self, latent_days=LATENT_DAYS, infectious_days=INFECTIOUS_DAYS, dt=DT, num_curves=NUM_CURVES):
        self.num_curves = num_curves
        self.dt = dt
        self.eta = dt / latent_days  # E -> I probability at each step
        self.mu = dt / infectious_days  # I -> R probability at each step
        self.backw_latent_period = round(SOT_BEFORE_ONSET_DAYS / dt)  # steps between E -> I and symptom onset
        self.num_steps = round(MAX_DURATION_DAYS / dt) + 1  # Maximum duration of transmissibility curves (in steps)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir", nargs="?", default=DEFAULT_OUT_DIR,
                        help="Output directory, created if needed (default: %(default)s).")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="RNG seed (default: %(default)s).")
    parser.add_argument("--latent-days", type=float, default=LATENT_DAYS,
                        help="Mean latent period in days (default: %(default)s).")
    parser.add_argument("--infectious-days", type=float, default=INFECTIOUS_DAYS,
                        help="Mean infectious period in days (default: %(default)s).")
    args = parser.parse_args()

    p = Params(args.latent_days, args.infectious_days)
    rng = np.random.default_rng(seed=args.seed)

    xt0 = time.time()
    # The symptomatic curves must be drawn first to match the baseline random stream.
    sympt_curves, sympt_sot = produce_valid_curves(rng, p)
    asymp_curves, _ = produce_valid_curves(rng, p)
    asymp_sot = - (1. / p.dt) * np.ones(p.num_curves, dtype=int)  # Convention for asymptomatic: onset at t = -1 day

    print("Exporting...")
    os.makedirs(args.out_dir, exist_ok=True)
    export_transm_curves(os.path.join(args.out_dir, "v01_sympt.csv"), sympt_curves, sympt_sot, p)
    export_transm_curves(os.path.join(args.out_dir, "v01_asymp.csv"), asymp_curves, asymp_sot, p)

    print("Total exec. time = {:0.4e} s".format(time.time() - xt0))


def produce_valid_curves(rng, p: Params):
    """Draws event times and redraws the invalid ones (curves that do not end within num_steps)."""

    te_array, ti_array, sot_array = generate_geometric_event_times(rng, p.num_curves, p.eta, p.mu, p.backw_latent_period)

    for i_try in range(MAX_ATTEMPTS):
        invalid_mask = np.logical_not(validate_times(te_array, ti_array, p.num_steps))
        num_invalid = invalid_mask.sum()
        if num_invalid == 0:
            break

        te_array[invalid_mask], ti_array[invalid_mask], sot_array[invalid_mask] = \
            generate_geometric_event_times(rng, num_invalid, p.eta, p.mu, p.backw_latent_period)

    else:
        raise RuntimeError(f"MAX_ATTEMPTS = {MAX_ATTEMPTS} reached. No outputs produced.")

    print("Number of retries   = {:d}".format(i_try))
    return generate_curves_from_times(te_array, ti_array, p.num_steps), sot_array


def generate_geometric_event_times(rng, num_curves: int, eta: float, mu: float, fwd_sot_steps: int):
    """
    Generates geometrically distributed latent (te) and infectious (ti) periods, in steps.
    Symptom onset (sot) is te + fwd_sot_steps.
    """
    te_array = rng.geometric(eta, num_curves)
    ti_array = 1 + rng.geometric(mu, num_curves)
    sot_array = te_array + fwd_sot_steps

    return te_array, ti_array, sot_array


def validate_times(te_array, ti_array, num_steps):
    """Returns a mask with True for the valid entries."""
    return te_array + ti_array < num_steps


def generate_curves_from_times(te_array, ti_array, num_steps):
    """Boolean array (num_steps, num_curves): True while the individual is infectious."""
    curves = np.zeros((num_steps, te_array.shape[0]), dtype=bool)
    for i_t in range(num_steps):
        np.logical_xor(i_t >= te_array, i_t >= te_array + ti_array, out=curves[i_t, :])

    return curves


def export_transm_curves(fname, curves, sot_array, p: Params):
    with open(fname, "w") as fp:

        # First line – file header
        fp.write("{:0.4f}, {:d}, {:d}\n".format(p.dt, p.num_steps, p.num_curves))

        # Second line – curve indexes
        fp.write(",")
        np.savetxt(fp, [1 + np.arange(p.num_curves)], fmt="%d", delimiter=',')

        # Third line – symptom onsets
        fp.write(",")
        np.savetxt(fp, [p.dt * sot_array], fmt="%0.5f", delimiter=',')

        # File content
        for i_step in range(p.num_steps):
            fp.write("{:f},".format(p.dt * i_step))
            np.savetxt(fp, [curves[i_step]], fmt="%d", delimiter=",")


if __name__ == "__main__":
    main()
