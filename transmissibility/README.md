# Transmissibility curve files

This folder contains the transmissibility (infectiousness) curves used by the simulations. Each file holds one curve per column. When a node is infected, the model assigns it one curve from the file (symptomatic or asymptomatic, depending on the node) and follows it over time.

Each folder has one file for symptomatic nodes and one for asymptomatic nodes. Their paths are set in the input file by `transm_curves_sympt_path` and `transm_curves_asymp_path` (see [`input_parameters.md`](../vibes-core/docs/input_parameters.md)).

## File format

Plain-text CSV (comma-separated, no quoting). A file with `num_steps` time steps and `num_curves` curves has `3 + num_steps` lines and `1 + num_curves` columns.

```
0.2500, 161, 10000                  <- line 1: header (dt, num_steps, num_curves)
,1,2,3,...,10000                    <- line 2: curve indexes
,9.75,12.5,...,3.25                 <- line 3: symptom onset time of each curve (days)
0.000000,0,0,0,...,0                <- line 4: step 0, t = 0 d
0.250000,0,0,1,...,0                <- line 5: step 1, t = 0.25 d
...
40.000000,0,0,0,...,0               <- last line: step num_steps - 1
```

The first column of lines 2 to 3 is empty (or holds a label such as `time`) and is ignored. Column 1 of the data lines holds the time stamps. Curve `j` (1-based) is in column `j + 1`.

### Line 1: header

Three comma-separated numbers:

| Field | Type | Meaning |
|---|---|---|
| `dt` | float | Length of one time step, in days (e.g. `0.25`). |
| `num_steps` | int | Number of time steps in each curve, i.e. the number of data lines. |
| `num_curves` | int | Number of curves, i.e. the number of data columns (not counting the time column). |

The simulation reads these values to size its arrays and to convert times in days into steps. The file is rejected if `dt` is not positive. It is also rejected if `num_steps` exceeds `MAX_TRANSM_CURVE_DURATION`, or if any line is shorter than the header announces.

### Line 2: curve indexes

The 1-based index of each curve (`1, 2, ..., num_curves`), with an empty first cell. This line is only a human-readable aid, and the importer skips it. Because it is skipped, nothing checks that the indexes are correct, so keep it consistent with the columns.

### Line 3: symptom onset times (SOT)

One value per curve, in **days after infection** (not in steps), with an empty first cell.

- A **negative** value (e.g. `-1`) means that the curve is asymptomatic and has no symptom onset.
- Values may exceed the duration of the curve. The code handles this case.
- In the symptomatic file, the SOT of each curve is its own. In the asymptomatic file, every entry is `-1`.

### Lines 4 onward: curve values

Line `4 + i` (for `i = 0, ..., num_steps - 1`) holds the values at step `i`.

- **Column 1 (time stamp):** the time of step `i` in days, equal to `i * dt`. The importer ignores it, so the real time axis is always `i * dt` from the header.
- **Columns 2 to `num_curves + 1`:** the transmissibility of each curve at that step.

Each curve starts at the moment of infection, which is `t = 0` in the first data line. So the value at step `i` is the transmissibility of a node `i * dt` days after it was infected. Values must be non-negative.

If `normalize_transm_curves` is true in the input file, the simulation divides each curve by its sum before use. Each curve then has a total weight of 1. This happens at load time, so the files themselves keep their original values
