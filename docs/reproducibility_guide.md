<a id="first-paper"></a>
# First paper: "VIBES: A multiscale modeling approach integrating within-host and between-hosts dynamics in epidemics"

This guide shows how to re-run and postprocess simulations published in the first paper (PNAS), [Ventura et al, 2026](https://doi.org/10.1073/pnas.2523055123).

## Running and processing main analysis simulations

In the paper, we ran a beta grid for the VIBES model without interventions and with 100% isolation of symptomatic infections, both with 50 independent executions for each value of beta. To reproduce these simulations, run the following commands on a terminal open in the root of this repository (this may take several minutes to hours and create ~50 gigabytes in files):
```bash
python run_beta_grid.py inputs/model_paper_in/vl_model_baseline.in outputs/model_paper_out/vl_model_baseline -n 50
python run_beta_grid.py inputs/model_paper_in/vl_model_symp-isol.in outputs/model_paper_out/vl_model_symp-isol -n 50
```

This should create the `outputs/model_paper_out` directory and populate it with results from the simulations. To postprocess these simulations and obtain the data files that were used for the main text Figures 2 and 3, run:

```bash
python analysis/postprocess_beta_grid.py outputs/model_paper_out/vl_model_baseline
python analysis/postprocess_beta_grid.py outputs/model_paper_out/vl_model_symp-isol
```

Which should result in additional files inside each of the output directories. The `seqsets_analysis.csv` contain most of the data needed to reproduce Figures 2 and 3 (see below for full details on the output files). Additionally, Supplementary Figures S3 and S4 use data from `outputs/model_paper_out/vl_model_baseline/sim_0001/run_00000`. 

> Notice: Since different computer systems implement different random number generation mechanisms, results should be similar but numerically different from the ones presented in our paper.

## Running and processing sensitivity analysis simulations

Run the following additional simulations (it may take several hours and require additional ~250 gigabytes of permanent storage):

```bash
# Distribution of individual infectiousness
python run_beta_grid.py inputs/model_paper_in/sensitivity_in/no-gamma_in/no-gamma.in outputs/model_paper_out/sensitivity/no-gamma -n 50

# Probability of developing symptoms
python run_beta_grid.py inputs/model_paper_in/sensitivity_in/probab-sympt_in/x0p5_no-isol.in outputs/model_paper_out/sensitivity/probab-sympt/x0p5_no-isol -n 50
python run_beta_grid.py inputs/model_paper_in/sensitivity_in/probab-sympt_in/x1p5_no-isol.in outputs/model_paper_out/sensitivity/probab-sympt/x1p5_no-isol -n 50

# Infectiousness of symptomatic individuals relative to asymptomatic ones
python run_beta_grid.py inputs/model_paper_in/sensitivity_in/sympt_transm_coef_in/x2p0_no-isol_readjusted.in outputs/model_paper_out/sensitivity/sympt_transm_coef/x2p0_no-isol_readjusted -n 50
python run_beta_grid.py inputs/model_paper_in/sensitivity_in/sympt_transm_coef_in/x3p0_no-isol_readjusted.in outputs/model_paper_out/sensitivity/sympt_transm_coef/x3p0_no-isol_readjusted -n 50
python run_beta_grid.py inputs/model_paper_in/sensitivity_in/sympt_transm_coef_in/x4p0_no-isol_readjusted.in outputs/model_paper_out/sensitivity/sympt_transm_coef/x4p0_no-isol_readjusted -n 50
python run_beta_grid.py inputs/model_paper_in/sensitivity_in/sympt_transm_coef_in/x5p0_no-isol_readjusted.in outputs/model_paper_out/sensitivity/sympt_transm_coef/x5p0_no-isol_readjusted -n 50
```

Postprocess the files by running `postprocess_beta_grid.py` with all outputs:
```bash
# Distribution of individual infectiousness
python analysis/postprocess_beta_grid.py outputs/model_paper_out/sensitivity/no-gamma

# Probability of developing symptoms
python analysis/postprocess_beta_grid.py outputs/model_paper_out/sensitivity/probab-sympt/x0p5_no-isol
python analysis/postprocess_beta_grid.py outputs/model_paper_out/sensitivity/probab-sympt/x1p5_no-isol

# Infectiousness of symptomatic individuals relative to asymptomatic ones
python analysis/postprocess_beta_grid.py outputs/model_paper_out/sensitivity/sympt_transm_coef/x2p0_no-isol_readjusted
python analysis/postprocess_beta_grid.py outputs/model_paper_out/sensitivity/sympt_transm_coef/x3p0_no-isol_readjusted
python analysis/postprocess_beta_grid.py outputs/model_paper_out/sensitivity/sympt_transm_coef/x4p0_no-isol_readjusted
python analysis/postprocess_beta_grid.py outputs/model_paper_out/sensitivity/sympt_transm_coef/x5p0_no-isol_readjusted
```

The files produced from the above postprocessing, along with the main analysis files, should provide all data needed to reproduce Supplementary Figures S6-S8. Finally, to calculate time-dependent metrics and reproduce Supplementary Figure S5, you must run:

```bash
python analysis/postprocess_time_series.py outputs/model_paper_out/vl_model_baseline/sim_0001/
python analysis/postprocess_time_series.py outputs/model_paper_out/vl_model_baseline/sim_0002/
python analysis/postprocess_time_series.py outputs/model_paper_out/vl_model_baseline/sim_0004/
```


---
<a id="second-paper"></a>
# Second paper: "The impact of the distributions of the latent and infectious periods on epidemic dynamics"

This guide shows how to re-run and postprocess simulations published in the second paper, currently available as a preprint ([Ventura et al, 2026](https://doi.org/10.21203/rs.3.rs-11154983/v1)) and submitted for publication. The paper compares two parameterizations of the same model: the viral load-driven parameterization (VLP), which is the VIBES model of the first paper, and the geometric parameterization (GP), a traditional model with memoryless (geometrically distributed) latent and infectious periods.

## Running and processing main analysis simulations

All main text results (Figures 2 to 5) come from two beta grids, one for each parameterization, both with 50 independent executions for each value of beta. To reproduce these simulations, run the following commands on a terminal open in the root of this repository (this may take several minutes or hours and create ~50 gigabytes in files):

```bash
# Viral load-driven parameterization (VLP)
python run_beta_grid.py inputs/comparison_paper_in/vl_model_baseline.in outputs/comparison_paper_out/vl_model_baseline -n 50

# Geometric parameterization (GP)
python run_beta_grid.py inputs/comparison_paper_in/geometric_model_baseline.in outputs/comparison_paper_out/geometric_model_baseline -n 50
```

This should create the `outputs/comparison_paper_out` directory and populate it with results from the simulations. To postprocess these simulations and obtain the data files that were used for the main text figures, run:

```bash
python analysis/postprocess_beta_grid.py outputs/comparison_paper_out/vl_model_baseline
python analysis/postprocess_beta_grid.py outputs/comparison_paper_out/geometric_model_baseline
```

Which should result in additional files inside each of the output directories (see [outputs.md](outputs.md) for full details on the output files).

> Notice: Since different computer systems implement different random number generation mechanisms, results should be similar but numerically different from the ones presented in our paper.

## Time-dependent metrics

Supplementary Figure S1 uses time-dependent metrics, which are calculated with:

```bash
python analysis/postprocess_time_series.py outputs/comparison_paper_out/vl_model_baseline/sim_0001/
python analysis/postprocess_time_series.py outputs/comparison_paper_out/geometric_model_baseline/sim_0001/
```

## Geometric infectiousness curves

The GP simulations read their transmissibility (infectiousness) curves from `transmissibility/geometric/`, which are already included in this repository. If you wish to re-create these files, run:

```bash
python generate_geometric_curves.py
```

By default, this writes to the same `transmissibility/geometric/` directory (overwriting the files there) and reproduces the curves used in the paper, since it uses the same random seed. The VLP simulations use the viral dynamics curves in `transmissibility/viral-dynamics/`, which are the same as in the first paper. The file format is documented in [transmissibility/README.md](../transmissibility/README.md).