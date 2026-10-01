# Uncertainty-Aware Near-Field Multi-Beam Management in XL-RIS Systems Under Mobility and Computational Constraints: reproducibility package

Code and saved data to reproduce the figures and tables of:

> I. G. Sarfo-Mainoo, J. K. Mensah, S. Amonovi, K. Ntiamoah-Sarpong, and E. A. Affum, "Uncertainty-Aware Near-Field Multi-Beam Management in XL-RIS Systems Under Mobility and Computational Constraints," pre-submission snapshot. DOI: to be added after publication

This is a **pre-submission snapshot**, not an archived record of a published paper. Update this README, `CITATION.cff` and the DOI after submission/acceptance.

## Contents
| Path | Description |
|---|---|
| `03_MATLAB/experiments/` | Experiment scripts, figure script (`plot_figs_v11b.m`) and saved `.mat` results |
| `03_MATLAB/src/` | Channel, solver and metric functions |
| `07_Manuscript/figures/` | Reference copies of the six MATLAB figures as used in the manuscript |
| `MANIFEST.csv` / `manifest/manifest.csv` | Every file, its purpose, and the figure/table it supports |
| `VERSIONS.md` | Software and hardware environment |

## Requirements
- MATLAB R2025b (the release used for the paper). Toolboxes: see `VERSIONS.md`.
- No Python or third-party packages are required.
- Platform of the reported timings: Intel Core i7-1165G7 @ 2.80 GHz, 32 GB RAM, 4 compute threads, double precision.

## Reproducing the results
Run from `03_MATLAB/experiments/`, **in a copy of this folder** (scripts save their outputs in the current folder and will overwrite the shipped `.mat` files):

```matlab
cd 03_MATLAB/experiments
addPaper2Paths;
```

| Goal | Command | Paper item | Approx. cost |
|---|---|---|---|
| Regenerate all six figures from the saved data | `plot_figs_v11b` | Figs. 3, 4, 6, 7, 8, 9 | seconds |
| Re-run the focusing-loss validation | `exp02a_focusing_loss_validation` | Figs. 3, 4 | under 1 min |
| Re-run the timing round | `exp05b_tableIV_consistent` | Table tab:budget | minutes |
| Re-run loss by direction | `exp03g_loss_by_direction` | Table tab:loss | not measured |
| Re-run closed-form vs MC | `exp06_mcpsi_baseline` | Table tab:mcpsi | not measured |
| Full rate/outage chain (long) | `exp04a_rate_outage` -> `exp12a_rate_outage_safe` -> `exp14a_rate_outage_bestiter` | Figs. 6, 7; Table tab:worstcase | not measured |
| Full sparsity chain (long) | `exp04b_feedback_sparsity` -> `exp12b_sparsity_safe` -> `exp12c_finish` -> `exp14b_sparsity_bestiter` | Fig. 8 | not measured |
| Rician check (long) | `exp14e_v5g48` | Fig. 9 | not measured |

`plot_figs_v11b.m` in the **Zenodo** archive is byte-identical to the author's file and contains a hard-coded root path on line 68 (`root = 'C:\Dev\...'`). Replace that line with
`root = fileparts(fileparts(fileparts(mfilename('fullpath'))));`
before running. The **GitHub** copy already contains this one-line fix (logged in `docs/decisions.md`).

## Randomness
All stochastic scripts set `rng(seed,'twister')` (seeds 11, 23, 2026, etc., listed in Table "Simulation parameters" of the paper). Iteration counts are reproducible across machines; wall-clock timings are not and vary by several percent.

## Known limitations
- Saved `.mat` files are provided so figures can be regenerated without re-running hours-long sweeps.
- Exact numerical values can differ slightly across MATLAB releases and BLAS libraries.
- Table tab:mcpsi timings were obtained on a different run than the shipped `exp06_mcpsi.mat` (see `docs/decisions.md`).

## Citation
See `CITATION.cff`. Cite this archive: https://doi.org/10.5281/zenodo.23071062

## License
Code: GPL-2.0-only. Data: CC-BY-4.0.
