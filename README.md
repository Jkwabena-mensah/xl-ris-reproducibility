# Uncertainty-Aware Near-Field Multi-Beam Management in XL-RIS Systems Under Mobility and Computational Constraints: reproducibility package

Code and saved data to reproduce the figures and tables of:

> I. G. Sarfo-Mainoo, J. K. Mensah, S. Amonovi, K. Ntiamoah-Sarpong, and E. A. Affum, "Uncertainty-Aware Near-Field Multi-Beam Management in XL-RIS Systems Under Mobility and Computational Constraints," manuscript submitted to Elsevier.

The DOI of the paper will be linked to the archived record through the Zenodo related-identifier field once it exists. This README is not edited after archiving.

Archived release (cite this): https://doi.org/10.5281/zenodo.23071062

The GitHub repository (tag `v1.0.0`) and the Zenodo archive contain identical files. `SHA256SUMS` lists the SHA-256 hash of every other file.

## Contents
| Path | Description |
|---|---|
| `03_MATLAB/experiments/` | Experiment scripts, the figure script `plot_figs_v11b.m`, and saved `.mat` results |
| `03_MATLAB/src/` | Channel, solver and metric functions |
| `07_Manuscript/figures/` | Reference copies of the six MATLAB-generated figures (Figs. 4-9). Fig. 3 of the paper is a schematic and is not produced by this code |
| `MANIFEST.csv` | Every file, its purpose, and the figure/table it supports |
| `SHA256SUMS` | SHA-256 hash of every file except itself |
| `VERSIONS.md` | Software and hardware environment, and the random seeds used |
| `CITATION.cff`, `LICENSE` | Citation metadata and licences |

## Requirements
- MATLAB R2025b (the release used for the paper). Toolboxes: see `VERSIONS.md`.
- No Python or third-party packages are required.
- Platform of the reported timings: Intel Core i7-1165G7 @ 2.80 GHz, 32 GB RAM, 4 compute threads, double precision.

## Reproducing the results
Run from `03_MATLAB/experiments/`, **in a copy of this folder**. Scripts save their outputs in the current folder and overwrite the shipped `.mat` files.

```matlab
cd 03_MATLAB/experiments
addPaper2Paths;
```

| Goal | Command | Paper item | Approx. cost |
|---|---|---|---|
| Regenerate all six figures from saved data | `plot_figs_v11b` | Figs. 4-9 | seconds |
| Re-run the focusing-loss validation | `exp02a_focusing_loss_validation` | Figs. 4, 5 | about 22 s |
| Re-run the timing round | `exp05b_tableIV_consistent` | Table 5 | minutes |
| Re-run loss by direction | `exp03g_loss_by_direction` | Table 7 | not measured |
| Re-run closed form vs Monte Carlo | `exp06_mcpsi_baseline` | Table 6 | not measured |
| Rate/outage chain (long) | `exp04a_rate_outage` -> `exp12a_rate_outage_safe` -> `exp14a_rate_outage_bestiter` | Figs. 6, 7; Table 8 | not measured |
| Sparsity chain (long) | `exp04b_feedback_sparsity` -> `exp12b_sparsity_safe` -> `exp12c_finish` -> `exp14b_sparsity_bestiter` | Fig. 8 | not measured |
| Rician check (long) | `exp14e_v5g48` | Fig. 9 | not measured |

Notes:
- `exp05b_tableIV_consistent` keeps an earlier table numbering in its file name; it produces Table 5 of the paper.
- `plot_figs_v11b` prints two lines beginning `LoS t:` and `Ric t:`. These are the paired t-statistics shown in Fig. 9(b), not debug output.
- `plot_figs_v11b.m` locates the package root relative to its own location. The line that set a hard-coded path in the authors' working copy was replaced for this release.
- `exp04b_feedback_sparsity_PUBLISHED_backup.mat` is a backup of earlier saved sparsity results. It is read by `exp12c_finish.m` for comparison and is not the source of Fig. 8.
- Saved results for the table experiments: `exp05b_tableIV.mat` (Table 5), `exp03g_loss.mat` (Table 7), `exp06_mcpsi.mat` and `exp06_mcpsi_RESULTS.mat` (Table 6).

## Verification performed before release
- `exp02a_focusing_loss_validation` was re-run on a scratch copy of the archive and matches the shipped results (21 numeric variables compared, maximum relative difference 0).
- The six figure PDFs were regenerated from the saved data.
- The two packages were checked to be byte-identical, and the zip was unpacked and re-hashed against `SHA256SUMS`.
- Not re-run during verification: the other experiment scripts. Their saved results are the `.mat` files listed above and in `03_MATLAB/experiments/`.

## Randomness
All stochastic scripts call `rng(seed,'twister')`. Seeds used: `exp04a`, `exp12a`, `exp14a`: 11; `exp04b`, `exp12b`, `exp12c`, `exp14b`: 23; `exp03g`: 31; `exp06`: 606; `exp02a`, `exp05b`: 2026; `exp14e`: 11 (first restarts) then 777 (remaining restarts); `plot_figs_v11b`: 9 (figure jitter only). Table 4 of the paper lists the main seeds; this list is complete. Iteration counts are reproducible across machines; wall-clock timings are not and vary by several percent.

## Known limitations
- Saved `.mat` files are provided so figures can be regenerated without re-running long sweeps.
- Exact numerical values can differ slightly across MATLAB releases and BLAS libraries.
- Timings in Table 5 are specific to the platform above and will not reproduce exactly.
- The Table 6 timings were obtained in a different run from the shipped `exp06_mcpsi.mat`.

## Citation
See `CITATION.cff`. Cite this archive: https://doi.org/10.5281/zenodo.23071062

## License
Code: GPL-2.0-only. Data: CC-BY-4.0.
