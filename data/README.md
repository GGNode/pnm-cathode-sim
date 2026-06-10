# Data

## Paper Reference

- `khan2021_pnm_lib_cathode.pdf` — Khan et al. (2021) original paper
- `paper.pdf` — Duplicate paper copy

## Paper Figures

`paper_figures/` — Cropped figure PNGs extracted from the paper for validation comparison:

- `figure4_crop_200dpi.png` — V-Q discharge curves at 0.2C, 0.5C, 1C, 3C
- `figure6_crop_200dpi.png` — Spatial SoL distribution at 1C and 3C (75% SoL)
- `figure7_crop_200dpi.png` — Electrolyte concentration profiles
- `page_08_200dpi.png` — `page_11_200dpi.png` — Full page renders

## Validation Results

`validation/` — Pre-computed validation outputs:

- `discharge_*.npz` — Discharge simulation results (0.2C, 0.5C, 1C, 3C)
- `spatial_*_75sol.npz` — Spatial snapshots at 75% state of lithiation
- `figure*_comparison.png` — Generated comparison plots
- `validation_metrics.json` — Quantitative validation metrics
- `parallel_results.json` — Parallel multi-scenario validation results

## Notes

- If the original Khan et al. (2021) NREL XCT data becomes available, it can be used to generate geometry-matched networks via `scripts/generate_network.py`.
- Current validation uses synthetic cubic networks; quantitative comparison with paper figures is qualitative only.
