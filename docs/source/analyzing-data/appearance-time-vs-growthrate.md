# Maximum growth rate

`plot_lagTime_vs_growthRate.m`

## Purpose

Plots scatter panels of lag time (h) versus maximum Gompertz growth rate μₘₐₓ (px/h)
for selected sample groups. Each point is coloured by log₁₀(local colony density)
using a viridis palette shared across all panels, so density effects are directly
comparable. A LOESS smooth (red line) is overlaid on each panel when n ≥ 10 fitted
colonies, with its degree-2 polynomial equation annotated in red. Spearman ρ, R², and
density-weighted Spearman ρ𝑤 are annotated per group. Gompertz fitting and all
statistics are pre-computed in STEP 6d of `preprocess_pipeline_data`; this script reads
and visualises only.

## Requirements

- MATLAB R2020a or later (`exportgraphics` required).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct,
  including `ht.groups(g).growth` from STEP 6d (Gompertz fitting).

## Usage

Pass the `ht` struct as the only argument. An interactive terminal prompt asks for
group selection before plotting.

```matlab
plot_lagTime_vs_growthRate(ht);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups` (with `ht.groups(g).growth` from STEP 6d), `ht.params`, `ht.labels`, and `ht.colors`. |

## Colony Filtering

Only colonies with a successful Gompertz fit are plotted. A colony is included if it
passes all of:

- `growth.fit_ok == true` (Gompertz converged)
- Finite `mu_max` (maximum growth rate from Gompertz fit)
- Finite `lag_time` for the same colony

## Figure Layout

One panel per selected group, arranged in a grid (columns = `ceil(√n)`, rows filled
top-to-bottom). A shared viridis colorbar (log₁₀ colony density) is placed to the right
of the grid.

| Element | Description |
|---|---|
| Scatter points | Filled circles, size 27 pt, 80% alpha. Colour encodes log₁₀(local colony density) via the shared viridis palette. |
| LOESS smooth | Red solid line (span 0.75) drawn when n ≥ 10. A degree-2 polynomial fitted to the smoothed curve is annotated in red below the statistics (format: `y = a x² + b x + c`). |
| Annotation (top-right) | Four lines in black: `n`, `ρ = x.xxx p = x.xxx (sig)`, `R² = x.xxx`, `ρ𝑤 = x.xxx`. Fifth line in red: LOESS polynomial equation (omitted when n < 10). Significance stars: `*** p<0.001`, `** p<0.01`, `* p<0.05`, `ns`. |
| Axes | x: Lag time (h), ticks every 8 h, shared across all panels. y: Maximum growth rate (px/h), ticks every 100 px/h, shared. Major + midpoint minor dashed grid. Mirror ticks on top and right axes. |

## Outputs

All outputs are written to `ht.params.out_dir` with a group-and-timestamp suffix:

| Item | Description |
|---|---|
| `LagTime_vs_GrowthRate_<groups>_<timestamp>.fig` | MATLAB figure — re-openable and editable |
| `LagTime_vs_GrowthRate_<groups>_<timestamp>.pdf` | Vector PDF for publication |
| `LagTime_vs_GrowthRate_<groups>_<timestamp>.png` | Raster PNG at 300 dpi |
| `plot_lagTime_vs_growthRate_<timestamp>.txt` | Log file — selected groups, per-group n / Spearman ρ / p / weighted ρ, and saved file paths |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session (must include STEP 6d)
preprocess_pipeline_data(data);

% Step 2 — plot lag time vs growth rate
plot_lagTime_vs_growthRate(ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| `ht.groups(g).growth` not found | Re-run `preprocess_pipeline_data` — STEP 6d (Gompertz fitting) must complete before calling this function. |
| "No successful Gompertz fits found" | No colony in the selected groups had a converged Gompertz fit with a finite `mu_max`. Check that imaging duration was sufficient and that colony growth was observed within the acquisition window. |
| LOESS line absent from a panel | Fewer than 10 colonies with successful fits in that group — the LOESS threshold. Stats annotation is still shown. |
| `exportgraphics` error | Requires MATLAB R2020a+. Use `print()` as a fallback on older versions. |
