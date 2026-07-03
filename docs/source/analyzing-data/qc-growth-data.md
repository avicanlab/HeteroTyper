# Growth data per plate

`plot_QC_growth_data.m`

## Purpose

Produces four per-plate QC figures for all plates that have `growth_quant == true`:
growth curves (smoothed size timecourses), colony size distributions, lag time
distributions, and early doubling time distributions. All per-plate data are read
directly from `ht.qc` — no access to the raw `data` struct is needed. Each figure is
arranged in the most square auto-computed grid that fits all plates, with unused slots
blanked. Colony counts are annotated in the top-left corner of every subplot.

## Requirements

- MATLAB R2020a or later (`exportgraphics` required).
- `preprocess_pipeline_data(data);` must be run first using the updated version of the
  script. `ht.qc.ix_growth` and the `ht.qc.plate_*` fields must be present — re-run the
  preprocessor if they are missing.

## Usage

Pass the `ht` struct as the only argument. No interactive prompts are shown. The
function runs automatically and saves all outputs to `ht.params.out_dir`.

```matlab
plot_QC_growth_data(ht);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.qc` (with `ht.qc.ix_growth` and the `ht.qc.plate_*` fields), `ht.params`, `ht.groups`, and `ht.colors`. If `ht.qc.ix_growth` is missing, re-run the updated `preprocess_pipeline_data`. |

## Plate selection

The function plots every plate index in `ht.qc.ix_growth`, which contains all plates
for which `data.processed{i}.growth_quant == true` during preprocessing — regardless of
whether the colony count falls within `[min_col, max_col]`. This is a broader set than
the QC-passing plates shown in `plot_QC_full_dataset` (which require both
`growth_quant` and colony count in range), so that growth kinetics can be inspected
even on plates that fail the count filter.

## Figures

Four figures are produced, each saved as `.fig`, `.pdf` (vector), and `.png` (300 dpi).
All use the same auto-square grid layout and group colouring from `ht.colors`.

| Figure | Description |
|---|---|
| **Plot 1 — Growth curves** | One subplot per plate. Each line is one valid colony's smoothed size timecourse (pre-masked by `flag_colony_ok`). Lines are drawn in the sample group colour at 45% alpha. x-axis runs from 0 to `max_lag` (h). |
| **Plot 2 — Colony size distributions** | Filled histogram of final colony sizes (last timepoint) per plate. Bin edges from `ht.xbins` (100 px width). All subplots share the global x-axis maximum across all plates. Red vertical line marks the median. |
| **Plot 3 — Lag time distributions** | Filled histogram of per-colony lag times (h, including `incTime`) per plate. Bin width 0.5 h; x-axis from `incTime` to `max_lag`. Red vertical line marks the median. |
| **Plot 4 — Early doubling time distributions** | Filled histogram of per-colony early doubling times (h) per plate. Bin width 0.5 h; x-axis from 0 to `max_lag`. Red vertical line marks the median. Colonies with non-positive or non-finite early DT are excluded. |

## Outputs

| Item | Description |
|---|---|
| `QC_growth_curves_<timestamp>.fig / .pdf / .png` | Plot 1 — Per-plate growth curve grids |
| `QC_size_dist_<timestamp>.fig / .pdf / .png` | Plot 2 — Per-plate colony size histogram grids |
| `QC_lag_time_dist_<timestamp>.fig / .pdf / .png` | Plot 3 — Per-plate lag time histogram grids |
| `QC_early_DT_dist_<timestamp>.fig / .pdf / .png` | Plot 4 — Per-plate early doubling time histogram grids |
| `plot_QC_growth_data_<timestamp>.txt` | Log file — plate count, grid dimensions, per-plate n_valid and distribution summary (min, max, median) for each of the four metrics. |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — growth QC plots
plot_QC_growth_data(ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| `ht.qc.ix_growth` not found | Re-run the updated `preprocess_pipeline_data`. The `ht.qc.ix_growth` and `ht.qc.plate_*` fields were added in the current version; older `ht` structs will not have them. |
| "No plates with growth data found" | No plates have `growth_quant == true` in `data.processed`. Check that the acquisition pipeline completed successfully for at least one plate. |
| Growth curve panel is empty for a plate | `ht.qc.plate_time{i}` or `ht.qc.plate_size_tc{i}` is empty for that plate. This can occur if `colonies.new.time_info` or `timecourse_size_smoothed` was not populated during acquisition. |
| Lag time or early DT panel is empty for a plate | `ht.qc.plate_lag{i}` or `ht.qc.plate_early_dt{i}` is empty. Re-run `preprocess_pipeline_data` to regenerate the cache, and verify that `colonies.new.lag_time` and `colonies.new.early_doublingtime` are present in `data.processed` for that plate. |
| `exportgraphics` error | Requires MATLAB R2020a+. Use `print()` as a fallback on older versions. |
