# Growth curves

`plot_combined_samples_growth_curves.m`

## Purpose

Plots colony growth curves for every sample group in the experiment. Each line
represents one colony, coloured distinctly. All panels share a common y-axis, so groups
are directly comparable.

## Requirements

- MATLAB R2020a or later (`exportgraphics` required).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct.

## Usage

Pass the `ht` struct as the only argument. No interactive prompts are shown.

```matlab
plot_combined_samples_growth_curves(ht);
```

The function runs automatically and saves all outputs to `ht.params.out_dir`.

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups`, `ht.params`, and `ht.labels`. |

## Colony Filtering

Only colonies that passed QC during preprocessing are plotted. A colony is included if
it meets all three criteria set in `preprocess_pipeline_data`:

- Final colony size `> size_threshold` (default: 100 px)
- Eccentricity `< ecc_threshold` (default: 0.70)
- Finite intensity-per-size ratio

:::{note}
`grp.size_timecourse` already contains only QC-passing colonies after preprocessing —
no further filtering is needed at plot time. The valid colony count per group is
printed to the terminal and the log file.
:::

## Outputs

Three figure files and one log are written to `ht.params.out_dir`, all with a
timestamp suffix:

| Parameter | Description |
|---|---|
| `GrowthCurves_<timestamp>.fig` | MATLAB figure, reopenable and editable |
| `GrowthCurves_<timestamp>.pdf` | Vector PDF for publication |
| `GrowthCurves_<timestamp>.png` | Raster PNG at 300 dpi |
| `plot_combined_samples_growth_curves_<timestamp>.txt` | Log file, valid/total colony counts and file paths |

The figure contains one panel per sample group, stacked vertically. All panels share
the same y-axis ceiling (global max final size + 100 px) and x-axis starting at
`p.incTime`.

## Adjustable Layout Parameters

Edit these constants near the top of the function to change figure geometry:

| Parameter | Description |
|---|---|
| `panel_w` | Panel width in inches (default: 5.0) |
| `panel_h` | Panel height in inches (default: 2.0) |
| `gap_y` | Vertical gap between panels in inches (default: 0.55) |
| `margin_l` / `margin_r` | Left and right figure margins in inches |
| `margin_b` / `margin_t` | Bottom and top figure margins in inches |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — plot growth curves
plot_combined_samples_growth_curves(ht);
```

## Troubleshooting

| Parameter | Description |
|---|---|
| Blank / "No valid colonies" panel | No colonies passed QC for that group. Review `size_threshold` and `ecc_threshold` in `preprocess_pipeline_data`. |
| `exportgraphics` error | Requires MATLAB R2020a+. Use `print()` as a fallback on older versions. |
| Log file not created | Check that `ht.params.out_dir` exists and is writable. |
| Y-axis too compressed | One group has very large colonies. Adjust `y_ceil` manually at the top of the function. |
