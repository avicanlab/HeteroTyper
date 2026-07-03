# Comparing lag time with other phenotypic features

`plot_lagTime_vs_morphology.m`

## Purpose

Plots scatter panels of lag time (h) versus selected morphological features. All
selected sample groups are overlaid in a single panel per feature, each in its own
colour. A linear LOESS-style regression line and a colour-matched annotation showing the
group label, Spearman rho (ρ), R², p-value, and polynomial equation are drawn for every
group.

## Requirements

- MATLAB R2020a or later (`exportgraphics` required).
- Statistics and Machine Learning Toolbox (`kstest2` in MATLAB).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct. It
  computes ρ, R², slope, and intercept in `ht.groups(g).morph_corr` for each lag-time
  vs morphology pair.

## Usage

Pass the `ht` struct as the only argument. An optional output struct captures per-group
survival statistics. The function opens an interactive terminal prompt for group and
morphology feature selection, then runs automatically.

```matlab
results = plot_lagTime_vs_morphology(ht);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups`, `ht.params`, `ht.labels`, and `ht.colors`. |

## Colony filtering

The script reads lag times directly from `ht.groups(g).lag_time`, which already
contains only QC-passing colonies from the preprocessor. No additional filtering is
applied here. A colony is included if it passed all three criteria during
preprocessing:

- Final colony size `> size_threshold` (default: 100 px)
- Eccentricity `< ecc_threshold` (default: 0.70)
- Finite intensity-per-size ratio

:::{note}
Colony filtering is identical across all pipeline scripts. All statistics in
`ht.groups(g).corr` are computed on the same filtered set, so ρ, R², p-value, and slope
equation shown in each panel are fully consistent with every other pipeline output.
:::

## Statistics – r, ρ, R², slope, intercept

All four values are computed once per group per feature in STEP 6c of
`preprocess_pipeline_data` and stored in `ht.groups(g).morph_corr.(field)`. The plot
script reads them directly without recomputing.

| Parameter | Feature |
|---|---|
| r (Pearson) / ρ (Spearman) | Pearson r: direction and strength of the linear relationship. Spearman ρ: rank correlation, robust to non-linear associations — this is the primary statistic shown in the figure annotation. Both range −1 to +1. |
| R² | Fraction of between-colony variance in the morphological feature explained by lag time. R² = r² for simple linear regression. Requested by reviewers as the standard effect-size complement to r. |
| Slope | Change in the morphological feature per hour of additional lag time. Used to draw the regression line. |
| Intercept | Predicted feature value at lag time = 0. Used together with slope to draw the regression line. |
| p-value | Two-tailed p-value for Spearman ρ (shown in the annotation) and Pearson r (stored in `morph_corr`). Significance marked as `*** p<0.001`, `** p<0.01`, `* p<0.05`, `ns`. |

## Interactive Selection

Two prompts are shown on launch:

- **Groups** → enter labels separated by spaces or commas (e.g. `3h 24h`), or type `ALL`
- **Features** → enter feature numbers separated by commas (e.g. `1 3 5`), or type `ALL`

The 9 available morphological features are:

| Parameter | Feature |
|---|---|
| 1 | Colony size (px) |
| 2 | Area (px) |
| 3 | Intensity |
| 4 | Mean intensity |
| 5 | Intensity / size |
| 6 | Perimeter (px) |
| 7 | Circularity |
| 8 | Eccentricity |
| 9 | Solidity |

Group labels are case-insensitive. Invalid entries trigger a warning and re-prompt
without exiting. Duplicates are silently removed.

## Figure Layout

One panel per selected feature, arranged in a single row. All selected groups are
overlaid in each panel.

| Parameter | Feature |
|---|---|
| Panel grid | 1 row × n_features columns. Figure width scales with feature count; height is fixed. |
| Scatter points | Filled circles, dot size 12 pt, alpha = 0.30, each group in its own colour. |
| Regression line | Linear fit (`polyfit` degree 1) spanning only the x-range of that group's own data — no extrapolation beyond the observed lag time range. Line colour = 90% of group colour (darker). |
| Grid | Major dashed y-grid lines (grey 0.72); one lighter dashed minor grid line exactly midway between each pair of major ticks. No x-grid. |
| Axis lines | Solid black, LineWidth 1.0. Outward ticks. No box. |
| Annotation | One line per group, stacked from the top-left corner downward (10% of y-range per step). Text colour matches the group's dot colour. Format: group label, then `ρ=x.xxx`, `R^2=x.xxx`, `p=x.xxx` (significance stars), and the polynomial equation of the LOESS smooth (shown in red on the line below). |
| x-axis | Lag Time (h), ticks every 4 h, label on every panel. |
| y-axis | Feature name label on every panel; tick values on every panel. |

## Outputs

Three figure files and one log file are written to `ht.params.out_dir`:

| Item | Description |
|---|---|
| `LagTime_vs_Morphology_<groups>_<timestamp>.fig` | MATLAB figure — re-openable and editable |
| `LagTime_vs_Morphology_<groups>_<timestamp>.pdf` | Vector PDF for publication |
| `LagTime_vs_Morphology_<groups>_<timestamp>.png` | Raster PNG at 300 dpi |
| `plot_lagTime_vs_morphology_<timestamp>.txt` | Log file — selected groups and features, full r / R² / p table for every group–feature combination, and saved file paths |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — plot scatter panels
plot_lagTime_vs_morphology(ht);
```

## Troubleshooting

| Error Message | Description |
|---|---|
| Annotation shows only group label, no ρ / R² / p | `ht.groups(g).morph_corr` not found — old preprocessor struct loaded. Re-run the updated `preprocess_pipeline_data`. A warning is printed to the terminal. |
| Regression line absent for a group | Same cause as above, or `morph_corr.slope` is `NaN` (fewer than 3 valid colonies). Re-run the preprocessor. |
| Regression line extends too far | Cannot happen — the line is clamped to `[min(lag_time), max(lag_time)]` of each group's own data. |
| Annotations overlap each other | Too many groups selected for the panel height. Reduce group count or increase `panel_h` at the top of the function. |
| y-axis range too compressed | `ht.global.*` is driven by outliers. Adjust `size_threshold` or `ecc_threshold` in `preprocess_pipeline_data`. |
| `exportgraphics` error | Requires MATLAB R2020a+. Use `print()` as a fallback on older versions. |
| Log file not created | Check that `ht.params.out_dir` exists and is writable. |
