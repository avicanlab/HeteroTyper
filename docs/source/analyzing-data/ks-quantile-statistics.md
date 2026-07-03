# Kolmogorov – Smirnov (KS) and Quantile Statistics

`plot_KS_quantile_statistics.m`

## Purpose

Visualises pairwise Kolmogorov–Smirnov (KS) D statistics and 90th-percentile quantile
differences (Q90) across selected sample groups and phenotypic features. Two horizontal
bar-chart figures are produced — one per metric — alongside an Excel workbook
containing all plotted values. The script performs no statistical computation itself;
all values are read directly from `ht.pairwise`, which is pre-computed by
`preprocess_pipeline_data`.

## Requirements

- MATLAB R2020a or later (`exportgraphics` required).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct,
  including `ht.pairwise`.

## Usage

Pass the `ht` struct as the only argument. An optional output struct captures per-group
survival statistics.

```matlab
plot_KS_quantile_statistics(ht);
```

The function opens an interactive terminal prompt for group selection, then runs
automatically.

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. |
| `ht.pairwise` | Required sub-struct containing pre-computed `KS_D`, `Q90`, `pair_labels`, `pair_idx`, and `feat_names`. Error is raised if absent. |

## Colony filtering

This script applies no colony filtering of its own. All KS and Q90 statistics in
`ht.pairwise` were computed during preprocessing on the same QC-filtered colony set
used by every other pipeline script. A colony was included if it passed all three
criteria in `preprocess_pipeline_data`:

- Final colony size `> size_threshold` (default: 100 px)
- Eccentricity `< ecc_threshold` (default: 0.70)
- Finite intensity-per-size ratio

:::{note}
Colony filtering is identical across `preprocess_pipeline_data`, `plot_combined_samples`,
`plot_combined_samples_growth_curves`, `plot_CDF_with_statistics`, and this script. All
scripts operate on the same validated colony set.
:::

## Interactive Selection

On launch, the full selection menu is printed in one screen showing all available
groups and features. Two selections are made in sequence:

- **Groups** → enter labels separated by spaces or commas (e.g. `3h 24h`), or type
  `ALL`. At least 2 groups must be selected for pairwise comparison.
- **Features** → enter feature numbers separated by spaces or commas (e.g. `1 2 3`), or
  type `ALL`

Group labels are case-insensitive. An unrecognised label or invalid feature number
triggers a warning and re-prompts without exiting. Duplicate entries are silently
removed.

Example session:

```text
Groups → enter labels e.g. [ 3h 24h ] or ALL: ALL
Features → enter numbers e.g. [ 1 2 3 ] or ALL: 1 2
```

## Metrics Explained

| Metric | Description |
|---|---|
| KS statistic D | Kolmogorov–Smirnov D statistic: the maximum absolute difference between the two empirical CDFs. Ranges 0–1. A value of 0 means identical distributions; 1 means completely non-overlapping. The dashed vertical line at `D = 0.5` provides a visual reference. |
| Quantile difference Q90 | Difference in the 90th-percentile values between two groups for a given feature. Positive values indicate the first-named group has a higher Q90; negative values the opposite. The x-axis is symmetric and auto-scaled per feature column. |

:::{note}
Both metrics are computed once during preprocessing across all group pairs and all 10
features, then cached in `ht.pairwise`. Re-running this script with different group or
feature selections does not recompute anything.
:::

## Figures Produced

Two horizontal bar-chart figures are generated, each with one column per selected
feature and one row per group pair:

| Item | Description |
|---|---|
| **Figure 1 — KS Statistics (D)** | Fixed x-axis 0 to 1 with a tick at 0.5. Dashed grey reference line at D = 0.5. Bars are greyscale, graduating light to dark across features. |
| **Figure 2 — Quantile Difference (Q90)** | Symmetric x-axis auto-scaled per feature column. Solid grey zero-line. Same greyscale palette as Figure 1. Pair labels shown on y-axis of the first column only. |

Figure width scales automatically with the number of selected features; figure height
scales with the number of pairs. Pair labels on the y-axis drive the left margin width
to prevent overflow.

## Saved Files

All outputs are written to `ht.params.out_dir` with a group-and-timestamp suffix:

| Item | Description |
|---|---|
| `KS_Statistics_<groups>_<timestamp>.fig / .pdf / .png` | KS D bar-chart figure |
| `Quantile_Q90_<groups>_<timestamp>.fig / .pdf / .png` | Q90 difference bar-chart figure |
| `KS_Quantile_Statistics_<groups>_<timestamp>.xlsx` | Excel workbook with two sheets: `KS_D` and `Q90_diff`. Rows = group pairs, columns = features. |
| `plot_KS_quantile_statistics_<timestamp>.txt` | Log file — selections made, full numeric table of KS D and Q90 values, and saved file paths. |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — plot survival curves and compute statistics
results = plot_KS_quantile_statistics(ht);
```

## Troubleshooting

| Error Message | Description |
|---|---|
| `ht.pairwise` not found | Re-run `preprocess_pipeline_data` — this field is only present in structs generated by the current version of the preprocessor. |
| No pairs found for selected groups | Fewer than 2 valid groups were selected. Select at least 2 groups to enable pairwise comparison. |
| Group label not recognised | Labels must match exactly as listed in the selection menu (case-insensitive). Check for trailing spaces or typos. |
| `xlsx` not saved | `writecell` requires MATLAB R2019a+. Check the log file warning for the exact error message. |
| `exportgraphics` error | Requires MATLAB R2020a+. Use `print()` as a fallback on older versions. |
