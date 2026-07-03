# Cumulative Distribution Function (CDF)

`plot_CDF_with_statistics.m`

## Purpose

Plots empirical survival functions (1 – CDF) of lag time across selected sample groups
and computes a slope metric that quantifies how rapidly colonies exit the lag phase.
Three figures are produced per run: a linear-scale survival plot, a log-scale survival
plot (same values, log y-axis), and a single bar chart of mean slope. Pairwise
statistical tests (KS and log-rank) are printed to the terminal and saved to a log file.

## Requirements

- MATLAB R2020a or later (`exportgraphics` required).
- Statistics and Machine Learning Toolbox (`kstest2` in MATLAB).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct.

## Usage

Pass the `ht` struct as the only argument. An optional output struct captures per-group
survival statistics. The function opens an interactive terminal prompt for group
selection, then runs automatically.

```matlab
results = plot_CDF_with_statistics(ht);
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
Colony filtering is identical across `preprocess_pipeline_data`, `plot_combined_samples`,
`plot_combined_samples_growth_curves`, and this script. All four scripts operate on the
same validated colony set stored in `ht.groups(g).lag_time` and related fields.
:::

## Interactive Selection

On launch, a table of available groups is printed showing each group's index, label,
and colony count. Enter group labels separated by spaces or commas, or type `ALL`:

- **Groups** → enter labels separated by spaces or commas (e.g. `3h 24h`), or type
  `ALL`

Labels are case-insensitive. Duplicate entries are silently removed. An unrecognised
label triggers a warning and re-prompts without exiting.

## Figure Outputs

Three figure files and one log are written to `ht.params.out_dir`, all with a timestamp
and feature-list suffix:

| Item | Description |
|---|---|
| **Figure 1 — Linear survival** | Survival `S = 1 – CDF` vs lag time (h) on a linear y-axis. One coloured line per group with circle markers. Legend shown top-right. |
| **Figure 2 — Log survival** | Identical data to Figure 1 displayed on a log y-axis (set via `YScale log`). Useful for identifying exponential exit kinetics. |
| **Figure 3 — Slope bar chart** | A single bar chart showing mean `|ds/dt|` (fraction h⁻¹) per group. Bars are colored by group. |

## Slope Metrics

One finite-difference slope metrics is computed from the empirical survival curve:

| Item | Description |
|---|---|
| Linear slope `\|dS/dt\|` | Absolute dropout rate (fraction h⁻¹). Measures how fast the surviving fraction falls per hour regardless of its current level. |

Both max and mean values are reported per group in the log file and in the returned
`results` struct.

## Statistical Tests

All pairwise group comparisons are computed and printed to the terminal and log file:

| Item | Description |
|---|---|
| KS test (`kstest2`) | Two-sample Kolmogorov–Smirnov test on lag time distributions. Reports p-value. |
| Log-rank test | Log-rank test on survival curves. Reports p-value. Returns `NaN` if the `logrank` function is unavailable. |

## Output Struct

The optional return value `results` is a struct array with one element per selected
group:

| Parameter | Feature |
|---|---|
| `results(g).group` | Group label string |
| `results(g).x` | Unique sorted lag time points (h) |
| `results(g).survival` | Survival values S at each time point |
| `results(g).n` | Number of colonies |
| `results(g).median` | Median lag time (h) |
| `results(g).q25 / q75` | 25th and 75th percentile lag time (h) |
| `results(g).min_lag / median_lag / max_lag` | Min, median, and max lag time (h) |
| `results(g).slope` | Struct with max and mean `\|dS/dt\|` |

## Saved Files

Six figure files and one log are written to `ht.params.out_dir` with a
group-and-timestamp suffix:

| Item | Description |
|---|---|
| `CDF_linear_<groups>_<timestamp>.fig / .pdf / .png` | Linear-scale survival figure |
| `CDF_log_<groups>_<timestamp>.fig / .pdf / .png` | Log-scale survival figure |
| `CDF_slope_<groups>_<timestamp>.fig / .pdf / .png` | Slope bar chart figure |
| `plot_CDF_with_statistics_<timestamp>.txt` | Log file — group selection, pairwise statistics, and slope summaries |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — plot survival curves and compute statistics
results = plot_CDF_with_statistics(ht);
```

## Troubleshooting

| Error Message | Description |
|---|---|
| Group label not recognised | Labels are case-insensitive and must match exactly as listed in the selection menu. |
| Log-rank p = `NaN` | The `logrank` function was not found. Ensure the Statistics and Machine Learning Toolbox is installed, or ignore this metric. |
| Empty survival curve for a group | `ht.groups(g).lag_time` is empty, meaning no colonies passed QC for that group. Review `size_threshold` and `ecc_threshold` in `preprocess_pipeline_data`. |
| `exportgraphics` error | Requires MATLAB R2020a+. Use `print()` as a fallback on older versions. |
