# Doubling time

`plot_doublingTime.m`

## Purpose

Plots normalised probability density histograms of colony doubling time (h) for
selected sample groups. Each panel shows one group with a red median line, Gini index
annotation, and group label. All panels share a common x-axis so groups are directly
comparable. The script reads pre-computed histogram, statistics, and Gini values from
the `ht` struct — no recomputation at plot time.

## Requirements

- MATLAB R2020a or later (`exportgraphics` required).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct,
  including `hist`, `gini`, `stats`, and `xbins`.

## Usage

Pass the `ht` struct as the only argument. The function opens an interactive terminal
prompt for group selection, then runs automatically.

```matlab
plot_doublingTime(ht);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by the updated `preprocess_pipeline_data`. Must contain `ht.groups`, `ht.params`, `ht.labels`, `ht.colors`, `ht.xbins.doublingTime`, and per-group `hist` / `stats` / `gini` fields for `doublingTime`. |

## Colony filtering

Doubling time values in `ht.groups(g).doublingTime` are computed inside the same valid
colony mask as every other feature in `preprocess_pipeline_data`. A colony is included
if it passes all three QC criteria:

- Final colony size `> size_threshold` (default: 100 px)
- Eccentricity `< ecc_threshold` (default: 0.70)
- Finite intensity-per-size ratio

Additionally, a doubling time value is set to `NaN` for any colony whose smoothed size
curve never reaches twice the threshold size within the observation window. `NaN`
values are excluded from all histogram, statistics, and Gini computations.

:::{note}
Colony filtering is identical across `preprocess_pipeline_data`, `plot_combined_samples`,
`plot_combined_samples_growth_curves`, `plot_CDF_with_statistics`,
`plot_KS_quantile_statistics`, and this script. All scripts operate on the same
validated colony set.
:::

## How Doubling Time Is Calculated

Doubling time is calculated once per colony during preprocessing from the smoothed size
timecourse:

- Find the first time point where colony size `>= size_threshold`. Call this `t_start`.
- Find the first subsequent time point where size `>= 2 x size(t_start)`. Call this
  `t_end`.
- Doubling time = `t_end - t_start` (h).
- If either time point cannot be found the value is `NaN` and the colony is excluded
  from all downstream analysis.

:::{note}
Doubling time reflects the time to first doubling after the colony crosses the size
threshold, not from time zero. Colonies that never double within the observation window
are excluded.
:::

## Interactive Selection

On launch, a table listing all available groups and their valid doubling time colony
counts is printed. Enter group labels separated by spaces or commas, or type `ALL`:

- **Groups** → enter labels separated by spaces or commas (e.g. `3h 24h`), or type
  `ALL`

Labels are case-insensitive. An unrecognised label triggers a warning and re-prompts
without exiting. Duplicate entries are silently removed.

## Figure Layout

One panel per selected group, stacked vertically, top group first. All panels share the
same x-axis range (from `ht.xbins.doublingTime`, 0.5 h bin width). Each panel shows:

- Filled probability density bar histogram in the group colour
- Black outline (stairs) trace over the bars
- Red vertical median line spanning the full panel height
- Smart-corner annotation: group label, n, median (h), and Gini index

The annotation is placed in the corner with the least data mass (left or right) to
avoid overlap with bars. The red median line is only drawn when the median is a finite
number; panels with no valid data show a grey "no valid data" message instead.

## Outputs

Three figure files, one log file, and one workspace variable are produced:

| Item | Description |
|---|---|
| `DoublingTime_<groups>_<timestamp>.fig` | MATLAB figure — re-openable and editable |
| `DoublingTime_<groups>_<timestamp>.pdf` | Vector PDF for publication |
| `DoublingTime_<groups>_<timestamp>.png` | Raster PNG at 300 dpi |
| `plot_doublingTime_<timestamp>.txt` | Log file — selected groups, per-group n / median / Gini, and saved file paths |
| `combined_doublingT` | Concatenated finite doubling time vector across all selected groups, written to the MATLAB base workspace |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — plot doubling time histograms
plot_doublingTime(ht);
```

## Troubleshooting

| Error Message | Description |
|---|---|
| Red median line missing | Median is `NaN`. This means the old preprocessor was used — `ht.groups(g).stats.doublingTime.median` was computed before `NaN` values were stripped. Re-run the updated `preprocess_pipeline_data` to fix. |
| `ht.xbins.doublingTime` not found | Re-run the updated `preprocess_pipeline_data`. The script falls back to on-the-fly bin computation with a warning, but pre-computed histogram, Gini, and stats will not be available. |
| Median shown as `NaN` h in annotation | Same root cause as above — old `ht` struct in workspace. Re-run the updated preprocessor. |
| No valid doubling time data for a group | All colonies either failed QC or never doubled within the observation window. Check `size_threshold`, `ecc_threshold`, and imaging duration in `preprocess_pipeline_data`. |
| `exportgraphics` error | Requires MATLAB R2020a+. Use `print()` as a fallback on older versions. |
| Log file not created | Check that `ht.params.out_dir` exists and is writable. |
