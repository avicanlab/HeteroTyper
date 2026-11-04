# Ridgeline plots with distribution intervals

`plot_ridgeplot.m`

## Purpose

Draws one ridgeline figure per selected colony feature (appearance time, size, area,
intensity, …), with one ridge per sample group stacked top to bottom. Each ridge is a
grey kernel density estimate with, on its baseline, three bars in the group colour
covering the central 95%, 80% and 50% of colonies and a black dot at the group median.
Per-group summary statistics (n, median, IQR, MAD_med, Gini_counts) are written next to
every ridge, and a Brown–Forsythe test of spread is reported below the axes for all
selected groups together and for any pairs chosen at the prompt. All statistics are
also written to the log and to an Excel workbook.

:::{note}
New in v1.1.0. This single script replaces the earlier experiment-specific appearance
time ridgeplot variants (`…_6s_stress`, `…_6srna_growth`) and extends them to every
colony feature in `ht.groups`.
:::

## Requirements

- MATLAB R2023b or later (`exportgraphics` required). No toolbox function is used — the
  density estimate, percentiles, F-distribution and BH adjustment are implemented
  locally.
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct,
  including `ht.groups(g).size_timecourse` and `ht.groups(g).time_h` (used for
  Gini_counts).

## Usage

Pass the `ht` struct as the only argument. Five interactive terminal prompts follow
before plotting.

```matlab
plot_ridgeplot(ht);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups`, `ht.params` (`out_dir`, `incTime`, `max_appearance`; `img_int` and `size_threshold` are used when present, otherwise 0.5 h and 100 px), `ht.labels`, and `ht.colors`. |

## Interactive Selection

Five selections are made in sequence:

- **Groups** → every group is listed with its colony count. Enter group labels
  (comma-separated if they contain spaces), index numbers, or `ALL`.
- **Display order** → enter the top-to-bottom order as numbers of the selected groups
  (e.g. `3 1 2`); groups not listed follow in their current order. Enter keeps the
  selection order.
- **Features** → every feature field present in `ht.groups` is listed with its unit.
  Enter numbers or names (comma-separated), or `ALL`. Enter alone selects appearance
  time only.
- **Statistical comparisons** → every pair of the selected groups is listed and
  numbered. Enter pair numbers (e.g. `1 4 6`), ranges (`1-3`), or `ALL`; Enter for
  none (the all-group test is always run). The same pairs are tested for every
  selected feature.
- **Appearance threshold** → size in px used for Gini_counts (Enter =
  `ht.params.size_threshold`, default 100 px).

Example session:

```text
  Groups — labels (comma-separated) / index numbers, or ALL: ALL
  New order as numbers from this list (e.g. 3 1 2; ...) [Enter = keep]: 2 1
  Features — numbers or names (comma-separated), or ALL [Enter = Appearance time]: 1 2 6
  Comparisons — pair numbers (e.g. 1 4 6, ranges 1-3), ALL [Enter = none]: ALL
  Appearance threshold in px for Gini_counts [Enter = 100]:
```

An unrecognised entry triggers a warning and re-prompts without exiting.

## Available Features

| # | Feature | Field | Unit |
|---|---|---|---|
| 1 | Appearance time | `appearance_time` | h |
| 2 | Size | `size` | px |
| 3 | Area | `area` | px |
| 4 | Intensity | `intensity` | a.u. |
| 5 | Mean intensity | `mean_intensity` | a.u. |
| 6 | Intensity per size | `int_per_size` | a.u./px |
| 7 | Perimeter | `perimeter` | px |
| 8 | Circularity | `circularity` | – |
| 9 | Eccentricity | `eccentricity` | – |
| 10 | Solidity | `solidity` | – |

Only features present in `ht.groups` are offered. Names, units and number formats are
defined in the local function `feature_table()` and can be edited there if the units in
your dataset differ.

## Colony Filtering

The script applies no colony filtering of its own. Feature values are read from
`ht.groups(g).<feature>` and only finite values are used, so the colony set is the same
QC-filtered set used by `plot_combined_samples` (size-filtered, censored colonies
excluded).

## Statistics

| Metric | Description |
|---|---|
| Median, Q25/Q75, IQR | Group median, quartiles and interquartile range (same percentile definition as MATLAB `prctile`). The darkest bar spans Q25–Q75. |
| P10/P90, P2.5/P97.5 | Percentiles bounding the 80% and 95% bars. |
| MAD_med | Mean absolute deviation from the group median — the per-group quantity compared by the Brown–Forsythe test. |
| Gini_counts | **Appearance time:** Gini of the appearance-count table — colonies whose smoothed size timecourse first reaches the appearance threshold at each imaging timepoint (newly appearing, not cumulative; zero rows included). **Other features:** Gini of the counts of the same appearance-counted colonies over 30 equal-width bins spanning the pooled range of the selected groups (shared edges, so groups are comparable). 0 = colonies spread evenly, 1 = maximally concentrated. This is a different quantity from the Gini in `ht.groups(g).gini`. |
| Brown–Forsythe test | Test of equal spread, robust to skew: z = \|x − group median\| per colony, one-way ANOVA F-test on z. Run for all selected groups together and for each chosen pair. Pair p-values are Benjamini–Hochberg (BH) adjusted across the chosen pairs within each feature. Effect sizes: MAD_med ratio and IQR ratio (first / second group; > 1 = first group more spread). Tests absolute spread in the feature's units, independent of a shift in the median. |

Significance stars: `*** p<0.001`, `** p<0.01`, `* p<0.05`, `ns`.

:::{note}
Colonies of all positions are pooled and treated as independent. With large n, small
differences in spread become significant — read the MAD and IQR ratios alongside the
p-values.
:::

## Figure Layout

One figure per selected feature. Height scales with the number of groups and with the
number of Brown–Forsythe lines below the axes.

| Element | Description |
|---|---|
| Ridges | Grey Gaussian kernel density estimate per group (Silverman bandwidth; for appearance time never narrower than the imaging interval, since appearance times sit on the imaging grid), scaled to the same height for every group. |
| Interval bars | Three stacked bars on each ridge baseline in shades of the group colour: 95% (lightest), 80%, and 50% / IQR (darkest). Black dot = group median. |
| Pooled median | Vertical dashed line at the median of all selected groups pooled, labelled `Median, all groups (…)`. |
| Group labels | Left of the axes, bold, in the group colour. |
| Info column | Right of the axes, per ridge: `n`, median with `[IQR]`, `MAD_med`, and `Gini`. |
| Legend | Top of the figure in neutral grey, explaining the density, median dot and bars. |
| Brown–Forsythe block | Below the x-axis: the all-group result `F(df1, df2), p`, then one line per chosen pair with `F`, `p_BH`, MAD ratio and IQR ratio. |
| x-axis | Appearance time: from `incTime` to `max_appearance` with ticks every 5 h. Other features: from the 0.5th to the 99.5th percentile of all selected groups pooled, with 5% padding and rounded ticks (never below 0 for non-negative features). Extreme outliers can lie off-axis but are still included in every statistic. |

## Outputs

All outputs are written to `ht.params.out_dir` with a timestamp suffix:

| Item | Description |
|---|---|
| `Ridgeplot_<feature>_<timestamp>.fig / .pdf / .png` | One ridgeline figure per selected feature (`<feature>` = field name, e.g. `appearance_time`, `size`). PDF is vector; PNG is 300 dpi. |
| `ridgeplot_stats_<timestamp>.xlsx` | Excel workbook with three sheets: `Stats` (per feature × group: n, median, quartiles, IQR, P10/P90, P2.5/P97.5, MAD_med, n_Gini, Gini_counts), `BrownForsythe` (all-group and pairwise F, df, p, p_BH, stars, MAD and IQR ratios), and `README` (column descriptions). |
| `plot_ridgeplot_<timestamp>.txt` | Log file — selections, per-feature statistics table, Brown–Forsythe results, and saved file paths. |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — ridgeline plots, e.g. appearance time and size for all groups
plot_ridgeplot(ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| `Input must be the ht struct produced by preprocess_pipeline_data` | Pass `ht`, not `data`. Run `preprocess_pipeline_data(data)` first. |
| `none of the feature fields is present in ht.groups` | The `ht` struct is incomplete or from an older version. Re-run `preprocess_pipeline_data`. |
| "no finite values in the selected groups — skipped" | The chosen feature has no finite values in any selected group; no figure is produced for it. |
| Gini shown as `NaN` / `n_Gini` = 0 | No colony in the group reached the appearance threshold, or `size_timecourse` / `time_h` is missing for that group. Lower the threshold or re-run the preprocessor. |
| Some colonies appear to be missing from the plot | The x-axis is clipped at the pooled 0.5th–99.5th percentile for features other than appearance time. Off-axis values are still included in all statistics. |
| `xlsx` not saved | The log shows `WARNING: could not write xlsx` with the exact error. Check that the file is not open in Excel. |
| `exportgraphics` error | Requires MATLAB R2023b or later. Use `print()` as a fallback on older versions. |
