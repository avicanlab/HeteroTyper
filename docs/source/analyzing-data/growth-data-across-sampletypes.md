# Exploratory growth data overview across sample types

`plot_growth_data_across_sampletypes.m`

## Purpose

An exploratory visualisation script that works directly on the raw `data` struct
returned by `heterotyper`. It splits plates into four
hard-coded time groups (3H, 7H, 24H, 48H) read from `data.metadata.original.Time`, then
generates a suite of per-plate and summary figures covering colony size, lag time, and
early doubling time distributions and pairwise scatter plots. Unlike the other pipeline
scripts, this function takes the raw `data` struct (not `ht`), uses hard-coded parameter
values in the script body, and does not write figure files to disk or produce a log
file.

## Requirements

- MATLAB (any version with basic Image Processing and Statistics toolboxes). The
  `dscatter` function must be on the MATLAB path.
- The raw `data` struct returned by `heterotyper` (the
  variable produced before running `preprocess_pipeline_data`).

## Usage

Pass the raw `data` struct. The function returns the same struct with a new
`data.across_sample_types` sub-struct appended.

```matlab
data = plot_growth_data_across_sampletypes(data);
```

## Input

| Parameter | Description |
|---|---|
| `data` | Raw data struct from `heterotyper`. Must contain `data.processed`, `data.metadata.original` (with a `Time` column), and per-plate colony fields. |

## Hard-coded Parameters

The following values are set directly in the script body and must be edited manually if
your experiment differs:

| Parameter | Default / Description |
|---|---|
| `incTime` | 20 h — room temperature incubation time added to all lag times |
| `max_lag` | 52 h — maximum lag time used for axis limits and bin ranges |
| `min_col` / `max_col` | 10 / 650 — colony count filter; plates outside this range are excluded |
| `n_groups` / `max_plate` | 4 / 6 — number of time groups and maximum plates per group for subplot grids |
| `max_val` | 4100 px — colony size axis ceiling; must be set after a first exploratory run |
| Time labels | `3H`, `7H`, `24H`, `48H` — matched case-insensitively against `data.metadata.original.Time` |

## Figures Produced

Ten interactive figure windows are opened (not saved to disk). They show: final colony
size vs intensity (scatter), per-plate size / lag time / early doubling time
histograms, pairwise density scatter plots (`dscatter`) for size–lag, size–earlyDT, and
lag–earlyDT, plain scatter equivalents of the same pairs, overlaid distribution line
plots with median markers, and a four-panel per-group median summary.

## Output Struct

The returned `data` struct gains a new sub-struct `data.across_sample_types`
containing:

| Field | Description |
|---|---|
| `median_size` | Per-plate median final colony size, one cell per time group |
| `median_lagtime` | Per-plate median lag time, one cell per time group |
| `median_earlyDT` | Per-plate median early doubling time, one cell per time group |
| `col_nr` | Per-plate colony count, one cell per time group |
| `sample_types` | Cell array of group labels: `{'3H', '7H', '24H', '48H'}` |
| `p_values` | 4×4 matrix of t-test p-values across groups (median values) |

## Notes

This is an exploratory legacy script intended for a quick visual overview. It does not
use the `ht` struct, writes no files to disk, and applies no `exportgraphics`
formatting. For publication-ready figures and statistical output, use the `ht`-based
scripts such as [`plot_combined_samples`](combined-samples.md),
[`plot_CDF_with_statistics`](cdf-with-statistics.md), and
[`plot_correlation_matrix`](correlation-matrix.md). The `max_val` parameter in
particular requires a first run without axis limits to determine the appropriate value
for your data.
