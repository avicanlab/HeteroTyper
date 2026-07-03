# Correlation matrices and non-parametric group comparison

`plot_correlation_matrix.m`

## Purpose

Produces two complementary statistical outputs for selected sample groups and
phenotypic features. Section 1 computes within-group Pearson correlation matrices and
renders them as heatmaps (blue–white–red). Section 2 runs a Kruskal–Wallis omnibus test
followed by Dunn's post-hoc test with Bonferroni correction for every selected feature,
then renders all pairwise adjusted p-values as a combined heatmap figure and exports
them to a CSV table. Data are treated as non-normally distributed throughout —
Kruskal–Wallis and Dunn's test are rank-based and valid for any continuous
distribution.

## Requirements

- MATLAB R2020a or later (`exportgraphics` required). Statistics and Machine Learning
  Toolbox (`kruskalwallis` required for Section 2).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct.

## Usage

Pass the `ht` struct as the only argument. An interactive terminal prompt asks for group
and feature selections before any analysis runs.

```matlab
plot_correlation_matrix(ht);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups`, `ht.params`, `ht.labels`, and `ht.colors`. |

## Interactive Selection

On launch, two prompts are shown in sequence:

- **Groups** → enter labels separated by spaces or commas (e.g. `3h 24h`), or type
  `ALL`.
- **Features** → enter feature numbers (1–10) separated by commas, or type `ALL`.
  Features: 1 Lag Time, 2 Final Size, 3 Area, 4 Intensity, 5 Mean Intensity, 6 Int /
  Size, 7 Perimeter, 8 Circularity, 9 Eccentricity, 10 Solidity.

## Colony Filtering

Colonies with any non-finite value for a selected feature are excluded on a per-feature
basis before computing Pearson r or Kruskal–Wallis. Groups with fewer than 3 colonies
with all-finite values for a given feature are skipped with a warning.

:::{note}
Colony filtering is identical across `preprocess_pipeline_data`, `plot_combined_samples`,
`plot_combined_samples_growth_curves`, `plot_CDF_with_statistics`, and this script. All
scripts operate on the same validated colony set.
:::

## Figures Produced

| Figure | Description |
|---|---|
| **Section 1 — Correlation heatmaps** | One heatmap per selected group. Custom blue–white–red colormap (−1 to +1). Cell labels show r to 2 decimal places. Groups arranged in up to 2 columns. |
| **Section 2 — KW/Dunn combined figure** | One panel per selected feature. Colour encodes −log₁₀(p_adj), saturating at 4 (≤ 0.0001), using a white-to-red palette. Cell text shows significance stars (`*** p<0.001`, `** p<0.01`, `* p<0.05`, `ns`). Panel title shows Kruskal–Wallis H and p. Shared colorbar on the right. |

## Outputs

| Item | Description |
|---|---|
| `CorrelationMatrix_<groups>_<timestamp>.fig / .pdf / .png` | Section 1 Pearson heatmap figure |
| `KWDunn_<groups>_<timestamp>.fig / .pdf / .png` | Section 2 Kruskal–Wallis / Dunn combined figure |
| `KWDunn_pairwise_<groups>_<timestamp>.csv` | CSV table with one row per feature: `Feature`, `KW_H_stat`, `KW_p_value`, then one column per group pair showing Bonferroni-adjusted Dunn p-values. |
| `plot_correlation_matrix_<timestamp>.txt` | Log file — selected groups and features, per-group n, Kruskal–Wallis H / p, and all pairwise Dunn z / p_adj values. |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — run correlation and group comparison
plot_correlation_matrix(ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| Group skipped: "insufficient data for correlation" | Fewer than 3 colonies had finite values across all selected features for that group. Reduce the feature selection or check QC thresholds. |
| KW / Dunn result is `NaN` for a feature | All selected groups had fewer than 3 finite values for that feature. The CSV row will contain `NA`. |
| `exportgraphics` error | Requires MATLAB R2020a+. Use `print()` as a fallback on older versions. |
| `kruskalwallis` not found | Requires the Statistics and Machine Learning Toolbox. |
