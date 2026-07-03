# Comparing counts and visualizing individual colony size distributions

`plot_QC_full_dataset.m`

## Purpose

Produces four diagnostic figures that verify segmentation quality, colony counts, and
size distributions across the full plate dataset of a HeteroTyper experiment. All
filtering parameters, per-plate colony counts, area distributions, and manual reference
counts are read directly from the `ht` struct produced by `preprocess_pipeline_data` —
no recomputation is performed at plot time. Raw plate images are only accessed for the
optional image montage (Plot 3).

## Requirements

- MATLAB R2020a or later (`exportgraphics` required).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct,
  including the `ht.qc` subfield.

## Usage

Pass the `ht` struct as the first argument. `data` is optional and only required for
Plot 3. No interactive prompts are shown. The function runs automatically and saves all
outputs to `ht.params.out_dir`.

```matlab
% Plots 1, 1A, and 2 only (no image montage)
plot_QC_full_dataset(ht);

% All four plots including the image montage
plot_QC_full_dataset(ht, data);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.qc`, `ht.params`, `ht.groups`, `ht.labels`, and `ht.colors`. If `ht.qc` is missing, re-run `preprocess_pipeline_data` with the current version of the script. |
| `data` | *(Optional)* Raw data struct from HeteroTyper acquisition. Required only for Plot 3 (image montage). Omit to run Plots 1, 1A, and 2 without loading images into memory. |

## Colony filtering

Colony counts and plate selection are fully consistent with `preprocess_pipeline_data`.
No separate thresholds are applied inside this function. Three distinct count metrics
are stored in `ht.qc` and used at different points in the figures:

| Parameter | Description |
|---|---|
| `colony_count_all` | Total objects returned by the region detector — length of `flag_colony_ok`. Includes objects that will later be rejected. |
| `colony_count_clean` | Objects that passed `flag_colony_ok` (non-zero). Used in Plot 1 and annotated on the image montage tiles. |
| `colony_count_final` | Objects from `colonies.new` whose area is `>= size_threshold`. Matches exactly the colony population analysed in `plot_combined_samples`. Used as the automated count in Plot 1A. |

A plate is included in Plot 2 (colony-size histograms) only when both of the following
are true:

- `data.processed{i}.growth_quant == 1`
- `colony_count_final` is within `[min_col, max_col]` from `ht.params`

:::{note}
Colony filtering is identical across `preprocess_pipeline_data`, `plot_combined_samples`,
`plot_combined_samples_growth_curves`, `plot_CDF_with_statistics`,
`plot_KS_quantile_statistics`, and this script. All scripts operate on the same
validated colony set.
:::

## Figures

Four figures are produced. All are saved to `ht.params.out_dir` as `.fig`, `.pdf`
(vector), and `.png` (300 dpi). Figure windows are sized to 85% of the screen.

| Figure | Description |
|---|---|
| **Plot 1 — Colony counts** | Semi-log scatter of all three count metrics (`colony_count_all`, `colony_count_clean`, `colony_count_final`) against plate index. Grey vertical lines connect final to all counts per plate. Dashed horizontal lines mark `min_col` and `max_col` thresholds. Covers all plates regardless of QC pass/fail. |
| **Plot 1A — Automated vs manual counts** | Only produced if `ht.qc.manual_counts` contains non-zero entries. Two panels: (left) semi-log overlay of automated and manual counts per plate; (right) scatter of automated vs manual on a linear scale, with 95% confidence band, regression line, and annotation showing Spearman ρ, p, R², and the linear equation. Dots are coloured by sample group. |
| **Plot 2 — Colony size distributions** | One subplot per QC-passing plate arranged in the most square grid possible. Each subplot shows a filled histogram (bin width 100 px, fill colour = sample group colour, step outline = 60% darkened) of final colony areas. All subplots share a common x-axis limit (global maximum across all passing plates, rounded up). Subplot title shows the plate filename from `ht.qc.fn` if available, otherwise the plate index. |
| **Plot 3 — Image montage** | Only produced when `data` is supplied. One subplot per plate arranged in the most square grid possible (`ceil(√nr_plates)` columns). Each tile is a 1001×1001 px crop centred at image position `[1500, 1500]`, blending the segmentation mask with a contrast-stretched greyscale background. The final colony count (`colony_count_final`) is annotated in the top-left corner of every tile using normalised axes coordinates. |

## Outputs

All outputs are written to `ht.params.out_dir` with a timestamp suffix. Each figure is
saved as three files:

| Item | Description |
|---|---|
| `QC_colony_counts_<timestamp>.fig / .pdf / .png` | Plot 1 — Colony count semi-log scatter (all plates) |
| `QC_count_comparison_<timestamp>.fig / .pdf / .png` | Plot 1A — Automated vs manual count comparison. Only written when manual counts are available. |
| `QC_size_distributions_<timestamp>.fig / .pdf / .png` | Plot 2 — Per-plate colony size histogram grid. Only written when at least one plate passes QC. |
| `QC_montage_<timestamp>.fig / .pdf / .png` | Plot 3 — Full-dataset image montage. Only written when `data` is supplied. |
| `plot_QC_full_dataset_<timestamp>.txt` | Log file — QC filter thresholds, plates passing QC, per-plate colony counts (all / clean / final), and automated vs manual count statistics (Pearson r, Spearman ρ, R², regression equation). |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — QC plots without image montage
plot_QC_full_dataset(ht);

% Step 3 — re-run with image montage if needed
plot_QC_full_dataset(ht, data);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| `ht.qc` not found | Re-run the current version of `preprocess_pipeline_data`. The `ht.qc` sub-struct was added in a later release; older `ht` structs will not have it. |
| Plot 1A not produced | `ht.qc.manual_counts` is empty or all zeros. Ensure the metadata column containing manual counts was selected during preprocessing, and that the column has at least one non-zero value. |
| Plot 2 not produced / "No plates pass QC filter" | All plates failed the colony count filter (`[min_col, max_col]`). Review the thresholds in `preprocess_pipeline_data` or check that `growth_quant == 1` for the relevant plates. |
| Plot 3 not produced / "data not supplied" | Call the function as `plot_QC_full_dataset(ht, data)` rather than `plot_QC_full_dataset(ht)`. |
| Image crop is entirely black | The `imcrop` region `[1500, 1500, 1000, 1000]` falls outside the plate image boundaries. The plate images may be smaller than expected. Adjust the crop coordinates in the script if needed. |
| `exportgraphics` error | Requires MATLAB R2020a+. Use `print()` as a fallback on older versions. |
