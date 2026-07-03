# Morphology (lag time and colony size) distributions per plate with Gini coefficient

`plot_morphology_colonies_Gini.m`

## Purpose

Plots violin and scatter distributions of two colony phenotypes, lag time (h) and final
colony size (px), per plate, grouped by sample time group. Each panel corresponds to
one sample time group and shows every plate in that group along the x-axis, labelled
with the biological/technical replicate identifiers chosen at the interactive prompt
(reused automatically if already cached by `plot_lagTime_colonies_Gini`). A Gini
coefficient is annotated above each plate's violin to quantify within-plate
heterogeneity, and both a Spearman and a Pearson correlation between colony count (N)
and Gini value are reported per group and pooled across groups. Four figures are
produced in total: lag time and colony size, each with a Spearman and a Pearson
variant.

## Requirements

- MATLAB R2020a or later (`exportgraphics` and "Statistics and Machine Learning
  Toolbox" required).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct.

## Usage

Pass the `ht` struct and your raw `data` struct as arguments. An optional output struct
captures per-group lag time and colony size statistics.

```matlab
stats = plot_morphology_colonies_Gini(ht, data);
```

On first call, the function prompts for the same three metadata columns as
`plot_lagTime_colonies_Gini` (Steps A–C below). Selections are cached in
`ht.plate_labels`, so subsequent calls to either function skip the prompt. To
re-select, run:

```matlab
ht = rmfield(ht, 'plate_labels');
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups`, `ht.params`, `ht.labels`, `ht.colors`, `ht.qc`, and `ht.metadata`. If `ht.per_plate_size` is present it is used directly; otherwise it is built automatically from `data`. |
| `data` | Raw data struct from HeteroTyper acquisition (the variable returned by `heterotyper`). Required for per-plate lag time and colony size extraction and metadata access. |

## Interactive label selection

All available metadata columns are listed with a preview of unique values. Three
selections are made in sequence:

- **Step A:** Biological replicate column (e.g. `Set` → `Set1`, `Set2` …). Enter `0` to
  use plate index numbers only.
- **Step B:** Technical replicate column (e.g. `Replicate` → `R1`, `R2`, `R3` …). Enter
  `0` to skip.
- **Step C:** Dilution column (e.g. `Dilution` → `-4(4)`, `-5` …). Enter `0` to skip.

Labels are combined as `BioRep_TechRep_Dilution` (e.g. `Set1_R2_-4(4)`). Any skipped
step is omitted from the label. A preview of the first eight plate labels is printed to
the Command Window before plotting begins.

## Colony Filtering

Per-plate lag vectors are extracted using the same QC thresholds as
`plot_lagTime_colonies_Gini`. A plate is included in the lag-time panels only when all
three criteria below are met. Per-plate colony size is read from `ht.per_plate_size`
when available; if that cache is absent or empty the function rebuilds it from
`raw_data`, applying the size threshold in `ht.params.size_threshold` and excluding
colonies with censored lag times:

- `data.processed{i}.growth_quant == 1`
- Number of valid (non-censored) lag times is within `[min_col, max_col]` from
  `ht.params`
- Right-censored lag times (`lag >= max_lag − 0.25 h`) are set to `NaN` and excluded
  from all lag-time statistics and violin rendering; colonies with censored lag times
  are also excluded from the colony-size distribution.

## Statistics — Gini Coefficient, Spearman and Pearson Correlation

| Metric | Description |
|---|---|
| Gini coefficient | Computed per plate, independently for lag time and colony size (0 = perfectly uniform, 1 = maximally heterogeneous). Plates with fewer than 3 valid colonies receive no annotation. Annotated above each violin at a shared per-group label baseline, formatted as `G=0.03 N=161`. |
| Spearman ρ (per group and pooled) | Rank correlation between colony count (N) and Gini value, computed within each time group and again pooling all plates across groups. Shown in the corresponding subplot title. Significance: `*** p<0.001`, `** p<0.01`, `* p<0.05`, `ns p≥0.05`. |
| Pearson r (per group and pooled) | Linear correlation between colony count (N) and Gini value, computed the same way as the Spearman rows but with Pearson's r. Comparing the Spearman and Pearson results checks whether any N–Gini relationship is driven by a few outlier plates (Pearson) or holds more generally across the rank order (Spearman). |

## Outputs

| Item | Description |
|---|---|
| `morphology_Gini_lagTime_Spearman/_Pearson_<timestamp>.fig / .png / .svg / .pdf` and `morphology_Gini_colonySize_Spearman/_Pearson_<timestamp>.fig / .png / .svg / .pdf` | Four multi-panel violin and scatter figures, one per metric–correlation-type combination. Each is saved as a MATLAB figure (`.fig`, re-openable), raster PNG at 300 dpi, vector SVG, and PDF. |
| `GiniTable_lag_<group>` and `GiniTable_size_<group>` (workspace) | Per-plate tables with columns `PlatePos`, `Label`, `N_colonies`, and `Gini_lag` or `Gini_size`. One pair of tables per time group (e.g. `GiniTable_lag_3h`, `GiniTable_size_3h`). |
| `GiniCorr_NvsGini_lag_Spearman`, `_lag_Pearson`, `_size_Spearman`, `_size_Pearson` (workspace) | Four correlation summary tables with columns `Group`, `N_plates`, `r`, and `P_value`. Each includes one row per time group plus a `POOLED` row combining all groups. |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — plot lag time and colony size violin figures with Gini coefficients
stats = plot_morphology_colonies_Gini(ht, data);

% Step 3 — inspect correlation results in workspace
disp(GiniCorr_NvsGini_lag_Spearman);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| All subplots empty, N=0 | The per-plate lag cache (`ht.per_plate`) or size cache (`ht.per_plate_size`) is stale. Run: `ht = rmfield(ht, 'per_plate');` and/or `ht = rmfield(ht, 'per_plate_size');` then re-run the function. The Command Window diagnostic lines show which plates are accepted or rejected and why. |
| X-axis shows `Plate1`, `Plate2` … | Default labels are cached. The function will auto-detect this and re-prompt for metadata columns on the next call. Alternatively: `ht = rmfield(ht, 'plate_labels');` then re-run. |
| No metadata found warning | `ht.metadata` is empty. Re-run `preprocess_pipeline_data` with the patched version of the script, which stores the metadata table in `ht` automatically. |
| `exportgraphics` error | Requires MATLAB R2020a or later. |
