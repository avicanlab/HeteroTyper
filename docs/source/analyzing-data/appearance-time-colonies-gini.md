# Appearance time distributions per plate with Gini coefficient

`plot_appearanceTime_colonies_Gini.m`

## Purpose

Plots violin and scatter distributions of appearance time (h) per plate, grouped by sample
time group. Each panel corresponds to one sample time group and shows every plate in
that group along the x-axis, labelled with the biological/technical replicate
identifiers chosen at the interactive prompt. A Gini coefficient is annotated above
each plate's violin to quantify within-plate heterogeneity in appearance time, and a Spearman
correlation between colony count (N) and Gini value is reported per group and pooled
across groups to check that Gini is not simply a density artefact.

## Requirements

- MATLAB R2023b or later (`exportgraphics` and "Statistics and Machine Learning
  Toolbox" required).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct.

## Usage

Pass the `ht` struct and your raw `data` struct as arguments. An optional output struct
captures per-group appearance time statistics.

```matlab
stats = plot_appearanceTime_colonies_Gini(ht, data);
```

On first call, the function prompts for three metadata columns (Steps A–C below).
Selections are cached in `ht.plate_labels` so subsequent calls skip the prompt. To
re-select, run:

```matlab
ht = rmfield(ht, 'plate_labels');
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups`, `ht.params`, `ht.labels`, `ht.colors`, `ht.qc`, and `ht.metadata`. |
| `data` | Raw data struct from HeteroTyper acquisition (the variable returned by `heterotyper`). Required for per-plate appearance extraction and metadata access. |

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

Per-plate appearance vectors are extracted from `raw_data.processed` using the same QC
thresholds as `preprocess_pipeline_data`. A plate is included only when all three
criteria are met:

- `data.processed{i}.growth_quant == 1`
- Number of valid (non-censored) appearance times is within `[min_col, max_col]` from
  `ht.params`
- Right-censored appearance times (`appearance >= max_appearance − 0.25 h`) are set to `NaN` and excluded
  from all statistics and violin rendering.

## Statistics — Gini Coefficient and Spearman Correlation

| Metric | Description |
|---|---|
| Gini coefficient | Computed per plate from the appearance time distribution (0 = perfectly uniform, 1 = maximally heterogeneous). Plates with fewer than 3 valid colonies receive no annotation. Annotated above each violin at a fixed y position (`max_appearance + 2 h`) at 30°, formatted as `G=0.03 N=161`. |
| Spearman r₀ (per group) | Rank correlation between colony count (N) and Gini value across all plates within a time group. Shown on the subplot title as "Spearman r₀=0.142, p=0.431 (ns)". Significance: `*** p<0.001`, `** p<0.01`, `* p<0.05`, `ns p≥0.05`. |
| Spearman r₀ (pooled) | Secondary check pooling all plates across all time groups. A non-significant pooled result alongside non-significant per-group results provides the strongest evidence that Gini is density-independent. Exported to workspace as `GiniCorr_NvsGini`. |

:::{note}
**Appearance-time Gini origin (v1.1.0).** The Gini coefficient is G = MD / (2·mean),
where MD is the mean absolute difference between two colonies. Adding a constant to
every colony leaves MD unchanged but inflates the mean. Appearance times in `ht` are
hours since inoculation (= time from imaging start + `incTime`), so the constant RT
incubation time would shrink G towards 0 for any spread. From v1.1.0 the
appearance-time Gini is therefore computed on `appearance_time − ht.params.gini_appearance_origin`,
which defaults to `incTime` (i.e. appearance measured from the start of imaging). Set
`cfg.gini_appearance_origin` in `preprocess_pipeline_data` to `'incTime'` (default),
`'none'` (v1.0.0 behaviour, hours since inoculation), or a number in hours. Gini values
for ratio-scale features (size, area, intensity, …) are unchanged. Expect appearance-time
Gini values to be higher than in v1.0.0 outputs.
:::

## Outputs

| Item | Description |
|---|---|
| `appearanceTime_colonies_Gini_<timestamp>.fig / .png / .svg` | Multi-panel violin and scatter figure. MATLAB figure (`.fig`, re-openable), raster PNG at 300 dpi, and vector SVG for publication. |
| `GiniTable_<group>` (workspace) | Per-plate table with columns `PlatePos`, `Label`, `N_colonies`, `Gini_appearance`. One table per time group (e.g. `GiniTable_3h`, `GiniTable_7h`). |
| `GiniCorr_NvsGini` (workspace) | Correlation summary table with columns `Group`, `N_plates`, `Spearman_rho`, `P_value`. Includes one row per time group plus a `POOLED` row combining all groups. |

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — plot appearance time violin figure with Gini coefficients
stats = plot_appearanceTime_colonies_Gini(ht, data);

% Step 3 — inspect correlation results in workspace
disp(GiniCorr_NvsGini);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| All subplots empty, N=0 | The per-plate cache (`ht.per_plate`) is stale. Run: `ht = rmfield(ht, 'per_plate');` then re-run the function. The Command Window diagnostic lines will show which plates are accepted or rejected and why. |
| X-axis shows `Plate1`, `Plate2` … | Default labels are cached. The function will auto-detect this and re-prompt for metadata columns on the next call. Alternatively: `ht = rmfield(ht, 'plate_labels');` then re-run. |
| No metadata found warning | `ht.metadata` is empty. Re-run `preprocess_pipeline_data` with the patched version of the script, which stores the metadata table in `ht` automatically. |
| `exportgraphics` error | Requires MATLAB R2023b or later. On older versions, the function falls back to `print()` automatically for PNG and SVG output. |
