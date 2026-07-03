# Lag time and other phenotypic features for each biological replicate

`plot_combined_samples_biological_replicates.m`

## Purpose

Produces the identical panel layout as [`plot_combined_samples`](combined-samples.md)
(rows = sample groups, columns = selected phenotypic features), but renders one
separate figure per biological replicate instead of a single figure pooling all groups
together. Each replicate's figure shows only the colonies from plates belonging to that
replicate, letting you check visually whether the group-level distributions and their
Gini heterogeneity hold up within each biological replicate individually, rather than
only in the pooled view.

## Requirements

- MATLAB R2020a or later (`exportgraphics` required).
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct.
- `ht.metadata` must be a non-empty table containing a column that identifies the
  biological replicate for each plate.

## Usage

Three interactive selections are made in sequence:

1. The metadata column identifying the biological replicate,
2. Sample groups,
3. Phenotypic features.

```matlab
plot_combined_samples_biological_replicates(ht);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups`, `ht.params`, `ht.labels`, `ht.colors`, `ht.xbins`, `ht.global`, `ht.metadata`, and `ht.qc.colony_count_final`. |

## Interactive selection

On launch, a menu is printed to the terminal showing all available groups and
features. Two selections are made in sequence:

- **Biological replicate column:** Every metadata column is listed with a preview of
  unique values; if `ht.params.bio_rep_col` was already set (e.g. from STEP 6f of
  `preprocess_pipeline_data`), it is offered as the default (press Enter to accept).
- **Groups** → enter labels separated by spaces or commas (e.g. `3h 24h`), or type
  `ALL` (same menu as `plot_combined_samples`)
- **Features** → enter feature numbers separated by commas (e.g. `1 3 5`), or type
  `ALL` (same menu as `plot_combined_samples`)

## Figure Layout

Identical panel geometry and styling to `plot_combined_samples`: filled bar histogram
(probability density, group colour), red median line, mirrored tick marks on top/right
axes. Differences: the figure title reads "Biological Replicate: `<label>`"; each panel
is annotated with "`<group label>` (n=`<count>`)" plus "Gini=`<value>`". If a group has
zero colonies within a given replicate, the panel displays "no data" instead of an
empty histogram.

## Outputs

One figure per biological replicate found among the selected groups, plus one shared
log file, written to `ht.params.out_dir`:

| Parameter | Description |
|---|---|
| `biorep_<label>_features_<ids>_<timestamp>.fig/.pdf/.png` | One set per biological replicate |
| `plot_biorep_samples_<timestamp>.txt` | Log of the chosen replicate column, groups/features, replicate labels found, and per-replicate/group/feature n, median, IQR, Gini |

## Statistics Used

- Median (red line) and Gini coefficient (annotated) are computed per panel,
  identically to `plot_combined_samples`.
- IQR is computed and logged to the text file only — not shown on the figure.
- No cross-replicate hypothesis testing (Pearson r, Spearman rho, KS test,
  Jensen-Shannon divergence) is performed here, and no Excel workbook is produced. That
  analysis remains the responsibility of `plot_combined_samples`'s `ht.biorep_corr`
  block. This function is a descriptive/visual complement, not a statistical test of
  replicate agreement.

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session
preprocess_pipeline_data(data);

% Step 2 — plot phenotypic feature histograms for replicates
plot_combined_samples_biological_replicates(ht);
```

At the prompt, enter groups and features. Example session:

- **Groups** → enter labels separated by spaces or commas (e.g. `3h 24h`), or type `ALL`
- **Features** → enter feature numbers separated by commas (e.g. `1 3 5`), or type `ALL`

## Troubleshooting

| Error Message | Description |
|---|---|
| `ht.metadata` is missing or empty | Re-run `preprocess_pipeline_data` so the metadata table is stored in `ht`. |
| No biological replicates found | The chosen metadata column had no populated values for the selected plates, please pick a different column or add relative information on metadata file. |
| "no data" in a panel | That group had zero valid colonies within that specific replicate, not an error. |
| Group label not recognised | Labels are case-insensitive; the label just has to match one shown in the menu. |
| `exportgraphics` error | Requires MATLAB R2020a+. |
