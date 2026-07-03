# Export data to directory from MATLAB — `export_data_to_output_directory.m`

## Purpose

Exports the complete morphology dataset to Excel (`.xlsx`) and `.mat` files. Two tables
are produced: Table 1 contains all colonies that passed the size filter (the same
population used by `plot_combined_samples`); Table 2 is a subset restricted to colonies
that also have a finite intensity-per-size value. All data are read directly from
`ht.groups` — no re-filtering from the raw `data` struct is performed.

## Requirements

- `preprocess_pipeline_data(data);` must have been run first.
- Raw `data` struct is needed only to resolve per-colony plate identifiers and
  filenames from `data.metadata`.

## Usage

No interactive prompts. Runs automatically and saves all outputs.

```matlab
export_data_to_output_directory(data, ht);
```

## Colony population

Colonies in Table 1 satisfy: `~isnan(grp.size)`, i.e. they passed the size threshold set
in `preprocess_pipeline_data` (default 100 px) and were not `NaN`-masked for any other
reason. This is identical to the population shown in `plot_combined_samples`. Censored
lag-time entries are included as `NaN` in the `LagTime` column. Table 2 additionally
requires `isfinite(grp.int_per_size)`.

## Output directory

All files are saved to a dedicated subfolder inside `ht.params.out_dir` so the main
output directory is not cluttered:

- `ht.params.out_dir/Exported Data/`

## Outputs

| Item | Description |
|---|---|
| `Morphology_AllParameters_<timestamp>.xlsx` | Table 1 — all size-filtered colonies. Columns: `SampleGroup`, `Name`, `Folder`, `PlateNo`, `LagTime`, `FinalSize`, `Area`, `Intensity`, `Perimeter`, `Circularity`, `Eccentricity`, `Solidity`, `Centroid_X`, `Centroid_Y`. |
| `Morphology_AllParameters_IntPerSize_<timestamp>.xlsx` | Table 2 — subset with finite `IntPerSize`. Same columns as Table 1 plus `IntPerSize`. |
| `allGroupsTable_<timestamp>.mat` | Single `.mat` file containing both tables as `allGroupsTable_Bright` and `allGroupsTable_Bright_2`. Also assigned to the MATLAB base workspace. |
| `export_data_to_output_directory_<timestamp>.txt` | Log file — per-group colony counts and output file paths. |

## Typical Workflow

```matlab
% Run once per session
preprocess_pipeline_data(data);

% Export morphology tables
export_data_to_output_directory(data, ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| Table 1 is empty | No colonies passed the size filter. Check `p.size_threshold` in `ht.params` and re-run `preprocess_pipeline_data` if necessary. |
| Table 2 is smaller than expected | Colonies without a finite intensity-per-size value are excluded. This is normal when intensity data is unavailable for some plates. |
| `PlateNo` column shows wrong values | The `Position` column in `data.metadata.original` is used for plate identifiers. Verify that the metadata table contains the correct `Position` values. |
