# Export data to directory from MATLAB — `export_colony_images.m`

## Purpose

Exports cropped colony images (blended raw greyscale + colour segmentation mask) for
colonies whose lag time exceeds a user-defined threshold. Images are annotated with
nine per-colony metrics. Uses the same colony population as `plot_combined_samples` —
filtering is applied consistently via `ht.groups`; no re-filtering from raw `data`.

## Requirements

- MATLAB Image Processing Toolbox (`imfuse`, `imadjust`, `insertText`).
- `preprocess_pipeline_data(data);` must have been run first.
- Raw `data` struct is needed for plate images and segmentation masks.

## Usage

One interactive prompt asks for the lag-time threshold before export begins.

```matlab
export_colony_images(data, ht);
```

## Interactive Prompt

| Prompt | Description |
|---|---|
| Lag-time threshold (h) | Only colonies with `lag_time > threshold` are exported. Enter a non-negative number (e.g. 38). The threshold is logged and saved with the output. |

## Colony population

Colonies must satisfy: `~isnan(grp.size)` (passed size filter, same as
`plot_combined_samples`) AND `isfinite(grp.lag_time)` AND `lag_time > threshold`. No
eccentricity filter is applied, consistent with `plot_combined_samples`.

## Image format

Each image is a 301×301 px (max) crop centred on the colony centroid, blending a
contrast-stretched greyscale background with the colour segmentation mask. Nine metrics
are printed in the top-right corner: `Lag`, `Size`, `Int`, `MeanInt`, `IntSize`, `Cir`,
`Ecc`, `Sol`, `Peri`.

## Output directory

Images are saved to a dedicated subfolder to keep the main output directory clean:

- `ht.params.out_dir/Colony Images/<SampleGroup>/`

## Outputs

| Item | Description |
|---|---|
| `<SampleGroup>_Plate<N>_Colony<N>.png` | One `.png` per exported colony, named by group, plate position, and colony index. |
| `export_colony_images_<timestamp>.txt` | Log file — chosen threshold, per-plate exported colony counts. |

## Typical Workflow

```matlab
% Export colonies with lag time above threshold
export_colony_images(data, ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| No images exported | No colonies exceed the lag-time threshold for any plate. Try a lower threshold. |
| `insertText` error | Requires Computer Vision Toolbox or Image Processing Toolbox with `insertText` support. |
| Crop is blank or misaligned | Centroid coordinates in `ht.groups` are `NaN` for size-rejected colonies. Only non-`NaN` centroid entries are exported — check that valid colonies exist for that plate. |
