# Export data to directory from MATLAB — `export_colony_images_img_intPerSize_threshold.m`

## Purpose

Exports colony images split into two categories — "Valid colonies" and "Filtered
colonies" — based on intensity-per-size (`IntPerSize`) thresholds chosen interactively.
The user selects whether to apply a lower bound only, an upper bound only, or both
bounds simultaneously. Images receive a green border around the selected colony and
nine annotated metrics. Uses the same colony population as `plot_combined_samples`; no
re-filtering from raw `data`.

## Requirements

- MATLAB Image Processing Toolbox (`imfuse`, `imadjust`, `imdilate`, `bwperim`,
  `bwlabel`, `label2rgb`, `insertText`).
- `preprocess_pipeline_data(data);` must have been run first.
- Raw `data` struct is needed for plate images, segmentation masks, and the binary
  colony mask.

## Usage

Two to three interactive prompts ask for the threshold mode and bound value(s).

```matlab
export_colony_img_intPerSize_threshold(data, ht);
```

## Interactive Prompts

| Prompt | Description |
|---|---|
| Threshold mode (1 / 2 / 3) | `1` = lower bound only (`IntPerSize < L`). `2` = upper bound only (`IntPerSize > H`). `3` = both bounds (`IntPerSize < L` or `IntPerSize > H`). |
| Lower bound `L` | Asked when mode is 1 or 3. Colonies with `IntPerSize < L` are placed in "Filtered colonies". |
| Upper bound `H` | Asked when mode is 2 or 3. Colonies with `IntPerSize > H` are placed in "Filtered colonies". When mode is 3, `H` must be greater than `L`. |

## Colony population and categorisation

The base population is the same as `plot_combined_samples`: non-`NaN` size AND finite
`IntPerSize`. Within this population, colonies that fall outside the chosen
threshold(s) are "Filtered colonies"; the remainder are "Valid colonies". Both
categories are exported.

## Image format

Same as `export_colony_images`: 301×301 px crop, contrast-stretched greyscale blended
with the jet-coloured segmentation mask. Additionally, a 2 px green border is drawn
around the specific colony of interest using `bwperim` and `imdilate`. Nine metrics are
annotated in the top-right corner.

## Output directory

Images are saved inside a threshold-specific subfolder so multiple runs with different
thresholds coexist without overwriting:

- `ht.params.out_dir/Colony Images - IntensitySize threshold/<threshold description>/<category>/<SampleGroup>/`

The threshold description encodes the chosen mode, for example: `IntPerSize_lt_25.00`,
`IntPerSize_gt_75.00`, or `IntPerSize_lt_25.00_or_gt_75.00`.

## Outputs

| Item | Description |
|---|---|
| `IntSize_<value>_<Group>_Plate<N>_Colony<N>.png` | One `.png` per exported colony. The `IntPerSize` value is encoded in the filename for quick sorting. |
| `export_colony_img_intPerSize_threshold_<timestamp>.txt` | Log file — chosen thresholds, per-plate valid and filtered colony counts. |

## Typical Workflow

```matlab
% Export colonies split by IntPerSize threshold
export_colony_img_intPerSize_threshold(data, ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| No images in "Filtered colonies" folder | No colonies fell outside the chosen threshold(s). Try widening the bounds. |
| "H must be greater than L" warning loops | In mode 3, the upper bound must exceed the lower bound. Re-enter valid values. |
| Green border missing on some colonies | The green border requires a valid `labeled_mask` entry for the colony. If `ok_idx` mapping fails (e.g. region_props mismatch), the image is still saved without the border. |
| Images from different threshold runs mixed together | Each run creates its own subfolder named by threshold description. Old runs are not overwritten; check that you are looking in the correct subfolder. |
