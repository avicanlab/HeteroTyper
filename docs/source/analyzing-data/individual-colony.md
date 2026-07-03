# Visualize individual colonies

`plot_individual_colony.m`

## Purpose

Interactively visualises a single colony from a selected plate. A cropped region
(300 × 300 px) centred on the chosen colony is shown with: a contrast-enhanced
greyscale background, a rainbow-coloured segmentation mask blended on top, and a bright
green outline (2 px dilated perimeter) highlighting the selected colony. Eight
morphological and growth metrics are annotated in the top-right corner of the image.

## Requirements

- MATLAB with Image Processing Toolbox (`bwlabel`, `label2rgb`, `bwperim`, `imdilate`,
  `insertText`).
- Raw `data` struct from `heterotyper`. The
  `colonies.new` struct for the selected plate must be populated (requires
  `growth_quant == 1`).

## Usage

Call the function with the raw `data` struct. Two GUI dialog boxes will appear
prompting for plate and colony indices.

```matlab
plot_individual_colony(data);
```

## Input

| Parameter | Description |
|---|---|
| `data` | Raw data struct from `heterotyper`. |

## Interactive Prompts

| Prompt | Description |
|---|---|
| Plate index (1–N)? | Index of the plate to inspect. |
| Colony number? | Index of the colony within that plate's `region_props` table. |

## Annotation

Eight metrics are printed in the top-right corner of the blended image:

| Label | Metric |
|---|---|
| Lag | Lag time (h) from `colonies.new.lag_time` |
| Size | Final colony size (px) from `timecourse_size_smoothed(end)` |
| Area | Colony area (px) from `region_props` |
| Int | Final colony intensity from `timecourse_intensity_smoothed(end)` |
| IntSize | Intensity / size ratio |
| Ecc | Eccentricity from `region_props` |
| Cir | Circularity from `region_props` |
| Sol | Solidity from `region_props` |

## Notes

No files are written to disk. The crop is 300 × 300 px centred on the colony centroid,
clamped to image boundaries. The greyscale background is contrast-stretched with
`imadjust([0.025 0.15], [])`. The green border is obtained by dilating the colony
perimeter with a disk structuring element of radius 1.
