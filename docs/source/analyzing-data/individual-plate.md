# Visualize individual plates

`plot_individual_plate.m`

## Purpose

Interactively visualises a single plate from the raw `data` struct. The plate image is
shown alongside its per-colony growth curves. Three colour modes are available:
rainbow segmentation (assigns a random jet colour to each colony mask), lag time
(parula palette from `incTime` to `max_lag`), and final colony size (parula palette
from 0 to max observed size). Colony indices are overlaid on the image. This script
works directly on the raw `data` struct and does not require `ht`.

## Requirements

- MATLAB (any version with Image Processing Toolbox).
- Raw `data` struct from `heterotyper`.

## Usage

Call the function with the raw `data` struct. Two GUI dialog boxes will appear
prompting for plate index and colour mode.

```matlab
plot_individual_plate(data);
```

## Input

| Parameter | Description |
|---|---|
| `data` | Raw data struct from `heterotyper`. |

## Interactive Prompts

| Prompt | Description |
|---|---|
| Plate position (1–N)? | Index of the plate to visualise. Must be a valid index in `data.processed`. |
| Color by (0) segmentation, (1) lag-time, (2) col size? | `0`: rainbow jet segmentation. `1`: parula palette mapped to lag time (`incTime` to `max_lag`). `2`: parula palette mapped to final colony size (0 to plate maximum). |

## Hard-coded Parameters

| Parameter | Default / Description |
|---|---|
| `incTime` | 20 h — added to raw lag times for display |
| `max_lag` | 52 h — upper bound of the lag time colour scale and growth curve x-axis |

## Figure Layout

A single figure window (1200 × 400 px) with three columns: the plate image occupies
columns 1–2 (with the colony mask overlay and colony indices labelled), and column 3
shows growth curves for all QC-passing colonies on that plate. If `growth_quant == 0`
for the selected plate, the growth curve panel is left empty.

## Notes

No files are written to disk. Figures are displayed interactively only. Colour-mode 0
(segmentation) labels colonies in red; modes 1 and 2 use grey labels for readability
over the coloured mask.
