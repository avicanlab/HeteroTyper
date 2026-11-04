# Manual vs automatic colony counts over time

`plot_timecourse_manual_auto_counts.m`

## Purpose

Validates the automatic (pipeline) colony counts against manual counts, plate by plate,
over the imaging time course. For every plate, the cumulative number of colonies that
have reached a chosen size threshold is computed at each imaging timepoint and combined
with manual counts read from an Excel sheet. The combined table is exported as an Excel
workbook (long and wide format), and a line-dot figure with one subplot per plate shows
the automatic and manual counts side by side.

:::{note}
New in v1.1.0. Use this script to check segmentation and appearance-time detection
against manual plate counts before interpreting appearance-time distributions.
:::

## Requirements

- MATLAB R2023b or later (`tiledlayout` and `sheetnames` required).
- Raw `data` struct returned by `heterotyper`, containing `data.processed` and
  `data.metadata`.
- An Excel file (`.xlsx`, `.xls` or `.xlsm`) with the manual counts in long format: one
  row per plate × timepoint, with a plate-name column, a time (h) column, and a total
  colony count column (e.g. `Position` / `Time` / `Total Counts`).
- `preprocess_pipeline_data` is **not** required. If `ht` is passed, its parameters are
  reused (see Input).

## Usage

Pass the raw `data` struct, and optionally `ht`. All further choices are made
interactively.

```matlab
plot_timecourse_manual_auto_counts(data);       % asks for incTime
plot_timecourse_manual_auto_counts(data, ht);   % reuses ht.params.incTime
```

## Input

| Parameter | Description |
|---|---|
| `data` | Raw data struct returned by `heterotyper`. Must contain `data.processed` (with `colonies.new.timecourse_size_smoothed` and `colonies.new.time_info.elapsed_time_h` per plate) and `data.metadata`. Plate names are taken from `data.metadata.original.Position` (e.g. `Pos0073`). |
| `ht` (optional) | Master struct produced by `preprocess_pipeline_data`. When given, `ht.params.size_threshold` becomes the default size threshold, `ht.params.incTime` is reused without asking, and `ht.params.out_dir` is the starting folder of the output dialog. |

## Interactive Prompts

| Prompt | Description |
|---|---|
| **1) Size threshold** | Count colonies of at least this many px. Enter keeps the default (`ht.params.size_threshold`, or 100 px). |
| **RT incubation time** | Only asked when `ht` is not passed. Added to the **automatic** time axis only (elapsed imaging time + `incTime`). Manual times are used exactly as given in Excel. Enter = 0 h. |
| **Output folder** | Folder dialog for the table, figure and log. Cancel to type a path instead; Enter keeps the default, and a missing folder can be created on the spot. |
| **3) Manual counts file** | File dialog for the Excel file (cancel to type the path). |
| **Sheet** | Numbered list of the workbook's sheets; enter the number of the sheet with the manual counts. |
| **Columns** | Numbered list of the sheet's columns; enter the plate name, time (h) and total count columns. Defaults are guessed from the column names (`Position` / `Pos…` / `plate`, `Time…`, `Total count` / `count`); Enter accepts the default. |
| **Plates to plot** | `[1]` only plates with manual counts (default), or `[2]` all plates with manual and/or automatic counts. |

Example session:

```text
  1) Count colonies of at least how many px? [default = 100]:
  RT incubation time (h) [default = 0]: 3
  Sheet with the manual counts (number): 1
  Plate name (Pos....) column (number) [default = 1, Position]:
  Time (h) column (number) [default = 2, Time]:
  Total colony count column (number) [default = 3, Total Counts]:
  Choice [1/2]: 1
```

## Counting Rules

| Item | Description |
|---|---|
| Automatic count | A colony is counted from the first timepoint at which its smoothed size timecourse reaches ≥ the threshold (same first-crossing rule as Step 1 of [`extract_cell_count_per_time`](cell-count-per-time.md)). Counts are **cumulative**: the total number of colonies on the plate at each timepoint. The `flag_colony_ok` mask is applied when available, as in the QC caching of `preprocess_pipeline_data`. |
| Automatic time | `time_info.elapsed_time_h + incTime` (hours since inoculation). |
| Plate matching | Plate names are normalised to `PosNNNN` on both sides (e.g. `pos73`, `Pos0073_48h_…` → `Pos0073`; numeric positions are zero-padded). |
| Time matching (wide table) | Each manual time is matched to the nearest automatic timepoint of the same plate within half an imaging interval. Unmatched manual times keep their own row with `NaN` for the automatic count. |

## Figure Layout

One figure with one subplot per plate, arranged in the most square grid that fits all
plates (sorted by plate number).

| Element | Description |
|---|---|
| Lines | Line with filled dots per source: Automatic (blue) and Manual (orange), no marker edge. Shared legend at the top. |
| Subplot title | Plate name (e.g. `Pos0073`). |
| Axes | x: Time (h), same range for all plates. y: Total colony count, from 0, scaled per plate. |
| Figure title | `Manual vs automatic colony counts — colonies ≥ <threshold> px`. |

## Outputs

All outputs are written to the chosen output folder:

| Item | Description |
|---|---|
| `plot_timecourse_manual_auto_counts_<thr>px_<timestamp>.xlsx` | Excel workbook with two sheets: `Combined_long` (`Plate`, `Time_h`, `Count`, `Source` = Automatic / Manual) and `Combined_wide` (`Plate`, `Time_h`, `Auto_count`, `Manual_count`, `Auto_minus_Manual`). |
| `plot_timecourse_manual_auto_counts_<thr>px_<timestamp>.fig / .png / .pdf` | Line-dot figure with all plates as subplots. PNG at 300 dpi. |
| `plot_timecourse_manual_auto_counts_log_<timestamp>.txt` | Log file — parameters, input file and sheet, plates combined, a per-plate summary at the last manual timepoint (automatic, manual, and difference), and saved file paths. |

## Typical Workflow

```matlab
% Option A — without preprocessing (asks for incTime)
plot_timecourse_manual_auto_counts(data);

% Option B — after preprocessing, reusing incTime and size_threshold
preprocess_pipeline_data(data);
plot_timecourse_manual_auto_counts(data, ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| `Input must be the pipeline data struct (with .processed and .metadata)` | Pass the raw `data` struct as the first argument, not `ht`. |
| `No plate in data.processed has colonies.new.timecourse_size_smoothed + time_info` | No plate has a size timecourse. Check that growth quantification ran for the dataset. |
| `No rows with a plate name, time and count found in sheet` | The wrong sheet or columns were chosen, or the cells are empty / non-numeric. Check the column numbers at the prompt. |
| `WARNING No automatic counts for manual plate(s)` | Plate names in Excel do not match `data.metadata.original.Position`, or those plates have no timecourse. Check the spelling (`Pos0073` vs `Pos73` is handled; other names are kept as given). |
| Automatic and manual curves are shifted in time | The RT incubation time is added to automatic times only. Make sure the manual times in Excel are also hours since inoculation, or set `incTime` to match. |
| No file dialog appears | MATLAB is running without a display. Cancel/skip — the script falls back to typing the path. |
| `exportgraphics` error | Requires MATLAB R2023b or later. Use `print()` as a fallback on older versions. |