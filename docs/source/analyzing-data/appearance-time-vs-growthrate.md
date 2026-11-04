# Maximum growth rate

`plot_appearanceTime_vs_growthRate.m`

## Purpose

Plots scatter panels of appearance time (h) versus maximum Gompertz growth rate μₘₐₓ
(px/h) for any set of sample groups, and optionally tests selected pairs of groups
against each other. Three figures are produced: colonies coloured by log₁₀(local colony
density) on a shared viridis palette (Figure 1), colonies coloured by sample group
(Figure 2), and one overlay panel per chosen pairwise comparison (Figure 3). A LOESS
smooth is overlaid on each panel when n ≥ 10 fitted colonies. Per-group statistics
(n, Spearman ρ, ρ², density-weighted Spearman ρ𝑤, LOESS polynomial) and all pairwise
comparisons are written to the log and to an Excel workbook. Gompertz fitting and the
per-group statistics are pre-computed in STEP 6d of `preprocess_pipeline_data`; this
script reads them and adds only the pairwise tests.

:::{note}
In v1.1.0 this single script replaces the earlier experiment-specific variants
(`…_6s_stress`, `…_6srna_growth`). Per-group statistics are no longer printed inside
the panels; they are in the log file and in the `Stats` sheet of the workbook.
:::

## Requirements

- MATLAB R2023b or later (`exportgraphics` required). No toolbox function is used.
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct,
  including `ht.groups(g).growth` from STEP 6d (Gompertz fitting).

## Usage

Pass the `ht` struct as the only argument. Three interactive terminal prompts follow
before plotting.

```matlab
plot_appearanceTime_vs_growthRate(ht);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups` (with `ht.groups(g).growth` from STEP 6d), `ht.params`, `ht.labels`, and `ht.colors`. |

## Interactive Selection

Three selections are made in sequence:

- **Groups** → a table lists every group with its total colony count (`n_col`) and the
  number of successful Gompertz fits (`n_fitted`). Enter group labels (comma-separated
  if they contain spaces, e.g. `WT CTRL, 6S CTRL`), index numbers, or `ALL`.
- **Display order** → the selected groups are listed with numbers. Enter a new order
  (e.g. `3 1 2`); groups not listed follow in their current order. Press Enter to keep
  the selection order. Skipped when only one group is selected.
- **Statistical comparisons** → every possible pair of the selected groups is listed
  and numbered. Enter pair numbers (e.g. `1 4 6`), ranges (`1-3`), or `ALL`. Press
  Enter (or type `none`) for no comparisons; Figure 3 is then skipped.

Example session:

```text
  Your selection: 3h, 24h, 48h
  New order as numbers from this list (e.g. 3 1 2; ...) [Enter = keep]:
  Comparisons — pair numbers (e.g. 1 4 6, ranges 1-3), ALL [Enter = none]: 1 3
```

An unrecognised label, index or pair number triggers a warning and re-prompts without
exiting.

## Colony Filtering

Only colonies with a successful Gompertz fit are plotted and tested. A colony is
included if it passes all of:

- `growth.fit_ok == true` (Gompertz converged)
- Finite `mu_max` (maximum growth rate from Gompertz fit)
- Finite `appearance_time` for the same colony

## Pairwise Comparisons

Run only for the pairs chosen at the prompt. Both tests are Benjamini–Hochberg (BH)
adjusted across the chosen pairs, separately per test.

| Test | Description |
|---|---|
| μₘₐₓ — Mann–Whitney U | Two-sided test of μₘₐₓ between the two groups (normal approximation with tie and continuity correction). Effect sizes: rank-biserial r (`> 0`: first group higher) and the median ratio (first / second). |
| ρ — Spearman comparison | Compares the two groups' Spearman ρ (appearance time vs μₘₐₓ) by Fisher z with the Bonett–Wright variance `(1 + ρ²/2) / (n − 3)`. Asks whether the appearance–growth relationship is equally strong in both groups. |

Significance stars: `*** p<0.001`, `** p<0.01`, `* p<0.05`, `ns`.

:::{note}
Colonies of all positions are pooled and treated as independent. With large n, small
differences become significant — read the effect sizes alongside the p-values.
:::

## Figure Layout

Figures 1 and 2 use the same panel grid: one panel per selected group (columns =
`ceil(√n)`, rows filled top-to-bottom, in the chosen display order). The group label is
shown top-left in the group colour.

| Element | Description |
|---|---|
| **Figure 1 — local colony density** | Filled circles, size 27 pt, 80% alpha, coloured by log₁₀(local colony density) on a viridis palette shared across all panels. A shared colorbar (`Colonies (log₁₀)`) is placed to the right of the grid. |
| **Figure 2 — colonies by sample group** | Same grid, axes and LOESS lines as Figure 1, but every colony is drawn in its sample-group colour (`ht.colors`) at 40% alpha. No colorbar. |
| **Figure 3 — pairwise overlay** | One panel per chosen comparison (up to 3 columns). Both groups overlaid in their group colours at 25% alpha, one LOESS line per group in a darker shade of the group colour, and a legend. Title = pair; subtitle = significance stars of the μₘₐₓ and ρ tests (BH-adjusted). A footnote below the grid explains the tests. Skipped when no comparison is chosen. |
| LOESS smooth | Span 0.75, drawn when n ≥ 10 (red in Figures 1–2). A degree-2 polynomial fitted to the smoothed curve (`y = a x² + b x + c`) is reported in the log and workbook. |
| Axes | x: Appearance time (h), ticks every 8 h, shared across all panels and figures. y: Maximum growth rate (px/h), ticks every 100 px/h, shared in Figures 1–2. In Figure 3 each panel has its own y-limit (≈15% headroom for the legend), so Figure 3 panels are **not** directly comparable by height. Major + midpoint minor dashed grid. Mirror ticks on top and right axes. |

## Outputs

All outputs are written to `ht.params.out_dir` with a group-and-timestamp suffix:

| Item | Description |
|---|---|
| `AppearanceTime_vs_GrowthRate_<groups>_<timestamp>.fig / .pdf / .png` | Figure 1 — colonies coloured by local density |
| `AppearanceTime_vs_GrowthRate_byGroup_<groups>_<timestamp>.fig / .pdf / .png` | Figure 2 — colonies coloured by sample group |
| `AppearanceTime_vs_GrowthRate_pairs_overlay_<groups>_<timestamp>.fig / .pdf / .png` | Figure 3 — pairwise overlays (only when comparisons were chosen) |
| `AppearanceTime_vs_GrowthRate_stats_<groups>_<timestamp>.xlsx` | Excel workbook: `Stats` (per group: n_fit, median μₘₐₓ, Spearman ρ / p / stars, ρ², weighted ρ, LOESS polynomial coefficients and equation), `Comparisons` (per chosen pair: n, medians, median ratio, U, rank-biserial r, raw and BH p-values with stars for both tests), and `README` (column descriptions). |
| `plot_appearanceTime_vs_growthRate_<timestamp>.txt` | Log file — selections, per-group statistics table, pairwise comparison table, Figure 3 y-axis limits, and saved file paths |

PDFs are vector; PNGs are raster at 300 dpi; `.fig` files are re-openable and editable.

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session (must include STEP 6d)
preprocess_pipeline_data(data);

% Step 2 — plot appearance time vs growth rate, choose comparisons at the prompt
plot_appearanceTime_vs_growthRate(ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| `ht.groups(g).growth` not found | Re-run `preprocess_pipeline_data` — STEP 6d (Gompertz fitting) must complete before calling this function. |
| "No successful Gompertz fits found" | No colony in the selected groups had a converged Gompertz fit with a finite `mu_max`. Check that imaging duration was sufficient and that colony growth was observed within the acquisition window. |
| LOESS line absent from a panel | Fewer than 10 colonies with successful fits in that group — the LOESS threshold. The polynomial columns in the workbook are empty/NaN for that group. |
| Figure 3 not produced | No pairwise comparison was chosen (Enter at the comparisons prompt). The log notes "Figure 3 skipped". |
| Group label not recognised | Labels containing spaces must be comma-separated (e.g. `WT CTRL, 6S CTRL`). Index numbers from the table also work. |
| `xlsx` not saved | The log shows `WARNING: could not write xlsx` with the exact error. Check that the file is not open in Excel and that the output folder is writable. |
| `exportgraphics` error | Requires MATLAB R2023b or later. Use `print()` as a fallback on older versions. |
