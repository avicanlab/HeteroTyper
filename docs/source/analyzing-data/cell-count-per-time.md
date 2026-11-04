# Colony appearance counts and early vs late appearing colonies

`extract_cell_count_per_time.m`

## Purpose

Runs one of two independent analyses of colony appearance, chosen at the start:

- **Step 1 — Colony appearance counts per timepoint.** For each selected sample group,
  counts how many colonies newly reach a pixel-size threshold at each imaging
  timepoint, i.e. the timepoint at which each colony's smoothed size timecourse first
  crosses the threshold. Each colony is counted once, so the counts summed over all
  timepoints reproduce the valid-colony count used by `plot_combined_samples`.
- **Step 2 — Early vs late appearing colonies.** Within one sample group, colonies are
  binned by the EARLY and LATE timepoints at which they first appeared, and the bins
  are compared by growth rate, doubling time, and size growth before and after
  appearance (seven figures, with statistics in the log).

Only the chosen step is executed.

:::{note}
New in v1.1.0.
:::

## Requirements

- MATLAB R2023b or later (`exportgraphics` required) and "Statistics and Machine
  Learning Toolbox" (`lillietest`, `ttest2`, `ranksum`, `kstest2`, `ttest`,
  `signrank`, `tcdf`) for Step 2.
- `preprocess_pipeline_data(data);` must be run first to generate the `ht` struct,
  including `ht.groups(g).size_timecourse` and `ht.groups(g).time_h`.
- Step 2 additionally requires `ht.groups(g).growth` (Gompertz fits, STEP 6d) and
  `ht.groups(g).doublingTime`.

## Usage

Pass the `ht` struct as the only argument. The function first asks which step to run,
then opens the prompts for that step.

```matlab
extract_cell_count_per_time(ht);
```

## Input

| Parameter | Description |
|---|---|
| `ht` | Master struct produced by `preprocess_pipeline_data`. Must contain `ht.groups` (with `size_timecourse`, `time_h`, `size`, `n_colonies`; for Step 2 also `growth.mu_max` and `doublingTime`), `ht.params` (`out_dir`, `size_threshold`), and `ht.labels`. |

## Interactive Selection

| Prompt | Step | Description |
|---|---|---|
| **Run which step?** | – | `1` = appearance counts per timepoint, `2` = EARLY vs LATE comparison. |
| **Groups** | 1 | Enter group labels separated by spaces (e.g. `3h 24h`) or `ALL`. |
| **Sample group** | 2 | Enter exactly **one** group label. The table shows `n_col` and `n_fitted` (successful Gompertz fits) per group. |
| **Colony size threshold** | 1, 2 | A colony "appears" the first time its smoothed size timecourse reaches this many px. Enter keeps the default (`ht.params.size_threshold`). |
| **Number of EARLY / LATE timepoints** | 2 | How many EARLY and how many LATE timepoints to compare (default 2 each; need not be equal, e.g. 1 vs 3). Before this prompt, a table of newly appearing colonies per timepoint index (`idx`, `time (h)`, `n_appear`) is printed to help choose. |
| **EARLY / LATE timepoint indices** | 2 | Timepoint indices into the group's time grid (e.g. EARLY `2 4`, LATE `10 14`). EARLY and LATE indices must not overlap; a warning is shown if EARLY timepoints are not strictly earlier than LATE ones. |

Example session (Step 2):

```text
  Run which step? [1/2]: 2
  Enter ONE group label: 24h
  Threshold in px [default = 100, press Enter to keep]:
  How many EARLY timepoints? [default = 2]:
  How many LATE timepoints? [default = 2]: 3
  EARLY timepoint indices (2 needed, e.g. 2 4): 6 7
  LATE timepoint indices  (3 needed, e.g. 10 14): 18 20 22
```

:::{warning}
If the appearance threshold is set **below** `ht.params.size_threshold`, colonies whose
final size was below the preprocessing filter were already NaN-masked across their
entire timecourse in `preprocess_pipeline_data` and cannot be recovered here. Re-run the
preprocessor with a lower `cfg.size_threshold` if those colonies are needed.
:::

## Step 1 — Appearance Counts

For each selected group, the log reports the total, valid and counted colonies, and a
table per timepoint with columns `Time (h)`, `New` (colonies first reaching the
threshold at that timepoint), `Cumulative`, and `% of valid`. Valid colonies that never
reach the threshold in their smoothed timecourse are reported in a note and excluded
from the counts. When all selected groups share the same time grid, a combined table
(new colonies per timepoint, one column per group) is added. Step 1 writes the log
only; no figures are produced.

## Step 2 — Early vs Late Appearing Colonies

Each selected timepoint defines a bin (`EARLY-1`, `EARLY-2`, …, `LATE-1`, …) holding
the colonies that **first** reached the threshold at that timepoint. The log starts
with a bin summary: `n_appear`, colonies with a growth rate (`n_mu`) and doubling time
(`n_dt`), and the plate composition of every bin.

:::{note}
Growth rate and doubling time are **not** recomputed. They are the per-colony values
from `preprocess_pipeline_data` (μₘₐₓ from the Gompertz fit over the full size
timecourse; doubling time from the log-linear fit over the post-appearance timecourse).
Binning only chooses *which* colonies are compared, so Figures 1–2 compare colonies
that appeared early or late by their whole-trajectory phenotype. Figures 3–7 use the
sizes measured at specific timepoints around each colony's appearance.
:::

| Figure | Description |
|---|---|
| **Figure 1 — Growth rate** | Jittered scatter of μₘₐₓ (px/h) per bin, x = timepoint (h) with `n` per bin, colour = EARLY (blue) / LATE (red), median marker per bin. Significance brackets for every pairwise bin comparison (narrowest first) and a bold bracket for the pooled EARLY vs pooled LATE contrast on top. |
| **Figure 2 — Doubling time** | Same layout and statistics for doubling time (h). |
| **Figure 3 — Size growth over +5 h** | For the same colonies in each bin: size at their appearance timepoint next to size ≈ 5 h later (paired points joined by a line per colony). One pooled OLS slope (px/h) for the EARLY bins and one for the LATE bins, each with equation, Pearson r, R² and p drawn on the plot, plus the EARLY vs LATE slope comparison. |
| **Figures 4, 5, 6 — Three timepoints per bin** | Size at X − L, X (appearance) and X + A for the same colonies (only colonies valid at all three): Figure 4 = X−5h, X, X+5h; Figure 5 = X−3h, X, X+5h; Figure 6 = X−2h, X, X+2h. Two pooled OLS slopes per cohort — pre-appearance (X−L → X, dotted) and post-appearance (X → X+A) — with side-panel bars comparing EARLY vs LATE steepness. |
| **Figure 7 — Pre-appearance windows** | X−3h → X and X−2h → X, each analysed on its own as two subpanels of one figure with a shared y-axis. Per window and cohort: pooled OLS slope, EARLY vs LATE slope comparison, and EARLY vs LATE comparison of per-colony slopes. Each subpanel has its own steepness bar panel (slope ± SE with p). |

A bin is excluded from Figures 3–7 when the required timepoints fall outside the imaging
window (noted in the log); a figure is skipped when no bin qualifies.

## Statistics

| Comparison | Test |
|---|---|
| Bin vs bin, pooled EARLY vs pooled LATE (Figures 1–2; +5 h size gain in Figure 3) | Independent samples. Per pair, Lilliefors normality test on both samples: both normal → Welch's t-test (`ttest2`, unequal variances); otherwise Wilcoxon rank-sum (`ranksum`). A two-sample KS test (`kstest2`) is reported alongside for distribution shape. Pairwise p-values (and separately the KS p-values) are BH-FDR corrected per feature. |
| Size at X vs X + 5 h (Figure 3) | Matched pairs (same colonies). Lilliefors on the paired differences → paired t-test (`ttest`), otherwise Wilcoxon signed-rank (`signrank`). |
| Pooled slope per cohort (Figures 3–7) | Ordinary least squares (`polyfit`, degree 1) of colony size vs time elapsed since each colony's own window start, pooled over all bins of the cohort; reported as equation, Pearson r, R², p. |
| EARLY vs LATE slope (Figures 3–7) | Independent colonies → t-test on the difference of the two slopes (`compare_slopes`: Δslope, SE, t, df, p). |
| Pre vs post slope within a cohort (Figures 4–6) | Same colonies → each colony's own pre and post slope compared with the paired test (Lilliefors → paired t-test or signed-rank). |
| EARLY vs LATE per-colony slopes (Figure 7) | Lilliefors-gated Welch t-test / rank-sum, plus KS. |

Significance stars: `*** p<0.001`, `** p<0.01`, `* p<0.05`, `ns`.

:::{warning}
Colonies in a bin can come from several plates (biological/technical replicates) and
are treated as independent draws; plate of origin is not modelled. Check the plate
composition in the bin summary before over-interpreting a significant p-value — if a
bin comes from one plate, or is dominated by one plate, the difference may be a batch
effect. Colonies are also binned by the timepoint at which they first reach the
threshold, so their size at X is ≥ threshold by construction: pre-appearance slopes
partly reflect this selection. Compare cohorts and windows rather than reading the
absolute px/h as an unbiased growth rate.
:::

## Outputs

All outputs are written to `ht.params.out_dir` with a timestamp suffix (`<group>` =
the selected group label):

| Item | Description |
|---|---|
| `extract_cell_count_per_time_Step<1\|2>_<timestamp>.txt` | Log file. Step 1: per-group appearance count tables and the combined table. Step 2: parameters, bin summary with plate composition, all statistics tables, and saved file paths. |
| `GrowthRate_EarlyLate_<group>_<timestamp>.fig / .pdf / .png` | Step 2, Figure 1 |
| `DoublingTime_EarlyLate_<group>_<timestamp>.fig / .pdf / .png` | Step 2, Figure 2 |
| `SizeGrowth5h_EarlyLate_<group>_<timestamp>.fig / .pdf / .png` | Step 2, Figure 3 |
| `SizeGrowth3pt_Xminus<L>h_Xplus<A>h_EarlyLate_<group>_<timestamp>.fig / .pdf / .png` | Step 2, Figures 4–6 (`Xminus5h_Xplus5h`, `Xminus3h_Xplus5h`, `Xminus2h_Xplus2h`) |
| `SizeGrowth_PreWindows_EarlyLate_<group>_<timestamp>.fig / .pdf / .png` | Step 2, Figure 7 |

PDFs are vector; PNGs are raster at 300 dpi.

## Typical Workflow

```matlab
% Step 1 — run preprocessor once per session (must include STEP 6d for Step 2)
preprocess_pipeline_data(data);

% Step 2 — appearance counts per timepoint (choose 1), then inspect the log
extract_cell_count_per_time(ht);

% Step 3 — compare early vs late appearing colonies in one group (choose 2),
%          using the count table from the previous run to pick timepoints
extract_cell_count_per_time(ht);
```

## Troubleshooting

| Error / Symptom | Resolution |
|---|---|
| `Input must be the ht struct produced by preprocess_pipeline_data` | Pass `ht`, not `data`. Run `preprocess_pipeline_data(data)` first. |
| `Step 2 requires ht.groups(g).growth` | Re-run `preprocess_pipeline_data` — STEP 6d (Gompertz fitting) must complete before Step 2. |
| `Selected group "…" has no size_timecourse data` | The group has no colonies with a size timecourse. Choose another group or check preprocessing. |
| "NOTE: N valid colonies never reached … px" | The threshold is above some colonies' smoothed maximum size. They are excluded from the counts; lower the threshold if needed. |
| `WARNING: … has n_mu=…, n_dt=…` | A bin has fewer than 3 colonies with that feature; its statistics are unreliable. Choose timepoints with more appearing colonies (see the `n_appear` table). |
| "no bin has a usable +5h timepoint" / Figure 3–7 skipped | The selected timepoints are too close to the start or end of the imaging window. Choose timepoints further from the edges. |
| `exportgraphics` error | Requires MATLAB R2023b or later. Use `print()` as a fallback on older versions. |
