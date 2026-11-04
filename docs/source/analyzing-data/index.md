# Analyzing Data

## How to use this chapter

This chapter documents every analysis and export script included in the HeteroTyper
pipeline. Each section corresponds to one script and follows the same format: Purpose,
Requirements, Inputs, Outputs, and common errors.

Before running any script in this chapter, run `preprocess_pipeline_data.m` on your
loaded dataset (see [Data preprocessing](../pipeline-setup.md#data-preprocessing)). This
step performs all shared calculations once and stores the results in the `ht` struct,
which almost every script in this chapter reads directly rather than recomputing.

:::{warning}
**Analysing more than one dataset in the same MATLAB session.**
We recommend preprocessing only one dataset per MATLAB session. `preprocess_pipeline_data`
always saves its result to the MATLAB workspace as a variable named `ht`, replacing any
existing `ht` without warning. `plot_appearanceTime_colonies_Gini` and
`plot_morphology_colonies_Gini` also save their plate-label selections and per-plate
caches back to the workspace as `ht` the first time they are run on a dataset.

If you do need to analyse a second dataset in the same session, follow this order for
each dataset before starting the next one:

1. Run the main script and `preprocess_pipeline_data` for the first dataset.
2. Run `plot_appearanceTime_colonies_Gini(ht, data)` and
   `plot_morphology_colonies_Gini(ht, data)` straight away, so that their selections
   and caches are stored in `ht`.
3. Rename the results, for example:

       ht_exp1   = ht;    clear ht
       data_exp1 = data;  clear data

4. Only then run the main script and `preprocess_pipeline_data` for the next dataset.

Afterwards, pass the renamed variables to the analysis scripts, e.g.
`plot_combined_samples(ht_exp1)`. If you later clear a cache from a renamed struct (for
example `ht_exp1 = rmfield(ht_exp1, 'plate_labels')`), the Gini scripts will write a
new `ht` to the workspace again. Rename that one as well, so it does not overwrite or
get mixed up with another dataset's `ht`.
:::

**Scripts that take `ht`.** Most analysis and plotting scripts take `ht` as their only
input, for example `plot_combined_samples(ht)`. A few also require the raw data variable
used during preprocessing, because they need access to raw images or plate metadata not
stored in `ht`: `export_colony_images(data, ht)`,
`export_colony_img_intPerSize_threshold(data, ht)`,
`export_data_to_output_directory(data, ht)`, and the optional image montage in
`plot_QC_full_dataset(ht, data)`. Two further scripts,
[`plot_appearanceTime_colonies_Gini`](appearance-time-colonies-gini.md) and
[`plot_morphology_colonies_Gini`](morphology-colonies-gini.md), take both `ht` and your
raw data variable to build plate-level replicate labels; on first use, each asks which
metadata columns correspond to the biological and technical replicate, then caches the
answer in `ht` so later calls skip the prompt.

**Scripts that take `data` directly.** A small number of exploratory and visualization
scripts, `plot_individual_colony.m`, `plot_individual_plate.m`, and
`plot_growth_data_across_sampletypes.m`, operate directly on the raw `data` variable
rather than `ht`, and prompt you interactively for a plate or colony index to inspect.
These do not require `preprocess_pipeline_data.m` to have been run first. The same
applies to [`plot_timecourse_manual_auto_counts(data)`](timecourse-manual-auto-counts.md),
which compares automatic and manual colony counts per plate; passing `ht` as an optional
second argument only reuses its parameters.

**Interactive selections.** Several scripts, including `plot_combined_samples.m`,
`plot_appearanceTime_vs_morphology.m`, `plot_correlation_matrix.m`,
`plot_KS_quantile_statistics.m`, `plot_ridgeplot.m`, and
`plot_appearanceTime_vs_growthRate.m`, open a dialog box asking which sample groups
and/or phenotypic features to include. `plot_ridgeplot.m` and
`plot_appearanceTime_vs_growthRate.m` additionally ask for the display order of the
groups and for which group pairs to test statistically. These selections apply only to
that function call and are not saved; run the script again to choose a different
combination.

**Logging and outputs.** Every script mirrors everything it prints to the Command
Window into a timestamped log file (`<script_name>_YYYYMMDD_HHMMSS.txt`), saved
automatically to the output directory defined in `ht.params.out_dir` during
preprocessing. Figures are saved there as `.fig`, `.pdf`, and `.png` files; some scripts
additionally export `.csv` or `.xlsx` tables.

**Suggested order.** There is no strict order for working through this chapter, but a
typical workflow is to: (1) start with
[Comparing counts and visualizing individual colony size distributions](qc-full-dataset.md)
to confirm segmentation quality and colony counts, and
[Manual vs automatic colony counts over time](timecourse-manual-auto-counts.md) to
validate them against manual plate counts; (2) move on to the distribution- and
group-comparison sections, [Comparing appearance time with other phenotypic
features](appearance-time-vs-morphology.md), [Correlation matrices and non-parametric group
comparison](correlation-matrix.md), [Kolmogorov – Smirnov (KS) and Quantile
Statistics](ks-quantile-statistics.md), and [Ridgeline plots with distribution
intervals](ridgeplot.md), to characterize and statistically compare
phenotypic heterogeneity; and (3) finish with the growth-kinetics sections, [Doubling
time](doubling-time.md), [Maximum growth rate](appearance-time-vs-growthrate.md), [Growth
curves](growth-curves.md), and [Colony appearance counts and early vs late appearing
colonies](cell-count-per-time.md), and the export and visualization sections, [Export data to
directory from MATLAB](export-data-to-output-directory.md), [Visualize individual
colonies](individual-colony.md), and [Visualize individual plates](individual-plate.md),
once you have identified the groups or colonies of interest.

:::{toctree}
:maxdepth: 1
:caption: Analyzing Data

qc-full-dataset
appearance-time-vs-morphology
correlation-matrix
cdf-with-statistics
doubling-time
growth-data-across-sampletypes
export-data-to-output-directory
export-colony-images
export-colony-images-intpersize-threshold
growth-curves
qc-growth-data
ks-quantile-statistics
combined-samples
morphology-colonies-gini
combined-samples-biological-replicates
appearance-time-colonies-gini
appearance-time-vs-growthrate
individual-colony
individual-plate
ridgeplot
timecourse-manual-auto-counts
cell-count-per-time
:::
