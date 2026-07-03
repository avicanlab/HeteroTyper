# Analyzing Data

## How to use this chapter

This chapter documents every analysis and export script included in the HeteroTyper
pipeline. Each section corresponds to one script and follows the same format: Purpose,
Requirements, Inputs, Outputs, and common errors.

Before running any script in this chapter, run `preprocess_pipeline_data.m` on your
loaded dataset (see [Data preprocessing](../pipeline-setup.md#data-preprocessing)). This
step performs all shared calculations once and stores the results in the `ht` struct,
which almost every script in this chapter reads directly rather than recomputing.

**Scripts that take `ht`.** Most analysis and plotting scripts take `ht` as their only
input, for example `plot_combined_samples(ht)`. A few also require the raw data variable
used during preprocessing, because they need access to raw images or plate metadata not
stored in `ht`: `export_colony_images(data, ht)`,
`export_colony_img_intPerSize_threshold(data, ht)`,
`export_data_to_output_directory(data, ht)`, and the optional image montage in
`plot_QC_full_dataset(ht, data)`. Two further scripts,
[`plot_lagTime_colonies_Gini`](lagtime-colonies-gini.md) and
[`plot_morphology_colonies_Gini`](morphology-colonies-gini.md), take both `ht` and your
raw data variable to build plate-level replicate labels; on first use, each asks which
metadata columns correspond to the biological and technical replicate, then caches the
answer in `ht` so later calls skip the prompt.

**Scripts that take `data` directly.** A small number of exploratory and visualization
scripts, `plot_individual_colony.m`, `plot_individual_plate.m`, and
`plot_growth_data_across_sampletypes.m`, operate directly on the raw `data` variable
rather than `ht`, and prompt you interactively for a plate or colony index to inspect.
These do not require `preprocess_pipeline_data.m` to have been run first.

**Interactive selections.** Several scripts, including `plot_combined_samples.m`,
`plot_lagTime_vs_morphology.m`, `plot_correlation_matrix.m`, and
`plot_KS_quantile_statistics.m`, open a dialog box asking which sample groups and/or
phenotypic features to include. These selections apply only to that function call and
are not saved; run the script again to choose a different combination.

**Logging and outputs.** Every script mirrors everything it prints to the Command
Window into a timestamped log file (`<script_name>_YYYYMMDD_HHMMSS.txt`), saved
automatically to the output directory defined in `ht.params.out_dir` during
preprocessing. Figures are saved there as `.fig`, `.pdf`, and `.png` files; some scripts
additionally export `.csv` or `.xlsx` tables.

**Suggested order.** There is no strict order for working through this chapter, but a
typical workflow is to: (1) start with
[Comparing counts and visualizing individual colony size distributions](qc-full-dataset.md)
to confirm segmentation quality and colony counts; (2) move on to the distribution- and
group-comparison sections, [Comparing lag time with other phenotypic
features](lagtime-vs-morphology.md), [Correlation matrices and non-parametric group
comparison](correlation-matrix.md), and [Kolmogorov – Smirnov (KS) and Quantile
Statistics](ks-quantile-statistics.md), to characterize and statistically compare
phenotypic heterogeneity; and (3) finish with the growth-kinetics sections, [Doubling
time](doubling-time.md), [Maximum growth rate](lagtime-vs-growthrate.md), and [Growth
curves](growth-curves.md), and the export and visualization sections, [Export data to
directory from MATLAB](export-data-to-output-directory.md), [Visualize individual
colonies](individual-colony.md), and [Visualize individual plates](individual-plate.md),
once you have identified the groups or colonies of interest.

:::{toctree}
:maxdepth: 1
:caption: Analyzing Data

qc-full-dataset
lagtime-vs-morphology
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
lagtime-colonies-gini
lagtime-vs-growthrate
individual-colony
individual-plate
:::
