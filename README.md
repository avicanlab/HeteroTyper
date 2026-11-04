# HeteroTyper

Automated, high-throughput system for imaging bacterial colonies and quantifying phenotypic heterogeneity, including appearance time, growth rate, doubling time, final colony size, and shape descriptors, at the single-colony level.

HeteroTyper consists of two components:

- **Imaging platform**: a repurposed laser cutter gantry fitted with a GoPro HERO 10 Black camera and an Arduino Nano controller, capable of capturing time-lapse images of up to 104 agar plates per run. The hardware architecture, lighting design, camera selection, and Arduino Nano controller circuit are documented in [Imaging Platform Information](https://heterotyper-imaging.readthedocs.io/).
- **Computational pipeline**: a MATLAB pipeline that segments individual colonies from the time-lapse images and extracts growth and morphology metrics for each colony.

## Universal plate-type pipeline

The computational pipeline runs a single, universal segmentation and analysis workflow across plate types, with no separate configuration required per agar type. It has been tested on:

- **Bright-background media**: Luria-Bertani agar (LA)
- **Dark-background media**: blood agar (BA), MacConkey agar (MC), *Yersinia* Selective Agar (YSA)

No programming background is required to run the pipeline, though basic familiarity with MATLAB is helpful. The only requirement is that colonies of the species of interest have a roughly circular appearance, as segmentation and morphology analysis are optimized for this geometry.

## Repository contents

- **`heterotyper.m`** — Main script for initializing and running the pipeline
- **`/Code`** — Downstream analysis scripts (preprocessing, quality control (QC), and analysis/plotting functions)
- **`/arduino`** — Arduino Nano script (`Ardunio_setupfile_v1.ino`) that coordinates image capture on the imaging platform

## Requirements

- MATLAB R2023b or later, on Windows or macOS
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

## Getting started

1. Install MATLAB (R2023b or later) with the toolboxes listed above.
2. Clone this repository and open the `HeteroTyper Pipeline` folder in MATLAB.
3. Initialize the pipeline from the Command Window:

   ```matlab
   NAME_OF_YOUR_DATA = heterotyper();
   ```

   When prompted, enter the full path to your image folder, metadata file, and output folder. Both Windows paths (e.g. `D:\MyExperiment\Images`) and macOS paths (e.g. `/Users/<name>/MyExperiment/Images`) are supported. On macOS, you can drag a folder from Finder into the Command Window to paste its full path.

4. Run preprocessing on your dataset:

   ```matlab
   preprocess_pipeline_data(NAME_OF_YOUR_DATA);
   ```

   Preprocessing saves its results to the MATLAB workspace as a struct named `ht`, which the downstream analysis scripts take as input.

> [!WARNING]
> We recommend preprocessing only one dataset per MATLAB session. `preprocess_pipeline_data` always saves its result as `ht` and replaces any existing `ht` without warning. To analyse a second dataset in the same session, first run `plot_appearanceTime_colonies_Gini` and `plot_morphology_colonies_Gini` on the first dataset, then rename `ht` (e.g. `ht_exp1 = ht; clear ht`) before preprocessing the next one. See the [Pipeline User Guide](https://heterotyper.readthedocs.io/) for details.

Full setup and usage instructions, including all preprocessing prompts and every downstream analysis script, are documented in the [Pipeline User Guide](https://heterotyper.readthedocs.io/).

## What's new in v1.1.0

- "Lag time" is renamed to **appearance time** throughout the pipeline (scripts, outputs, and documentation). Data and `ht` structs saved with v1.0.0 must be regenerated.
- New analyses: ridgeline plots with tests of spread (`plot_ridgeplot`), validation of automatic against manual colony counts (`plot_timecourse_manual_auto_counts`), and colony appearance counts with early vs late appearing colony comparisons (`extract_cell_count_per_time`).
- Minimum MATLAB version raised to R2023b.

## Applications

Potential applications include basic research into bacterial population heterogeneity and antibiotic sensitivity testing to optimize treatment of bacterial infections.
