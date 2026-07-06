# HeteroTyper

Automated, high-throughput system for imaging bacterial colonies and quantifying phenotypic heterogeneity, including lag time, growth rate, doubling time, final colony size, and shape descriptors, at the single-colony level.

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

## Getting started

1. Install MATLAB (R2020a or later).
2. Clone this repository and open the `HeteroTyper Pipeline` folder in MATLAB.
3. Initialize the pipeline from the Command Window:

   ```matlab
   NAME_OF_YOUR_DATA = heterotyper();
   ```

4. Run preprocessing on your dataset:

   ```matlab
   preprocess_pipeline_data(NAME_OF_YOUR_DATA);
   ```

Full setup and usage instructions, including all preprocessing prompts and every downstream analysis script, are documented in the [Pipeline User Guide](https://heterotyper.readthedocs.io/).

## Applications

Potential applications include basic research into bacterial population heterogeneity and antibiotic sensitivity testing to optimize treatment of bacterial infections.