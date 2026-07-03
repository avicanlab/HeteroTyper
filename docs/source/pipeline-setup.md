# HeteroTyper Pipeline Setup

## Setup

To setup the HeteroTyper pipeline, MATLAB must be installed on your computer. Once
installed, please download the required pipeline from GitHub using
[https://github.com/avicanlab/HeteroTyper](https://github.com/avicanlab/HeteroTyper).

1. Run MATLAB on your computer.
2. Using the folder panel, navigate inside the HeteroTyper Pipeline folder.

:::{note}
Before running the pipeline, make sure the folder panel shows the HeteroTyper Pipeline
folder as the current directory, with the `Code` subfolder and the main
`heterotyper.m` script visible.
:::

<!-- TODO: insert screenshot of the MATLAB Current Folder panel showing the
     HeteroTyper Pipeline folder, matching the original manual's Setup screenshot. -->

## Initialize pipeline

Using the Command Window, run the following to initialize the pipeline.

```matlab
NAME_OF_YOUR_DATA = heterotyper();
```

HeteroTyper is an interactive pipeline, and following questions will be asked to user
before it's initialized:

0. Path to the target image folder
1. Number of plates in the experiment
2. Minimum colony numbers on plates
3. Maximum colony numbers on plates
4. File directory and name for metadata (Excel file: `*.xlsx`)
5. Sheet name used in metadata
6. Column name to use in plate center summary figure
7. Select a number to enable/disable figure saving (`0`: disable | `1`: enable and enter
   output directory path)

Once you're done with input parameters, a final outlook of parameters will be shown and
pipeline will start running.

## Terminate pipeline

To terminate the pipeline, use the following combinations on your keyboard:

- **Windows:** `Ctrl + C`
- **MacOS:** `Ctrl + C` or `Command + .`

## Data preprocessing

Navigate to the "Code" folder inside the pipeline using the Current Folder panel on your
MATLAB.

<!-- TODO: insert screenshot of the Code folder and preprocess_pipeline_data.m open in
     the MATLAB editor, matching the original manual's Data preprocessing screenshot. -->

To initialize data preprocessing script, run the following command in your MATLAB
Command Window.

```matlab
preprocess_pipeline_data(NAME_OF_YOUR_DATA);
```

The user will be asked to provide the following information prior to initializing the
preprocessing:

0. **RT incubation time** — total time the plates spent at room temperature before
   imaging begins (e.g. 20 or 20.5 h).
1. **Imaging interval** — first select the time unit (1 = minutes, 2 = hours), then
   enter the interval value (e.g. 30 min or 0.5 h).
2. **Maximum possible lag time** — auto-calculated from image count as
   `incTime + (n_images − 1) × interval`. Press ENTER to accept the calculated value, or
   type an override.
3. **Output directory** — full path to the folder where figures and log files are
   saved. Press ENTER for the current MATLAB directory. If the path does not exist, the
   script offers to create it.
4. **Metadata column for sample grouping** — shows a numbered list of all metadata
   columns with value previews; enter the column number that identifies the sample
   group (e.g. `Time`, `Condition`, `Strain`).
5. **Biological replicate column** — identifies the biological sample or strain used
   as the x-axis label in violin plots (e.g. `Set`, `BioRep`, `SampleName`). Enter `0`
   to skip and use plate position numbers only.
6. **Technical replicate column** — identifies replicate plates of the same
   biological sample; appended to the bio-rep label as `BioRep_TechRep` (e.g.
   `Replicate`, `TechRep`, `R1`/`R2`/`R3`). Enter `0` to skip.
7. **Manual count column for QC Plot 1A** — asked separately after preprocessing
   starts; select the metadata column containing manually counted colony numbers to
   enable the automated vs manual comparison plot. Enter `0` to skip.

The pipeline will generate a log file (`preprocess_log_YYYYMMDD_HHMMSS.txt`) that will
save everything printed on the Command Window onto the output directory defined by the
user in Step 3.
