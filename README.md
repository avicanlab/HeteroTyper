# HeteroTyper
This repository contains the computational pipeline in MATLAB within HeteroTyper to quantify bacterial  heterogeneity in given population from time-lapse agar plate images by extracting lag time, morphology-related parameters with robust segmentation on bright and dark media.

## Repository Organization

### 1) Computational Pipeline
The comptutational pipeline in ```MATLAB``` includes different optimizations for bright and dark agar plates. The implementations are available under:
- _HeteroTyper (Bright)/_ for bright plates,
- _HeteroTyper (Dark)/_ for dark plates


### 2) Bright & Dark Plate Images
Sample images to test the computational pipeline are also included in this repository. Images can be found under: 
- _Images/Bright/_ for bright plates (Sample Groups: 3h, 7h, 24h, 48h),
- _Images/Dark/_ for dark plates (Sample Groups: MLN, Spleen (SPL), Liver)


### 3) Sample Information
Information regarding all the samples is provided in ```Metadata.xlsx```. This file contains the experimental metadata associated with the samples images. It includes two sheets: 
- ```Meta_Bright```, which provides information for bright plates (columns: ```Position```, ```Time```, ```Replicate```, ```Dilution```, ```Count```),
- ```Meta_Dark```, which provides information for dark plates (columns: ```Position```, ```Mouse```, ```Organ```, ```Dilution```, ```Count```)

### 4) Implementation Manual
The full user manual — setup, initialization, preprocessing, and a reference page for
every analysis/export script — is published at [ReadTheDocs](https://heterotyper.readthedocs.io)
and built from the [`docs/`](docs/) folder in this repository. See
[Documentation](#documentation) below for how to build it locally.

## Documentation

The documentation lives in [`docs/`](docs/) as MyST Markdown and is built with
[Sphinx](https://www.sphinx-doc.org/) using the
[PyData Sphinx Theme](https://pydata-sphinx-theme.readthedocs.io/). It is published
automatically on [ReadTheDocs](https://readthedocs.org/) from the `main` branch.

To build it locally:

```bash
cd docs
pip install -r requirements.txt
make html        # On Windows: .\make.bat html
```

The rendered site is written to `docs/build/html/index.html`.

## Contacts
For queries on the implementation and data, please contact:
- kemal.avican@umu.se
- sena.gizem@umu.se


## Funding
This work was supported by Swedish Research Council (No. 2021-02466), Kempestiftelserna (JCK22-0017), and the Medical Faculty at Umeå University (FS 2.1.6-281-22) to K. Avican, by Swedish Research Council Excellence Center grant (No. 2022-06543) for the Center for Modeling Adaptive Mechanisms in Living Systems Under Stress to K. Avican. We acknowledge µNordic Single Cell Hub (µNiSCH) and Small Animal Research Imaging Facility (SARIF) at Umeå University. K. Kochanowski is supported by the Spanish Ministry of Research and Innovation (Grant number: RYC2021-033035-I /AEI/10.13039/501100011033).
