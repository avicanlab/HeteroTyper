# HeteroTyper Pipeline User Manual

This is the user reference for the **HeteroTyper computational pipeline**, the analysis
software that complements the HeteroTyper imaging system. The HeteroTyper system
captures time-lapse images of bacterial colonies growing on agar plates and records how
each colony develops over time. The primary aim of the pipeline is to quantify **lag time**
at the single-colony level — the time each individual colony takes before it adjusts to
the new environment and initiates growth.

By measuring lag time across thousands of colonies simultaneously, the pipeline enables
analysis of how variable this behaviour is within a bacterial population. In addition to
lag time, the pipeline extracts a range of complementary growth and morphological
features for each colony, including final colony size, intensity per size, growth rate,
doubling time, and shape descriptors such as circularity.

This handbook is intended for researchers who wish to quantify bacterial colony growth
from time-lapse imaging data. While the pipeline was developed alongside the HeteroTyper
imaging system, it is not dependent on it, and any researcher with time-lapse images of
agar plates may use it. The only requirement is that the species of interest forms
colonies with a roughly circular appearance, as the segmentation and morphology analysis
are optimized for this colony geometry. No programming background is required, though
basic familiarity with MATLAB will help when following the instructions.

The handbook guides the user through installing and initializing the pipeline,
preprocessing raw imaging data, and running all analysis scripts to explore, compare, and
export results. Each section describes one script or pipeline stage: what it does, what
it requires to run, what inputs it accepts, what outputs it produces, and how to resolve
the most common errors. Definitions of all quantitative metrics are provided in the
sections where they are first introduced.

:::{toctree}
:maxdepth: 2
:caption: Contents

introduction
pipeline-setup
analyzing-data/index
:::
