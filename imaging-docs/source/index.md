# HeteroTyper Imaging Platform

This is the system reference for the **HeteroTyper imaging platform**, the hardware
that complements the [HeteroTyper computational pipeline](https://heterotyper.readthedocs.io/).
HeteroTyper is an automated, high-throughput system for analysing phenotypic
heterogeneity among bacterial colonies growing on agar plates. The system is composed
of an imaging platform and a computational pipeline: the imaging platform records
time-lapse images of up to 104 agar plates, while the computational pipeline analyses
properties of individual bacterial colonies such as lag time and growth rate.

This document describes the core of the HeteroTyper imaging platform, including its
architecture, lighting design, camera selection, and the Arduino-based controller
circuit that coordinates image capture. It is intended for anyone operating,
maintaining, or extending the physical imaging hardware.

:::{toctree}
:maxdepth: 2
:caption: Contents

introduction
implementation
future-improvements
appendix
:::
