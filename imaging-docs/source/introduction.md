# Introduction

A bacterial colony is a population of bacteria originating from a single bacterium.
Despite this clonal origin, colonies can diverge substantially in growth dynamics and
stress tolerance, a phenomenon known as phenotypic heterogeneity. Characterizing this
heterogeneity at the level of individual colonies, rather than as a population average,
requires imaging and analysis tools capable of tracking many colonies in parallel over
time.

HeteroTyper addresses this need as an automated, high-throughput platform combining
time-lapse imaging with computational analysis of colony-level growth parameters,
including lag time and growth rate, across up to 104 agar plates simultaneously.
Potential applications include basic research into the origins of phenotypic
variability and antibiotic sensitivity testing to optimize treatment of bacterial
infections.

The imaging platform rests on four legs and has a rectangular shape (see below).
Plates are placed on a plane with a sloping surface beneath. The machine is a
repurposed laser cutter whose 2D gantry arm moves the camera over each plate position
in sequence.

```{figure} _static/machine-3d-model.png
:alt: 3D model of the HeteroTyper imaging platform rendered in Autodesk Fusion
:width: 70%
:align: center

3D model of HeteroTyper 1.0 rendered in Autodesk Fusion.
```

## Light Source

The imaging platform is illuminated by LED strips attached to the inner walls of the
machine, on the walls that extrude above the plate grid. Although all plates receive
approximately the same light intensity on average, the plates closest to the edges are
unevenly lit.

## Choice of Camera

One objective was to select an appropriate camera for detailed imaging of colonies on
agar plates. To evaluate cameras, the following parameters were considered: megapixel
(MP) value, sensor size, lens type, and focal length.

The MP value determines resolution and image detail. Higher MP values allow more
detailed imaging, which is important for heterogeneity analysis. Most common camera
sensors have 10–20 MP and professional cameras approximately 12–50 MP.

The sensor size is equally important: a larger sensor gathers more light and produces
higher quality images. Sensor size can be described by its dimensions (A × B mm),
area, and crop factor. The crop factor is the ratio of the diagonal of a 35 mm film
sensor to that of the camera sensor. A larger crop factor means a more cropped image
relative to a full-frame 35 mm image.

Focal length, expressed in mm, is the optical distance between the sensor and the
point of convergence of light in the lens. It determines angle of view, magnification,
and depth of field. A smaller focal length gives a wider angle of view with lower
magnification. Prime lenses have a fixed focal length and are typically lighter and
less expensive; zoom lenses offer variable focal lengths and greater versatility.

The selected camera for HeteroTyper is the GoPro HERO 10 Black. It provides 23 MP in
photo mode with a 1/2.3" CMOS sensor (4.55 × 6.17 mm). The lens is a prime lens with a
16 mm focal length (35 mm equivalent) and an aperture of f/2.8. The crop factor of the
1/2.3" CMOS sensor is 5.6.
