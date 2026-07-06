# Appendix

## Summarised User Manual for Imaging System

The following is a step-by-step procedure for running an imaging session:

1. Connect power to the machine, GoPro camera, and Arduino controller.
2. Set the desired imaging cycle interval using the timing controls on the interface
   board.
3. Insert a memory card into the GoPro, preferably an empty one.
4. Power on the laser cutter using the key on the side panel.
5. Select the appropriate LightBurn file for the plate configuration to be used and
   start the run.
6. The machine will begin the scan cycle automatically. The Arduino controller will
   trigger the GoPro shutter at each plate position.
7. The machine restarts automatically after the configured inter-cycle interval. It is
   safe to retrieve the memory card between imaging cycles to preview images and check
   for alignment issues.
8. To stop, turn off the laser cutter using the key and disconnect power to the
   Arduino and camera.

## Arduino Nano Circuit Diagram

The diagram below shows a simplified overview of how all electronic modules are
attached. It only shows signals, not power connections, and it does not show
individual cables, voltage levels, or the communication protocols used. Its purpose is
to give an overview of the setup, while leaving the details to visual inspection on
site.

```{figure} _static/circuit-diagram.png
:alt: Simplified circuit diagram of the HeteroTyper imaging system
:width: 90%
:align: center

Simplified circuit diagram showing connections between modules. Signals to/from the
laser cutter are marked with arrows, hardware components in rectangular boxes, and the
camera and Arduino in a rounded box. The laser cutter gantry is shown as a dashed
container.
```

## List of Components

| Component                          | Purpose                                                        |
|-------------------------------------|------------------------------------------------------------------|
| GoPro HERO 10 Black                 | Imaging of agar plates (23 MP, prime 16 mm lens, f/2.8).         |
| Arduino Nano                        | Microcontroller (coordinates image capture timing).             |
| BSS138 Logic Level Shifter (3.3V / 5V) | Converts 5V Arduino signals to 3.3V GoPro input.              |
| Potentiometer                       | Sets the 30-minute imaging interval on the Arduino Nano.         |

*List of electronic and optical components used in HeteroTyper.*
