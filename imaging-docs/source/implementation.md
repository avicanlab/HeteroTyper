# Implementation

## Machine Structure and Space

HeteroTyper is built on a repurposed laser cutter frame, which provides the 2D gantry
arm used to move the camera over the plate grid. The machine has a rectangular shape
and stands on four legs connected by an H-shaped support on either side. Beneath the
existing plate grid there is a sloping metal pane. This pane had no existing circuitry
connected to it, meaning its removal does not affect the core functionality of the
system.

## Current Lighting

The imaging platform is illuminated by LED strips mounted along the inner walls of the
enclosure, above the plate grid. This wall-mounted arrangement provides approximately
uniform average light intensity across the grid but produces uneven illumination on
plates positioned at the edges, which receive different intensities than centrally
located plates.

## Camera and Mounting

The GoPro HERO 10 Black is selected as the imaging camera. It provides 23 MP
resolution in photo mode, a 1/2.3" CMOS sensor with a crop factor of 5.6, and a fixed
16 mm prime lens (35 mm equivalent) with f/2.8 aperture. Its compact form factor is
well suited to the existing camera mount on the machine arm.

The GoPro is fitted onto the existing camera mount on the machine's moving arm. Power
cables connect to the GoPro, and a signal cable runs from the laser output of the
machine to the Arduino controller, which triggers image capture at the correct time.

The software controlling image capture is designed to require minimal user
interaction. The user need only ensure a memory card is present while the machine is
running. Image sequences are saved automatically, with each sequence in its own
folder. For example, if sequences are captured every hour for 24 hours, 24 folders are
created on the storage medium.

## Arduino Nano Controller and Camera Circuit

An Arduino Nano microcontroller serves as the central controller for coordinating
image capture. It interfaces with the GoPro HERO 10 Black and the laser cutter's
existing signal output to automate the imaging process.

### Circuit Overview

The Arduino Nano controller circuit consists of the following functional blocks:

- **Arduino Nano**: the main microcontroller, running at 5V via USB or a regulated 5V
  supply. It reads the trigger signal from the laser cutter and sends a shutter signal
  to the GoPro.
- **Logic Level Shifter (3.3V ↔ 5V)**: the GoPro HERO 10 Black operates its I/O at
  3.3V. A bidirectional logic level shifter (e.g., BSS138-based) is placed between the
  Arduino Nano (5V logic) and the GoPro trigger input to prevent damage to the camera.
- **Laser Trigger Input**: the laser cutter outputs a signal when a scan cycle begins.
  This signal is routed through a resistor voltage divider (if necessary to bring it
  within 0–5V range) to a digital input pin on the Arduino Nano (e.g., D2, configured
  as an interrupt pin).
- **GoPro Shutter Control**: a digital output pin on the Arduino (e.g., D3) is
  connected through the logic level shifter to the GoPro's wired remote/shutter input.
  A brief HIGH pulse triggers the camera to capture an image.
- **Timing Control**: the Arduino firmware monitors the laser trigger input and, after
  a configurable delay (to allow the gantry to settle above a plate), fires the
  shutter pulse. The delay and inter-sequence timing are configurable in firmware.
- **Power Supply**: the Arduino Nano is powered from a 5V regulated supply shared with
  the rest of the low-power electronics. The GoPro is powered from its own dedicated
  supply or battery.

### Pin Assignments

| Arduino Pin  | Connected To                | Function                                                                     |
|--------------|------------------------------|-------------------------------------------------------------------------------|
| D2 (INT0)    | Laser cutter signal output   | Interrupt-driven trigger: fires on rising edge when laser scan cycle starts   |
| D3           | Logic level shifter (LV side)| Shutter pulse output to GoPro (3.3V logic via shifter)                       |
| 5V / GND     | Level shifter HV/LV power    | Supplies HV (5V) and LV (3.3V) rails to the logic level shifter               |
| VIN / GND    | Regulated 5V supply          | Power input to the Arduino Nano                                              |

*Arduino Nano pin assignments for the HeteroTyper imaging system.*

### Circuit Operation

On power-up, the Arduino Nano initialises all pins and attaches an interrupt to D2
(rising edge). The laser cutter's signal line goes HIGH at the start of each scan
cycle, firing the interrupt service routine. The firmware then waits a configurable
settling delay (default: 200 ms) to allow the gantry arm to position and stabilise
above the plate, then asserts D3 HIGH for 100 ms to trigger the GoPro shutter. D3 is
then returned LOW. This sequence repeats for every scan position until the full plate
grid has been imaged.

Between full imaging cycles, the Arduino monitors an elapsed-time counter. When the
configured inter-cycle interval has passed, it issues a reset signal to the laser
cutter controller via the existing laser output line to begin the next scan. This
allows fully automated, unattended multi-hour imaging sessions.
