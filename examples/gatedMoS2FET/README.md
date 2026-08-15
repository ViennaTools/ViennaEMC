# Gated MoS₂ FET

Self-consistent (EMC + Poisson) simulation of a back-gated
monolayer-MoS₂ field-effect transistor in a 3D device domain: `x` is the
transport direction (source to drain), `y` the width, and `z` the gate
axis.

**Executable:** `build/examples/gatedMoS2FET/gatedMoS2FET`

## What it models

The MoS₂ monolayer sits in a single `z`-layer with `k_z = 0` confinement
(`electron2DDevice.hpp`); the gate with its HfO₂ oxide is built into the
gate contact at the bottom of the domain, and source and drain are n+
ohmic contacts at the two `x`-boundaries. Metal (Schottky) and
carrier-reservoir contacts are also available, for contacting an undoped
channel the way real 2D-material FETs are. The electrostatic potential is
solved self-consistently with an SOR Poisson solver and a cloud-in-cell
particle-mesh scheme.

The monolayer transport model (valleys and scattering stack) is shared
with the [singleLayerMoS2](../singleLayerMoS2/) example, using the
Kaasbjerg parameter set.

## What you can study

- Transfer characteristics: drain current vs. gate voltage.
- Channel carrier-density and potential profiles under gate control.

## Running the example

```bash
cd build/examples/gatedMoS2FET
./gatedMoS2FET
```

Device geometry, doping, contact configuration, and bias points are set
in `gatedMoS2FET.cpp`.
