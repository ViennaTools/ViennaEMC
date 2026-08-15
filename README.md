# ViennaEMC

ViennaEMC is a header-only C++17 library for semiconductor carrier transport simulation using the **Multi-Valley Ensemble Monte Carlo (EMC)** method. It is developed at the [Institute for Microelectronics](http://www.iue.tuwien.ac.at/) at [TU Wien](https://www.tuwien.ac.at/).

---

## Table of Contents

1. [What it does](#what-it-does)
2. [Library architecture](#library-architecture)
3. [System requirements](#system-requirements)
4. [Building](#building)
   - [Step 1 — Clone and configure](#step-1--clone-and-configure)
   - [Step 2 — Build dependencies](#step-2--build-dependencies)
   - [Step 3 — Build the examples](#step-3--build-the-examples)
   - [Step 4 — Run an example](#step-4--run-an-example)
5. [Examples](#examples)
6. [Building and running tests](#building-and-running-tests)
7. [Installing as a library](#installing-as-a-library)
8. [Integrating into your own CMake project](#integrating-into-your-own-cmake-project)
9. [Plotting results](#plotting-results)
10. [Generating documentation](#generating-documentation)
11. [Troubleshooting](#troubleshooting)
12. [Authors and license](#authors-and-license)

---

## What it does

The Ensemble Monte Carlo method simulates the quantum-mechanical transport of charge carriers (electrons, holes) through a semiconductor by tracking a statistical ensemble of representative particles. Each particle drifts under the applied electric field and is scattered stochastically by phonons, impurities, and other mechanisms. The distribution of carrier momenta and energies evolves from these microscopic events, giving macroscopic observables such as drift velocity, energy, and valley occupancy.

ViennaEMC extends the classical EMC method with:

- **Multi-valley analytical band structures** — parabolic and non-parabolic, isotropic and anisotropic effective masses, with the Herring–Vogt transformation for anisotropic valleys.
- **Pluggable scatter mechanisms** — acoustic phonon, optical phonon (zero- and first-order intervalley), Fröhlich (polar-optical) in 3D and 2D, Coulomb scattering, grain boundary scattering, and more. New mechanisms can be added by subclassing `emcScatterMechanism`.
- **Self-scattering (null-scattering)** — the standard technique that lets every particle share the same time step regardless of its local scatter rate.
- **Particle–particle interactions via FMM** — long-range Coulomb forces between moving carriers and fixed donors, computed with [ScalFMM](https://gitlab.inria.fr/solverstack/ScalFMM).
- **Poisson solver coupling** — for full device simulations (MOSFETs, resistors) the electrostatic potential is solved self-consistently on a grid.
- **Hot-carrier dynamics in metal-halide perovskites** — a dynamic phonon bath (`emcPhononBath`) tracks the out-of-equilibrium LO phonon occupation N_q per wave-vector bin, coupled to an acoustic reservoir in a three-temperature carrier/LO/acoustic model. Fröhlich scatter rates are rebuilt every time step from the evolving occupation, capturing the hot-phonon bottleneck and the acoustic cascade. The Fröhlich vertex can be dynamically screened by the free-carrier gas (`emcPlasmonScreening` + `emcScreenedFroehlichInteraction`), with the rate and final-state angle drawn from the coupling-weighted occupation over the accessible wavevector window. On top of this the `hotCarrierMHP` example builds a full **bipolar device model** — electrons and holes sharing the bath, carrier–carrier and electron–hole scattering, radiative and Auger recombination, Pauli band filling, and energy-selective extraction from which open-circuit voltage and power conversion efficiency are computed. The same machinery drives the `hotPhononGa2O3` high-field transport example for wide-bandgap β-Ga₂O₃.
- **Monolayer MoS₂ transport and gas sensing** — a first-principles-parameterised model for 2D transition-metal dichalcogenides, with the full intrinsic and extrinsic scattering stack (screened polar-optical and piezoelectric phonons, remote surface-optical substrate phonons, charged-impurity and surface-roughness scattering), screened multivalley (K+Q) transport with per-valley Pauli exclusion, a self-consistent gated-FET device, and an ambient/adsorbate extension for chemical sensing (charge-transfer doping, adsorbate Coulomb and 1/f noise, humidity screening).

---

## Library architecture

ViennaEMC is **header-only**: all classes live under `include/`. No library file needs to be compiled — you only link against ScalFMM (for FMM examples) and OpenMP.

```
include/
├── emcDevice.hpp              — simulation domain (geometry, doping, temperature)
├── emcMaterial.hpp            — material constants (dielectric, mass density, …)
├── emcSimulation.hpp          — top-level device simulation loop
├── emcScatterHandler.hpp      — manages scatter table construction and selection
├── emcConstants.hpp           — physical constants (SI units)
│                                # hot-carrier modules (examples/hotCarrierMHP)
├── emcPhononBath.hpp          — three-temperature LO/acoustic phonon bath (HPB + acoustic cascade), wavevector-resolved N_q with optional per-bin lifetime profile
├── emcPlasmonScreening.hpp    — dynamic Debye screening state (q_s from live density and carrier temperature)
├── emcCarrierCarrierScatter.hpp   — intra-species carrier–carrier scattering
├── emcInterCarrierScatter.hpp     — electron–hole scattering (reduced-mass frame)
├── emcRecombination.hpp       — radiative (Bnp) + Auger (CCCH/CHHS) recombination
├── emcEnergySelectiveContact.hpp  — energy-selective carrier extraction
├── emcPauliExclusion.hpp      — Lugli–Ferry band filling on a k-space occupancy grid
├── emcHotCarrierOutput.hpp    — Fermi–Dirac fit → carrier temperature, V_OC, PCE
│                                # monolayer-MoS₂ modules (examples/singleLayerMoS2, gatedMoS2FET)
├── emcMultiValleyPauliExclusion.hpp — per-valley/subvalley band filling (K+Q)
├── Ambient/
│   ├── emcAdsorbateChargeTransfer.hpp — surface charge-transfer doping (donor/acceptor)
│   ├── emcAdsorbateNoise.hpp          — adsorption/desorption KMC (1/f + RTN)
│   └── emcHumidityDielectric.hpp      — humidity-dependent surface permittivity
├── ParticleType/
│   ├── emcParticleType.hpp    — base class; holds valleys + scatter mechanisms
│   ├── emcElectron.hpp        — concrete electron particle type
│   ├── emcHole.hpp            — concrete hole particle type (bipolar transport)
│   └── emcAcceptor.hpp        — acceptor particle type
├── ValleyTypes/
│   ├── emcParabolicIsotropValley.hpp
│   ├── emcNonParabolicIsotropValley.hpp
│   ├── emcParabolicAnisotropValley.hpp
│   ├── emcNonParabolicAnistropValley.hpp
│   └── …SingleLayer variants  — 2D transport (e.g. MoS₂)
└── ScatterMechanisms/
    ├── emcAcousticScatterMechanism.hpp
    ├── emcFroehlichInteraction.hpp          — equilibrium Fröhlich (3D)
    ├── emcHotPhononFroehlichMechanism.hpp   — HPB Fröhlich with dynamic N_q
    ├── emcScreenedFroehlichInteraction.hpp  — screened Fröhlich (equilibrium + HPB, q-resolved rate and angle)
    ├── emcZeroOrderInterValleyScatterMechanism.hpp
    ├── emcFirstOrderInterValleyScatterMechanism.hpp
    ├── emcCoulombScatterMechanism.hpp
    │                                        # 2D-material (MoS₂) mechanisms
    ├── emcFroehlichInteractionSingleLayer.hpp      — 2D polar-optical (Fröhlich)
    ├── emcScreenedIntravalleyOpticalMechanism.hpp  — screened intravalley optical
    ├── emcPiezoelectricSingleLayerScatterMechanism.hpp — piezoelectric acoustic
    ├── emcRemoteSurfaceOpticalPhononMechanism.hpp  — remote substrate SO phonons
    ├── emc2DChargedImpurityScatterMechanism.hpp    — charged-impurity sheet
    ├── emcSurfaceRoughnessScatterMechanism.hpp     — interface roughness
    └── emc2DScreening.hpp                          — free-carrier static screening
```

The `examples/` directory contains standalone programs that use the library. Each example has its own `CMakeLists.txt` and can be studied as a template for new simulations.

---

## System requirements

| Requirement | Minimum version | Notes |
|---|---|---|
| Linux | any recent | Tested on Ubuntu 22.04+ |
| C++ compiler | GCC 9 / Clang 10 | Must support C++17 and OpenMP |
| CMake | 3.12 | |
| Git | any | Required to fetch ScalFMM |
| Internet access | — | Only during `make buildDependencies` |
| Intel MKL | optional | Used by ScalFMM for BLAS/LAPACK; LAPACK/OpenBLAS can substitute |

ScalFMM is the only external dependency and is **downloaded and built automatically** by CMake the first time you run `make buildDependencies`.

---

## Building

### Step 1 — Clone and configure

```bash
git clone https://github.com/ViennaTools/ViennaEMC.git
cd ViennaEMC
mkdir build && cd build

# Configure (release build, examples enabled)
cmake .. -DCMAKE_BUILD_TYPE=Release -DVIENNAEMC_BUILD_EXAMPLES=ON
```

If you want to install the library to a custom prefix instead of `/usr/local`:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DVIENNAEMC_BUILD_EXAMPLES=ON \
         -DCMAKE_INSTALL_PREFIX=/path/to/install
```

### Step 2 — Build dependencies

This downloads ScalFMM from its Git repository and compiles it. It only needs to run once and may take several minutes depending on your machine and internet connection.

```bash
make buildDependencies
```

After it completes, re-run CMake so it can find the newly built ScalFMM:

```bash
cmake ..
```

> **Note:** If you already have ScalFMM installed, skip this step and pass `-Dscalfmm_DIR=/path/to/scalfmm/install` to the `cmake` command instead.

### Step 3 — Build the examples

```bash
make -j$(nproc)
```

This compiles all examples in parallel. Each executable is placed under `build/examples/<name>/`.

### Step 4 — Run an example

The executables write their output to the directory from which they are invoked, so change into the example's build folder first:

```bash
cd build/examples/bulkSimulation
./bulkSimulation
```

Output files (`.txt`) appear in the same directory and can then be plotted with the provided Python scripts (see [Plotting results](#plotting-results)).

---

## Examples

Each example is a standalone program with its own `CMakeLists.txt` and its
own README describing the model, the parameters, and the output files.
Executables build to `build/examples/<name>/` and write their output to
the directory they are run from.

| Example | Description |
|---|---|
| [Bulk Silicon](examples/bulkSimulation/) | Electron transport in bulk silicon under a constant field: six non-parabolic anisotropic X-valleys, acoustic and zero-order intervalley phonon scattering. Velocity–field curves, mean energy, valley occupation. The recommended starting point. |
| [Bulk Silicon with FMM](examples/bulkSimulationFMM/) | The bulk example plus real-space Coulomb forces between carriers and fixed donors, computed with the Fast Multipole Method (ScalFMM). |
| [Single-layer MoS₂](examples/singleLayerMoS2/) | 2D transport in monolayer MoS₂: deformation-potential, polar-optical, and piezoelectric phonons, remote substrate phonons, charged impurities, surface roughness, free-carrier screening, multivalley (K+Q) transport with per-valley Pauli exclusion. Kaasbjerg, Li, and Pilotto parameter sets. |
| [Ambient gas sensing in MoS₂](examples/singleLayerMoS2/) | The `ambientMoS2` driver in the same directory: adsorbate charge-transfer doping, adsorbate Coulomb scattering, humidity screening, and adsorption/desorption noise (1/f, RTN) for chemical sensing. |
| [Gated MoS₂ FET](examples/gatedMoS2FET/) | Self-consistent (EMC + Poisson) back-gated monolayer-MoS₂ transistor in 3D, with ohmic, Schottky, and carrier-reservoir contacts. Transfer characteristics and channel profiles. |
| [2D MOSFET](examples/mosfet2D/) | n-channel silicon MOSFET with self-consistent Poisson coupling; validated against the CEMC simulator from ViennaWD. |
| [2D Resistor](examples/resistor2D/) | Uniformly doped silicon resistor — the simplest self-consistent device, a stepping stone to the MOSFET. |
| [Hot-Carrier Dynamics in Metal-Halide Perovskites](examples/hotCarrierMHP/) | Photo-excited cooling and hot-carrier device operation in bulk MHPs: wavevector-resolved hot-phonon bottleneck with dynamic screening, acoustic cascade, bipolar transport, recombination, band filling, and energy-selective extraction (V_OC, PCE). Presets for MAPbI₃, CsSnI₃, FAPbI₃, CsPbI₃, CsPbBr₃. |
| [Hot-Phonon Transport in β-Ga₂O₃](examples/hotPhononGa2O3/) | High-field electron transport in bulk β-Ga₂O₃ with equilibrium vs. wavevector-resolved non-equilibrium phonon occupations: velocity–field curves, five-mode polar-phonon set, plasmon screening, impurity scattering. |

---

## Building and running tests

```bash
cd build
cmake .. -DVIENNAEMC_BUILD_TESTS=ON
make buildTests
ctest --output-on-failure
```

---

## Installing as a library

To install the headers and CMake config files to a system or custom path:

```bash
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install
make install
```

This copies:
- Headers to `<prefix>/include/ViennaEMC/`
- CMake package files to `<prefix>/lib/cmake/ViennaEMC/`

To uninstall:

```bash
make uninstall
```

---

## Integrating into your own CMake project

After installing, add the following to your project's `CMakeLists.txt`:

```cmake
find_package(ViennaEMC REQUIRED
             PATHS "/path/to/install/lib/cmake/ViennaEMC")

add_executable(myApp main.cpp)
target_include_directories(myApp PUBLIC ${VIENNAEMC_INCLUDE_DIRS})
target_link_libraries(myApp PRIVATE ${VIENNAEMC_LIBRARIES})
```

A minimal simulation program has this structure:

```cpp
#include <emcDevice.hpp>
#include <ParticleType/emcElectron.hpp>
#include <ValleyTypes/emcNonParabolicIsotropValley.hpp>
#include <ScatterMechanisms/emcFroehlichInteraction.hpp>
// ... include your particle handler

int main() {
    // 1. Define material and device geometry
    emcMaterial<double> material(eps_r, mass_density, ni, v_sound, bandgap);
    emcDevice<double, 3> device{material, maxPos, spacing, tempK};
    device.addConstantDopingRegion({0,0,0}, maxPos, dopingDensity);

    // 2. Define carrier type, valley(s), and scatter mechanism(s)
    std::map<int, std::unique_ptr<emcParticleType<...>>> particleTypes;
    particleTypes[0] = std::make_unique<emcElectron<double, DeviceType>>(nrEnergies, maxEnergy);
    particleTypes[0]->addValley(std::make_unique<emcNonParabolicIsotropValley<double>>(
        relEffMass, constants::me, degFactor, alpha));
    particleTypes[0]->addScatterMechanism({0},
        std::make_unique<emcFroehlichEmission3D<double>>(
            0, phononEnergy, relEffMass, eps_hi, eps_lo, tempK, "myMaterial"));

    // 3. Create particle handler, generate particles, run loop
    MyParticleHandler handler(device, particleTypes, field, dt);
    handler.generateInitialParticles();
    for (SizeType step = 0; step < nrSteps; step++)
        handler.moveParticles(dt);
}
```

The `examples/bulkSimulation/` example is the recommended starting point.

---

## Plotting results

A Python helper package `emcPlottingFiles` is provided in `helper/emcPlottingFiles/`. Install it once:

```bash
cd helper
pip install emcPlottingFiles
```

Required Python packages (installed automatically): `numpy`, `matplotlib`, `pandas`, `scipy`.

Most examples contain a Python plotting script (see each example's README). Run it from the directory where the output `.txt` files were written:

```bash
cd build/examples/bulkSimulation
python3 ../../../examples/bulkSimulation/plotBulkSimulationResults.py
```

---

## Generating documentation

The API documentation is generated with [Doxygen](https://www.doxygen.nl/):

```bash
cd docs/doxygen
./make_doxygen.sh
firefox html/index.html
```

---

## Troubleshooting

### ScalFMM build fails with missing Intel MKL include path

If CMake reports a path like `/opt/intel/oneapi/compiler/2025.3/include` does not exist, this is a mismatch between different oneAPI component versions. Two workarounds:

```bash
# Option A: create the missing directory
sudo mkdir -p /opt/intel/oneapi/compiler/2025.3/include

# Option B: remove the conflicting version from pkg-config search path
export PKG_CONFIG_PATH=$(echo $PKG_CONFIG_PATH | tr ':' '\n' \
    | grep -v '2025.3' | tr '\n' ':')
cmake ..
```

### `cannot find -lscalfmm` linker error

This means the CMake config is passing a plain string `scalfmm` instead of the imported target `scalfmm::scalfmm`. Check that `cmake/ViennaEMCConfig.cmake.in` uses:

```cmake
list(APPEND VIENNAEMC_LIBRARIES scalfmm::scalfmm)
```

### Hot phonon bath: N_q does not change from N_0

This is almost always a sign that the phonon mode density `Dph` has `Vsim` in the denominator instead of the numerator. The correct formula in `emcPhononBath.hpp` is:

```cpp
T Dph = q * q * dq * Vsim / (T(2) * constants::pi * constants::pi);
```

`Dph` should be a large number (order 10⁵ for a 100 nm box with dq = 2×10⁹ m⁻¹).

### Very noisy N_q oscillations

This occurs when there are too few phonon modes per bin (small `Dph`), so each scatter event changes N_q by a large fraction. Solutions:

- Use fewer, wider bins: reduce `nrPhononBins` and increase `dq`.
- For the single-mode approximation, use `nrPhononBins = 1` and `dq = 2e9` (covers the full Fröhlich q-range).

### Compilation errors: stray `@` or `*/` in scatter mechanism headers

C++ block comments are terminated by the first `*/`. Symbols like `m*/m_e` inside `/** … */` comments will prematurely close the comment block. Use `m_star/m_e` (ASCII only) in block comments. Similarly, avoid Unicode characters such as `ε` — use ASCII descriptions (`eps_inf`) instead.

---

## Authors and license

**Current maintainer:** Lado Filipovic  
**Former contributors:** Laura Gollner, Robin Steiner, Anna Benzer

This library was developed at the [Institute for Microelectronics](http://www.iue.tuwien.ac.at/), [TU Wien](https://www.tuwien.at/).

See file `LICENSE` in the base directory for the license terms.
