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
   - [Bulk Silicon](#bulk-silicon)
   - [Bulk Silicon with FMM particle–particle interactions](#bulk-silicon-with-fmm-particleparticle-interactions)
   - [Single-layer MoS₂](#single-layer-mos)
   - [2D MOSFET](#2d-mosfet)
   - [2D Resistor](#2d-resistor)
   - [Hot Carrier Cooling in Metal-Halide Perovskites (HPB)](#hot-carrier-cooling-in-metal-halide-perovskites-hpb)
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
- **Hot Phonon Bottleneck (HPB)** — a dynamic phonon bath (`emcPhononBath`) tracks the out-of-equilibrium LO phonon occupation N_q per wave-vector bin. Fröhlich scatter rates are updated every time step to reflect the evolving phonon temperature. This models the hot-phonon effect in metal-halide perovskites (MHP) and other polar semiconductors.

---

## Library architecture

ViennaEMC is **header-only**: all classes live under `include/`. No library file needs to be compiled — you only link against ScalFMM (for FMM examples) and OpenMP.

```
include/
├── emcDevice.hpp              — simulation domain (geometry, doping, temperature)
├── emcMaterial.hpp            — material constants (dielectric, mass density, …)
├── emcSimulation.hpp          — top-level device simulation loop
├── emcPhononBath.hpp          — dynamic LO phonon bath (for Hot Phonon Bottleneck)
├── emcScatterHandler.hpp      — manages scatter table construction and selection
├── emcConstants.hpp           — physical constants (SI units)
├── ParticleType/
│   ├── emcParticleType.hpp    — base class; holds valleys + scatter mechanisms
│   └── emcElectron.hpp        — concrete electron particle type
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
    ├── emcZeroOrderInterValleyScatterMechanism.hpp
    ├── emcFirstOrderInterValleyScatterMechanism.hpp
    ├── emcCoulombScatterMechanism.hpp
    └── …SingleLayer variants
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

### Bulk Silicon

**Path:** `examples/bulkSimulation/`
**Executable:** `build/examples/bulkSimulation/bulkSimulation`

Simulates electron transport in bulk silicon under a constant applied electric field. The conduction band is modelled with six equivalent X-valleys (non-parabolic, anisotropic). Scatter mechanisms: acoustic phonon deformation potential, zero-order intervalley phonon scattering.

**What you can study:**
- Steady-state drift velocity vs. field strength (velocity–field characteristic)
- Average carrier energy vs. field
- Valley occupancy as electrons scatter between the six X-valleys

**Key parameters to tune in `bulkSimulation.cpp`:**

| Parameter | Default | Description |
|---|---|---|
| `appliedFieldStrength` | 1 kV/cm | Magnitude of the background electric field |
| `appliedFieldDirection` | `{1,0,0}` | Direction (automatically normalised) |
| `totalTime` | 3 ps | Total simulation time |
| `dt` | 1 fs | Time step |

**Output files:** `avgEnergy<valley>.txt`, `avgVelocity<valley>.txt`, `valleyOccupation.txt`

---

### Bulk Silicon with FMM particle–particle interactions

**Path:** `examples/bulkSimulationFMM/`
**Executable:** `build/examples/bulkSimulationFMM/bulkSimulationFMM`

Same as the bulk silicon example but adds real-space Coulomb forces between electrons and between electrons and fixed ionised donors, computed with the Fast Multipole Method (ScalFMM). Useful for studying screening and carrier–carrier scattering effects.

**Additional requirement:** ScalFMM must be built (handled by `make buildDependencies`).

---

### Single-layer MoS₂

**Path:** `examples/singleLayerMoS2/`
**Executable:** `build/examples/singleLayerMoS2/singleLayerMoS2`

Simulates 2D electron transport in a monolayer of MoS₂. The simulation is 2D: valleys and scatter mechanisms use the `SingleLayer` variants in `include/`. Scatter mechanisms: acoustic and optical phonons (2D Fröhlich interactions available).

**What you can study:**
- 2D drift velocity vs. field
- Effects of different valley parameterisations (several literature models are available via a compile-time flag)

---

### 2D MOSFET

**Path:** `examples/mosfet2D/`
**Executable:** `build/examples/mosfet2D/mosfet2D`

Simulates an n-channel silicon MOSFET in 2D. The device geometry, doping profile, contacts, and gate potential are defined in the source file. The Poisson equation is solved self-consistently at each time step.

**What you can study:**
- I–V characteristics (drain current vs. drain–source voltage at different gate voltages)
- Carrier concentration and velocity profiles along the channel

Results in this example were validated against the classical EMC simulator CEMC from [ViennaWD](https://www.iue.tuwien.ac.at/software/viennawd/). See the example's `README.md` for a detailed description of the device structure and parameter definitions.

---

### 2D Resistor

**Path:** `examples/resistor2D/`
**Executable:** `build/examples/resistor2D/resistor2D`

Simulates a uniformly doped silicon resistor in 2D. A simpler device than the MOSFET — useful as a starting point for understanding the device simulation workflow before tackling the MOSFET.

---

### Hot Carrier Cooling in Metal-Halide Perovskites (HPB)

**Path:** `examples/hotCarrierMHP/`
**Executable:** `build/examples/hotCarrierMHP/hotCarrierMHP`

Simulates the cooling of photo-excited hot carriers in a bulk metal-halide perovskite (MHP) semiconductor, with and without the **Hot Phonon Bottleneck (HPB)** effect. This example implements the single-mode Fröhlich + phonon bath model of Lugli (1987), as applied to MHP in Faber, Filipovic, Koster, *J. Phys. Chem. Lett.* **2024**, 15, 12601.

#### Physics

When carriers relax they emit LO phonons (Fröhlich scattering). In the HPB regime the emission rate is so fast that the LO phonon population N_q rises well above its thermal equilibrium value N_0 (Bose–Einstein). The hot phonon bath in turn stimulates re-absorption, slowing further carrier cooling. ViennaEMC models this with:

1. **`emcPhononBath`** — tracks N_q in configurable wave-vector bins. After each time step it updates N_q according to:

   ```
   Nq[i] += (nEmitted[i] - nAbsorbed[i]) / Dph  -  (dt / tauLO) * (Nq[i] - N0)
   ```

   where `Dph = q² dq Vsim / (2π²)` is the number of phonon modes in bin `i`, and `tauLO` is the Klemens phonon lifetime.

2. **`emcHotPhononFroehlichAbsorption3D` / `emcHotPhononFroehlichEmission3D`** — Fröhlich scatter mechanisms that read N_q from the shared phonon bath instead of using the fixed N_0. After each scatter event the transferred wave vector |q| = |k_new − k_old| is recorded in the bath.

3. **`reinitScatterTables()`** — called after each phonon bath update to rebuild the self-scattering table with the new N_q-dependent rates.

#### Running the example

The example contains a compile-time switch:

```cpp
static constexpr bool USE_HOT_PHONON = true;  // false = equilibrium phonons
```

To compare HPB vs. equilibrium, compile and run the executable twice, changing this flag:

```bash
cd build/examples/hotCarrierMHP

# Run 1: HPB enabled (default)
./hotCarrierMHP
# Produces: avgEnergyHPB.txt, phononOccupationHPB.txt

# Edit hotCarrierMHP.cpp: set USE_HOT_PHONON = false, then recompile
cd /path/to/ViennaEMC
# edit examples/hotCarrierMHP/hotCarrierMHP.cpp
cd build && make hotCarrierMHP

# Run 2: equilibrium phonons
cd examples/hotCarrierMHP
./hotCarrierMHP
# Produces: avgEnergy.txt, phononOccupation.txt
```

Each run takes approximately 5–15 minutes on a single core for 10 ps of simulation time with 1000 carriers.

#### Output files

| File | Columns | Description |
|---|---|---|
| `avgEnergyHPB.txt` | `t [s]   <E> [eV]` | Average carrier energy vs. time, HPB run |
| `avgEnergy.txt` | `t [s]   <E> [eV]` | Average carrier energy vs. time, equilibrium run |
| `phononOccupationHPB.txt` | `t [s]   N_q` | Mean phonon occupation vs. time, HPB run |
| `phononOccupation.txt` | `t [s]   N_q` | Mean phonon occupation vs. time, equilibrium run (constant N_0) |

#### Expected results

With the default MAPbI₃-like parameters (m* = 0.2 mₑ, ε∞ = 10, εₛ = 30, ℏω₀ = 33 meV, τ_LO = 2 ps, n = 10²⁴ m⁻³, T_init ≈ 3000 K):

- N_q rises from N_0 = 0.387 to a peak of ≈ 0.463 (+19 %) at t ≈ 1.2 ps, then decays on the τ_LO = 2 ps timescale.
- During 0–2 ps, HPB and equilibrium carriers cool at similar rates (increased N_q raises both emission and absorption).
- After t ≈ 3 ps, HPB carriers plateau at a **higher energy** than the equilibrium run: hot phonons near threshold (E ≈ ℏω₀ = 33 meV) continually re-excite carriers, setting a higher effective carrier temperature floor.

#### Key parameters in `hotCarrierMHP.cpp`

| Parameter | Default | Description |
|---|---|---|
| `USE_HOT_PHONON` | `true` | Enable/disable Hot Phonon Bottleneck |
| `relEffMass` | 0.2 | Relative effective mass m*/mₑ |
| `eps_hi` | 10 | Optical (high-frequency) relative permittivity ε∞ |
| `eps_lo` | 30 | Static (low-frequency) relative permittivity εₛ |
| `phononEnergy` | 0.033 eV | LO phonon energy ℏω₀ |
| `tauLO` | 2×10⁻¹² s | LO phonon Klemens lifetime |
| `latticeTempK` | 300 K | Lattice temperature |
| `initCarrierTempK` | 3000 K | Initial hot carrier temperature |
| `boxSide` | 100 nm | Cubic simulation box side length |
| `carrierDensity` | 10²⁴ m⁻³ | Carrier density (≈ 10¹⁸ cm⁻³) |
| `dt` | 5×10⁻¹⁵ s | EMC time step (5 fs) |
| `totalTime` | 10×10⁻¹² s | Total simulation time (10 ps) |
| `nrPhononBins` | 1 | Number of phonon wave-vector bins (1 = single-mode) |
| `dq` | 2×10⁹ m⁻¹ | Bin width in wave-vector space |

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

Each example in `examples/` contains a Python plotting script. Run it from the directory where the output `.txt` files were written:

```bash
cd build/examples/bulkSimulation
python3 ../../../examples/bulkSimulation/plotBulkResults.py
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
