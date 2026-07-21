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
   - [Gated MoS₂ FET](#gated-mos-fet)
   - [Ambient gas sensing in monolayer MoS₂](#ambient-gas-sensing-in-monolayer-mos)
   - [2D MOSFET](#2d-mosfet)
   - [2D Resistor](#2d-resistor)
   - [Hot-Carrier Dynamics in Metal-Halide Perovskites](#hot-carrier-dynamics-in-metal-halide-perovskites)
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
- **Hot-carrier dynamics in metal-halide perovskites** — a dynamic phonon bath (`emcPhononBath`) tracks the out-of-equilibrium LO phonon occupation N_q per wave-vector bin, coupled to an acoustic reservoir in a three-temperature carrier/LO/acoustic model. Fröhlich scatter rates are rebuilt every time step from the evolving occupation, capturing the hot-phonon bottleneck and the acoustic cascade. On top of this the `hotCarrierMHP` example builds a full **bipolar device model** — electrons and holes sharing the bath, carrier–carrier and electron–hole scattering, radiative and Auger recombination, Pauli band filling, and energy-selective extraction from which open-circuit voltage and power conversion efficiency are computed.
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
├── emcPhononBath.hpp          — three-temperature LO/acoustic phonon bath (HPB + acoustic cascade)
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

Simulates 2D electron transport in a monolayer of MoS₂ (a transition-metal dichalcogenide). Valleys and scatter mechanisms use the `SingleLayer` and 2D variants in `include/`. The scattering stack covers deformation-potential acoustic and intervalley phonons, polar-optical (Fröhlich) and piezoelectric coupling, and free-carrier screening, plus, for supported devices, remote surface-optical substrate phonons and charged-impurity and surface-roughness scattering. Screened multivalley (K+Q) transport with per-valley Pauli exclusion is supported. Parameter sets from Kaasbjerg, Li, and Pilotto are provided (see the example's `README.md`).

**What you can study:**
- 2D drift velocity vs. field, per valley
- Phonon-limited and substrate-limited mobility and its temperature dependence
- The effect of different literature parameterisations and scatter mechanisms

---

### Gated MoS₂ FET

**Path:** `examples/gatedMoS2FET/`
**Executable:** `build/examples/gatedMoS2FET/gatedMoS2FET`

Self-consistent (EMC + Poisson) simulation of a back-gated monolayer-MoS₂ field-effect transistor in 3D: `x` is the transport direction (source to drain), `y` the width, and `z` the gate axis. The monolayer sits in one `z`-layer (`electron2DDevice`), the gate with its HfO₂ oxide is at the bottom, and the source and drain are n+ ohmic contacts. Metal (Schottky) and carrier-reservoir contacts are available for contacting an undoped channel the way real 2D-material FETs are.

**What you can study:**
- Transfer characteristics (drain current vs. gate voltage)
- Channel carrier-density and potential profiles under gate control

---

### Ambient gas sensing in monolayer MoS₂

**Path:** `examples/singleLayerMoS2/ambientMoS2.cpp`
**Executable:** `build/examples/singleLayerMoS2/ambientMoS2`

Runs the ambient/adsorbate extension for monolayer MoS₂ used for chemical sensing. Three demonstrations go through the same adsorbate model: an O₂ partial-pressure sweep giving the conductivity response σ(P)/σ₀ (charge-transfer depletion plus adsorbate Coulomb scattering), a relative-humidity sweep giving μ(RH) (dielectric screening), and an adsorption/desorption noise series N_ads(t) (1/f and random-telegraph noise). Mobility is obtained from a short low-field EMC drift-velocity run, and σ = n e μ.

**What you can study:**
- Gas-sensing transduction: conductivity and mobility response to O₂, humidity, and adsorbate coverage (donors such as NH₃ vs. acceptors such as NO₂)
- Low-frequency noise from dynamic adsorption and desorption

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

### Hot-Carrier Dynamics in Metal-Halide Perovskites

**Path:** `examples/hotCarrierMHP/`
**Executable:** `build/examples/hotCarrierMHP/hotCarrierMHP`

Simulates photo-excited hot-carrier cooling and device operation in a bulk metal-halide perovskite (MHP). It extends the single-carrier hot-phonon-bottleneck model of Faber, Filipovic, Koster, *J. Phys. Chem. Lett.* **2024**, 15, 12601 into a **bipolar, three-temperature device simulation**.

#### What it models

- **Hot-phonon bottleneck (HPB)** — carriers emit LO phonons (Fröhlich); the LO occupation N_q rises above its equilibrium value N_0 and stimulated re-absorption slows further cooling. `emcPhononBath` tracks N_q per wave-vector bin and updates it each step as

  ```
  Nq[i] += (nEmitted[i] - nAbsorbed[i]) / Dph  -  (dt / tauLO) * (Nq[i] - Nstar)
  ```

  where `Dph = q² dq Vsim / (2π²)` is the number of phonon modes in the bin and `tauLO` is the Klemens lifetime. The Fröhlich rates are rebuilt every step from the current N_q.
- **Acoustic cascade** — a three-temperature carrier/LO/acoustic model in which the LO decay target `Nstar` follows a finite-lifetime acoustic reservoir (lifetime `tau_ac`) that the LO modes back-heat, so the bottleneck persists longer. An optional **Ridley** channel (LO → LA + TO) competes with Klemens decay.
- **Bipolar transport** — electrons and holes share the phonon bath, with intra-species carrier–carrier and inter-species electron–hole scattering.
- **Recombination** — radiative (Bnp) and Auger (CCCH/CHHS); Auger returns the gap energy to a surviving carrier (Auger reheating).
- **Band filling** — Lugli–Ferry Pauli blocking on a k-space occupancy grid.
- **Energy-selective contact (ESC)** — extraction through a window [E_ex ± ΔE/2], from which the working open-circuit voltage and PCE are computed.

Optionally, `--use_multimode 1` resolves the two IR-active MAPbI₃ LO branches instead of a single effective mode.

#### Running the example

Every mechanism and material parameter is a command-line flag — there are no compile-time switches. Run from the build folder:

```bash
cd build/examples/hotCarrierMHP

# Full default model (MAPbI3, all mechanisms on)
./hotCarrierMHP

# Equilibrium-phonon reference (HPB off)
./hotCarrierMHP --use_hot_phonon 0

# CsSnI3 presets with a selective device contact
./hotCarrierMHP --material sn --E_ex 0.2 --delta_E 0.05

# Acoustic-lifetime study, reproducible, larger box for low shot noise
./hotCarrierMHP --tau_ac 10e-12 --box 200e-9 --seed 1
```

At the end of a run the console prints the working `V_OC`, extraction `yield`, and `PCE`. A default run (100 nm box, 10 ps) finishes in well under a minute; the larger boxes used for low-noise or band-filling runs take longer.

#### Key command-line flags

Defaults are for MAPbI₃; `--material sn` switches to CsSnI₃ presets. Times are in seconds, energies in eV, densities in m⁻³. Boolean flags take `0`/`1`.

| Flag | Default (pb / sn) | Description |
|---|---|---|
| `--material` | `pb` | Presets: `pb` = MAPbI₃, `sn` = CsSnI₃ |
| `--use_hot_phonon` | `1` | Hot-phonon bottleneck (dynamic N_q) |
| `--use_acoustic` | `1` | Acoustic reservoir (three-temperature model) |
| `--use_holes` | `1` | Bipolar transport (electrons + holes) |
| `--use_recomb` | `1` | Radiative + Auger recombination |
| `--use_cc` | `1` | Carrier–carrier scattering |
| `--use_esc` | `1` | Energy-selective extraction |
| `--use_bf` | `0` | Pauli band filling (selects the large box) |
| `--use_multimode` | `0` | Two-branch (multi-mode) Fröhlich coupling |
| `--ridley_w` | `0` | Ridley LO→LA+TO branching fraction (0–1) |
| `--E_photon` | `3.1` | Pump photon energy |
| `--E_ex` / `--delta_E` | `0.4` / `0.05` | ESC window centre / width |
| `--density` | `1e24` | Carrier density (10²⁴ m⁻³ ≈ 10¹⁸ cm⁻³) |
| `--tau_ac` | `30e-12` | Acoustic-reservoir lifetime |
| `--tau_lo` | `2e-12` / `0.1e-12` | LO (Klemens) phonon lifetime |
| `--box` | auto | Cubic box side [m]; `0` auto-selects (100 nm, or 464 nm with `--use_bf`) |
| `--total_time` | `10e-12` | Simulation time |
| `--seed` | `0` | RNG seed; `0` = wall-clock, nonzero = reproducible |

Further overrides: `--eps_hi`, `--eps_lo`, `--mass_e`, `--mass_h`, `--omega_lo`, `--omega_to`, `--tau_to`, `--tau_ex`, `--E_gap`, `--T_lat`.

#### Output files

The driver writes five time series. Each filename carries a suffix encoding the enabled mechanisms — `HPB` or `EQ`, then `_AC` (acoustic), `_MM` (multi-mode), `_RID` (Ridley), `_BF` (band filling), `_1C` (single-carrier), `_ESC` (extraction). The full default model, for example, writes `carrierTempHPB_AC_ESC.txt`.

| File (before suffix) | Columns | Description |
|---|---|---|
| `avgEnergyElectrons` / `avgEnergyHoles` | `t   <E>[eV]` | Mean carrier energy vs time |
| `phononOccupation` | `t   N_LO   [N_ac  T_ac]   [N_TO  T_TO]` | LO occupation, plus acoustic (and Ridley TO) reservoir when enabled |
| `nrCarriers` | `t   N_e  N_h  N_esc_e  N_esc_h` | Surviving and extracted carrier counts |
| `carrierTemp` | `t   T_MB_e  T_MB_h  T_FD_e  μ_e  T_FD_h  μ_h` | Maxwell–Boltzmann and Fermi–Dirac carrier temperatures / chemical potentials |

> **Reproducibility:** particle moves run over OpenMP threads, each with its own RNG stream, so the seed **and** the thread count together define a trajectory. Pin `OMP_NUM_THREADS` and pass a nonzero `--seed` for repeatable runs.

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
