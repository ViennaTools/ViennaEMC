/**
 * Hot Carrier Cooling in Metal Halide Perovskites with Hot Phonon Bottleneck.
 *
 * Simulates hot carrier cooling in a perovskite-like bulk semiconductor
 * using the ensemble Monte Carlo method with the hot phonon bottleneck (HPB)
 * effect, following the methodology of:
 *
 *   Faber, Filipovic, Koster, J. Phys. Chem. Lett. 2024, 15, 12601.
 *
 * Material parameters (MAPbI₃-like, single-mode Fröhlich approximation):
 *   m*   = 0.2 m_e
 *   ε_∞  = 10  (optical dielectric constant)
 *   ε_s  = 30  (static dielectric constant)
 *   ħω₀  ≈ 33 meV  (LO phonon, ω₀ = 50 ps^{-1})
 *   τ_LO = 0.1 ps  (LO phonon Klemens lifetime)
 *
 * The simulation outputs:
 *   avgEnergyHPB.txt  –  average carrier energy [eV] vs time [s]  (with HPB)
 *   avgEnergy.txt     –  same without HPB (equilibrium phonon occupation)
 *
 * Run this example twice, toggling the USE_HOT_PHONON define, to compare.
 */

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>

// Reuse the basicBulkParticleHandler from the bulk simulation example.
// Copy or symlink that file, or adjust the path as needed.
#include "../bulkSimulation/basicBulkParticleHandler.hpp"

#include <ParticleType/emcElectron.hpp>
#include <ScatterMechanisms/emcFroehlichInteraction.hpp>
#include <ScatterMechanisms/emcHotPhononFroehlichMechanism.hpp>
#include <ValleyTypes/emcNonParabolicIsotropValley.hpp>
#include <emcDevice.hpp>
#include <emcPhononBath.hpp>

// Set to true to enable the hot phonon bottleneck; false for equilibrium.
static constexpr bool USE_HOT_PHONON = true;

const SizeType Dim = 3;
using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using ParticleHandler = basicBulkParticleHandler<NumType, DeviceType>;
using MapIdxTypeToPartType = ParticleHandler::MapIdxToParticleTypes;

// ---- Material parameters (MAPbI₃-like) --------------------------------
const NumType eps_hi = 10.;    // optical dielectric constant ε_∞
const NumType eps_lo = 30.;    // static dielectric constant ε_s
const NumType relEffMass = 0.2; // relative effective mass (m*/m_e)
const NumType alpha_np = 0.;   // non-parabolicity factor [1/eV] (0 = parabolic)
const NumType phononEnergy = 0.033; // LO phonon energy [eV]
const NumType tauLO = 2e-12;   // LO phonon Klemens lifetime [s] (MHP: 1-10 ps)
const NumType latticeTempK = 300.; // lattice temperature [K]

// ---- Simulation space --------------------------------------------------
// Box of (100 nm)³ → carrier density 10^18 cm^{-3} gives ~1000 particles.
const NumType boxSide = 100e-9;  // [m]
const std::array<NumType, 3> maxPos = {boxSide, boxSide, boxSide};
const std::array<NumType, 3> spacing = {10e-9, 10e-9, 10e-9};
const NumType Vsim = boxSide * boxSide * boxSide; // [m^3]

// Carrier density 10^18 cm^{-3} = 10^24 m^{-3}
const NumType carrierDensity = 1e24; // [m^{-3}]

// ---- Simulation parameters ---------------------------------------------
// Hot carriers are initialised at elevated temperature (~3000 K) using the
// device temperature trick; the phonon bath stays at the lattice temperature.
const NumType initCarrierTempK = 3000.; // [K]  initial "hot" carrier temperature

// Time step: must resolve the fastest phonon scattering rate (~10 fs).
const NumType dt = 5e-15;       // [s]  5 fs per step
const NumType totalTime = 10e-12; // [s]  10 ps cooling window
const SizeType nrStepsBetweenOutput = 10;  // save every 0.05 ps

// ---- Phonon bath parameters --------------------------------------------
// Single-mode LO phonon approximation: one bin spanning 0..2e9 m^{-1}.
// All Fröhlich scatter events (regardless of |q|) are pooled into this
// one mode, matching the single-mode HPB model of Lugli (1987) and
// Faber et al. JPCL 2024. With ~100 000 modes in the 100 nm³ box this
// gives smooth, low-noise N_q evolution (~0.003 change per timestep).
const SizeType nrPhononBins = 1;
const NumType dq = 2e9;               // bin width [m^{-1}] covers full Fröhlich q-range

// -----------------------------------------------------------------------

int main() {
  std::cout << "=== Hot Carrier Cooling in MHP"
            << (USE_HOT_PHONON ? " (HPB enabled)" : " (equilibrium phonons)")
            << " ===\n\n";

  // MHP-like material: use eps_lo as the dielectric constant of the device
  // (used only for Debye length; irrelevant for bulk / periodic BC).
  MaterialType material(eps_lo,    // relative dielectric constant
                        4000.,     // mass density [kg/m^3] (perovskite ~4 g/cm^3)
                        1e10,      // intrinsic carrier density [m^{-3}] (wide gap)
                        2000.,     // sound velocity [m/s]
                        1.6);      // band gap [eV] (MAPbI₃ ~ 1.6 eV)

  // Device temperature is used to initialize carrier momenta via Maxwellian.
  // We use an elevated temperature to mimic photoexcited hot carriers.
  DeviceType device{material, maxPos, spacing, initCarrierTempK};
  device.addConstantDopingRegion({0, 0, 0}, maxPos, carrierDensity);

  // ---- Particle type (conduction band electrons) ----------------------
  MapIdxTypeToPartType particleTypes;
  // 1000 energy levels up to 4 eV covers the hot-carrier regime
  particleTypes[0] = std::make_unique<emcElectron<NumType, DeviceType>>(
      1000, 4., /*usePotentialForInit=*/false);

  // Single non-parabolic isotropic valley (Gamma point of MHP)
  particleTypes[0]->addValley(
      std::make_unique<emcNonParabolicIsotropValley<NumType>>(
          relEffMass, constants::me, /*degFactor=*/1, alpha_np));

  // ---- Phonon bath (shared between abs + em scatter mechanisms) -------
  auto phononBath = std::make_shared<emcPhononBath<NumType>>(
      nrPhononBins, dq, tauLO, phononEnergy, latticeTempK, Vsim);

  std::cout << "Phonon bath: N_0 = " << phononBath->getN0()
            << "  (equilibrium at " << latticeTempK << " K)\n";

  // ---- Scatter mechanisms: Fröhlich absorption + emission -------------
  if (USE_HOT_PHONON) {
    // HPB version: N_q read from shared phonon bath
    particleTypes[0]->addScatterMechanism(
        {0},
        std::make_unique<emcHotPhononFroehlichAbsorption3D<NumType>>(
            0, phononEnergy, relEffMass, eps_hi, eps_lo, phononBath, "MHP"));
    particleTypes[0]->addScatterMechanism(
        {0},
        std::make_unique<emcHotPhononFroehlichEmission3D<NumType>>(
            0, phononEnergy, relEffMass, eps_hi, eps_lo, phononBath, "MHP"));
  } else {
    // Equilibrium version: N_q = N_0 (Bose-Einstein at lattice temperature)
    particleTypes[0]->addScatterMechanism(
        {0},
        std::make_unique<emcFroehlichAbsorption3D<NumType>>(
            0, phononEnergy, relEffMass, eps_hi, eps_lo, latticeTempK, "MHP"));
    particleTypes[0]->addScatterMechanism(
        {0},
        std::make_unique<emcFroehlichEmission3D<NumType>>(
            0, phononEnergy, relEffMass, eps_hi, eps_lo, latticeTempK, "MHP"));
  }

  // ---- Particle handler (bulk, periodic BC, no applied field) ---------
  std::array<NumType, 3> noField = {0., 0., 0.};
  ParticleHandler handler(device, particleTypes, noField, 0.);

  // Generate initial hot carrier distribution
  std::cout << "Creating hot carrier ensemble...\n";
  handler.generateInitialParticles();
  handler.printNrParticles();

  // ---- Main simulation loop -------------------------------------------
  const SizeType nrSteps = static_cast<SizeType>(std::ceil(totalTime / dt));
  std::cout << "Simulation: " << nrSteps << " steps × " << dt * 1e15
            << " fs = " << totalTime * 1e12 << " ps\n\n";

  std::string suffix = USE_HOT_PHONON ? "HPB" : "";
  std::ofstream energyFile("avgEnergy" + suffix + ".txt");
  std::ofstream phononFile("phononOccupation" + suffix + ".txt");

  // Record initial state
  auto avgE0 = handler.getAvgEnergy(0);
  energyFile << 0.0 << " " << avgE0[0] << "\n";
  phononFile << 0.0 << " " << phononBath->getMeanNq() << "\n";

  auto tStart = std::chrono::high_resolution_clock::now();

  for (SizeType step = 1; step <= nrSteps; step++) {
    // 1. Move all particles (drift + scatter), recording phonon events in bath
    handler.moveParticles(dt);

    // 2. Evolve phonon bath (only needed for HPB run)
    if (USE_HOT_PHONON) {
      phononBath->update(dt);
      // 3. Rebuild scatter tables with updated N_q
      particleTypes[0]->reinitScatterTables();
    }

    // Write output at regular intervals
    if (step % nrStepsBetweenOutput == 0) {
      auto avgE = handler.getAvgEnergy(0);
      NumType t = step * dt;
      energyFile << t << " " << avgE[0] << "\n";
      phononFile << t << " " << phononBath->getMeanNq() << "\n";

      std::cout << "  t = " << t * 1e12 << " ps | <E> = " << avgE[0]
                << " eV | N_q = " << phononBath->getMeanNq() << "\n";
    }
  }

  auto tEnd = std::chrono::high_resolution_clock::now();
  std::cout << "\nCPU time: "
            << std::chrono::duration_cast<std::chrono::seconds>(tEnd - tStart)
                   .count()
            << " s\n";

  energyFile.close();
  phononFile.close();

  std::cout << "\nOutput written to avgEnergy" << suffix << ".txt and"
            << " phononOccupation" << suffix << ".txt\n";
  return 0;
}
