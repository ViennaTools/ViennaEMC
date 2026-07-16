// Self-consistent gated monolayer-MoS2 FET (EMC + Poisson), 3D device:
//   x = transport (source->drain), y = width, z = perpendicular (gate axis).
// The MoS2 monolayer sits in one z-layer (k_z = 0 confinement, electron2DDevice);
// the gate sits at ZMIN with the HfO2 oxide built into the gate contact; source &
// drain are n+ ohmic contacts at XMIN/XMAX. Uses the 3D-capable CIC PM scheme.
#include <chrono>
#include <iostream>
#include <memory>
#include <numeric>
#include <fstream>

#include <ParticleHandler/emcBasicParticleHandler.hpp>
#include <PMSchemes/emcCICScheme.hpp>
#include <PoissonSolver/emcSORSolver.hpp>
#include <emcDevice.hpp>
#include <emcSimulation.hpp>

#include "../singleLayerMoS2/parameterKaasbjerg.hpp"
#include "electron2DDevice.hpp"

const SizeType Dim = 3;
const std::string fileNamePrefix = "gatedMoS2FET";
using namespace std::chrono;
using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using PMScheme = emcCICScheme<NumType, DeviceType>;
using ParticleHandler = emcBasicParticleHandler<NumType, DeviceType, PMScheme>;
using PoissonSolver = emcSORSolver<NumType, DeviceType, ParticleHandler>;
using SimulationType = emcSimulation<NumType, DeviceType, PoissonSolver,
                                     ParticleHandler, PMScheme>;
using ValueVec = DeviceType::ValueVec;

const NumType temperature = 300;

NumType adaptPotential(const NumType &pot, const DeviceType &device) {
  return device.getMaterial().getBandGap() / 2 - pot * device.getThermalVoltage();
}

// MoS2 as an effective 3D medium for the electrostatics: eps_r, rho[kg/m3],
// Ni[1/m3], vSound, bandGap. (Ni moderate for numerical stability; the FET
// carriers come from doping + gate, not the intrinsic density.)
MaterialType getMoS2Material() { return MaterialType{4.0, 5060., 1e16, 6700., 1.8}; }
// sheet -> volume doping over the 0.65 nm monolayer: 1e13 cm^-2 -> 1.5e26 m^-3
const NumType tMono = 0.65e-9;

// HfO2 gate-dielectric substrate parameters for the B1-B3 channels under the gate
// (Ma & Jena, PRX 4, 011043 (2014), Tab. I): high-/low-freq permittivity + the two
// remote surface-optical (SO) phonon mode energies [eV].
const NumType hfo2EpsInf = 5.03, hfo2Eps0 = 23.0;
const NumType hfo2wSO1 = 0.0124, hfo2wSO2 = 0.0484;
// B2: fixed charge / trapped charge at the MoS2/HfO2 interface (representative,
// HfO2 is trap-rich); B3 screening uses the static MoS2/HfO2 average permittivity.
const NumType interfaceChargeDensity = 1e16; // N_it ~ 1e12 cm^-2
const NumType interfaceEpsAvg = 12.0;        // (eps_MoS2 + eps_HfO2,static)/2

// Real monolayer-MoS2 FET: a UNIFORM flake with only a low intrinsic
// (unintentional n) doping -- NO chemical S/D doping. The mobile channel charge
// is induced electrostatically by the GATE; carriers enter/leave through metal
// Ti/Au Schottky contacts. The screening + substrate scattering use the
// representative ON-state gate-induced sheet density (the operating point).
const NumType intrinsicDoping = 5e22;         // uniform unintentional n (~3e10 cm^-2)
const NumType schottkyBarrier = 0.05;         // Ti/MoS2 electron barrier [V] (good contact)
const NumType channelOnStateSheet = 9.75e16;  // ~1e13 cm^-2 (on-state screening)

const std::vector<NumType> Vgs = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0}; // transfer sweep (wide)
const std::vector<NumType> Vds = {1.0};

int main() {
#ifdef _OPENMP
  omp_set_num_threads(16);
  std::cout << ">> Parallel: " << omp_get_max_threads() << " threads\n";
#endif
  std::ofstream idvg("gatedMoS2FET_IdVg.txt");
  idvg << "# Vg[V]  Vd[V]  Id[A]\n";
  for (auto Vg : Vgs) {
    for (auto Vd : Vds) {
      // ---- geometry (monolayer channel); coarse in-plane cells so each thin
      //      monolayer cell holds ~10 carriers ----
      const ValueVec maxPos = {200e-9, 1500e-9, 0.65e-9}; // x,y,z (m)
      const ValueVec spacing = {5e-9, 10e-9, 0.65e-9}; // finer along transport (x)
      DeviceType device{getMoS2Material(), maxPos, spacing};

      // ---- UNIFORM intrinsic (unintentional n) doping over the whole flake ----
      //      (NO chemical S/D doping -- this is a 2D-material FET, not a MOSFET)
      device.addConstantDopingRegion({0, 0, 0}, {200e-9, 1500e-9, tMono}, intrinsicDoping);

      // ---- contacts: Ti/Au Schottky source/drain ON TOP (ZMAX) over the S/D
      //      footprints (real 2D-FET geometry: metal sits on top of the flake,
      //      not on its atomically-thin edge); HfO2 back-gate on the bottom
      //      (ZMIN) under the channel. Carriers inject from the metal over the
      //      Schottky barrier, flow laterally, and the back-gate accumulates the
      //      channel between the two top contacts.
      // top contacts at the ends; the gate is pulled IN to 60..140nm, leaving
      // 20nm ungated intrinsic SPACERS (40..60, 140..160) between each contact and
      // the gated channel so the contacts cannot refill the channel when the gate
      // depletes it -> a real off-state.
      device.addSchottkyContact(emcBoundaryPos::ZMAX, 0.0, {0, 0}, {40e-9, 1500e-9}, schottkyBarrier);    // source (top)
      device.addSchottkyContact(emcBoundaryPos::ZMAX, Vd, {160e-9, 0}, {200e-9, 1500e-9}, schottkyBarrier); // drain (top)
      device.addGateContact(emcBoundaryPos::ZMIN, Vg, {60e-9, 0}, {140e-9, 1500e-9},
                            20.0, 2e-9, device.getMaterial().getBandGap() / 2.); // HfO2 back-gate, 2nm oxide

      std::cout << "gated MoS2 FET  Vg=" << Vg << " V  Vd=" << Vd << " V\n";

      PoissonSolver solver(device, 1e-3, 1.8);
      PMScheme pmScheme;

      // (c) quantum capacitance: degenerate 2D Fermi-Dirac electron statistics in
      //     BOTH the Poisson closure and the particle generation (same beta).
      //     beta = Ni*tMono / N2D,  N2D = g*m*kT/(2 pi hbar^2) the 2D DOS. The
      //     density-vs-potential derivative saturates at 1/beta = the quantum
      //     capacitance, instead of growing exponentially (classical Boltzmann).
      const NumType dosDegeneracy = 2.0; // spin-valley locked K,K' -> g = 2
      const NumType mEff = MoS2Kaasbjerg::relEffMassK * constants::me;
      const NumType N2D = dosDegeneracy * mEff * (constants::kB * temperature) /
                          (2 * constants::pi * constants::hbar * constants::hbar);
      const NumType degeneracyBeta = getMoS2Material().getNi() * tMono / N2D;
      solver.setDegeneracyBeta(degeneracyBeta);

      emcSimulationParameter<NumType, DeviceType> param;
      param.setTimes(3.7e-12, 1e-16, 0.7e-12); // simTime 3.7ps = 0.7ps settle + 3.0ps (30k-step) avg
      param.setAdaptPotentialForWriteFunction(adaptPotential);
      param.setNamePrefix(fileNamePrefix + "Vg" + std::to_string((int)(Vg * 1e3)));
      param.setNrStepsBetweenShowProgress(1000);
      param.setNrStepsForFinalAvg(30000);

      // ---- MoS2 particle + intrinsic scattering (calibrated Kaasbjerg set) ----
      auto electron = std::make_unique<electron2DDevice<NumType, DeviceType>>();
      electron->setDegeneracyBeta(degeneracyBeta); // same statistics as the solver
      MoS2Kaasbjerg::addValleys(electron);
      std::vector<int> idxRegions(device.getDopingProfile().getNrDopingRegions());
      std::iota(idxRegions.begin(), idxRegions.end(), 0);
      MoS2Kaasbjerg::addAcousticScatterMechanisms(electron, idxRegions, temperature);
      MoS2Kaasbjerg::addZeroOrderIntervalleyScatterMechanisms(electron, idxRegions, temperature);
      MoS2Kaasbjerg::addFirstOrderIntervalleyScatterMechanisms(electron, idxRegions, temperature);
      // (a) free-carrier screening of the polar (Froehlich) + piezoelectric
      //     interactions, PER REGION, by the local sheet density n = doping*tMono
      //     (n+ S/D are strongly screened, the channel weakly). envPermittivity=1
      //     matches the bulk Kaasbjerg calibration convention.
      // Single uniform flake (region 0); screening/substrate use the ON-state
      // gate-induced density (the operating point).
      const int channelRegion = 0;
      const auto &dopProf = device.getDopingProfile();
      for (int r : idxRegions) {
        NumType nSheet = (r == channelRegion) ? channelOnStateSheet
                                              : dopProf.getDoping(r) * tMono;
        MoS2Kaasbjerg::addFroehlichScatterMechanisms(electron, {r}, temperature, nSheet, 1.0);
        MoS2Kaasbjerg::addPiezoelectricScatterMechanisms(electron, {r}, temperature, nSheet, 1.0);
      }
      // (b) substrate channels B1-B3 in the channel (region 1), which sits on the
      //     HfO2 gate dielectric: B1 remote surface-optical phonon, B2 charged
      //     interface impurities, B3 interface roughness. These pull the channel
      //     mobility down from ~150 (intrinsic) to the realistic supported value.
      const NumType nSheetCh = channelOnStateSheet;
      MoS2Kaasbjerg::addRemoteSurfaceOpticalPhonon(
          electron, {channelRegion}, temperature, nSheetCh, hfo2EpsInf, hfo2Eps0,
          hfo2wSO1, hfo2wSO2, /*epsEnv=*/1.0);
      MoS2Kaasbjerg::addChargedImpurityScatterMechanism(
          electron, {channelRegion}, temperature, interfaceChargeDensity, nSheetCh,
          interfaceEpsAvg, /*remoteDistance=*/0.0, /*chargeNumber=*/1.0);
      MoS2Kaasbjerg::addSurfaceRoughnessScatterMechanism(
          electron, {channelRegion}, temperature, nSheetCh, interfaceEpsAvg);
      param.addParticleType(std::move(electron));

      SimulationType simulation(param, device, solver, pmScheme);
      // Ramo drift current over the CENTRAL gated channel (70..130nm) only -- away
      // from the resistive ungated spacers (40..60,140..160) where hot carriers in
      // the high access field would otherwise spike the drift sum.
      simulation.setChannelCurrentRegion(60e-9, 140e-9, 80e-9); // full gated channel
      simulation.setPoissonInterval(5); // sub-cycle for a smooth self-consistent field
      auto start = high_resolution_clock::now();
      simulation.execute();
      auto end = high_resolution_clock::now();
      NumType Id = simulation.getAvgDriftCurrent();
      std::cout << "Vg=" << Vg << " Vd=" << Vd << "  I_drift = " << Id << " A"
                << "  (CPU " << duration_cast<seconds>(end - start).count()
                << " s)\n";
      idvg << Vg << " " << Vd << " " << Id << "\n";
      idvg.flush();
    }
  }
  return 0;
}
