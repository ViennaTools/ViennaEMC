/*! \brief Phase 4 item 4.1c: the ANALYTIC control for the full-band n+/n/n+.
 *
 * The same device as examples/fullBandResistor with a profile - same box,
 * grid, contacts, Poisson, cloud-in-cell, carriers per particle, times - run
 * with ViennaEMC's own analytic silicon model (X valleys, acoustic, f/g
 * intervalley, Brooks-Herring Coulomb from SiliconFunctions.hpp). Only the
 * carrier model differs, so whatever differs in the result is the carrier
 * model. Same command line as fullBandResistor minus the package:
 *   analyticResistorNNN [doping cm^-3] [V] [t_ps] [transient_ps] [dt_fs]
 *                       [threads] [width_um] [poisson_every]
 *                       [carriers_per_particle] [profile]
 */
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>

#include <PMSchemes/emcCICScheme.hpp>
#include <ParticleHandler/emcBasicParticleHandler.hpp>
#include <ParticleType/emcElectron.hpp>
#include <ParticleType/emcHole.hpp>
#include <ValleyTypes/emcParabolicIsotropValley.hpp>
#include <PoissonSolver/emcSORSolver.hpp>
#include <emcSimulation.hpp>

#include "../SiliconFunctions.hpp"

const SizeType Dim = 2;
using NumType = double;
using DeviceType = emcDevice<NumType, Dim>;
using PMScheme = emcCICScheme<NumType, DeviceType>;
using ParticleHandler = emcBasicParticleHandler<NumType, DeviceType, PMScheme>;
using PoissonSolver = emcSORSolver<NumType, DeviceType, ParticleHandler>;
using SimulationType = emcSimulation<NumType, DeviceType, PoissonSolver, ParticleHandler, PMScheme>;
using ValueVec = DeviceType::ValueVec;

int main(int argc, char **argv) {
  const double dopingCm3 = argc > 1 ? std::atof(argv[1]) : 1e16;
  const double voltage = argc > 2 ? std::atof(argv[2]) : 0.05;
  const double tTotal = (argc > 3 ? std::atof(argv[3]) : 50.0) * 1e-12;
  const double tTrans = (argc > 4 ? std::atof(argv[4]) : 20.0) * 1e-12;
  const double dT = (argc > 5 ? std::atof(argv[5]) : 1.0) * 1e-15;
  const int nThreads = argc > 6 ? std::atoi(argv[6]) : 4;
  const double LzUm = argc > 7 ? std::atof(argv[7]) : 1.0;
  const int poissonEvery = argc > 8 ? std::atoi(argv[8]) : 1;
  const SizeType carriersPerPart = argc > 9 ? std::atoi(argv[9]) : 1;
  const std::string profile = argc > 10 ? argv[10] : "";
  const bool holes = argc > 11 && std::string(argv[11]) == "hole";
  const double dopSign = holes ? -1.0 : 1.0;
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#endif
  const double Lx = 1e-6, Ly = 1e-6, Lz = LzUm * 1e-6, hx = 1e-8, hy = 5e-8;
  auto material = Silicon::getSiliconMaterial<NumType>();
  ValueVec maxPos = {Lx, Ly}, spacing = {hx, hy}, origin = {0, 0};
  DeviceType device{material, maxPos, spacing};
  device.setDeviceWidth(Lz);
  int nRegions = 0;
  if (profile.empty()) {
    device.addConstantDopingRegion(origin, maxPos, dopSign * dopingCm3 * 1e6); nRegions = 1;
  } else {
    double x0 = 0; std::string rest = profile;
    while (!rest.empty()) {
      const auto comma = rest.find(','); std::string tok = rest.substr(0, comma);
      rest = comma == std::string::npos ? "" : rest.substr(comma + 1);
      const auto at = tok.find('@');
      const double N = std::atof(tok.substr(0, at).c_str());
      const double x1 = at == std::string::npos ? Lx : std::atof(tok.substr(at + 1).c_str()) * 1e-6;
      ValueVec lo = {x0, 0}, hi = {x1, Ly};
      device.addConstantDopingRegion(lo, hi, dopSign * N * 1e6); nRegions++;
      std::printf("# doping region [%.2f, %.2f) um: %.3g cm^-3\n", x0 * 1e6, x1 * 1e6, N);
      x0 = x1;
    }
  }
  device.addOhmicContact(emcBoundaryPos::XMAX, 0, {origin[1]}, {maxPos[1]});
  device.addOhmicContact(emcBoundaryPos::XMIN, voltage, {origin[1]}, {maxPos[1]});
  PoissonSolver solver(device, 1e-4, 1.8);
  solver.setMobileSpeciesHoles(holes);   // p-type: the simulated species is holes
  PMScheme pmScheme;
  emcSimulationParameter<NumType, DeviceType> param;
  param.setTimes(tTotal, dT, tTrans);
  param.setNrCarriersPerPart(carriersPerPart);
  param.setNamePrefix("anres");
  param.setNrStepsBetweenShowProgress(5000);
  param.setNrStepsForFinalAvg(static_cast<SizeType>((tTotal - tTrans) / dT));
  std::vector<int> regions; for (int r = 0; r < nRegions; r++) regions.push_back(r);
  if (holes) {
    // EMPIRICAL-BAND SILICON HOLES, the standard fitted model (Jacoboni &
    // Reggiani, Rev. Mod. Phys. 55, 645 (1983), simplified to isotropic
    // parabolic bands): heavy hole m* = 0.49, light hole m* = 0.16, both at
    // the valence-band edge, degeneracy 1 each; acoustic deformation potential
    // 5.0 eV; non-polar optical D_tK = 6.0e8 eV/cm, hbar*omega = 61.2 meV,
    // intra-band and inter-band (HH <-> LH); Brooks-Herring Coulomb per
    // region. Every number here is a FIT to silicon data - this is the model
    // the full-band package replaces, run as the control.
    //
    // RESULT (2026-09-04): NOT a usable control as parameterised. In a 1 um
    // resistor at 500 V/cm it reads mu_h = 2775 (1e16) and 1706 (1e17)
    // against the measured ~450 / ~330 and the full-band 447 / 344 - it
    // under-scatters ~5x. Not an injection artefact (the snapshot is 80/20
    // heavy/light, thermal is 84/16): the literature parameters were fitted
    // together with warped bands and overlap factors this isotropic
    // implementation does not have, and a hole has ONE optical partner band
    // where a silicon electron has five intervalley partners. Making it read
    // 450 would mean tuning D and Xi to the answer - a fit, not a control.
    // Kept as a demonstration of why no empirical silicon hole model ships
    // with the framework; the full-band hole device is compared against the
    // full-band bulk and experiment instead.
    auto h = std::make_unique<emcHole<NumType, DeviceType>>(1000, 4, false);
    using HoleValley = emcParabolicIsotropValley<NumType>;
    h->addValley(std::make_unique<HoleValley>(0.49, h->getMass(), 1));   // 0: heavy
    h->addValley(std::make_unique<HoleValley>(0.16, h->getMass(), 1));   // 1: light
    const NumType sigmaAcH = 5.0, dOptH = 6.0e10, hwOptH = 0.0612;
    // the intervalley mechanism's map is over SUB-valleys (degeneracy 1 here:
    // {0 -> {0}}); the final VALLEY is a separate argument, so intra-band
    // (HH->HH, LH->LH) and inter-band (HH<->LH) optical are two instances
    const std::map<SizeType, std::vector<SizeType>> sub0 = {{0, {0}}};
    for (SizeType v = 0; v < 2; v++) {
      h->addScatterMechanism(regions, std::make_unique<Acoustic>(v, sigmaAcH, device));
      h->addScatterMechanism(regions, std::make_unique<ZeroIvAb>("Oi", v, sub0, dOptH, hwOptH, device));
      h->addScatterMechanism(regions, std::make_unique<ZeroIvEm>("Oi", v, sub0, dOptH, hwOptH, device));
      h->addScatterMechanism(regions, std::make_unique<ZeroIvAb>("Oe", v, 1 - v, sub0, dOptH, hwOptH, device));
      h->addScatterMechanism(regions, std::make_unique<ZeroIvEm>("Oe", v, 1 - v, sub0, dOptH, hwOptH, device));
      h->addScatterMechanism(regions, std::make_unique<emcCoulombScatterMechanism<NumType, DeviceType>>(v, Silicon::epsR, device));
    }
    param.addParticleType(std::move(h));
  } else {
    auto electrons = std::make_unique<emcElectron<NumType, DeviceType>>(1000, 4, false);
    Silicon::addXValley(electrons);
    Silicon::addAcousticScattering(0, electrons, device, regions);
    Silicon::addZeroOrderInterValleyScattering(0, electrons, device, regions);
    Silicon::addFirstOrderInterValleyScattering(0, electrons, device, regions);
    Silicon::addCoulombScattering(0, electrons, device, regions);
    param.addParticleType(std::move(electrons));
  }
  SimulationType simulation(param, device, solver, pmScheme);
  simulation.setPoissonInterval(poissonEvery);
  const auto t0 = std::chrono::high_resolution_clock::now();
  simulation.execute();
  const auto t1 = std::chrono::high_resolution_clock::now();
  std::printf("# analytic resistor (%s): V=%.3f V over %.2f um -> F=%.0f V/cm, profile '%s', %d region(s), CIC, %zu carriers/particle\n",
              holes ? "HOLES: HH 0.49 + LH 0.16, ac 5 eV, opt 6e8 eV/cm @ 61 meV, BH" : "electrons: Si X valleys",
              voltage, Lx * 1e6, voltage / Lx * 1e-2, profile.c_str(), nRegions, (size_t)carriersPerPart);
  std::printf("#   wall %.0f s\n", std::chrono::duration<double>(t1 - t0).count());
  return 0;
}
