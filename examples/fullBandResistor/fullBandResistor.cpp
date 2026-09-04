/*! \brief Phase 4 milestone 1: "bulk in a box".
 *
 * The resistor2D example - a 1 um x 1 um uniformly doped box with two ohmic
 * contacts - driven by FULL-BAND carriers from a .matpkg through
 * emcFullBandParticleHandler, with Poisson, the particle-mesh scheme,
 * contacts and the step loop unchanged from ViennaEMC. At 0.05 V over 1 um
 * the field is 500 V/cm, the field every bulk number in this project was
 * measured at, so the acceptance test is direct: the time-weighted drift
 * velocity divided by V/L must reproduce the bulk mobility for the same
 * package and doping level (Si electrons, 1e16: 1220 +/- 38 cm^2/Vs).
 *
 *   fullBandResistor <pkg> [doping cm^-3 = 1e16] [V = 0.05] [t_ps = 50]
 *                    [transient_ps = 20] [dt_fs = 1] [threads = 4] [width_um = 1]
 *                    [poisson_every = 1] [carriers_per_particle = 1]
 *                    [profile = "" | "N1@x1,N2@x2,N3"]
 * profile: doping regions along x, e.g. "1e17@0.3,1e16@0.7,1e17" = 1e17 for
 * x < 0.3 um, 1e16 for 0.3..0.7 um, 1e17 beyond - an n+/n/n+ resistor. Each
 * region scatters at its own impurity-ladder level (milestone 2).
 */
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>

#include <FullBand/emcFullBandParticleHandler.hpp>
#include <PMSchemes/emcCICScheme.hpp>
#include <PMSchemes/emcNGPScheme.hpp>
#include <ParticleType/emcElectron.hpp>
#include <ParticleType/emcHole.hpp>
#include <PoissonSolver/emcSORSolver.hpp>
#include <emcSimulation.hpp>

#include "../SiliconFunctions.hpp"

const SizeType Dim = 2;
using NumType = double;
using DeviceType = emcDevice<NumType, Dim>;
// CLOUD-IN-CELL BY DEFAULT. With nearest-grid-point the interior mobility
// read 1080-1105 against the bulk 1170 +/- 10 (same package, same rates,
// same <E>) - and 1080 again with 10x the particles per cell, so it is NGP's
// systematic self-force, not statistical mesh noise. Cloud-in-cell reads
// 1211. -DPM_NGP keeps the comparison available.
#ifdef PM_NGP
using PMScheme = emcNGPScheme<NumType, DeviceType>;
#else
using PMScheme = emcCICScheme<NumType, DeviceType>;
#endif
using ParticleHandler = emcFullBandParticleHandler<NumType, DeviceType, PMScheme>;
using PoissonSolver = emcSORSolver<NumType, DeviceType, ParticleHandler>;
using SimulationType = emcSimulation<NumType, DeviceType, PoissonSolver, ParticleHandler, PMScheme>;
using ValueVec = DeviceType::ValueVec;

int main(int argc, char **argv) {
  if (argc < 2) { std::fprintf(stderr, "usage: fullBandResistor <pkg> [doping_cm3] [V] [t_ps] [transient_ps] [dt_fs] [threads]\n"); return 1; }
  const std::string pkg = argv[1];
  const double dopingCm3 = argc > 2 ? std::atof(argv[2]) : 1e16;
  const double voltage = argc > 3 ? std::atof(argv[3]) : 0.05;
  const double tTotal = (argc > 4 ? std::atof(argv[4]) : 50.0) * 1e-12;
  const double tTrans = (argc > 5 ? std::atof(argv[5]) : 20.0) * 1e-12;
  const double dT = (argc > 6 ? std::atof(argv[6]) : 1.0) * 1e-15;
  const int nThreads = argc > 7 ? std::atoi(argv[7]) : 4;
  const double LzUm = argc > 8 ? std::atof(argv[8]) : 1.0;   // device width: particles per cell scale with it
  const int poissonEvery = argc > 9 ? std::atoi(argv[9]) : 1; // solve Poisson every n steps (frozen field test)
  const SizeType carriersPerPart = argc > 10 ? std::atoi(argv[10]) : 1;
  const std::string profile = argc > 11 ? argv[11] : "";
  const bool holes = argc > 12 && std::string(argv[12]) == "hole";   // p-type: emcHole + negative doping
  const double dopSign = holes ? -1.0 : 1.0;
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#endif
  const double Lx = 1e-6, Ly = 1e-6, Lz = LzUm * 1e-6;
  const double hx = 1e-8, hy = 5e-8;

  // full-band handler configuration (static: emcSimulation builds the handler)
  ParticleHandler::packagePath = pkg;
  ParticleHandler::temperature = 300.0;
  ParticleHandler::dopingCm3 = dopingCm3;
  ParticleHandler::avgFromStep = static_cast<SizeType>(tTrans / dT);
  if (std::getenv("ENGINE_NBANDS")) ParticleHandler::nrBands = std::atoi(std::getenv("ENGINE_NBANDS"));

  auto material = Silicon::getSiliconMaterial<NumType>();
  ValueVec maxPos = {Lx, Ly}, spacing = {hx, hy}, origin = {0, 0};
  DeviceType device{material, maxPos, spacing};
  device.setDeviceWidth(Lz);
  if (profile.empty()) {
    device.addConstantDopingRegion(origin, maxPos, dopSign * dopingCm3 * 1e6);   // m^-3, negative = acceptors
  } else {
    // "N1@x1,N2@x2,N3": region i spans [x_{i-1}, x_i), the last to Lx
    double x0 = 0; std::string rest = profile;
    while (!rest.empty()) {
      const auto comma = rest.find(','); std::string tok = rest.substr(0, comma);
      rest = comma == std::string::npos ? "" : rest.substr(comma + 1);
      const auto at = tok.find('@');
      const double N = std::atof(tok.substr(0, at).c_str());
      const double x1 = at == std::string::npos ? Lx : std::atof(tok.substr(at + 1).c_str()) * 1e-6;
      ValueVec lo = {x0, 0}, hi = {x1, Ly};
      device.addConstantDopingRegion(lo, hi, dopSign * N * 1e6);
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
  param.setNamePrefix("fbres");
  param.setNrStepsBetweenShowProgress(5000);
  param.setNrStepsForFinalAvg(static_cast<SizeType>((tTotal - tTrans) / dT));
  // analytic electron type for BOOKKEEPING ONLY (charge, contact reservoir
  // counts, initial carriers per cell). emcParticleType::check() refuses a
  // moved type with no valley, so it gets the analytic X valley - which the
  // full-band handler never calls - and NO scattering mechanisms, so its
  // scatter tables are inert zeros.
  if (holes) {
    auto h = std::make_unique<emcHole<NumType, DeviceType>>(1000, 4, false);
    Silicon::addXValley(h);   // dummy valley for check(); never called
    param.addParticleType(std::move(h));
  } else {
    auto electrons = std::make_unique<emcElectron<NumType, DeviceType>>(1000, 4, false);
    Silicon::addXValley(electrons);
    param.addParticleType(std::move(electrons));
  }

  SimulationType simulation(param, device, solver, pmScheme);
  simulation.setPoissonInterval(poissonEvery);
  const auto t0 = std::chrono::high_resolution_clock::now();
  simulation.execute();
  const auto t1 = std::chrono::high_resolution_clock::now();

  const auto &h = simulation.getParticleHandler();
  const double F = voltage / Lx;                      // V/m, along -x for electrons? sign below
  const double vx = h.meanVelocity(0);
  const double mu = -h.carrier() * vx / F * 1e4;      // electrons drift against E (mu = -v/F), holes with it
  std::printf("# full-band resistor: V=%.3f V over %.2f um -> F=%.0f V/cm, doping %.3g cm^-3, width %.1f um, Poisson every %d, PM=%s\n",
              voltage, Lx * 1e6, F * 1e-2, dopingCm3, Lz * 1e6, poissonEvery,
#ifdef PM_NGP
              "NGP");
#else
              "CIC");
#endif
  std::printf("#   <vx> = %.4e m/s   mu = %.1f cm2/Vs   <E>-CBM = %.4f eV   bands:", vx, mu, h.meanEnergy());
  for (int b = 0; b < h.bands() && b < 8; b++) std::printf(" b%d=%.3f", b, h.bandOccupancy(b));
  std::printf("\n#   events: %lld real, %lld self (%.1f%% self), %lld removed at contacts, %lld reflections\n",
              h.realEvents(), h.selfEvents(),
              100.0 * h.selfEvents() / std::max(1LL, h.realEvents() + h.selfEvents()),
              h.removedAtContacts(), h.reflections());
  std::printf("#   wall %.0f s\n", std::chrono::duration<double>(t1 - t0).count());
  // x-profiles: interior (0.2..0.8 um) vs contact regions
  const SizeType nx = h.profileSize();
  double vIn = 0, wIn = 0, eIn = 0;
  // per-cell fields (velocity, energy, band share) are written by the
  // framework as <prefix>FB*Avg.txt; only the interior summary is printed here
  for (SizeType i = 0; i < nx; i++) {
    const double x = i * hx;
    if (x >= 0.2e-6 && x <= 0.8e-6) { vIn += h.profileVx(i) * h.profileWeight(i); eIn += h.profileE(i) * h.profileWeight(i); wIn += h.profileWeight(i); }
  }
  if (wIn > 0)
    std::printf("#   INTERIOR 0.2-0.8 um: <vx> = %.4e m/s  <E>-CBM = %.4f eV  -> mu(nominal F) = %.1f cm2/Vs\n",
                vIn / wIn, eIn / wIn, -h.carrier() * (vIn / wIn) / F * 1e4);
  return 0;
}
