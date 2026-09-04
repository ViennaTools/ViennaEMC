/*! \brief Phase 4 item 4.1e: the mosfet2D geometry with FULL-BAND carriers.
 *
 * Same device as examples/mosfet2D (125 x 100 nm, 1 nm grid, p substrate
 * 5e17, source/drain n+ 5e19, channel region 5e18 acceptors under a 24 nm
 * gate with 1.2 nm oxide, substrate/source/drain ohmic contacts), driven by
 * emcFullBandParticleHandler on a .matpkg. Differences from the analytic
 * example, both deliberate: cloud-in-cell instead of the example's local
 * NEC-VWD scheme, and DOS-weighted thermal injection instead of its
 * velocity-weighted variant. Each doping region scatters at its own
 * impurity-ladder level. Oxide interface: specular reflection.
 *
 *   fullBandMosfet <pkg> [Vd = 1] [Vg = 1] [t_ps = 10] [transient_ps = 5]
 *                  [dt_fs = 0.15] [carriers_per_particle = 1] [threads = 6]
 */
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <numeric>

#include <FullBand/emcFullBandParticleHandler.hpp>
#include <PMSchemes/emcCICScheme.hpp>
#include <ParticleType/emcElectron.hpp>
#include <PoissonSolver/emcSORSolver.hpp>
#include <emcSimulation.hpp>

#include "../SiliconFunctions.hpp"

const SizeType Dim = 2;
using NumType = double;
using DeviceType = emcDevice<NumType, Dim>;
using PMScheme = emcCICScheme<NumType, DeviceType>;
using ParticleHandler = emcFullBandParticleHandler<NumType, DeviceType, PMScheme>;
using PoissonSolver = emcSORSolver<NumType, DeviceType, ParticleHandler>;
using SimulationType = emcSimulation<NumType, DeviceType, PoissonSolver, ParticleHandler, PMScheme>;
using ValueVec = DeviceType::ValueVec;

NumType adaptPotential(const NumType &pot, const DeviceType &device) {
  return device.getMaterial().getBandGap() / 2 - pot * device.getThermalVoltage();
}

int main(int argc, char **argv) {
  if (argc < 2) { std::fprintf(stderr, "usage: fullBandMosfet <pkg> [Vd] [Vg] [t_ps] [transient_ps] [dt_fs] [carriers/particle] [threads]\n"); return 1; }
  const std::string pkg = argv[1];
  const double Vd = argc > 2 ? std::atof(argv[2]) : 1.0;
  const double Vg = argc > 3 ? std::atof(argv[3]) : 1.0;
  const double tTotal = (argc > 4 ? std::atof(argv[4]) : 10.0) * 1e-12;
  const double tTrans = (argc > 5 ? std::atof(argv[5]) : 5.0) * 1e-12;
  const double dT = (argc > 6 ? std::atof(argv[6]) : 0.15) * 1e-15;
  const SizeType carriersPerPart = argc > 7 ? std::atoi(argv[7]) : 1;
  const int nThreads = argc > 8 ? std::atoi(argv[8]) : 6;
#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#endif
  ParticleHandler::packagePath = pkg;
  ParticleHandler::temperature = 300.0;
  ParticleHandler::dopingCm3 = 0;             // per-region ladder levels from the device
  ParticleHandler::avgFromStep = static_cast<SizeType>(tTrans / dT);

  const ValueVec deviceMaxPos = {125e-9, 100e-9}, spacing = {1e-9, 1e-9};
  DeviceType device{Silicon::getSiliconMaterial<NumType>(), deviceMaxPos, spacing};
  device.setDeviceWidth(1e-6);
  device.addConstantDopingRegion({0, 30e-9}, {125e-9, 100e-9}, -5e23);      // substrate p 5e17
  device.addConstantDopingRegion({0, 0}, {51e-9, 30e-9}, 5e25);             // source n+ 5e19
  device.addConstantDopingRegion({51e-9, 0}, {75e-9, 30e-9}, -5e24);        // channel p 5e18
  device.addConstantDopingRegion({75e-9, 0}, {125e-9, 30e-9}, 5e25);        // drain n+ 5e19
  device.addOhmicContact(emcBoundaryPos::YMAX, 0, {0}, {125e-9});           // substrate
  device.addOhmicContact(emcBoundaryPos::YMIN, 0, {0}, {51e-9});            // source
  device.addGateContact(emcBoundaryPos::YMIN, Vg, {51e-9}, {75e-9}, 3.9, 1.2e-9,
                        device.getMaterial().getBandGap() / 2.);            // gate
  device.addOhmicContact(emcBoundaryPos::YMIN, Vd, {75e-9}, {125e-9});      // drain

  PoissonSolver solver(device, 1e-4, 1.8);
  PMScheme pmScheme;
  emcSimulationParameter<NumType, DeviceType> param;
  param.setTimes(tTotal, dT, tTrans);
  param.setNrCarriersPerPart(carriersPerPart);
  param.setAdaptPotentialForWriteFunction(adaptPotential);
  char pre[64]; std::snprintf(pre, sizeof pre, "fbmosVd%dVg%d", (int)(Vd * 1e3), (int)(Vg * 1e3));
  param.setNamePrefix(pre);
  param.setNrStepsBetweenShowProgress(5000);
  param.setNrStepsForFinalAvg(static_cast<SizeType>((tTotal - tTrans) / dT));
  auto electrons = std::make_unique<emcElectron<NumType, DeviceType>>(1000, 4, false);
  Silicon::addXValley(electrons);            // bookkeeping only; never called
  param.addParticleType(std::move(electrons));

  SimulationType simulation(param, device, solver, pmScheme);
  const auto t0 = std::chrono::high_resolution_clock::now();
  simulation.execute();
  const auto t1 = std::chrono::high_resolution_clock::now();
  const auto &h = simulation.getParticleHandler();
  std::printf("# full-band MOSFET: Vd=%.2f V Vg=%.2f V  dt=%.2f fs  %s\n", Vd, Vg, dT * 1e15, pkg.c_str());
  std::printf("#   <E>-CBM (whole device) = %.4f eV   bands:", h.meanEnergy());
  for (int b = 0; b < h.bands() && b < 8; b++) std::printf(" b%d=%.3f", b, h.bandOccupancy(b));
  std::printf("\n#   events: %lld real, %lld self (%.1f%% self), %lld removed at contacts, %lld reflections\n",
              h.realEvents(), h.selfEvents(), 100.0 * h.selfEvents() / std::max(1LL, h.realEvents() + h.selfEvents()),
              h.removedAtContacts(), h.reflections());
  std::printf("#   wall %.0f s\n", std::chrono::duration<double>(t1 - t0).count());
  return 0;
}
