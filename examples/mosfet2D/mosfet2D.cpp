#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>

#include <ParticleHandler/emcBasicParticleHandler.hpp>
#include <PoissonSolver/emcSORSolver.hpp>
#include <SurfaceScatterMechanisms/emcConstantSurfaceScatterMechanism.hpp>
#include <emcDevice.hpp>
#include <emcSimulation.hpp>

#include "../SiliconFunctions.hpp"
#include "NECSchemeVWD.hpp"
#include "electronVWD.hpp"

const SizeType Dim = 2;
const std::string fileNamePrefix = "mosfet";

using namespace std::chrono;
using NumType = double;

using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using PMScheme = emcNECSchemeVWD<NumType, DeviceType>;
using ParticleHandler = emcBasicParticleHandler<NumType, DeviceType, PMScheme>;
using PoissonSolver = emcSORSolver<NumType, DeviceType, ParticleHandler>;
using SimulationType = emcSimulation<NumType, DeviceType, PoissonSolver,
                                     ParticleHandler, PMScheme>;

using ValueVec = DeviceType::ValueVec;

/// adapts potential to the ViennaWD-Output
NumType adaptPotential(const NumType &pot, const DeviceType &device) {
  return device.getMaterial().getBandGap() / 2 -
         pot * device.getThermalVoltage();
}

const std::vector<NumType> Vds = {1.};
const std::vector<NumType> Vgs = {1.};

bool includeParameterInFileName = true;

/// parameters to include surface roughness at the Si-SiO2 interface
bool addRoughnessToOxideInterface = false;
/// probability of a specular scattering (in the other times diffusive
/// scattering is performed)
NumType probSpecularScattering = 0.5;

//! number of used threads (in parallel region)
const SizeType nrThreads = 4;

/// creates an nChannel-MOSFET using the same parameter as ViennaWD Example 1
int main() {
#ifdef _OPENMP
  omp_set_num_threads(nrThreads);
  std::cout << ">> Parallel version, using " << nrThreads << " threads.\n\n";
#else
  std::cout << "\n>> Sequential version.\n\n";
#endif
  for (auto &Vg : Vgs) {
    for (auto &Vd : Vds) {
      // create Device (n-channel MOSFET)
      const ValueVec deviceMaxPos = {125e-9, 100e-9}; // m
      const ValueVec spacing = {1e-9, 1e-9};          // m
      DeviceType device{Silicon::getSiliconMaterial<NumType>(), deviceMaxPos,
                        spacing};
      device.setDeviceWidth(1e-6);
      // add doping regions
      device.addConstantDopingRegion({0, 30e-9}, {125e-9, 100e-9},
                                     -5e23);                        // bulk
      device.addConstantDopingRegion({0, 0}, {51e-9, 30e-9}, 5e25); // source
      device.addConstantDopingRegion({51e-9, 0}, {75e-9, 30e-9}, -5e24); // gate
      device.addConstantDopingRegion({75e-9, 0}, {125e-9, 30e-9},
                                     5e25); // drain

      // add Contacts
      device.addOhmicContact(emcBoundaryPos::YMAX, 0, {0},
                             {125e-9}); // substrate
      device.addOhmicContact(emcBoundaryPos::YMIN, 0, {0}, {51e-9}); // source
      device.addGateContact(emcBoundaryPos::YMIN, Vg, {51e-9}, {75e-9}, 3.9,
                            1.2e-9,
                            device.getMaterial().getBandGap() / 2.); // gate
      device.addOhmicContact(emcBoundaryPos::YMIN, Vd, {75e-9},
                             {125e-9}); // drain

      std::cout << "Applied Voltages:\n";
      std::cout << "\tVd = " << Vd << " V\n";
      std::cout << "\tVg = " << Vg << " V\n";

      // create Poisson Solver + PMScheme
      PoissonSolver solver(device, 1e-4, 1.8);
      PMScheme pmScheme;

      // set simulation parameter + create simulation
      emcSimulationParameter<NumType, DeviceType> param;
      param.setTimes(10e-12, 1.5e-16, 5e-12);
      param.setAdaptPotentialForWriteFunction(adaptPotential);
      std::string parameter = "";
      if (includeParameterInFileName) {
        parameter += "Vd" + std::to_string((int)(Vd * 1e3));
        parameter += "Vg" + std::to_string((int)(Vg * 1e3));
      }
      param.setNamePrefix(fileNamePrefix + parameter);
      param.setNrStepsBetweenShowProgress(1000);
      param.setNrStepsForFinalAvg(6667);

      // add particle types that should be simulated
      auto electron = std::make_unique<electronVWD<NumType, DeviceType>>();
      Silicon::addXValley(electron);

      std::vector<int> idxRegions(
          device.getDopingProfile().getNrDopingRegions());
      std::iota(idxRegions.begin(), idxRegions.end(), 0);
      Silicon::addAcousticScattering(0, electron, device, idxRegions);
      Silicon::addCoulombScattering(0, electron, device, idxRegions);
      Silicon::addZeroOrderInterValleyScattering(0, electron, device,
                                                 idxRegions);

      if (addRoughnessToOxideInterface) {
        electron->setSurfaceScatterMechanism(
            emcBoundaryPos::YMIN,
            std::make_unique<
                emcConstantSurfaceScatterMechanism<NumType, DeviceType>>(
                probSpecularScattering, device.getMaxPos()));
      }

      param.addParticleType(std::move(electron));

      SimulationType simulation(param, device, solver, pmScheme);

      // execute + time simulation
      auto start = high_resolution_clock::now();
      simulation.execute();
      auto end = high_resolution_clock::now();
      std::cout << "CPU time: " << duration_cast<seconds>(end - start).count()
                << " s\n";
    }
  }
  return 0;
}