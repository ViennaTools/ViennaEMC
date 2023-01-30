#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>

#include "../SiliconFunctions.hpp"
#include "basicBulkParticleHandler.hpp"

#include <ParticleType/emcElectron.hpp>
#include <emcDevice.hpp>

/**
 * Simulate electron transport in bulk silicon with applied background
 * field.
 *
 * This example simulates bulk silicon, it uses periodic boundary condition.
 * During the simulation the ensemble average of the drift velocity and the
 * energy is tracked and written to files in the end of the simulation.
 * The resulting files can be plotted with plotBulkSimulationResults.py
 * in case the path to the files and the parameters are adapted in the file.
 */

const SizeType Dim = 3;
const std::string fileNamePrefix = "bulkSimulation";

using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using ParticleHandler = basicBulkParticleHandler<NumType, DeviceType>;
using MapIdxTypeToPartType = ParticleHandler::MapIdxToParticleTypes;

const NumType temperature = 300; // in K
const NumType doping = 1e23;     // in 1. / m3

/// bulk extent + grid spacing which represent the simulation
/// space (with periodic boundary conditions).
/// careful: if doping is adapted, also adapt extent of bulk
/// to get a specific number of particles!
const std::array<NumType, 3> maxPos = {5e-7, 5e-7, 5e-7};
const std::array<NumType, 3> spacing = {1e-7, 1e-7, 1e-7};

//! characteristics of the applied field. Defined by the
//! direction of the field (is normalized in simulation)
//! and the strength of the field!
const std::vector<NumType> appliedFieldStrength = {10e5}; // in V / m
const std::array<NumType, 3> appliedFieldDirection = {-1, 0, 0};

//! boolean that determines if the parameters of the simulation should be
//! included in the filename of the resulting files (needed if multiple applied
//! fields are tested)
bool includeParameterInFileName = false;

// simulation parameter
const NumType dt = 1e-16;        // time of a step [s]
const NumType totalTime = 4e-12; // total simulation time [s]
const SizeType nrStepsBetweenOutput = 10000;

//! Prints the expected number of particles that will be created!
//! If this number is too high / low, either adapt doping or
//! maxPos or spacing!
template <class DerivedParticleType>
void printExpectedNrParticles(
    std::unique_ptr<DerivedParticleType> &particleType,
    const DeviceType &device) {
  DeviceType::SizeVec coord;
  emcGrid<NumType, Dim> pot(device.getGridExtent(), 0);
  SizeType nrPart = 0;
  for (coord.fill(0); !device.isEndCoord(coord); device.advanceCoord(coord))
    nrPart += particleType->getInitialNrParticles(coord, device, pot);

  std::cout << "Expected Nr. of Created Particles: ~ " << nrPart << " "
            << particleType->getName() << "\n";
}

//! number of used threads (in parallel region)
const SizeType nrThreads = 4;

int main() {
#ifdef _OPENMP
  omp_set_num_threads(nrThreads);
  std::cout << ">> Parallel version, using " << nrThreads << " threads.\n\n";
#else
  std::cout << "\n>> Sequential version.\n\n";
#endif

  // create geometry (only extent, spacing + dielectric constant important)
  DeviceType device{Silicon::getSiliconMaterial<NumType>(), maxPos, spacing,
                    temperature};
  device.addConstantDopingRegion({0, 0, 0}, maxPos, doping);

  // add particle types of interest
  MapIdxTypeToPartType particleTypes;
  particleTypes[0] =
      std::make_unique<emcElectron<NumType, DeviceType>>(1000, 1., false);

  // add valleys to particle type
  Silicon::addXValley(particleTypes[0]);

  // add scattering mechanisms to particleType
  Silicon::addAcousticScattering(0, particleTypes[0], device, {0});
  // Silicon::addCoulombScattering(0, particleTypes[0], device, {0});
  Silicon::addZeroOrderInterValleyScattering(0, particleTypes[0], device, {0});
  Silicon::addFirstOrderInterValleyScattering(0, particleTypes[0], device, {0});
  // Silicon::addGrainScattering(particleTypes[0], 0.5, 1);

  // create particle handler
  ParticleHandler handler(device, particleTypes, appliedFieldDirection);
  // handler.setSeed(1);

  // write the simulation parameter
  const SizeType nrSteps = std::ceil(totalTime / dt);
  std::cout << "Simulation Parameter ...\n";
  std::cout << "\tSimulation Time = " << totalTime << " s\n";
  std::cout << "\tStep Time = " << dt << " s\n";
  std::cout << "\tNr. Steps = " << nrSteps << "\n";

  // do simulation for each applied field strength
  auto start = std::chrono::high_resolution_clock::now();
  for (auto &field : appliedFieldStrength) {
    handler.resetAppliedFieldStrength(field);
    std::cout << "\tApplied El. Field = { " << handler.getAppliedField()
              << " } in V / m\n";

    // create initial particles + write their characteristics
    std::cout << "Creating Particles...\n";
    handler.generateInitialParticles();
    handler.printNrParticles();
    handler.print(fileNamePrefix, "Eq");

    // write parameter string, add used electric field, temperature and number
    // of particles
    std::string parameter = "";
    if (includeParameterInFileName) {
      parameter += "E" + std::to_string((int)field);
      parameter += "T" + std::to_string((int)temperature);
      parameter += "N" + std::to_string(handler.getNrParticles(0));
    }

    std::vector<std::vector<NumType>> avgEnergy(nrSteps + 1);
    std::vector<std::vector<NumType>> avgDriftVel(nrSteps + 1);
    std::vector<std::vector<NumType>> valleyOcc(nrSteps + 1);

    // store initial characteristics
    avgEnergy[0] = handler.getAvgEnergy(0);
    avgDriftVel[0] = handler.getAvgDriftVelocity(0);
    valleyOcc[0] = handler.getValleyOccupationProbability(0);

    // perform simulation
    std::cout << "Starting Simulation ...\n";
    for (SizeType idxStep = 1; idxStep <= nrSteps; idxStep++) {
      handler.moveParticles(dt);

      // get average particle characteristics
      avgEnergy[idxStep] = handler.getAvgEnergy(0);
      avgDriftVel[idxStep] = handler.getAvgDriftVelocity(0);
      valleyOcc[idxStep] = handler.getValleyOccupationProbability(0);

      if (idxStep % nrStepsBetweenOutput == 0) {
        std::cout << "\tStep Nr. " << std::to_string(idxStep) << " / "
                  << nrSteps << "\n";
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout
        << "CPU time: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " s\n";

    // write output for current results
    std::ofstream energyFile, occFile, driftVelocityFile;
    energyFile.open(fileNamePrefix + "AvgEnergy" + parameter + ".txt");
    driftVelocityFile.open(fileNamePrefix + "AvgDriftVelocity" + parameter +
                           ".txt");
    occFile.open(fileNamePrefix + "valleyOccupation" + parameter + ".txt");
    for (SizeType idxStep = 0; idxStep < avgEnergy.size(); idxStep++) {
      energyFile << idxStep * dt << " " << avgEnergy[idxStep] << "\n";
      driftVelocityFile << idxStep * dt << " " << avgDriftVel[idxStep] << "\n";
      occFile << idxStep * dt << " " << valleyOcc[idxStep] << "\n";
    }
    energyFile.close();
    driftVelocityFile.close();
    occFile.close();

    handler.deleteParticles();
  }

  return 0;
}