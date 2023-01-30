#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>

/**
 * Simulate electron transport in bulk silicon with applied background
 * field and includes the real-space particle-particle interactions with the
 * help of the FMM.
 *
 * This example simulates bulk silicon, it uses periodic boundary condition.
 * During the simulation the ensemble average of the drift velocity and the
 * energy is tracked and written to files in the end of the simulation.
 * The resulting files can be plotted with @file plotBulkSimulationFMMResults.py
 * in case the path to the files and the parameters are adapted in the file.
 *
 * Additionally, real-space particle-particle interactions are included with the
 * help of scalFMM.
 */

//! if this is defined a cutoff radius is used for
//! the particle-particle interaction
#define USE_CUTOFF_KERNEL

#include "../SiliconFunctions.hpp"
#include "FMMBulkParticleHandler.hpp"
#include "customHotElectron.hpp"
#include <ParticleType/emcDonor.hpp>
#include <ParticleType/emcElectron.hpp>

const SizeType Dim = 3;
const std::string fileNamePrefix = "bulkSimulationFMM";

using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using ParticleHandler = FMMBulkParticleHandler<NumType, DeviceType>;

using MapIdxTypeToPartType = ParticleHandler::MapIdxToParticleTypes;

const NumType temperature = 300; // in K
const NumType doping = 1e23;     // in 1. / m3

// bulk extent + grid spacing
// careful: if doping is adapted, also adapt extent of bulk
// to get a specific number of particles!
const std::array<NumType, 3> maxPos = {13e-8, 13e-8, 13e-8};
const std::array<NumType, 3> spacing = {1e-8, 1e-8, 1e-8};

//! characteristics of the applied field. Defined by the
//! direction of the field (is normalized in simulation)
//! and the strength of the field!
const NumType appliedFieldStrength = 0; // in V / m
const std::array<NumType, 3> appliedFieldDirection = {0, -1, 0};

//! boolean that determines if the parameters of the simulation should be
//! included in the filename of the resulting files (needed if multiple applied
//! fields are tested)
bool includeParameterInFileName = false;

//! simulation parameter
const NumType dt = 1e-16;        // time of a step [s]
const NumType totalTime = 2e-12; // total simulation time [s]
const SizeType nrStepsBetweenOutput = 100;

//! determines how often the simulation box with all the particles
//! is repeated for potential and force calculation with FMM.
//! number should be in (-1, 0, 1, ...).
//! same as parameter inUpperLeaf in
//! scalFMM/include/Core/FFMMAlgotihmPeriodic.hpp.
const SizeType periodicityParam = 1;

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

  std::cout << "Expected Nr. of Created Particles: " << nrPart << " "
            << particleType->getName() << "\n";
}

int main() {
  // create device (only extent, spacing + dielectric constant important)
  DeviceType device{Silicon::getSiliconMaterial<NumType>(), maxPos, spacing,
                    temperature};
  device.addConstantDopingRegion({0, 0, 0}, maxPos, doping);

  // create particles that should be simulated
  MapIdxTypeToPartType particleTypes;
  // particleTypes[0] =
  //     std::make_unique<emcElectron<NumType, DeviceType>>(1000, 1, false);
  particleTypes[0] =
      std::make_unique<customHotElectron<NumType, DeviceType>>(1, 0.4);
  // add valley + scattering mechanisms to particle type
  Silicon::addXValley(particleTypes[0]);
  Silicon::addAcousticScattering(0, particleTypes[0], device, {0});
  Silicon::addZeroOrderInterValleyScattering(0, particleTypes[0], device, {0});
  Silicon::addFirstOrderInterValleyScattering(0, particleTypes[0], device, {0});

  particleTypes[1] = std::make_unique<emcDonor<NumType, DeviceType>>();

  ParticleHandler handler(device, particleTypes, appliedFieldDirection,
                          appliedFieldStrength, periodicityParam);

  std::cout << "Creating Particles...\n";
  handler.generateInitialParticles();
  handler.printNrParticles();
  handler.print(fileNamePrefix, "Eq");

  std::cout << "Nr repeated Simulation Boxes: " << handler.getNrRepeatedBoxes()
            << "\n";

  SizeType nrSteps = std::ceil(totalTime / dt);
  std::cout << "Simulation Parameter ...\n";
  std::cout << "\tSimulation Time = " << totalTime << " s\n";
  std::cout << "\tStep Time = " << dt << " s\n";
  std::cout << "\tNr. Steps = " << nrSteps << "\n";
  std::cout << "\tApplied El. Field = { " << handler.getAppliedField()
            << " } in V / m\n";

  auto start = std::chrono::high_resolution_clock::now();
  std::vector<std::vector<NumType>> avgEnergy(nrSteps + 1);
  std::vector<std::vector<NumType>> avgDriftVel(nrSteps + 1);
  std::vector<std::vector<NumType>> valleyOcc(nrSteps + 1);

  // store initial characteristics
  avgEnergy[0] = handler.getAvgEnergy(0);
  avgDriftVel[0] = handler.getAvgDriftVelocity(0);
  valleyOcc[0] = handler.getValleyOccupationProbability(0);

  std::cout << "Starting Simulation ..." << std::endl;
  for (SizeType idxStep = 1; idxStep <= nrSteps; idxStep++) {

    handler.executeFMM();
    handler.moveParticles(dt);

    // get average particle characteristics
    avgEnergy[idxStep] = handler.getAvgEnergy(0);
    avgDriftVel[idxStep] = handler.getAvgDriftVelocity(0);
    valleyOcc[idxStep] = handler.getValleyOccupationProbability(0);

    if (idxStep % nrStepsBetweenOutput == 0) {
      std::cout << "\tStep Nr. " << std::to_string(idxStep) << " / " << nrSteps
                << std::endl;
      // handler.print(fileNamePrefix, std::to_string(idxStep));
    }

    handler.resetForcesAndPotential();
    handler.rearrangeTree();
  }
  handler.print(fileNamePrefix, "Final");

  auto end = std::chrono::high_resolution_clock::now();
  std::cout
      << "CPU time: "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " s\n";

  // write parameter string, add used electric field, temperature and number
  // of particles
  std::string parameter = "";
  if (includeParameterInFileName) {
    parameter += "E" + std::to_string((int)appliedFieldStrength);
    parameter += "T" + std::to_string((int)temperature);
    parameter += "N" + std::to_string(handler.getNrParticles(0));
  }

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

  return 0;
}