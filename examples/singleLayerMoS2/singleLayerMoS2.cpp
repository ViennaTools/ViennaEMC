#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <emcDevice.hpp>
#include <emcUtil.hpp>

#include "../bulkSimulation/basicBulkParticleHandler.hpp"
#include "electron2D.hpp"
#include "parameterKaasbjerg.hpp"
#include "parameterLi.hpp"
#include "parameterPilotto.hpp"

const SizeType Dim = 3;
const std::string fileNamePrefix = "singleLayerMoS2";

using namespace std::chrono;
using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using ParticleType = emcParticleType<NumType, DeviceType>;
using ParticleHandler = basicBulkParticleHandler<NumType, DeviceType>;

using MapIdxTypeToPartType = ParticleHandler::MapIdxToParticleTypes;
using ParticleType = ParticleHandler::ParticleType;
using ValueVec = DeviceType::ValueVec;
using SizeVec = DeviceType::SizeVec;

// options for different applied electric field strengths (in V / m)
enum struct AppliedFieldsType : SizeType { NO, LOW, HIGH, CUSTOM };
const std::vector<NumType> noField = {0};
const std::vector<NumType> lowFields = {
    1e4,  2e4,  3e4,  4e4,  5e4,  6e4,  7e4,  8e4,  9e4,  10e4,
    11e4, 12e4, 13e4, 14e4, 15e4, 16e4, 17e4, 18e4, 19e4, 20e4};
const std::vector<NumType> highFields = {1e5,   2e5,   5e5,   10e5,  20e5,
                                         40e5,  60e5,  80e5,  100e5, 150e5,
                                         200e5, 250e5, 300e5, 350e5, 400e5};

// options to select the parameter for the valleys and scatter mechanisms
// given in different papers
enum struct PaperType : SizeType { KAASBJERG, LI, PILOTTO };

// customizable simulation parameter ----------------------------------------
// bulk parameter (size of simulated region + spacing)
const ValueVec maxPos = {5e-7, 5e-7, 0.65e-9};  // in m
const ValueVec spacing = {1e-8, 1e-8, 0.65e-9}; // in m

const NumType temperature = 300; // in K

const NumType dt = 1e-16;        // time step in s
const NumType totalTime = 2e-12; // total sim. time in s
const SizeType nrStepsBetweenOutput = 10000;
// const SizeType nrStepsBetweenOutput = 100;

// set the direction and strengths (in V / m) of the tested electric fields
const ValueVec appliedFieldDir = {1, 0, 0};
AppliedFieldsType appliedFields = AppliedFieldsType::HIGH;
// those applied fields are used if CUSTOM is selected as applied fields type
const std::vector<NumType> customFields = {40e5};

// sets the paper from which the parameter for the conduction band valley
// and scatter mechanisms are taken from (see more in the included files)
PaperType selectedPaperForParameter = PaperType::PILOTTO;

// boolean that determines if the parameters of the simulation should be
// included in the filename of the resulting files (needed if multiple applied
// fields are tested)
bool includeParameterInFileName = true;
// boolean that determines if the velocity of all particles is sampled every
// time the number of steps is outputted. Needed only for the calculation of the
// velocity auto correlation function.
bool sampleVelocityOfParticles = false;

//! number of used threads (in parallel region)
const SizeType nrThreads = 4;
// -------------------------------------------------------------------

template <class DerivedParticleType>
void setKaasbjergParameter(std::unique_ptr<DerivedParticleType> &particleType) {
  MoS2Kaasbjerg::addValleys(particleType);
  MoS2Kaasbjerg::addAcousticScatterMechanisms(particleType, {0}, temperature);
  MoS2Kaasbjerg::addZeroOrderIntervalleyScatterMechanisms(particleType, {0},
                                                          temperature);
  MoS2Kaasbjerg::addFirstOrderIntervalleyScatterMechanisms(particleType, {0},
                                                           temperature);
  MoS2Kaasbjerg::addFroehlichScatterMechanisms(particleType, {0}, temperature);
}

template <class DerivedParticleType>
void setLiParameter(std::unique_ptr<DerivedParticleType> &particleType) {
  MoS2Li::addValleys(particleType);
  MoS2Li::addAcousticScatterMechanisms(particleType, {0}, temperature);
  MoS2Li::addZeroOrderIntervalleyScatterMechanisms(particleType, {0},
                                                   temperature);
}

template <class DerivedParticleType>
void setPilottoParameter(std::unique_ptr<DerivedParticleType> &particleType) {
  MoS2Pilotto::addValleys(particleType);
  MoS2Pilotto::addAcousticScatterMechanisms(particleType, {0}, temperature);
  MoS2Pilotto::addZeroOrderIntervalleyScatterMechanisms(particleType, {0},
                                                        temperature);
}

int main() {
#ifdef _OPENMP
  omp_set_num_threads(nrThreads);
  std::cout << ">> Parallel version, using " << nrThreads << " threads.\n\n";
#else
  std::cout << "\n>> Sequential version.\n\n";
#endif

  // print simulation parameter
  SizeType nrSteps = std::ceil(totalTime / dt);
  std::cout << "Simulation Parameter ...\n";
  std::cout << "\tSimulation Time = " << totalTime << " s\n";
  std::cout << "\tStep Time = " << dt << " s\n";
  std::cout << "\tNr. Steps = " << nrSteps << "\n";
  std::cout << "\tTemperature = " << temperature << " K\n";

  // set applied electric fields
  std::vector<NumType> appliedFieldStrengths;
  switch (appliedFields) {
  case AppliedFieldsType::NO:
    appliedFieldStrengths = noField;
    std::cout << "\tTested Fields = No Field\n";
    break;
  case AppliedFieldsType::LOW:
    appliedFieldStrengths = lowFields;
    std::cout << "\tTested Fields = Low-Fields\n";
    break;
  case AppliedFieldsType::HIGH:
    appliedFieldStrengths = highFields;
    std::cout << "\tTested Fields = High-Fields\n";
    break;
  default: // CUSTOM case
    appliedFieldStrengths = customFields;
    std::cout << "\tTested Fields = Custom Fields\n";
  }

  // set device characteristics: material class + amount of doping in doping
  // region not important. only maxPos + spacing important as it determines the
  // number of used electrons.
  MaterialType MoS2{1, 1, 1, 1, 1};
  DeviceType device{MoS2, maxPos, spacing};
  device.addConstantDopingRegion({0, 0, 0}, maxPos, 1);

  // set particle type: electron with adapted valleys and scatter mechanisms
  MapIdxTypeToPartType particleTypes;
  particleTypes[0] = std::make_unique<electron2D<NumType, DeviceType>>();
  switch (selectedPaperForParameter) {
  case PaperType::KAASBJERG:
    setKaasbjergParameter(particleTypes[0]);
    std::cout << "\tUsed Parameter from Paper = Kaasbjerg\n";
    break;
  case PaperType::LI:
    setLiParameter(particleTypes[0]);
    std::cout << "\tUsed Parameter from Paper = Li\n";
    break;
  case PaperType::PILOTTO:
    setPilottoParameter(particleTypes[0]);
    std::cout << "\tUsed Parameter from Paper = Pilotto\n";
    break;
  }

  ParticleHandler handler(device, particleTypes, appliedFieldDir);

  auto start = std::chrono::high_resolution_clock::now();
  for (auto &appliedField : appliedFieldStrengths) {
    std::cout << "Current Field : " << appliedField / 1e5 << " kV / cm\n";

    handler.resetAppliedFieldStrength(appliedField);
    std::cout << "\tCreating Particles...\n";
    handler.generateInitialParticles();
    handler.printNrParticles();
    // handler.print(fileNamePrefix, "Eq" + parameter);

    // write parameter string, add used electric field, temperature and number
    // of particles
    std::string parameter = "";
    if (includeParameterInFileName) {
      parameter += "E" + std::to_string((int)appliedField);
      parameter += "T" + std::to_string((int)temperature);
      parameter += "N" + std::to_string(handler.getNrParticles(0));
    }

    std::vector<std::vector<NumType>> avgEnergy(nrSteps + 1);
    std::vector<std::vector<NumType>> avgDriftVel(nrSteps + 1);
    std::vector<std::vector<NumType>> valleyOcc(nrSteps + 1);

    std::ofstream osVel;
    if (sampleVelocityOfParticles) {
      osVel.open(fileNamePrefix + "velocityAtDiffTimes" + parameter + ".txt");
      osVel << dt * nrStepsBetweenOutput << " " << handler.getNrParticles(0)
            << "\n";
    }

    // store initial characteristics
    avgEnergy[0] = handler.getAvgEnergy(0);
    avgDriftVel[0] = handler.getAvgDriftVelocity(0);
    valleyOcc[0] = handler.getValleyOccupationProbability(0);

    // perform simulation
    std::cout << "\tStarting Simulation ...\n";
    for (SizeType idxStep = 1; idxStep <= nrSteps; idxStep++) {
      handler.moveParticles(dt);

      avgEnergy[idxStep] = handler.getAvgEnergy(0);
      avgDriftVel[idxStep] = handler.getAvgDriftVelocity(0);
      valleyOcc[idxStep] = handler.getValleyOccupationProbability(0);

      if (idxStep % nrStepsBetweenOutput == 0) {
        std::cout << "\t\tStep Nr. " << std::to_string(idxStep) << " / "
                  << nrSteps << "\n";
        if (sampleVelocityOfParticles)
          handler.printVelocities(osVel);
      }
    }

    if (sampleVelocityOfParticles)
      osVel.close();

    // handler.print(fileNamePrefix, "Final" + parameter);

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

    auto end = std::chrono::high_resolution_clock::now();
    std::cout
        << "\tCurrent time passed: "
        << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
        << " s\n";
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::cout
      << "Total time passed: "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " s\n";

  return 0;
}