#include <chrono>
#include <memory>

#include <PMSchemes/emcNGPScheme.hpp>
#include <ParticleHandler/emcBasicParticleHandler.hpp>
#include <ParticleType/emcElectron.hpp>
#include <emcGrainScatterMechanism.hpp>
#include <emcSimulation.hpp>

//Different Solvers
#include <PoissonSolver/emcSORSolver.hpp>
#include <PoissonSolver/emcJacobiSolver.hpp>
#include <PoissonSolver/emcBicgSTABSolver.hpp>

#include "../SiliconFunctions.hpp"

using namespace std::chrono;

const SizeType Dim = 2;
const std::string fileNamePrefix = "resistor";

using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using PMScheme = emcNGPScheme<NumType, DeviceType>;
using ParticleHandler = emcBasicParticleHandler<NumType, DeviceType, PMScheme>;
using PoissonSolver = emcSORSolver<NumType, DeviceType, ParticleHandler>; //!< adapt, if change the solver
using SimulationType = emcSimulation<NumType, DeviceType, PoissonSolver,
                                     ParticleHandler, PMScheme>;

using GridType = emcGrid<NumType, Dim>;
using ValueVec = DeviceType::ValueVec;
using SizeVec = DeviceType::SizeVec;

// boolean that determines if the parameters of the simulation should be
// included in the filename of the resulting files
bool includeParameterInFileName = true;

emcMaterial<NumType> material = Silicon::getSiliconMaterial<NumType>();
NumType device_extent_x = 1e-6; // m
NumType device_extent_y = 1e-6; // m
NumType device_extent_z = 1e-6; // m

NumType spacing_x = 1e-8; // m
NumType spacing_y = 5e-8; // m

std::vector<NumType> appliedVoltages = {0.05}; // [V]

NumType dT = 1e-15;            // s
NumType totalSimTime = 5e-11;  // s
NumType transientTime = 2e-11; // s

// add grain scattering
bool addGrainScattering = false;
// mean rate at which grain scattering appears
NumType grainScatterRate = 1e13; // 1 / s
// probability of transmission through grain boundary
NumType transmissionProb = 1;

//! number of used threads (in parallel region)
const SizeType nrThreads = 10;

/// example device: Resistor
int main() {
#ifdef _OPENMP
  omp_set_num_threads(nrThreads);
  std::cout << ">> Parallel version, using " << nrThreads << " threads.\n\n";
#else
  std::cout << "\n>> Sequential version.\n\n";
#endif
  for (auto voltage : appliedVoltages) {
    std::cout << "Current applied voltage = " << voltage << "\n";
    auto startSetup = high_resolution_clock::now();

    // create device
    ValueVec deviceMaxPos = {device_extent_x, device_extent_y};
    ValueVec spacing = {spacing_x, spacing_y};
    DeviceType device{material, deviceMaxPos, spacing};
    device.setDeviceWidth(device_extent_z);

    // add n-doped region
    ValueVec origin = {0, 0};
    device.addConstantDopingRegion(origin, deviceMaxPos, 1e22);

    // add Contacts
    device.addOhmicContact(emcBoundaryPos::XMAX, 0, {origin[1]},
                           {deviceMaxPos[1]});
    device.addOhmicContact(emcBoundaryPos::XMIN, voltage, {origin[1]},
                           {deviceMaxPos[1]});

    /*! \brief Set poisson solver
     * for BICGSTAB and Jacobi the third parameter give the max. allowed Iterations
     * for SOR represent the third parameter omega (appropriate values in (0,2))
     */
    PoissonSolver solver(device, 1e-5, 1.8); //!< adapt if needed
    PMScheme pmScheme;

    // set Simulation Parameter
    emcSimulationParameter<NumType, DeviceType> param;
    param.setTimes(totalSimTime, dT, transientTime);

    std::string parameter = "";
    if (includeParameterInFileName) {
      // applied voltage in mV appended to filename
      parameter += "V" + std::to_string((int)(voltage * 1000));
      // time in attoseconds appended to filename
      parameter += "as" + std::to_string((int)(dT * 1e18));
    }

    param.setNamePrefix(fileNamePrefix + parameter);
    param.setNrStepsBetweenShowProgress(5000);
    param.setNrStepsForFinalAvg(20000);

    // add simulated Particle Type, Valleys and Scatter Mechanisms
    auto electrons =
        std::make_unique<emcElectron<NumType, DeviceType>>(1000, 4, false);
    Silicon::addXValley(electrons);
    Silicon::addAcousticScattering(0, electrons, device, {0});
    Silicon::addZeroOrderInterValleyScattering(0, electrons, device, {0});
    Silicon::addFirstOrderInterValleyScattering(0, electrons, device, {0});
    Silicon::addCoulombScattering(0, electrons, device, {0});

    if (addGrainScattering) {
      electrons->setGrainScatterMechanism(
          std::make_unique<emcGrainScatterMechanism<NumType>>(
              transmissionProb, grainScatterRate));
    }

    param.addParticleType(std::move(electrons));

    SimulationType simulation(param, device, solver, pmScheme);

    // execute + time simulation
    auto start = high_resolution_clock::now();
    simulation.execute();
    auto end = high_resolution_clock::now();
    std::cout << "CPU time: " << duration_cast<seconds>(end - start).count() << " s\n";
  }

  return 0;
}