#ifndef EMC_SIMULATION_HPP
#define EMC_SIMULATION_HPP

#include <PMSchemes/emcAbstractPMScheme.hpp>
#include <ParticleHandler/emcAbstractParticleHandler.hpp>
#include <PoissonSolver/emcAbstractSolver.hpp>
#include <emcDevice.hpp>
#include <emcOutput.hpp>
#include <emcParticleInitialization.hpp>
#include <emcSimulationParameter.hpp>
#include <emcSimulationResults.hpp>
#include <emcUtil.hpp>

/*! \brief class performing the multi-valley ensemble Monte Carlo simulation for
 * devices.
 *
 * @param param the parameter for the simulation.
 * @param device the simulated device structure.
 * @param solver the poisson-solver (for the calculation of the
 * potential within the device).
 * @param pmScheme the particle-mesh (PM) scheme which relates the continuous
 * particle position and the discrete representation of the device.
 * @param particleHandler the particle handler (stores and handles the
 * simulated particles of all particle types)
 * @param results class that stores the simulation results.
 */
template <class T, class DeviceType, class PoissonSolver, class ParticleHandler,
          class PMScheme>
class emcSimulation {
  static const SizeType Dim = DeviceType::Dimension;

  static_assert(std::is_same<DeviceType, emcDevice<T, Dim>>::value,
                "DeviceType in emcSimulation is not of the right type.");

  static_assert(
      std::is_base_of<emcAbstractSolver<T, DeviceType, ParticleHandler>,
                      PoissonSolver>::value,
      "PoissonSolver in emcSimulation is not of the right type.");

  static_assert(
      std::is_base_of<emcAbstractPMScheme<T, DeviceType>, PMScheme>::value,
      "PMScheme in emcSimulation is not of the right type.");

  static_assert(
      std::is_base_of<emcAbstractParticleHandler<T, DeviceType, PMScheme, Dim>,
                      ParticleHandler>::value,
      "ParticleHandler in emcSimulation is not of the right type.");

  emcSimulationParameter<T, DeviceType> &param;
  const DeviceType &device;
  PoissonSolver &solver;
  PMScheme pmScheme;
  ParticleHandler particleHandler;
  emcSimulationResults<T, DeviceType> results;

public:
  emcSimulation() = delete;

  emcSimulation(emcSimulationParameter<T, DeviceType> &inParam,
                DeviceType &inDevice, PoissonSolver &inSolver,
                PMScheme inPMScheme)
      : param(inParam), device(inDevice), solver(inSolver),
        pmScheme(inPMScheme),
        particleHandler(device, pmScheme, param.particleTypes,
                        param.nrCarriersPerPart, param.seedRNG),
        results(device, param) {
    param.check();
    param.print();
    checkDopingProfile();
  }

  /**
   * @brief executes the EMC simulation.
   *
   * First equilibrium conditions of the device are calculated
   * and the particles are initialized. Next, the particles are
   * moved step by step. Until the simulation time is over and
   * the final results are calculated and written to files.
   */
  void execute() {
    bool resetBC = true;
    auto totalSteps = param.getNrSteps();

    std::cout << "Equilibrium Characteristics ..." << std::endl;
    calcEquilibriumCharacteristics();
    writeCurrentResultsToFiles("Eq");

    std::cout << "Monte Carlo Procedure ..." << std::endl;
    for (SizeType nrStep = 0; nrStep < totalSteps; nrStep++) {
      performEMCStep(resetBC, param.isTransientStep(nrStep));
      if (nrStep % param.nrStepsBetweenShowProgress == 0) {
        resetBC = false;
        showProgress(nrStep);
      }
      if (nrStep >= totalSteps - param.nrStepsForFinalAvg)
        results.updateAverageCharacteristics();
    }
    showProgress(totalSteps);
    writeFinalResults();
  }

private:
  /**
   * @brief calculates the equilibrium conditions in the device and
   * initializes simulated particles based on them.
   *
   * First, equilibrium potential and el. field are calculated. Next,
   * particles are generated. Finally, particles are assigned to grid
   * points to calculate equilibrium particle concentration at
   * each grid point.
   */
  void calcEquilibriumCharacteristics() {
    solver.calcEquilibriumPotential(results.currPot, device);
    pmScheme.calcEField(results.eField, results.currPot, device);
    particleHandler.generateInitialParticles(results.currPot);
    assignParticlesToGrid();
    results.updateCurrentParticleConcentrations(device);
    particleHandler.printNrParticles();
  }

  /**
   * @brief performs one step of the EMC workflow.
   *
   * Distinguishes the case when particle-particle interactions
   * are included and the case where they are not included. These two
   * cases differ in the way the force at each particle is calculated.
   *
   * In general this step calculates the electric field in the device, then
   * the particles are moved based on the calculated electric field and the
   * equilibrium condition at Ohmic contacts is handled.
   *
   * @param resetBC boolean that determines if the boundary conditions for
   * the poisson equation should be re-evaluated.
   * @param isTransient boolean that determines if the current step is a
   * transient one.
   */
  void performEMCStep(bool resetBC = true, bool isTransient = false) {
    if (particleHandler.calcsPartPartInteraction()) {
      solver.calcBackgroundPotential(results.currPot, device, particleHandler,
                                     true);
      pmScheme.calcEField(results.eField, results.currPot, device);
      auto nrRemPart =
          particleHandler.driftScatterParticles(param.stepTime, results.eField);
      auto nrInjPart = particleHandler.handleOhmicContacts();
      if (!isTransient)
        results.updateCurrent(nrRemPart, nrInjPart);
    } else {
      solver.calcNonEquilibriumPotential(results.currPot, device,
                                         results.currConc[0], resetBC);
      pmScheme.calcEField(results.eField, results.currPot, device);
      auto nrRemPart =
          particleHandler.driftScatterParticles(param.stepTime, results.eField);
      auto nrInjPart = particleHandler.handleOhmicContacts();
      assignParticlesToGrid();
      results.updateCurrentParticleConcentrations(device);
      if (!isTransient)
        results.updateCurrent(nrRemPart, nrInjPart);
    }
  }

  /// calculate the carrier concentration from the simulated particles
  /// at each grid point.
  void assignParticlesToGrid() {
    for (SizeType idxType = 0; idxType < param.particleTypes.size();
         idxType++) {
      results.nrPart[idxType].fill(0);
      particleHandler.assignParticlesToMesh(idxType, results.nrPart[idxType]);
    }
  }

  /// outputs the current steps and writes current results to files
  void showProgress(SizeType nrStep) {
    std::cout << "\tNr. Iteration: \t\t" << nrStep << " / "
              << param.getNrSteps() << std::endl;
    // writeCurrentResultsToFiles(std::to_string(nrStep));
  }

  void writeCurrentResultsToFiles(std::string nameSuffix) {
    results.writeCurrentResults(nameSuffix, device);
    particleHandler.print(param.namePrefix, nameSuffix);
  }

  void writeFinalResults() {
    results.writeFinalResults(device);
    particleHandler.print(param.namePrefix, "Final");
  }

  /// @brief function checks if doping region is added to each discrete grid
  /// point in device. -
  void checkDopingProfile() const {
    auto idxDopingRegions = device.getDopingProfile().getDopingRegionIdx();
    std::array<SizeType, Dim> coord;
    for (coord.fill(0); !idxDopingRegions.isEndCoord(coord);
         idxDopingRegions.advanceCoord(coord)) {
      if (idxDopingRegions[coord] == -1)
        break;
    }
    if (!device.isEndCoord(coord)) {
      std::string strCoord;
      std::for_each(coord.begin(), coord.end(),
                    [this, &strCoord](SizeType val) {
                      strCoord += std::to_string(val) + " ";
                    });
      emcMessage::getInstance()
          .addError(
              "Doping Regions must be added to every discrete grid point in "
              "the simulated device. Found missing doping region at coordinate "
              "(" +
              strCoord + ").")
          .print();
    }
  }
};

#endif // EMC_SIMULATION_HPP