#ifndef EMC_SIMULATION_PARAMETER_HPP
#define EMC_SIMULATION_PARAMETER_HPP

#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <string>

#include <ParticleType/emcParticleType.hpp>
#include <emcMessage.hpp>
#include <emcUtil.hpp>

/**
 * @brief Class that stores all parameters for a simulation with
 * emcSimulation.
 *
 * @param simTime total simulation time [s]
 * @param stepTime time of one EMC step (dt) [s]
 * @param transientTime time until steady state is reached [s]
 * @param nrCarriersPerPart number of carriers that are represented per
 * simulated particle (if > 1 superparticles are used, testing required)
 * @param namePrefix prefix for the names of the output files
 * @param adaptPotentialForWrite function pointer that is called when the
 * potential is written to a file, allows adaptation of potential for write
 * process
 * @param nrStepsBetweenShowProgress number of steps before the progress is
 * shown again
 * @param nrStepsForFinalAvg number of steps that are used at the end to average
 * the final potential and final particle concentration
 * @param particleTypes indexed particle types that should be used in the
 * simulation
 * @param seedRNG seed for the random number generator(s), if the seed is set,
 * the simulation results can be reproduced by setting the same seed again
 */
template <class T, class DeviceType> class emcSimulationParameter {
  typedef emcParticleType<T, DeviceType> ParticleType;
  typedef std::map<SizeType, std::unique_ptr<ParticleType>>
      MapIdxToParticleTypes;

  static const SizeType Dim = DeviceType::Dimension;

  T simTime, stepTime, transientTime;
  SizeType nrCarriersPerPart{1}, nrStepsBetweenShowProgress{100},
      nrStepsForFinalAvg{100};
  std::string namePrefix = "";
  T (*adaptPotentialForWrite)(const T &, const DeviceType &) = nullptr;
  MapIdxToParticleTypes particleTypes;
  SizeType seedRNG =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();

public:
  emcSimulationParameter()
      : emcSimulationParameter(10e-12, 1.5e-16, 3e-12, 1) {}

  /**
   * @brief Construct a new emcSimulationParameter object.
   *
   * @param inSimTime total simulation time [s]
   * @param inStepTime time of one step of the EMC procedure [s]
   * @param inTransientTime time until steady state is reached [s]
   */
  emcSimulationParameter(T inSimTime, T inStepTime, T inTransientTime)
      : emcSimulationParameter(inSimTime, inStepTime, inTransientTime, 1) {}

  /**
   * @brief Construct a new emcSimulationParameter object.
   *
   * @param inSimTime total simulation time [s]
   * @param inStepTime time of one step of the EMC procedure [s]
   * @param inTransientTime time until steady state is reached [s]
   * @param inNrCarriersPerPart number of carriers that are represented per
   * simulated particle (if > 1 superparticles are used, testing required)
   */
  emcSimulationParameter(T inSimTime, T inStepTime, T inTransientTime,
                         SizeType inNrCarriersPerPart)
      : simTime(inSimTime), stepTime(inStepTime),
        transientTime(inTransientTime), nrCarriersPerPart(inNrCarriersPerPart),
        namePrefix(""), adaptPotentialForWrite(nullptr),
        nrStepsBetweenShowProgress(100), nrStepsForFinalAvg(100),
        particleTypes() {
    checkTimes();
    checkNrCarriersPerParticles();
  }

  void setTimes(T inSimTime, T inStepTime, T inTransientTime) {
    simTime = inSimTime;
    stepTime = inStepTime;
    transientTime = inTransientTime;
    checkTimes();
  }

  void setSeed(SizeType inSeed) { seedRNG = inSeed; }

  void setSimTime(T inSimTime) { simTime = inSimTime; }

  void setStepTime(T inStepTime) { stepTime = inStepTime; }

  void setTransientTime(T inTransientTime) { transientTime = inTransientTime; }

  void setNrCarriersPerPart(SizeType inNrCarriersPerParticle) {
    nrCarriersPerPart = inNrCarriersPerParticle;
    checkNrCarriersPerParticles();
  }

  void setNamePrefix(std::string inNamePrefix) { namePrefix = inNamePrefix; }

  void setAdaptPotentialForWriteFunction(
      T (*inAdaptPotentialForWrite)(const T &, const DeviceType &)) {
    adaptPotentialForWrite = inAdaptPotentialForWrite;
  }

  /// sets number of steps that are performed between two progress updates
  void setNrStepsBetweenShowProgress(SizeType nrSteps) {
    nrStepsBetweenShowProgress = nrSteps;
  }

  /// set nr of steps at the end that are used for averaging the potential and
  /// the electron concentration
  void setNrStepsForFinalAvg(SizeType nrSteps) { nrStepsForFinalAvg = nrSteps; }

  /// add a particle type that should be simulated
  template <class DerivedParticleType>
  typename std::enable_if<
      std::is_base_of<ParticleType, DerivedParticleType>::value>::type
  addParticleType(std::unique_ptr<DerivedParticleType> &&newParticleType) {
    newParticleType->check();
    particleTypes[particleTypes.size()] = std::move(newParticleType);
  }

  SizeType getNrParticleTypes() const { return particleTypes.size(); }

  SizeType getNrSteps() const { return std::ceil(simTime / stepTime); }

  SizeType getNrNonTransientSteps() const {
    return getNrSteps() - getNrTransientSteps();
  }
  SizeType getNrTransientSteps() const {
    return std::ceil(transientTime / stepTime);
  }

  bool isTransientStep(SizeType nrStep) const {
    return nrStep < getNrTransientSteps();
  }

  void print() const {
    std::cout << "Simulation Parameter ...\n";
    std::cout << "\tTotal Simulation Time:\t" << simTime << " s\n";
    std::cout << "\tStep Time:\t\t" << stepTime << " s\n";
    std::cout << "\tTransient Time:\t\t" << transientTime << " s\n";
    std::cout << "\tNr. of Steps:\t\t" << getNrSteps() << "\n";
    std::cout << "\tNr. of Steps for Avg:\t" << getNrNonTransientSteps()
              << "\n";
  }

  void check() const {
    checkTimes();
    checkNrSteps();

    if (particleTypes.size() == 0) {
      emcMessage::getInstance()
          .addError("Add at least one particleType to simulation parameter!")
          .print();
    }
    SizeType countMovingTypes = 0;
    for (const auto &type : particleTypes) {
      if (type.second->isMoved())
        countMovingTypes++;
    }
    if (countMovingTypes == 0) {
      emcMessage::getInstance()
          .addError(
              "Add at least one moving particleType to simulation parameter!")
          .print();
    }
  }

private:
  void checkTimes() const {
    if (stepTime < 0 || simTime < 0 || transientTime < 0) {
      emcMessage::getInstance()
          .addError("The time parameter can't be negative!")
          .print();
    }
    if (stepTime == 0) {
      emcMessage::getInstance().addError("Step Time can't be zero!").print();
    }
    if (stepTime > simTime || transientTime > simTime) {
      emcMessage::getInstance()
          .addError("stepTime, avgTime and transientTime have to be smaller "
                    "than simTime!")
          .print();
    }
  }

  void checkNrCarriersPerParticles() const {
    if (nrCarriersPerPart < 1) {
      emcMessage::getInstance()
          .addError(
              "Nr. of Carriers per simulated particle has to be at least 1!")
          .print();
    }
  }

  void checkNrSteps() const {
    if (nrStepsForFinalAvg > getNrSteps()) {
      emcMessage::getInstance()
          .addError("nrStepsForFinalAvg has to be smaller than total step nr.!")
          .print();
    }
  }

  template <class, class> friend class emcSimulationResults;
  template <class, class, class, class, class> friend class emcSimulation;
};

#endif // EMC_SIMULATION_PARAMETER_HPP