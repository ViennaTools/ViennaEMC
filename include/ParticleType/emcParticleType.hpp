#ifndef EMC_PARTICLE_TYPE_HPP
#define EMC_PARTICLE_TYPE_HPP

#include <cmath>
#include <memory>
#include <random>
#include <string>

#include <SurfaceScatterMechanisms/emcSurfaceScatterMechanism.hpp>
#include <ValleyTypes/emcAbstractValley.hpp>
#include <emcBoundaryPos.hpp>
#include <emcGrainScatterMechanism.hpp>
#include <emcGrid.hpp>
#include <emcMessage.hpp>
#include <emcParticle.hpp>
#include <emcScatterHandler.hpp>
#include <emcUtil.hpp>

/// abstract class that holds characteristics and behaviour
/// of specific particleType
/// @param valleys vector of valleys that describe the dispersion relation of a
/// particle (only needed if particle is moving)
/// @param scatterHandler vector of scatterHandler ( each handler handles
/// scattering events of one valley) (only needed if particle is moving)
template <class T, class DeviceType> struct emcParticleType {
  static const SizeType Dim = DeviceType::Dimension;

  typedef typename DeviceType::ValueVec ValueVec;
  typedef typename DeviceType::SizeVec SizeVec;
  typedef emcAbstractValley<T> AbstractValley;
  typedef emcScatterMechanism<T> AbstractScatterMechanism;
  typedef emcSurfaceScatterMechanism<T, DeviceType>
      AbstractSurfaceScatterMechanism;

  std::vector<std::unique_ptr<AbstractValley>> valleys;
  emcScatterHandler<T, DeviceType> scatterHandler;
  mutable std::uniform_real_distribution<T> uniDistLog{1e-6, 1.};

  emcParticleType(SizeType inHandlerNrEnergyLevels = 1000,
                  T inHandlerMaxEnergy = 4.)
      : valleys(), scatterHandler(emcScatterHandler<T, DeviceType>(
                       inHandlerNrEnergyLevels, inHandlerMaxEnergy)) {}

  virtual ~emcParticleType() = default;

  /// returns name of particle type, e.g. "electrons"
  virtual std::string getName() const = 0;

  /// returns mass of particle type
  virtual T getMass() const {
    addError("getMass", "isMoved");
    return 0;
  }

  /// returns charge of particle type
  virtual T getCharge() const = 0;

  /// tells if particle is moved during the simulation
  virtual bool isMoved() const = 0;

  /// tells if particles at contact should be injected during simulation
  virtual bool isInjected() const = 0;

  /// calculates the initial nr of particles near each coordinate.
  virtual T getInitialNrParticles(const SizeVec &coord,
                                  const DeviceType &device,
                                  const emcGrid<T, Dim> &potential) = 0;

  /// generates a particle with initial random characteristics.
  /// is only called if isMoved() returns true
  virtual emcParticle<T> generateInitialParticle(const SizeVec &coord,
                                                 const DeviceType &device,
                                                 emcRNG &rng) {
    addError("generateIntialParticle", "isMoved");
    return emcParticle<T>();
  }

  /// calculates the number of particles that should be near an ohmic contact.
  /// is only called if isInjected() returns true
  virtual T getExpectedNrParticlesAtContact(const SizeVec &coord,
                                            const DeviceType &device) {
    addError("getExpectedNrParticlesAtContact", "isInjected");
    return 0;
  }

  /// generates an injected particle with random characteristics.
  /// is only called if isInjected() returns true
  virtual emcParticle<T> generateInjectedParticle(const SizeVec &coord,
                                                  const DeviceType &device,
                                                  emcRNG &rng) {
    addError("generateInjectedParticle", "isInjected");
    return emcParticle<T>();
  }

  /// scatters the particle (used during free-flight-scatter-function)
  /// is only called if isMoved() returns true
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    scatterHandler.scatterParticle(particle, rng);
  }

  /// scatters the particle at boundary
  /// is only called if isMoved() returns true
  void scatterParticleAtBoundary(emcBoundaryPos boundary,
                                 emcParticle<T> &particle,
                                 const DeviceType &device,
                                 std::array<T, Dim> &pos, emcRNG &rng) const {
    scatterHandler.scatterParticleAtBoundary(boundary, particle, device, pos,
                                             rng);
  }

  /// scatters the particle at grain
  /// is only called if isMoved() returns true
  void scatterParticleAtGrain(emcParticle<T> &particle, emcRNG &rng) const {
    scatterHandler.scatterParticleAtGrain(particle, rng);
  }

  SizeType getNrValleys() const { return valleys.size(); }

  auto getValley(SizeType idxValley) const {
    checkIdxValley(idxValley);
    return valleys[idxValley].get();
  }

  /// add a valley to particleType
  template <class DerivedValley>
  typename std::enable_if<
      std::is_base_of<AbstractValley, DerivedValley>::value>::type
  addValley(std::unique_ptr<DerivedValley> &&newValleyType) {
    newValleyType->check();
    // TODO check if particle mass is right
    valleys.push_back(std::move(newValleyType));
  }

  void setGrainScatterMechanism(
      std::unique_ptr<emcGrainScatterMechanism<T>> &&newMechanism) {
    scatterHandler.setGrainScatterMechanism(
        std::forward<std::unique_ptr<emcGrainScatterMechanism<T>> &&>(
            newMechanism));
  }

  template <class DerivedScatterMechanism>
  typename std::enable_if<std::is_base_of<AbstractScatterMechanism,
                                          DerivedScatterMechanism>::value>::type
  addScatterMechanism(const std::vector<int> &regions,
                      std::unique_ptr<DerivedScatterMechanism> &&newMechanism) {
    auto idxValley = newMechanism->getIdxValley();
    checkIdxValley(idxValley);
    newMechanism->setPtrValley(valleys);
    newMechanism->check();
    scatterHandler.addScatterMechanism(
        std::forward<std::unique_ptr<DerivedScatterMechanism> &&>(newMechanism),
        regions);
  }

  template <class DerivedSurfaceScatterMechanism>
  typename std::enable_if<
      std::is_base_of<AbstractSurfaceScatterMechanism,
                      DerivedSurfaceScatterMechanism>::value>::type
  setSurfaceScatterMechanism(
      emcBoundaryPos boundaryPosition,
      std::unique_ptr<DerivedSurfaceScatterMechanism> &&newMechanism) {
    scatterHandler.setSurfaceScatterMechanism(
        std::forward<std::unique_ptr<DerivedSurfaceScatterMechanism> &&>(
            newMechanism),
        boundaryPosition);
  }

  void initScatterTables() { scatterHandler.initScatterTables(); }

  T getGrainTau() const { return scatterHandler.getGrainTau(); }

  T getTau(SizeType idxValley, SizeType idxRegion) const {
    return scatterHandler.getTau(idxRegion, idxValley);
  }

  /**
   * @brief Helper function that returns new random remaining free-flight time
   * (tau) for a given particle
   * @param idxValley index of the current valley of the particle
   * @param idxRegion index of the current region of the partice
   * @param rng random number generator
   */
  T getNewTau(SizeType idxValley, SizeType idxRegion, emcRNG &rng) const {
    return -std::log(uniDistLog(rng)) * getTau(idxValley, idxRegion);
  }

  /**
   * @brief Helper function that returns new random remaining free-flight time
   * (tau) until the next grain scattering event for a given particle
   * @param idxValley index of the current valley of the particle
   * @param idxRegion index of the current region of the partice
   * @param rng random number generator
   */
  T getNewGrainTau(emcRNG &rng) const {
    return -std::log(uniDistLog(rng)) * getGrainTau();
  }

  void check() const {
    if (isMoved()) {
      if (valleys.empty()) {
        emcMessage::getInstance()
            .addError("Moving Particle Type " + getName() +
                      " has to at least have one added valley.")
            .print();
      }
    }
  }

private:
  void checkIdxValley(SizeType idxValley) const {
    if (idxValley >= valleys.size()) {
      emcMessage::getInstance()
          .addError("Used index for Valley for" + getName() + "is invalid.")
          .print();
    }
  }

  void addError(std::string nameFunc, std::string dependentFunc) const {
    emcMessage::getInstance()
        .addError("Function " + nameFunc + "() is not implemented for " +
                  getName() + ". Either set the return value of " +
                  dependentFunc +
                  "() to false or implement required Function for given "
                  "ParticleType.")
        .print();
  }
};

#endif // EMC_PARTICLE_TYPE_HPP