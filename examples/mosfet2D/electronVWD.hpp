#ifndef ELECTRON_VWD_HPP
#define ELECTRON_VWD_HPP

#include <math.h>

#include <ParticleType/emcParticleType.hpp>
#include <emcConstants.hpp>
#include <emcParticleDrift.hpp>
#include <emcParticleInitialization.hpp>
#include <emcUtil.hpp>

/**
 * @brief Electrons that behave as implemented in ViennaWD.
 */
template <class T, class DeviceType>
struct electronVWD : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::ValueVec ValueVec;
  typedef typename DeviceType::SizeVec SizeVec;

  static const SizeType Dim = DeviceType::Dimension;

  std::uniform_real_distribution<T> dist{0, 1};

  electronVWD(SizeType inHandlerNrEnergyLevels = 1000,
              T inHandlerMaxEnergy = 4.)
      : emcParticleType<T, DeviceType>(inHandlerNrEnergyLevels,
                                       inHandlerMaxEnergy){};

  std::string getName() const { return "Electrons"; }

  T getMass() const { return constants::me; }

  T getCharge() const { return -constants::q; }

  bool isMoved() const { return true; }

  bool isInjected() const { return true; }

  //! Diff to emcElectron: nr of initial particles always based on potential.
  //! Diff to emcElectron: round number of expected electrons.
  T getInitialNrParticles(const SizeVec &coord, const DeviceType &device,
                          const emcGrid<T, Dim> &potential) {
    T eDensity = std::exp(potential[coord]) * device.getMaterial().getNi();
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 || coord[idxDim] == potential.getSize(idxDim) - 1)
        eDensity *= 0.5;
    }
    return std::round(eDensity * device.getCellVolume());
  }

  //! Diff to emcElectron: round number of expected electrons.
  T getExpectedNrParticlesAtContact(const SizeVec &coord,
                                    const DeviceType &device) {
    T doping = device.getDopingProfile().getDoping(coord);
    T expectedNrPart = device.getCellVolume() * doping;
    const auto &extent = device.getGridExtent();
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 || coord[idxDim] == extent[idxDim] - 1)
        expectedNrPart *= 0.5;
    }
    return std::round(expectedNrPart);
  }

  emcParticle<T> generateInitialParticle(const SizeVec &coord,
                                         const DeviceType &device,
                                         emcRNG &rng) {
    emcParticle<T> part;
    part.region = device.getDopingProfile().getDopingRegionIdx(coord);
    part.valley = std::floor(this->getNrValleys() * dist(rng));
    auto valley = this->getValley(part.valley);
    part.subValley = std::floor(valley->getDegeneracyFactor() * dist(rng));
    initParticleKSpaceMaxwellian(part, coord, device, valley, rng);
    part.tau = this->getNewTau(part.valley, part.valley, rng);
    part.grainTau = this->getNewGrainTau(rng);
    return part;
  }

  emcParticle<T> generateInjectedParticle(const SizeVec &coord,
                                          const DeviceType &device,
                                          emcRNG &rng) {
    emcParticle<T> part;
    part.region = device.getDopingProfile().getDopingRegionIdx(coord);
    part.valley = std::floor(this->getNrValleys() * dist(rng));
    auto valley = this->getValley(part.valley);
    part.subValley = std::floor(valley->getDegeneracyFactor() * dist(rng));
    initParticleKSpaceMaxwellian(part, coord, device, valley, rng);
    part.tau = this->getNewTau(part.valley, part.valley, rng);
    part.grainTau = this->getNewGrainTau(rng);
    return part;
  }
};

#endif // ELECTRON_VWD_HPP