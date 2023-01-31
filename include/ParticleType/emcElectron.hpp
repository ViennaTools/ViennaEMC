#ifndef EMC_ELECTRON_HPP
#define EMC_ELECTRON_HPP

#include <ParticleType/emcParticleType.hpp>
#include <emcConstants.hpp>
#include <emcParticleInitialization.hpp>
#include <emcUtil.hpp>

/**
 * @brief Electron Particle Type with random characteristics in thermal
 * equilibrium in beginning.
 *
 * @tparam T Numeric Type
 * @tparam DeviceType Device Type
 * @param usePotentialForInit boolean that determines how the number of initial
 * particles is calculated (if true potential is used, else the amount of doping
 * is used for this)
 */
template <class T, class DeviceType>
struct emcElectron : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::ValueVec ValueVec;
  typedef typename DeviceType::SizeVec SizeVec;

  static const SizeType Dim = DeviceType::Dimension;

  std::uniform_real_distribution<T> dist{1e-6, 1.};
  bool usePotentialForInit;

  emcElectron(SizeType inHandlerNrEnergyLevels = 1000,
              T inHandlerMaxEnergy = 4., bool inUsePotentialForInit = true)
      : emcParticleType<T, DeviceType>(inHandlerNrEnergyLevels,
                                       inHandlerMaxEnergy),
        usePotentialForInit(inUsePotentialForInit){};

  std::string getName() const { return "Electrons"; }

  T getMass() const { return constants::me; }

  T getCharge() const { return -constants::q; }

  bool isMoved() const { return true; }

  bool isInjected() const { return true; }

  T getInitialNrParticles(const SizeVec &coord, const DeviceType &device,
                          const emcGrid<T, Dim> &potential) {
    T eDensity;
    if (usePotentialForInit) // init nr particles based on potential
      eDensity = std::exp(potential[coord]) * device.getMaterial().getNi();
    else // init nr particles based on doping
      eDensity = device.getDopingProfile().getDoping(coord);

    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 || coord[idxDim] == potential.getSize(idxDim) - 1)
        eDensity *= 0.5;
    }
    return eDensity * device.getCellVolume();
  }

  T getExpectedNrParticlesAtContact(const SizeVec &coord,
                                    const DeviceType &device) {
    T expectedNrPart =
        device.getCellVolume() * device.getDopingProfile().getDoping(coord);
    const auto &extent = device.getGridExtent();
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 || coord[idxDim] == extent[idxDim] - 1)
        expectedNrPart *= 0.5;
    }
    return expectedNrPart;
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
    part.tau = this->getNewTau(part.valley, part.region, rng);
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
    part.tau = this->getNewTau(part.valley, part.region, rng);
    part.grainTau = this->getNewGrainTau(rng);
    return part;
  }
};

#endif // EMC_ELECTRON_HPP