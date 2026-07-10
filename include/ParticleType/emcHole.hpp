#ifndef EMC_HOLE_HPP
#define EMC_HOLE_HPP

#include <ParticleType/emcParticleType.hpp>
#include <emcConstants.hpp>
#include <emcParticleInitialization.hpp>
#include <emcUtil.hpp>

/**
 * @brief Mobile hole particle type initialised in thermal equilibrium.
 *
 * Symmetric counterpart to emcElectron: charge = +q, moved = true.
 * Valence band structure is supplied via addValley() using any of the
 * existing valley types with the appropriate hole effective mass.
 *
 * Initialisation density:
 *  - usePotentialForInit = true  (device simulations):
 *      p = ni * exp(-phi/Vt)  — inverse of electron response to potential
 *  - usePotentialForInit = false (bulk / photo-excited simulations):
 *      p = |doping|  — use |getDoping()| so both positive (photo-generated)
 *                       and negative (acceptor) doping values give the
 *                       correct hole count
 *
 * @tparam T          numeric type (float or double)
 * @tparam DeviceType device type
 */
template <class T, class DeviceType>
struct emcHole : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::ValueVec ValueVec;
  typedef typename DeviceType::SizeVec SizeVec;

  static const SizeType Dim = DeviceType::Dimension;

  std::uniform_real_distribution<T> dist{1e-6, 1.};
  bool usePotentialForInit;
  T initEnergyEV; // > 0: fixed photoexcitation energy; 0: Maxwellian at device T

  emcHole(SizeType inHandlerNrEnergyLevels = 1000, T inHandlerMaxEnergy = 4.,
          bool inUsePotentialForInit = false, T inInitEnergyEV = T(0))
      : emcParticleType<T, DeviceType>(inHandlerNrEnergyLevels,
                                       inHandlerMaxEnergy),
        usePotentialForInit(inUsePotentialForInit),
        initEnergyEV(inInitEnergyEV) {}

  std::string getName() const { return "Holes"; }

  T getMass() const { return constants::me; }

  T getCharge() const { return +constants::q; }

  bool isMoved() const { return true; }

  bool isInjected() const { return true; }

  T getInitialNrParticles(const SizeVec &coord, const DeviceType &device,
                          const emcGrid<T, Dim> &potential) {
    T hDensity;
    if (usePotentialForInit)
      // potential is normalised by Vt; holes are anti-correlated with electrons
      hDensity = std::exp(-potential[coord]) * device.getMaterial().getNi();
    else
      hDensity = std::fabs(device.getDopingProfile().getDoping(coord));

    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 ||
          coord[idxDim] == potential.getSize(idxDim) - 1)
        hDensity *= T(0.5);
    }
    return hDensity * device.getCellVolume();
  }

  T getExpectedNrParticlesAtContact(const SizeVec &coord,
                                    const DeviceType &device) {
    T expectedNrPart = device.getCellVolume() *
                       std::fabs(device.getDopingProfile().getDoping(coord));
    const auto &extent = device.getGridExtent();
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 || coord[idxDim] == extent[idxDim] - 1)
        expectedNrPart *= T(0.5);
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
    if (initEnergyEV > T(0))
      initParticleKSpaceFixed(part, initEnergyEV, coord, device, valley, rng);
    else
      initParticleKSpaceMaxwellian(part, coord, device, valley, rng);
    part.tau = this->getNewTau(part.valley, part.region, rng);
    part.grainTau = this->getNewGrainTau(rng);
    return part;
  }

  emcParticle<T> generateInjectedParticle(const SizeVec &coord,
                                          const DeviceType &device,
                                          emcRNG &rng) {
    return generateInitialParticle(coord, device, rng);
  }
};

#endif // EMC_HOLE_HPP
