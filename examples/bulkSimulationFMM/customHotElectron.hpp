#include <ParticleType/emcParticleType.hpp>
#include <emcConstants.hpp>
#include <emcParticleInitialization.hpp>

/** Ensemble of hot and cold electrons.
 *
 * A portion (initProb) of the particles is initialized as hot
 * with the given energy (initEnergy).
 *
 * @param initEnergy initial energy of hot electrons.
 * @param initProb probability of electron being hot (number between 0 and 1).
 */
template <class T, class DeviceType>
struct customHotElectron : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::ValueVec ValueVec;
  typedef typename DeviceType::SizeVec SizeVec;

  static const SizeType Dim = DeviceType::Dimension;

  std::uniform_real_distribution<T> dist{0, 1};

  T hotInitEnergy;
  T hotProb;

  customHotElectron(T inProbability, T inEnergy)
      : emcParticleType<T, DeviceType>(1000, 1.), hotProb(inProbability),
        hotInitEnergy(inEnergy) {
    std::cout << "\tHot Electron probability =\t" << hotProb << "\n";
    std::cout << "\tHot Electron energy =\t" << hotInitEnergy << "\n";
  }

  std::string getName() const { return "HotElectrons"; }

  T getMass() const { return constants::me; }

  T getCharge() const { return -constants::q; }

  bool isMoved() const { return true; }

  bool isInjected() const { return false; }

  // init electrons based on doping (assume all dopants released electron)
  T getInitialNrParticles(const SizeVec &coord, const DeviceType &device,
                          const emcGrid<T, Dim> & /*potential*/) {
    T eDensity = device.getDopingProfile().getDoping(coord);
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 ||
          coord[idxDim] == device.getGridExtent()[idxDim] - 1)
        eDensity *= 0.5;
    }
    return eDensity * device.getCellVolume();
  }

  emcParticle<T> generateInitialParticle(const SizeVec &coord,
                                         const DeviceType &device,
                                         emcRNG &rng) {
    emcParticle<T> part;
    part.valley = std::floor(this->getNrValleys() * dist(rng));
    auto valley = this->getValley(part.valley);
    part.subValley = std::floor(valley->getDegeneracyFactor() * dist(rng));
    part.region = 0;

    part.energy =
        -1.5 * device.getThermalVoltage() * std::log(this->uniDistLog(rng));
    if (dist(rng) < hotProb)
      part.energy = hotInitEnergy;

    part.k = initRandomDirection(valley->getNormWaveVec(part.energy), dist(rng),
                                 dist(rng));
    part.tau = this->getNewTau(part.valley, part.subValley, rng);
    part.grainTau = this->getNewGrainTau(rng);
    return part;
  }
};
