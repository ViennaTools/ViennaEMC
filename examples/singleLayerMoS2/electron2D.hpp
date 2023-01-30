#include <ParticleType/emcParticleType.hpp>
#include <emcConstants.hpp>

/**
 * @brief Particle Type that describes electrons in 2-dimensional material.
 *
 * Nr of initialized particles near grid point is fixed. All particles are
 * initialized in random subvalley of first added valley (K-valley). Energy is
 * assigned randomly around mean (2D) energy of boltzmann distribution and
 * wave-vector is assigned randomly in (x,y)-direction, z-component is assumed
 * to be 0.
 */
template <class T, class DeviceType>
struct electron2D : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::ValueVec ValueVec;
  typedef typename DeviceType::SizeVec SizeVec;

  static const SizeType Dim = DeviceType::Dimension;

  std::uniform_real_distribution<T> dist;
  std::uniform_real_distribution<T> distForLog;

  electron2D()
      : dist(0., 1.),
        distForLog(1e-6, 1.), emcParticleType<T, DeviceType>(5000, 0.5) {}

  std::string getName() const { return "Electrons"; }

  T getMass() const { return constants::me; }

  T getCharge() const { return -constants::q; }

  bool isMoved() const { return true; }

  // no contacts are used (injection not necessary)
  bool isInjected() const { return false; }

  // fix number of particles near grid point
  T getInitialNrParticles(const SizeVec &coord, const DeviceType & /*device*/,
                          const emcGrid<T, Dim> & /*potential*/) {
    return 4;
  }

  emcParticle<T> generateInitialParticle(const SizeVec & /*coord*/,
                                         const DeviceType &device,
                                         emcRNG &rng) {
    emcParticle<T> part;

    // initialize all electrons in random subvalley of K-valley
    part.valley = 0;
    auto &valley = this->valleys[0];
    part.subValley = std::floor(valley->getDegeneracyFactor() * dist(rng));

    // for bulk simulation only one region is used
    part.region = 0;

    part.energy = -device.getThermalVoltage() * std::log(distForLog(rng));

    // randomly assign wave-vector direction in x- and y-direction
    auto normK = valley->getNormWaveVec(part.energy);
    T angle = 2 * constants::pi * dist(rng);
    part.k[0] = normK * std::cos(angle);
    part.k[1] = normK * std::sin(angle);
    part.k[2] = 0;

    part.tau = this->getNewTau(part.valley, part.region, rng);
    part.grainTau = this->getNewGrainTau(rng);
    return part;
  }
};