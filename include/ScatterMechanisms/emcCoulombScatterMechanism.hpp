#ifndef EMC_COULOMB_SCATTER_MECHANISM_HPP
#define EMC_COULOMB_SCATTER_MECHANISM_HPP

#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>

/*! \brief Coulomb Scattering (Brooks-Herring Approach).
 *
 * Assumes that scattering is elastic and anisotrop.
 */
template <class T, class DeviceType>
class emcCoulombScatterMechanism : public emcScatterMechanism<T> {

private:
  DeviceType &device;
  T scatterConst1;
  T scatterConst2;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcCoulombScatterMechanism() = delete;

  emcCoulombScatterMechanism(SizeType inIdxValley, T epsR, DeviceType &inDevice)
      : emcScatterMechanism<T>(inIdxValley), dist(0., 1.), device(inDevice) {
    T epsMat = constants::eps0 * epsR;
    T temperature = device.getTemperature();
    T Vt = device.getThermalVoltage();
    scatterConst1 = std::sqrt(2 * constants::q) *
                    std::pow(constants::kB * temperature, 2) /
                    (constants::pi * pow(constants::hbar, 4));
    scatterConst2 = 8 * epsMat * Vt / (constants::hbar * constants::hbar);
  }

  std::string getName() const { return "Coulomb"; }

  T getScatterRate(T energy, SizeType idxRegion) const {
    T regionDoping = std::fabs(device.getDopingProfile().getDoping(idxRegion));
    auto currValley = this->ptrValley[this->idxValley];
    T md = currValley->getEffMassDOS();
    T mc = currValley->getEffMassCond();
    T alpha = currValley->getNonParabolicity();
    T gamma = currValley->getGamma(energy);
    return scatterConst1 * pow(md, 3. / 2.) / regionDoping * std::sqrt(gamma) *
           (2 * alpha * energy + 1.0) /
           (1 + (scatterConst2 * mc * gamma / regionDoping));
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto &valley = this->ptrValley[this->idxValley];
    T gamma = valley->getGamma(particle.energy);
    T regionDoping =
        std::fabs(device.getDopingProfile().getDoping(particle.region));
    T mc = valley->getEffMassCond();
    T rand = dist(rng);
    T debyeEnergy = regionDoping / (scatterConst2 * mc);
    T cosTheta = 1.0 - rand * 2.0 / ((1 - rand) * gamma / debyeEnergy + 1.0);
    particle.k = initRandomDirectionWithRespectToCurrentK(particle.k, cosTheta,
                                                          dist(rng));
  }
};

#endif