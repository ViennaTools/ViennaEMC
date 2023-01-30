#ifndef ACOUSTIC_SCATTER_MECHANISM_HPP
#define ACOUSTIC_SCATTER_MECHANISM_HPP

#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>

/*! \brief Intravalley acoustic phonon scattering
 *
 * Assumes that scattering is elastic and isotrop.
 *
 * @param scatterConst constant that is used for the calculation of
 * the scatter rates
 */
template <class T>
class emcAcousticScatterMechanism : public emcScatterMechanism<T> {
private:
  T scatterConst;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcAcousticScatterMechanism() = delete;

  /*! \brief Constructor that gets all information about material from given
   * device.
   *
   * @param inIdxValley initial and final idx of the corresponding valley
   * @param sigma deformation potential (in eV)
   * @param device current used device
   */
  template <class DeviceType>
  emcAcousticScatterMechanism(SizeType inIdxValley, T sigma,
                              const DeviceType &device)
      : emcAcousticScatterMechanism(
            inIdxValley, sigma, device.getMaterial().getRho(),
            device.getMaterial().getVelSound(), device) {}

  /*! \brief Constructor that needs information about material as input
   * parameter.
   *
   * @param inIdxValley initial and final idx of the corresponding valley
   * @param sigma deformation potential (in eV)
   * @param materialDensity density of the used material (in kg / m³)
   * @param velSound velocity of sound in material [m / s]
   * @param device current used device
   */
  template <class DeviceType>
  emcAcousticScatterMechanism(SizeType inIdxValley, T sigma, T materialDensity,
                              T velSound, const DeviceType &device)
      : emcScatterMechanism<T>(inIdxValley), dist(0., 1.) {
    T cL = materialDensity * std::pow(velSound, 2); // elastic constant

    scatterConst = std::sqrt(2.0 * constants::q) *
                   std::pow(sigma * constants::q, 2) * constants::kB *
                   device.getTemperature() /
                   (constants::pi * cL * pow(constants::hbar, 4));
  }

  std::string getName() const { return "Acoustic"; }

  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    auto currValley = this->ptrValley[this->idxValley];
    T md = currValley->getEffMassDOS();
    T alpha = currValley->getNonParabolicity();
    T gamma = currValley->getGamma(energy);
    return scatterConst * pow(md, 3. / 2.) * std::sqrt(gamma) *
           (2 * alpha * energy + 1.0);
  }

  /// scatter particle elastically in random direction
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    particle.k = initRandomDirection(norm(particle.k), dist(rng), dist(rng));
  }
};

#endif // ACOUSTIC_SCATTER_MECHANISM_HPP