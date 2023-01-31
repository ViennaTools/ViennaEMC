#ifndef ACOUSTIC_SINGLE_LAYER_SCATTER_MECHANISM_HPP
#define ACOUSTIC_SINGLE_LAYER_SCATTER_MECHANISM_HPP

#include <math.h>

#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Intravalley Acoustic phonon scattering for single-layer of
 * Material.
 *
 * For single layer it is assumed that there is no movement in z-direction.
 * Formula for scattering from Paper of Kaasbjerg et al. with given link:
 * (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.115317)
 *
 * @param scatterConst pre-computed constant needed for calculation of
 * scatter rate
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcAcousticSingleLayerMechanism : public emcScatterMechanism<T> {
private:
  T scatterConst;
  mutable std::uniform_real_distribution<T> dist;
  std::string nameSuffix;

public:
  emcAcousticSingleLayerMechanism() = delete;

  /*! \brief Constructor.
   * @param inValley index of initial valley (valley to which scatter mechanism
   * is assigned)
   * @param sigma deformation potential [eV]
   * @param densityMaterial density of material [1 / m^2]
   * @param velSound sound velocity in material [m / s]
   * @param temperature temperature of material (device)
   * @param inNameSuffix suffix for filename
   */
  emcAcousticSingleLayerMechanism(SizeType inValley, T sigma, T densityMaterial,
                                  T velSound, T temperature,
                                  std::string inNameSuffix = "")
      : dist(0., 1.),
        nameSuffix(inNameSuffix), emcScatterMechanism<T>(inValley) {

    T cL = densityMaterial * std::pow(velSound, 2); // elastic constant
    scatterConst = std::pow(sigma * constants::q, 2) * constants::kB *
                   temperature / (cL * std::pow(constants::hbar, 3));
  }

  std::string getName() const { return "AcousticSL" + nameSuffix; }

  /// \brief scatter-rate calculation can be used for parabolic /
  /// nonparabolic and istrop / anistrop valleys.
  T getScatterRate(T energy, SizeType /*region*/) const {
    T md = this->ptrValley[this->idxValley]->getEffMassDOS();
    T alpha = this->ptrValley[this->idxValley]->getNonParabolicity();
    T rate = md * scatterConst * (1 + 2 * alpha * energy);
    return rate;
  }

  /// \brief Elastic, Isotrop Intravalley Scattering.
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    T angle = 2 * constants::pi * dist(rng);
    auto currValley = this->ptrValley[particle.valley];
    auto vogtFactor = currValley->getVogtTransformationFactor();

    // random direction weighted by vogtFactor (see Esseni p. 276)
    std::array<T, 3> k{std::cos(angle) / vogtFactor[0],
                       std::sin(angle) / vogtFactor[1], 0};
    T factor = 1. / (std::sqrt(k[0] * k[0] + k[1] * k[1]));
    k[0] *= factor;
    k[1] *= factor;

    particle.k =
        currValley->transformToDeviceCoord(particle.subValley, particle.k);
    T normK = currValley->getNormWaveVec(particle.energy);
    k[0] *= normK;
    k[1] *= normK;
    particle.k = k;
  }
};

#endif // ACOUSTIC_SINGLE_LAYER_SCATTER_MECHANISM_HPP