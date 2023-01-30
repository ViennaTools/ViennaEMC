#ifndef EMC_NONPARABOLIC_ISOTROP_SINGLELAYER_VALLEY_HPP
#define EMC_NONPARABOLIC_ISOTROP_SINGLELAYER_VALLEY_HPP

#include <math.h>

#include <ValleyTypes/emcAbstractValley.hpp>
#include <emcConstants.hpp>
#include <emcMessage.hpp>
#include <emcUtil.hpp>

/*! \brief Nonparabolic Isotropical Valley for Single Layer (assumes
 * 2D-extension in x and y-direction!)
 *
 * @param alpha non - parabolicity (in 1 / eV)
 * @param degFactor number of subvalleys in that valley
 * @param effMass effective mass (in each direction)
 * @param vogtFactor factor for herring-vogt transformation, (1,1,0) for isotrop
 * 2d valley
 * @param bottomEnergy energy of the bottom of the valley (in eV)
 */
template <class T>
class emcNonParabolicIsotropSingleLayerValley : public emcAbstractValley<T> {
private:
  T alpha;
  SizeType degFactor;
  T effMass;
  std::array<T, 3> vogtFactor;
  T bottomEnergy;

public:
  emcNonParabolicIsotropSingleLayerValley() = delete;

  /*! \brief Constructor.
   *
   * @param inRelEffMass relative effective mass of valley
   * @param inParticleMass mass of particles
   * @param inDegFactor nr of equivalent subvalleys
   * @param inAlpha non - parabolicity (in 1 / eV)
   * @param inBottomEnergy energy of the bottom of the valley (in eV)
   */
  emcNonParabolicIsotropSingleLayerValley(T inRelEffMass, T inParticleMass,
                                          SizeType inDegFactor, T inAlpha,
                                          T inBottomValleyEnergy = 0.)
      : degFactor(inDegFactor), effMass(inRelEffMass * inParticleMass),
        vogtFactor({1, 1, 0}), alpha(inAlpha),
        bottomEnergy(inBottomValleyEnergy) {}

  T getBottomEnergy() const { return bottomEnergy; }

  T getEffMassDOS(T energy = 0) const {
    return effMass * std::pow(1 + 2 * alpha * energy, 3.);
  }

  T getEffMassCond(T energy = 0) const {
    return effMass * (1 + 2 * energy * alpha);
  }

  T getNonParabolicity() const { return alpha; }

  SizeType getDegeneracyFactor() const { return degFactor; }

  T getNormWaveVec(T energy) const {
    return std::sqrt(2 * effMass * constants::q * getGamma(energy)) /
           constants::hbar;
  }

  T getEnergy(const std::array<T, 3> &k) const {
    T gamma = constants::hbar * constants::hbar * (k[0] * k[0] + k[1] * k[1]) /
              (effMass * constants::q);
    return gamma / (1 + std::sqrt(1 + 2 * alpha * gamma));
  }

  const std::array<T, 3> &getVogtTransformationFactor() const {
    return vogtFactor;
  }

  std::array<T, 3> getVelocity(const std::array<T, 3> &k, T energy,
                               SizeType /*idxSubValley*/) const {
    T factor = constants::hbar /
               (effMass * std::sqrt(1 + 4 * alpha * getGamma(energy)));
    return scale(k, factor);
  }

  T getGamma(T energy) const { return energy * (1 + alpha * energy); }

  std::array<T, 3> transformToEllipseCoord(SizeType /*idxSubValley*/,
                                           const std::array<T, 3> &vec) const {
    return vec;
  }

  std::array<T, 3> transformToDeviceCoord(SizeType /*idxSubValley*/,
                                          const std::array<T, 3> &vec) const {
    return vec;
  }
};

#endif // EMC_NONPARABOLIC_ISOTROP_SINGLELAYER_VALLEY_HPP