#ifndef EMC_PARABOLIC_ISOTROP_SINGLELAYER_VALLEY_HPP
#define EMC_PARABOLIC_ISOTROP_SINGLELAYER_VALLEY_HPP

#include <math.h>

#include <ValleyTypes/emcAbstractValley.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Parabolic Isotrop Valley for Single Layer (assumes 2D-extension in x
 * and y-direction!)
 *
 * @param degFactor number of subvalleys in that valley
 * @param effMass effective mass (in each direction)
 * @param vogtFactor factor for herring-vogt transformation, (1,1,0) for isotrop
 * 2d valley
 * @param bottomEnergy energy of the bottom of the valley (in eV)
 */
template <class T>
class emcParabolicIsotropSingleLayerValley : public emcAbstractValley<T> {
private:
  SizeType degFactor;
  T effMass;
  std::array<T, 3> vogtFactor;
  T bottomEnergy;

public:
  emcParabolicIsotropSingleLayerValley() = delete;

  /*! \brief Constructor.
   *
   * @param inRelEffMass relative effective mass of valley
   * @param inParticleMass mass of particles
   * @param inDegFactor nr of equivalent subvalleys
   * @param inBottomEnergy energy of the bottom of the valley (in eV)
   */
  emcParabolicIsotropSingleLayerValley(T inRelEffMass, T inParticleMass,
                                       SizeType inDegFactor,
                                       T inBottomValleyEnergy = 0.)
      : degFactor(inDegFactor), effMass(inRelEffMass * inParticleMass),
        vogtFactor({1, 1, 0}), bottomEnergy(inBottomValleyEnergy) {}

  T getBottomEnergy() const { return bottomEnergy; }

  T getEffMassDOS(T /*energy*/ = 0) const { return effMass; }

  T getEffMassCond(T /*energy*/ = 0) const { return effMass; }

  T getNonParabolicity() const { return 0; }

  SizeType getDegeneracyFactor() const { return degFactor; }

  T getNormWaveVec(T energy) const {
    return std::sqrt(2 * effMass * constants::q * energy) / constants::hbar;
  }

  T getEnergy(const std::array<T, 3> &k) const {
    return constants::hbar * constants::hbar * (k[0] * k[0] + k[1] * k[1]) /
           (2 * effMass * constants::q);
  }

  const std::array<T, 3> &getVogtTransformationFactor() const {
    return vogtFactor;
  }

  std::array<T, 3> getVelocity(const std::array<T, 3> &k, T /*energy*/,
                               SizeType /*idxSubValley*/) const {
    return scale(k, constants::hbar / effMass);
  }

  T getGamma(T energy) const { return energy; }

  std::array<T, 3> transformToEllipseCoord(SizeType /*idxSubValley*/,
                                           const std::array<T, 3> &vec) const {
    return vec;
  }

  std::array<T, 3> transformToDeviceCoord(SizeType /*idxSubValley*/,
                                          const std::array<T, 3> &vec) const {
    return vec;
  }
};

#endif // EMC_PARABOLIC_ISOTROP_SINGLELAYER_VALLEY_HPP