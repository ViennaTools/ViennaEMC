#ifndef EMC_PARABOLIC_ISOTROP_VALLEY_HPP
#define EMC_PARABOLIC_ISOTROP_VALLEY_HPP

#include <math.h>

#include <ValleyTypes/emcAbstractValley.hpp>
#include <emcConstants.hpp>

/** \brief class for parabolic valley of a semiconductor
 * with isotropical effective mass
 * @param relEffMass relative effective mass of particle
 * @param particleMass mass of one simulated particle
 * @param effMass effective mass of particle
 * @param degFactor degeneracy factor of energy-valley / nr of equal valleys
 * @param bottomEnergy energy of the bottom of the valley (in eV)
 */
template <class T>
class emcParabolicIsotropValley : public emcAbstractValley<T> {
private:
  T relEffMass;
  T particleMass;
  T effMass;
  SizeType degFactor;
  std::array<T, 3> vogtFactor;
  T bottomEnergy;

public:
  emcParabolicIsotropValley() = delete;

  /*! \brief Constructor.
   *
   * @param inRelEffMass relative effective mass of valley
   * @param inParticleMass mass of particles
   * @param inDegFactor nr of equivalent subvalleys
   * @param inBottomEnergy energy of the bottom of the valley (in eV)
   */
  emcParabolicIsotropValley(T inRelEffMass, T inParticleMass,
                            SizeType inDegFactor, T inBottomEnergy = 0.)
      : relEffMass(inRelEffMass), particleMass(inParticleMass),
        degFactor(inDegFactor), effMass(inRelEffMass * inParticleMass),
        vogtFactor({1, 1, 1}), bottomEnergy(inBottomEnergy) {}

  /// returns DOS effMass (not energy dependent)
  T getEffMassDOS(T /*energy*/ = 0) const { return effMass; }

  /// returns counductive effMass (not energy dependent)
  T getEffMassCond(T /*energy*/ = 0) const { return effMass; }

  T getBottomEnergy() const { return bottomEnergy; }

  /// nonParabolicity is 0 for parabolic valleys
  T getNonParabolicity() const { return 0; }

  SizeType getDegeneracyFactor() const { return degFactor; }

  /// returns norm of wave-vector: k = sqrt(2 * m * E) / h
  T getNormWaveVec(T energy) const {
    return std::sqrt(2 * effMass * constants::q * energy) / constants::hbar;
  }

  /// returns energy (in eV): E = k^2 * h^2 / (2 * m)
  T getEnergy(const std::array<T, 3> &k) const {
    return constants::hbar * constants::hbar * square(k) /
           (2 * effMass * constants::q);
  }

  /// Vogt-Transformation Factor not needed for isotrop valleys
  /// is set to (1,1,1) in beginning
  const std::array<T, 3> &getVogtTransformationFactor() const {
    return vogtFactor;
  }

  /// returns velocity = h * k / m
  std::array<T, 3> getVelocity(const std::array<T, 3> &k, T /*energy*/,
                               SizeType /*idxSubValley*/) const {
    return scale(k, constants::hbar / effMass);
  }

  /// return gamma = energy (in eV)
  T getGamma(T energy) const { return energy; }

  /// isotrop valley characteristics are independent of coordinate system,
  /// just returns vector in given coordinate system
  std::array<T, 3> transformToEllipseCoord(SizeType /*idxSubValley*/,
                                           const std::array<T, 3> &vec) const {
    return vec;
  }

  /// isotrop valley characteristics are independent of coordinate system,
  /// just returns vector in given coordinate system
  std::array<T, 3> transformToDeviceCoord(SizeType /*idxSubValley*/,
                                          const std::array<T, 3> &vec) const {
    return vec;
  }
};

#endif // EMC_PARABOLIC_ISOTROP_VALLEY_HPP