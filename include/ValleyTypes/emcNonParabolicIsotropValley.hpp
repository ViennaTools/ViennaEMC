#ifndef EMC_NONPARABOLIC_ISOTROP_VALLEY_HPP
#define EMC_NONPARABOLIC_ISOTROP_VALLEY_HPP

#include <math.h>

#include <ValleyTypes/emcAbstractValley.hpp>
#include <emcConstants.hpp>

/** \brief class for non-parabolic valley of a semiconductor (Kane Modell)
 * and isotropical effective masses
 * @param relEffMass relative effective mass
 * @param particleMass mass of one simulated particle
 * @param degFactor degeneracy factor of energy-valley / nr of equal valleys
 * @param alpha non-parabolicity factor [1 / eV]
 * @param bottomEnergy energy of the bottom of the valley (in eV)
 * @param effMass eff mass of particle
 */
template <class T>
class emcNonParabolicIsotropValley : public emcAbstractValley<T> {
private:
  T relEffMass;
  T particleMass;
  SizeType degFactor;
  T alpha;
  T bottomEnergy;

  T effMass;
  std::array<T, 3> vogtFactor;

public:
  emcNonParabolicIsotropValley() = delete;

  /*! \brief Constructor.
   *
   * @param inRelEffMass relative effective mass of valley
   * @param inParticleMass mass of particles
   * @param inDegFactor nr of equivalent subvalleys
   * @param inAlpha non-parabolicity factor (in 1 / eV)
   * @param inBottomEnergy energy of the bottom of the valley (in eV)
   */
  emcNonParabolicIsotropValley(T inRelEffMass, T inParticleMass,
                               SizeType inDegFactor, T inAlpha,
                               T inBottomEnergy = 0.)
      : relEffMass(inRelEffMass), particleMass(inParticleMass),
        effMass(inRelEffMass * inParticleMass), degFactor(inDegFactor),
        alpha(inAlpha), vogtFactor({1, 1, 1}), bottomEnergy(inBottomEnergy) {}

  /// returns DOS effMass (energy-dependent)
  T getEffMassDOS(T energy = 0) const {
    return effMass * std::pow(1 + 2 * alpha * energy, 3.);
  }

  /// returns conductive effMass (energy-dependent)
  T getEffMassCond(T energy = 0) const {
    return effMass * (1 + 2 * energy * alpha);
  }

  T getBottomEnergy() const { return bottomEnergy; }

  /// returns nonParabolicity alpha (in 1 / eV)
  T getNonParabolicity() const { return alpha; }

  SizeType getDegeneracyFactor() const { return degFactor; }

  /// returns norm of wave-vector: k = sqrt(2 * m * Gamma) / h
  T getNormWaveVec(T energy) const {
    return std::sqrt(2 * effMass * constants::q * getGamma(energy)) /
           constants::hbar;
  }

  /// returns energy (in eV): E = 2 * gamma / (1 + 4 * alpha * gamma)
  T getEnergy(const std::array<T, 3> &k) const {
    T gamma = constants::hbar * constants::hbar * square(k) /
              (effMass * constants::q);
    return gamma / (1 + std::sqrt(1 + 2 * alpha * gamma));
  }

  /// returns velocity = h * k / (m * sqrt(1 + 4 * alpha * Gamma))
  std::array<T, 3> getVelocity(const std::array<T, 3> &k, T energy,
                               SizeType /*idxSubValley*/) const {
    T factor = constants::hbar /
               (effMass * std::sqrt(1 + 4 * alpha * getGamma(energy)));
    return scale(k, factor);
  }

  /// Vogt-Transformation Factor not needed for isotrop valleys
  /// is set to (1,1,1) in beginning
  const std::array<T, 3> &getVogtTransformationFactor() const {
    return vogtFactor;
  }

  /// returns gamma(k) = h^2 * k^2 / 2 * m (in eV)
  /// @param energy energy of particle (in eV)
  T getGamma(T energy) const { return energy * (1 + alpha * energy); }

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

#endif // EMC_NONPARABOLIC_ISOTROP_VALLEY_HPP