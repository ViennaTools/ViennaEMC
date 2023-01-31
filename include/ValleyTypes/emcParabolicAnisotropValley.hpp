#ifndef EMC_PARABOLIC_ANISOTROP_VALLEY_HPP
#define EMC_PARABOLIC_ANISOTROP_VALLEY_HPP

#include <algorithm>
#include <math.h>

#include <ValleyTypes/emcAbstractValley.hpp>
#include <emcConstants.hpp>
#include <emcMessage.hpp>

/** \brief class for parabolic valley of a semiconductor
 * and anisotropical effective masses
 * @param relEffMass vector holding the relative effective mass in each
 * direction
 * @param particleMass mass of one simulated particle
 * @param degFactor degeneracy factor of energy-valley / nr of equal valleys
 * @param rotMatrices rotation matrix that relates device and ellipse coord
 * system for each subvalley
 * @param bottomEnergy energy of the bottom of the valley (in eV)
 */
template <class T>
class emcParabolicAnisotropValley : public emcAbstractValley<T> {
private:
  std::array<T, 3> relEffMass;
  T particleMass;
  SizeType degFactor;
  std::vector<std::array<T, 9>> rotMatrices;
  T bottomEnergy;

  T effMassCond;
  T effMassDOS;
  std::array<T, 3> vogtFactor;

public:
  emcParabolicAnisotropValley() = delete;

  /*! \brief Constructor.
   *
   * @param inRelEffMass vector holding the relative effective mass in each
   * direction
   * @param inParticleMass mass of particles
   * @param inDegFactor nr of equivalent subvalleys
   * @param inBottomEnergy energy of the bottom of the valley (in eV)
   */
  emcParabolicAnisotropValley(std::array<T, 3> inRelEffMass, T inParticleMass,
                              SizeType inDegFactor, T inBottomEnergy = 0.)
      : relEffMass(inRelEffMass), particleMass(inParticleMass),
        degFactor(inDegFactor), bottomEnergy(inBottomEnergy) {
    initRotationMatrices();
    initEffMasses();
    initVogtFactor();
  }

  /*! \brief Sets the ECS system of a specific subvalley.
   *
   * Note: The 3 input-directions have to be orthogonal to each other.
   * They don't have to be normalized.
   *
   * @param idxSubValley index of the subvalley
   * @param dir1 direction of m_1
   * @param dir2 direction of m_2
   * @param dir3 direction of m_3
   */
  void setSubValleyEllipseCoordSystem(SizeType idxSubValley,
                                      std::array<T, 3> dir1,
                                      std::array<T, 3> dir2,
                                      std::array<T, 3> dir3) {
    checkIdxSubValley(idxSubValley);
    checkOrthogonalBasis(dir1, dir2, dir3);
    normalize(dir1);
    std::copy(dir1.begin(), dir1.end(), rotMatrices[idxSubValley].begin());
    normalize(dir2);
    std::copy(dir2.begin(), dir2.end(), rotMatrices[idxSubValley].begin() + 3);
    normalize(dir3);
    std::copy(dir3.begin(), dir3.end(), rotMatrices[idxSubValley].begin() + 6);
  }

  /// returns DOS effMass (not energy dependent).
  /// geometric mean of effective masses in each direction
  T getEffMassDOS(T /*energy*/ = 0) const { return effMassDOS; }

  /// returns conductive effMass (not energy dependent)
  /// harmonic mean of effective masses in each direction
  T getEffMassCond(T /*energy*/ = 0) const { return effMassCond; }

  T getBottomEnergy() const { return bottomEnergy; }

  /// nonParabolicity is 0 for parabolic valleys
  T getNonParabolicity() const { return 0; }

  /// returns degeneracy Factor of valley
  SizeType getDegeneracyFactor() const { return degFactor; }

  /// returns norm of wave-vector: k = sqrt(2 * effMassCond * E) / h
  T getNormWaveVec(T energy) const {
    return std::sqrt(2 * effMassCond * constants::q * energy) / constants::hbar;
  }

  /// returns energy (in eV): E = k^2 * h^2 / (2 * effMassCond)
  T getEnergy(const std::array<T, 3> &k) const {
    return constants::hbar * constants::hbar * square(k) /
           (2 * effMassCond * constants::q);
  }

  /// returns velocity = h * k * vogt / effMassCond  (in ellipse coords)
  std::array<T, 3> getVelocity(const std::array<T, 3> &k, T /*energy*/,
                               SizeType idxSubValley) const {
    std::array<T, 3> vel = transformToEllipseCoord(idxSubValley, k);
    std::transform(vel.begin(), vel.end(), vogtFactor.begin(), vel.begin(),
                   [this](const T &waveVec, const T &vogt) {
                     return constants::hbar * vogt * waveVec / effMassCond;
                   });
    return transformToDeviceCoord(idxSubValley, vel);
  }

  /// returns Vogt-Transformation Factor (in ellipse coords)
  const std::array<T, 3> &getVogtTransformationFactor() const {
    return vogtFactor;
  }

  /// return gamma = energy (in eV)
  T getGamma(T energy) const { return energy; }

  /// helper that transforms vec from device to ellipse coordinates
  /// assumes: rows of rotMatrix are normalized
  std::array<T, 3> transformToEllipseCoord(SizeType idxSubValley,
                                           const std::array<T, 3> &vec) const {
    std::array<T, 3> res;
    res[0] = vec[0] * rotMatrices[idxSubValley][0] +
             vec[1] * rotMatrices[idxSubValley][1] +
             vec[2] * rotMatrices[idxSubValley][2];
    res[1] = vec[0] * rotMatrices[idxSubValley][3] +
             vec[1] * rotMatrices[idxSubValley][4] +
             vec[2] * rotMatrices[idxSubValley][5];
    res[2] = vec[0] * rotMatrices[idxSubValley][6] +
             vec[1] * rotMatrices[idxSubValley][7] +
             vec[2] * rotMatrices[idxSubValley][8];
    return res;
  }

  /// helper that transforms vec from ellipse to device coordinates
  /// assumes: rows of rotMatrix are normalized
  std::array<T, 3> transformToDeviceCoord(SizeType idxSubValley,
                                          const std::array<T, 3> &vec) const {
    std::array<T, 3> res;
    res[0] = vec[0] * rotMatrices[idxSubValley][0] +
             vec[1] * rotMatrices[idxSubValley][3] +
             vec[2] * rotMatrices[idxSubValley][6];
    res[1] = vec[0] * rotMatrices[idxSubValley][1] +
             vec[1] * rotMatrices[idxSubValley][4] +
             vec[2] * rotMatrices[idxSubValley][7];
    res[2] = vec[0] * rotMatrices[idxSubValley][2] +
             vec[1] * rotMatrices[idxSubValley][5] +
             vec[2] * rotMatrices[idxSubValley][8];
    return res;
  }

private:
  /// rotation matrices for each subValley are initialized as identity matrices
  void initRotationMatrices() {
    std::array<T, 9> rotMatrix;
    std::fill(rotMatrix.begin(), rotMatrix.end(), 0);
    rotMatrix[0] = 1;
    rotMatrix[4] = 1;
    rotMatrix[8] = 1;
    for (SizeType idxSubValley = 0; idxSubValley < degFactor; idxSubValley++)
      rotMatrices.push_back(rotMatrix);
  }

  /// helper that calculates conductive and DOS effective mass
  void initEffMasses() {
    effMassDOS = std::accumulate(relEffMass.begin(), relEffMass.end(), 1.,
                                 std::multiplies<T>());
    effMassDOS = std::pow(effMassDOS, 1. / 3.) * particleMass;

    effMassCond =
        std::accumulate(relEffMass.begin(), relEffMass.end(), 0.,
                        [](T &acc, const T &inc) { return acc += 1. / inc; });
    effMassCond = 3. * particleMass / effMassCond;
  }

  /// helper that calculates the Vogt-factor = sqrt(EffMassCond / effMass)
  void initVogtFactor() {
    std::transform(relEffMass.begin(), relEffMass.end(), vogtFactor.begin(),
                   [&, this](const T &effMass) {
                     return std::sqrt(getEffMassCond() /
                                      (effMass * particleMass));
                   });
  }

  /// helper function that checks if the three directions are all orthogonal
  /// that is needed for them to build a coordinate system
  void checkOrthogonalBasis(const std::array<T, 3> &dir1,
                            const std::array<T, 3> &dir2,
                            const std::array<T, 3> &dir3) {
    if (innerProduct(dir1, dir2) != 0 || innerProduct(dir1, dir3) != 0 ||
        innerProduct(dir2, dir3) != 0) {
      emcMessage::getInstance()
          .addError(
              "The given coordinate system for a subvalley is not orthogonal, "
              "adapt that.")
          .print();
    }
  }

  void checkIdxSubValley(SizeType idxSubValley) {
    if (idxSubValley >= degFactor) {
      emcMessage::getInstance()
          .addError("Using invalid subvalley index.")
          .print();
    }
  }
};

#endif // EMC_PARABOLIC_ANISOTROP_VALLEY_HPP