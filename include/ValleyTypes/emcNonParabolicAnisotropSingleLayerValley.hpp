#ifndef EMC_NONPARABOLIC_ANISOTROP_SINGLELAYER_VALLEY_HPP
#define EMC_NONPARABOLIC_ANISOTROP_SINGLELAYER_VALLEY_HPP

#include <math.h>

#include <ValleyTypes/emcAbstractValley.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Non-Parabolic Anisotropical Valley for Single Layer (assumes
 * 2D-extension in x and y-direction!)
 *
 * @param alpha non-parabolicity (in 1 / eV)
 * @param degFactor number of subvalleys in that valley
 * @param effMassLong longitudinal effective mass
 * @param effMassTrans transversal effective mass
 * @param effMassCond conductive effective mass
 * @param effMassDOS density of states effective mass
 * @param vogtFactor factor for herring-vogt transformation, (1,1,0) for isotrop
 * 2d valley
 * @param bottomEnergy energy of the bottom of the valley (in eV)
 * @param rotationAngle angles which describe the orientation of the subvalleys
 * w.r.t. the x-axis (and the long. effective mass) in radians between 0 and 2 *
 * pi
 */
template <class T>
class emcNonParabolicAnisotropSingleLayerValley : public emcAbstractValley<T> {
private:
  T alpha;
  SizeType degFactor;
  T effMassLong, effMassTrans;
  T effMassCond, effMassDOS;
  T bottomEnergy;
  std::array<T, 3> vogtFactor;
  std::vector<T> rotationAngle;
  std::vector<T> sinAngle, cosAngle;

public:
  emcNonParabolicAnisotropSingleLayerValley() = delete;

  /*! \brief Constructor.
   *
   * @param relEffMassLongitudinal longitudinal rel. effective mass
   * @param relEffMassTransversal transversal rel. effective mass
   * @param inParticleMass mass of particles
   * @param inDegFactor nr of equivalent subvalleys
   * @param inAlpha non-parabolicity factor (in 1 / eV)
   * @param inBottomEnergy energy of the bottom of the valley (in eV)
   */
  emcNonParabolicAnisotropSingleLayerValley(T relEffMassLongitudinal,
                                            T relEffMassTransversal,
                                            T inParticleMass,
                                            SizeType inDegFactor, T inAlpha,
                                            T inBottomEnergy = 0.)
      : emcNonParabolicAnisotropSingleLayerValley(
            relEffMassLongitudinal, relEffMassTransversal, inParticleMass,
            inDegFactor, inAlpha, std::vector<T>(degFactor, 0.),
            inBottomEnergy) {}

  /*! \brief Constructor.
   *
   * @param relEffMassLongitudinal longitudinal rel. effective mass
   * @param relEffMassTransversal transversal rel. effective mass
   * @param inParticleMass mass of particles
   * @param inDegFactor nr of equivalent subvalleys
   * @param inAlpha non-parabolicity factor (in 1 / eV)
   * @param inRotationAngles angles which describe the orientation of each
   * subvalley w.r.t. the x-axis (and the long. effective mass) in radians
   * between 0 and 2 * pi
   * @param inBottomEnergy energy of the bottom of the valley (in eV)
   */
  emcNonParabolicAnisotropSingleLayerValley(T relEffMassLongitudinal,
                                            T relEffMassTransversal,
                                            T inParticleMass,
                                            SizeType inDegFactor, T inAlpha,
                                            std::vector<T> inRotationAngles,
                                            T inBottomEnergy = 0.)
      : degFactor(inDegFactor),
        effMassLong(relEffMassLongitudinal * inParticleMass),
        effMassTrans(relEffMassTransversal * inParticleMass),
        effMassCond(2. / (1. / effMassLong + 1. / effMassTrans)),
        effMassDOS(std::sqrt(effMassLong * effMassTrans)),
        vogtFactor({std::sqrt(effMassDOS / effMassLong),
                    std::sqrt(effMassDOS / effMassTrans), 0}),
        alpha(inAlpha), bottomEnergy(inBottomEnergy),
        rotationAngle(inRotationAngles), sinAngle(rotationAngle.size(), 0.),
        cosAngle(rotationAngle.size(), 0.) {
    if (rotationAngle.size() != degFactor) {
      emcMessage::getInstance()
          .addError("Wrong size of rotationAngle vector!")
          .print();
    }
    for (SizeType idxSubValley = 0; idxSubValley < degFactor; idxSubValley++)
      setSubValleyEllipseCoordSystem(idxSubValley, rotationAngle[idxSubValley]);
  }

  void setSubValleyEllipseCoordSystem(SizeType idxSubValley,
                                      T newRotationAngle) {
    if (newRotationAngle < 0 || newRotationAngle > 2 * constants::pi) {
      emcMessage::getInstance()
          .addError("Rotation angle has to between 0 and 2 * PI.")
          .print();
    }
    rotationAngle[idxSubValley] = newRotationAngle;
    sinAngle[idxSubValley] = std::sin(newRotationAngle);
    cosAngle[idxSubValley] = std::cos(newRotationAngle);
  }

  T getBottomEnergy() const { return bottomEnergy; }

  T getEffMassDOS(T energy = 0) const {
    return effMassDOS * std::pow(1 + 2 * alpha * energy, 3.);
  }

  T getEffMassCond(T energy = 0) const {
    return effMassCond * (1 + 2 * energy * alpha);
  }

  T getNonParabolicity() const { return alpha; }

  SizeType getDegeneracyFactor() const { return degFactor; }

  T getNormWaveVec(T energy) const {
    return std::sqrt(2 * effMassDOS * constants::q * getGamma(energy)) /
           constants::hbar;
  }

  /// returns energy (in eV): E (1 + alpha E) = (k_x^2 + k_y^2) * h^2 / (2 * m)
  T getEnergy(const std::array<T, 3> &k) const {
    T gamma = constants::hbar * constants::hbar * (k[0] * k[0] + k[1] * k[1]) /
              (effMassDOS * constants::q);
    return gamma / (1 + std::sqrt(1 + 2 * alpha * gamma));
  }

  const std::array<T, 3> &getVogtTransformationFactor() const {
    return vogtFactor;
  }

  std::array<T, 3> getVelocity(const std::array<T, 3> &k, T energy,
                               SizeType idxSubValley) const {
    std::array<T, 3> vel = transformToEllipseCoord(idxSubValley, k);
    auto nonParaFactor = std::sqrt(1 + 4 * alpha * getGamma(energy));
    std::transform(vel.begin(), vel.end() - 1, vogtFactor.begin(), vel.begin(),
                   [this, &nonParaFactor](const T &waveVec, const T &vogt) {
                     return constants::hbar * vogt * waveVec /
                            (effMassDOS * nonParaFactor);
                   });
    return transformToDeviceCoord(idxSubValley, vel);
  }

  T getGamma(T energy) const { return energy * (1 + alpha * energy); }

  // return rotated coordinate system (use predefined roationAngle)
  std::array<T, 3> transformToEllipseCoord(SizeType idxSubValley,
                                           const std::array<T, 3> &vec) const {
    auto res = vec;
    res[0] = vec[0] * cosAngle[idxSubValley] - vec[1] * sinAngle[idxSubValley];
    res[1] = vec[0] * sinAngle[idxSubValley] + vec[1] * cosAngle[idxSubValley];
    res[2] = 0;
    return res;
  }

  // rotate coordinate system back (-rotationAngle)
  std::array<T, 3> transformToDeviceCoord(SizeType idxSubValley,
                                          const std::array<T, 3> &vec) const {
    auto res = vec;
    res[0] = vec[0] * cosAngle[idxSubValley] + vec[1] * sinAngle[idxSubValley];
    res[1] = -vec[0] * sinAngle[idxSubValley] + vec[1] * cosAngle[idxSubValley];
    res[2] = 0;
    return res;
  }
};

#endif // EMC_NONPARABOLIC_ANISOTROP_SINGLELAYER_VALLEY_HPP