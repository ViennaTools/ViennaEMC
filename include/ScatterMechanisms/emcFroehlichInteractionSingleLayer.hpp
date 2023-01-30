#ifndef EMC_FROEHLICH_INTERACTION_SINGLE_LAYER_HPP
#define EMC_FROEHLICH_INTERACTION_SINGLE_LAYER_HPP

#include <math.h>

#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Froehlich Interaction Absorption for single-layer of
 * Material.
 *
 * For single layer it is assumed that there is no movement in z-direction.
 * Formula for scattering from Paper of Kaasbjerg et al. with given link:
 * (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.115317)
 * Assumes parabolic bands.
 *
 * TODO: adapt current scattering for 2D case!
 *
 * @param scatterConst pre-computed constant needed for calculation of
 * scatter rate
 * @param phononEnergy energy of involved phonon [eV]
 * @param effWidth effective layer thickness [m]
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcFroehlichInteractionAbsorptionSL : public emcScatterMechanism<T> {
private:
  T phononEnergy;
  T effWidth;
  T scatterConst;
  mutable std::uniform_real_distribution<T> dist;
  std::string nameSuffix;

public:
  emcFroehlichInteractionAbsorptionSL() = delete;

  /*! \brief Constructor.
   *
   * @param inValley index of initial valley (valley to which scatter mechanism
   * is assigned)
   * @param inPhononEnergy energy of the involved phonon [eV]
   * @param couplingConst coupling constant [eV / m]
   * @param temperature temperature of material (device)
   * @param inNameSuffix suffix for filename
   */
  emcFroehlichInteractionAbsorptionSL(SizeType inValley, T inPhononEnergy,
                                      T couplingConstant, T effectiveWidth,
                                      T temperature,
                                      std::string inNameSuffix = "")
      : dist(0., 1.), nameSuffix(inNameSuffix),
        phononEnergy(inPhononEnergy), emcScatterMechanism<T>(inValley),
        effWidth(effectiveWidth) {
    T exponent = constants::q * phononEnergy / (constants::kB * temperature);
    T nrPhonons = 1. / (std::exp(exponent) - 1.);
    scatterConst = std::pow(couplingConstant * constants::q, 2) * nrPhonons /
                   (2 * constants::pi * std::pow(constants::hbar, 3));
  }

  std::string getName() const { return "froehlichAbsorptionSL" + nameSuffix; }

  /*! \brief Calculation of scatter rate uses numerical integration, assumes
   * parabolic bands. */
  T getScatterRate(T energy, SizeType /*region*/) const {
    T md = this->ptrValley[this->idxValley]->getEffMassDOS();
    T k = this->ptrValley[this->idxValley]->getNormWaveVec(energy);
    T eFactor = phononEnergy / energy;
    T integral =
        approximateIntegralMidpoint(0., 2 * constants::pi, 10000, k, eFactor);
    return integral * scatterConst * md;
  }

  /// \brief Inelastic, Anisotrop Scattering.
  /// Formula found in Semiconductor Transport - Ferry page 220
  /// TODO find out if formula really usable for 2D case?
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    // adapt + store energy
    T initEnergy = particle.energy;
    T finalEnergy = particle.energy + phononEnergy;
    particle.energy = finalEnergy;
    assert(particle.energy > 0);

    // calculate the initial angle to the x-axis
    T phi = std::atan(particle.k[1] / particle.k[0]);

    // calculate angle between final and initial state
    T f = 2 * std::sqrt(initEnergy * finalEnergy) /
          std::pow(std::sqrt(initEnergy) - std::sqrt(finalEnergy), 2);
    T theta = std::acos((1 + f - std::pow(1 + 2 * f, dist(rng))) / f);

    T normK = this->ptrValley[this->idxValley]->getNormWaveVec(particle.energy);
    particle.k[0] = normK * std::cos(phi - theta);
    particle.k[1] = normK * std::sin(phi - theta);
  }

private:
  /// function whose integral needs to be approximated to calculate the
  /// scatter rate.
  T integrand(T theta, T k, T eFactor) const {
    T cosTheta = std::cos(theta);
    T root = std::sqrt(cosTheta * cosTheta + eFactor);
    return (-cosTheta + root) / root *
           std::pow(std::erfc(effWidth * k / 2. * (-cosTheta + root)), 2);
  }

  /// approximates the required integral using the midpoint rule
  template <typename... Targs>
  T approximateIntegralMidpoint(T a, T b, SizeType nrIntervals,
                                Targs... Fargs) const {
    T dx = (b - a) / (T)nrIntervals;
    T result = 0;
    for (T currPt = a + dx / 2.; currPt <= b - dx / 2.; currPt += dx)
      result += integrand(currPt, Fargs...);
    return result * dx;
  }

  // /// approximates the required integral using the trapezoidal rule
  // template <typename... Targs>
  // T approximateIntegralTrapez(T a, T b, SizeType nrIntervals,
  //                             Targs... Fargs) const {
  //   T dx = (b - a) / (T)nrIntervals;
  //   T result = 0;
  //   for (T currPt = a; currPt < b; currPt += dx)
  //     result += integrand(currPt, Fargs...) + integrand(currPt + dx,
  //     Fargs...);
  //   return result * dx / 2.;
  // }
};

/*! \brief Froehlich Interaction Emission for single-layer of
 * Material.
 *
 * For single layer it is assumed that there is no movement in z-direction.
 * Formula for scattering from Paper of Kaasbjerg et al. with given link:
 * (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.115317)
 * Assumes parabolic bands.
 *
 * @param scatterConst pre-computed constant needed for calculation of
 * scatter rate
 * @param phononEnergy energy of involved phonon [eV]
 * @param effWidth effective layer thickness [m]
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcFroehlichInteractionEmissionSL : public emcScatterMechanism<T> {
private:
  T phononEnergy;
  T effWidth;
  T scatterConst;
  std::string nameSuffix;

  mutable std::uniform_real_distribution<T> dist;

public:
  emcFroehlichInteractionEmissionSL() = delete;

  /*! \brief Constructor.
   *
   * @param inValley index of initial valley (valley to which scatter mechanism
   * is assigned)
   * @param inPhononEnergy energy of the involved phonon [eV]
   * @param couplingConst coupling constant [eV / m]
   * @param temperature temperature of material (device)
   * @param inNameSuffix suffix for filename
   */
  emcFroehlichInteractionEmissionSL(SizeType inValley, T inPhononEnergy,
                                    T couplingConstant, T effectiveWidth,
                                    T temperature,
                                    std::string inNameSuffix = "")
      : dist(0., 1.), nameSuffix(inNameSuffix),
        phononEnergy(inPhononEnergy), emcScatterMechanism<T>(inValley),
        effWidth(effectiveWidth) {
    T exponent = constants::q * phononEnergy / (constants::kB * temperature);
    T nrPhonons = std::exp(exponent) / (std::exp(exponent) - 1.);
    scatterConst = std::pow(couplingConstant * constants::q, 2) * nrPhonons /
                   (2 * constants::pi * std::pow(constants::hbar, 3));
  }

  std::string getName() const { return "froehlichEmissionSL" + nameSuffix; }

  /*! \brief Calculation of scatter rate uses numerical integration, assumes
   * parabolic bands. */
  T getScatterRate(T energy, SizeType /*region*/) const {
    if (energy > phononEnergy) {
      T md = this->ptrValley[this->idxValley]->getEffMassDOS();
      T k = this->ptrValley[this->idxValley]->getNormWaveVec(energy);
      T eFactor = phononEnergy / energy;
      T thetaMax = std::acos(std::sqrt(eFactor));
      T integral =
          approximateIntegralMidpoint(-thetaMax, thetaMax, 10000, k, eFactor);
      return integral * scatterConst * md;
    } else {
      return 0;
    }
  }

  /// \brief Inelastic, Anisotrop Scattering.
  /// Formula found in Semiconductor Transport - Ferry page 220
  /// TODO find out if formula really usable for 2D case?
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    // adapt + store energy
    T initEnergy = particle.energy;
    T finalEnergy = particle.energy - phononEnergy;
    particle.energy = finalEnergy;
    assert(particle.energy > 0);

    // calculate initial ange between k and x-axis
    T phi = std::atan(particle.k[1] / particle.k[0]);

    // calculate angle between final and initial state
    T f = 2 * std::sqrt(initEnergy * finalEnergy) /
          std::pow(std::sqrt(initEnergy) - std::sqrt(finalEnergy), 2);
    T theta = std::acos((1 + f - std::pow(1 + 2 * f, dist(rng))) / f);

    T normK = this->ptrValley[this->idxValley]->getNormWaveVec(particle.energy);
    particle.k[0] = normK * std::cos(phi + theta);
    particle.k[1] = normK * std::sin(phi + theta);
  }

private:
  /// function whose integral needs to be approximated to calculate the
  /// scatter rate
  T integrand(T theta, T k, T eFactor) const {
    T cosTheta = std::cos(theta);
    T root = std::sqrt(cosTheta * cosTheta - eFactor);

    T partPlus = (cosTheta + root);
    partPlus *= std::pow(std::erfc(effWidth * k / 2. * (cosTheta + root)), 2);

    T partMinus = (cosTheta - root);
    partMinus *= std::pow(std::erfc(effWidth * k / 2. * (cosTheta - root)), 2);
    return (partPlus + partMinus) / root;
  }

  /// approximates the required integral using the trapezoidal rule
  template <typename... Targs>
  T approximateIntegralMidpoint(T a, T b, SizeType nrIntervals,
                                Targs... Fargs) const {
    T dx = (b - a) / (T)nrIntervals;
    T result = 0;
    for (T currPt = a + dx / 2.; currPt <= b - dx / 2.; currPt += dx)
      result += integrand(currPt, Fargs...);
    return result * dx;
  }

  // /// approximates the required integral using the trapezoidal rule
  // template <typename... Targs>
  // T approximateIntegralTrapez(T a, T b, SizeType nrIntervals, Targs... Fargs)
  // const {
  //   T dx = (b - a) / (T)nrIntervals;
  //   T result = 0;
  //   for (T currPt = a + dx; currPt < b; currPt += dx)
  //     result += integrand(currPt, Fargs...) + integrand(currPt + dx,
  //     Fargs...);
  //   return result * dx / 2.;
  // }
};

#endif // EMC_FROEHLICH_INTERACTION_SINGLE_LAYER_HPP