#ifndef EMC_FROEHLICH_INTERACTION_SINGLE_LAYER_HPP
#define EMC_FROEHLICH_INTERACTION_SINGLE_LAYER_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include <ScatterMechanisms/emc2DScreening.hpp>
#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Samples the in-plane deflection angle psi between the initial and
 * final wave vector for 2D (single-layer) Froehlich (polar-LO) scattering.
 *
 * In a strictly 2D system the differential scattering rate into a final state
 * at deflection angle psi is proportional to the squared 2D Froehlich matrix
 * element with a finite-width form factor (Kaasbjerg et al., PRB 85, 115317;
 * Sohier, Calandra & Mauri, PRB 94, 085415) and, optionally, free-carrier
 * screening (see emc2DScreening.hpp):
 *
 *     dS/dpsi  ~  |g(q)|^2  =  erfc(q * w / 2)^2 / ( q * eps(q)^2 ) ,
 *
 * where the magnitude of the momentum transfer is fixed by momentum and energy
 * conservation,
 *
 *     q(psi) = sqrt( k^2 + kPrime^2 - 2 k kPrime cos(psi) ) ,
 *
 * and kPrime = |k'| is set by energy conservation (E' = E +- hbar*omega). This
 * replaces the previous, physically inconsistent use of the 3D polar-optical
 * angular distribution. The (unscreened) distribution is forward-peaked, but
 * screening shifts the peak to a finite angle, so the sample is drawn by
 * numerical inversion of the cumulative distribution rather than by rejection.
 * The weight depends only on cos(psi), so the magnitude is sampled on [0, pi]
 * and the sign is then chosen at random.
 *
 * @param k        norm of the initial wave vector [1 / m]
 * @param kPrime   norm of the final wave vector [1 / m]
 * @param effWidth effective layer thickness (form-factor width) [m]
 * @param screeningWavevector 2D free-carrier screening wavevector q_s [1 / m]
 *        (0 = unscreened)
 */
template <class T>
T sampleSingleLayerFroehlichDeflectionAngle(
    T k, T kPrime, T effWidth, T screeningWavevector, emcRNG &rng,
    std::uniform_real_distribution<T> &dist) {
  auto weight = [&](T psi) -> T {
    T q2 = k * k + kPrime * kPrime - 2 * k * kPrime * std::cos(psi);
    T q = std::sqrt(std::max(T(0), q2));
    if (q <= T(0))
      return T(0);
    T ff = std::erfc(effWidth * q / 2);
    return ff * ff * twoDScreeningFactor(q, screeningWavevector) / q;
  };
  // Sample the deflection magnitude in [0, pi] by numerical inversion of the
  // (unnormalized) cumulative distribution built with the midpoint rule. This
  // is robust for both the unscreened (forward-peaked) and screened
  // (interior-peaked) distributions and needs no rejection envelope.
  constexpr SizeType N = 128;
  std::array<T, N + 1> cdf;
  T dpsi = constants::pi / N;
  cdf[0] = T(0);
  for (SizeType i = 1; i <= N; ++i)
    cdf[i] = cdf[i - 1] + weight((i - T(0.5)) * dpsi);
  T total = cdf[N];
  if (!(total > T(0)))
    return 2 * constants::pi * dist(rng); // degenerate: fall back to isotropic
  T target = dist(rng) * total;
  SizeType lo = 1;
  while (lo < N && cdf[lo] < target)
    ++lo;
  T frac = (target - cdf[lo - 1]) / (cdf[lo] - cdf[lo - 1]);
  T psiMag = (T(lo) - 1 + frac) * dpsi; // deflection magnitude in [0, pi]
  return (dist(rng) < T(0.5)) ? psiMag : (2 * constants::pi - psiMag);
}

/*! \brief Froehlich Interaction Absorption for single-layer of
 * Material.
 *
 * For single layer it is assumed that there is no movement in z-direction.
 * Formula for scattering from Paper of Kaasbjerg et al. with given link:
 * (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.115317)
 * Assumes parabolic bands.
 *
 * The scatterParticle() method samples the final-state deflection angle from
 * the proper 2D polar-optical distribution via
 * sampleSingleLayerFroehlichDeflectionAngle().
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
  T screeningWavevector;

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
   * @param inScreeningWavevector 2D free-carrier screening wavevector q_s
   *        [1 / m]; 0 (default) leaves the interaction unscreened
   */
  emcFroehlichInteractionAbsorptionSL(SizeType inValley, T inPhononEnergy,
                                      T couplingConstant, T effectiveWidth,
                                      T temperature,
                                      std::string inNameSuffix = "",
                                      T inScreeningWavevector = 0)
      : dist(0., 1.), nameSuffix(inNameSuffix),
        phononEnergy(inPhononEnergy), emcScatterMechanism<T>(inValley),
        effWidth(effectiveWidth), screeningWavevector(inScreeningWavevector) {
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

  /// \brief Inelastic, Anisotrop 2D Scattering.
  /// The deflection angle is sampled from the 2D polar-optical distribution
  /// consistent with the 2D Froehlich matrix element (see
  /// sampleSingleLayerFroehlichDeflectionAngle()).
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto valley = this->ptrValley[this->idxValley];
    T initEnergy = particle.energy;
    particle.energy = initEnergy + phononEnergy;
    assert(particle.energy > 0);

    // wave-vector magnitudes before and after the (inelastic) scattering
    T k = valley->getNormWaveVec(initEnergy);
    T kPrime = valley->getNormWaveVec(particle.energy);

    // initial in-plane angle to x-axis; atan2 handles k[0]==0 correctly
    T phi = std::atan2(particle.k[1], particle.k[0]);
    T psi = sampleSingleLayerFroehlichDeflectionAngle(
        k, kPrime, effWidth, screeningWavevector, rng, dist);

    particle.k[0] = kPrime * std::cos(phi + psi);
    particle.k[1] = kPrime * std::sin(phi + psi);
    particle.k[2] = 0;
  }

private:
  /// function whose integral needs to be approximated to calculate the
  /// scatter rate. Optionally includes 2D free-carrier screening 1/eps(q)^2.
  T integrand(T theta, T k, T eFactor) const {
    T cosTheta = std::cos(theta);
    T root = std::sqrt(cosTheta * cosTheta + eFactor);
    T q = k * (-cosTheta + root); // momentum transfer at this angle
    return (-cosTheta + root) / root *
           std::pow(std::erfc(effWidth * q / 2.), 2) *
           twoDScreeningFactor(q, screeningWavevector);
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
 * The scatterParticle() method samples the final-state deflection angle from
 * the proper 2D polar-optical distribution via
 * sampleSingleLayerFroehlichDeflectionAngle().
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
  T screeningWavevector;

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
   * @param inScreeningWavevector 2D free-carrier screening wavevector q_s
   *        [1 / m]; 0 (default) leaves the interaction unscreened
   */
  emcFroehlichInteractionEmissionSL(SizeType inValley, T inPhononEnergy,
                                    T couplingConstant, T effectiveWidth,
                                    T temperature,
                                    std::string inNameSuffix = "",
                                    T inScreeningWavevector = 0)
      : dist(0., 1.), nameSuffix(inNameSuffix),
        phononEnergy(inPhononEnergy), emcScatterMechanism<T>(inValley),
        effWidth(effectiveWidth), screeningWavevector(inScreeningWavevector) {
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

  /// \brief Inelastic, Anisotrop 2D Scattering.
  /// The deflection angle is sampled from the 2D polar-optical distribution
  /// consistent with the 2D Froehlich matrix element (see
  /// sampleSingleLayerFroehlichDeflectionAngle()).
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto valley = this->ptrValley[this->idxValley];
    T initEnergy = particle.energy;
    particle.energy = initEnergy - phononEnergy;
    assert(particle.energy > 0);

    // wave-vector magnitudes before and after the (inelastic) scattering
    T k = valley->getNormWaveVec(initEnergy);
    T kPrime = valley->getNormWaveVec(particle.energy);

    // initial in-plane angle to x-axis; atan2 handles k[0]==0 correctly
    T phi = std::atan2(particle.k[1], particle.k[0]);
    T psi = sampleSingleLayerFroehlichDeflectionAngle(
        k, kPrime, effWidth, screeningWavevector, rng, dist);

    particle.k[0] = kPrime * std::cos(phi + psi);
    particle.k[1] = kPrime * std::sin(phi + psi);
    particle.k[2] = 0;
  }

private:
  /// function whose integral needs to be approximated to calculate the
  /// scatter rate. Optionally includes 2D free-carrier screening 1/eps(q)^2.
  T integrand(T theta, T k, T eFactor) const {
    T cosTheta = std::cos(theta);
    T root = std::sqrt(cosTheta * cosTheta - eFactor);

    T qPlus = k * (cosTheta + root); // momentum transfer, "+" branch
    T partPlus = (cosTheta + root);
    partPlus *= std::pow(std::erfc(effWidth * qPlus / 2.), 2) *
                twoDScreeningFactor(qPlus, screeningWavevector);

    T qMinus = k * (cosTheta - root); // momentum transfer, "-" branch
    T partMinus = (cosTheta - root);
    partMinus *= std::pow(std::erfc(effWidth * qMinus / 2.), 2) *
                 twoDScreeningFactor(qMinus, screeningWavevector);
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