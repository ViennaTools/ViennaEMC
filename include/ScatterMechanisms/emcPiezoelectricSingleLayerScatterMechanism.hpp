#ifndef PIEZOELECTRIC_SINGLE_LAYER_SCATTER_MECHANISM_HPP
#define PIEZOELECTRIC_SINGLE_LAYER_SCATTER_MECHANISM_HPP

#include <array>
#include <cmath>

#include <ScatterMechanisms/emc2DScreening.hpp>
#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Intravalley piezoelectric acoustic-phonon scattering for a single
 * layer of a piezoelectric 2D semiconductor (e.g. monolayer MoS2).
 *
 * Monolayer MoS2 lacks an inversion center and is therefore piezoelectric:
 * acoustic phonons produce a long-range piezoelectric potential that scatters
 * carriers. Following Kaasbjerg, Thygesen & Jauho, PRB 87, 235312 (2013), the
 * piezoelectric electron-phonon matrix element (their Eq. 13) is
 *
 *     |M^PE_{q,lambda}|^2 = (e11 e / eps0)^2 * q^2 * erfc(q sigma/2)^2
 *                           * A_lambda(q_hat)^2 ,
 *
 * with the anisotropy factor A_TA = -sin(3 theta_q), A_LA = cos(3 theta_q),
 * angle-averaged to <A^2> = 1/2 (their Eq. 15). Unlike the short-range
 * deformation-potential coupling (|M^DP|^2 ~ q^2), the piezoelectric coupling
 * is long-range and hence screened by the free carriers; it dominates at low
 * temperature / long wavelength.
 *
 * Treated as a quasi-elastic process in the equipartition (high-temperature)
 * limit, the total out-scattering rate is (cf. their Eqs. 12 & 15)
 *
 *     1/tau(E) = m_d (e11 e/eps0)^2 kB T / (2 hbar^3 rho c_lambda^2)
 *                * (1 + 2 alpha E)
 *                * (1/pi) \int_0^pi erfc(q sigma/2)^2 / eps(q)^2 dtheta ,
 *
 * where q = 2 k sin(theta/2) is the elastic momentum transfer, eps(q) is the
 * optional 2D free-carrier dielectric function (emc2DScreening.hpp), and the
 * angular integral tends to 1 in the long-wavelength limit, reproducing the
 * analytic result Ksi^2 -> (1/2)(e11 e/eps0)^2. The final-state deflection
 * angle is sampled from the same erfc(q sigma/2)^2 / eps(q)^2 distribution.
 *
 * @param piezoConst piezoelectric constant e11 [C / m]
 * @param effWidth effective width of the electronic wave functions sigma [m]
 * @param densityMaterial 2D mass density rho [kg / m^2]
 * @param velSound acoustic sound velocity of the branch [m / s]
 * @param screeningWavevector 2D free-carrier screening wavevector q_s [1 / m]
 *        (0 = unscreened)
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcPiezoelectricSingleLayerMechanism : public emcScatterMechanism<T> {
private:
  T effWidth;
  T scatterConst;
  T screeningWavevector;
  mutable std::uniform_real_distribution<T> dist;
  std::string nameSuffix;

  static constexpr SizeType nrAngleSteps = 128;

  /// (unnormalized) angular weight erfc(q sigma/2)^2 / eps(q)^2 at scattering
  /// angle theta for the elastic momentum transfer q = 2 k sin(theta/2).
  T angularWeight(T theta, T k) const {
    T q = 2 * k * std::sin(theta / 2);
    T ff = std::erfc(effWidth * q / 2);
    return ff * ff * twoDScreeningFactor(q, screeningWavevector);
  }

public:
  emcPiezoelectricSingleLayerMechanism() = delete;

  emcPiezoelectricSingleLayerMechanism(SizeType inValley, T piezoConst,
                                       T inEffWidth, T densityMaterial,
                                       T velSound, T temperature,
                                       std::string inNameSuffix = "",
                                       T inScreeningWavevector = 0)
      : emcScatterMechanism<T>(inValley), effWidth(inEffWidth),
        screeningWavevector(inScreeningWavevector), dist(0., 1.),
        nameSuffix(inNameSuffix) {
    // effective piezoelectric coupling energy U = e11 e / eps0 [J]; the factor
    // 1/2 = <A_lambda^2> is the angular mean of the piezoelectric anisotropy.
    T couplingEnergy = piezoConst * constants::q / constants::eps0;
    scatterConst = 0.5 * couplingEnergy * couplingEnergy * constants::kB *
                   temperature /
                   (densityMaterial * velSound * velSound *
                    std::pow(constants::hbar, 3));
  }

  std::string getName() const { return "PiezoelectricSL" + nameSuffix; }

  /// \brief total out-scattering rate; valid for parabolic / nonparabolic and
  /// (through the angular integral) the finite-width, screened coupling.
  T getScatterRate(T energy, SizeType /*region*/) const {
    auto valley = this->ptrValley[this->idxValley];
    T md = valley->getEffMassDOS();
    T alpha = valley->getNonParabolicity();
    T k = valley->getNormWaveVec(energy);
    // angular integral I(E) = (1/pi) \int_0^pi weight dtheta (midpoint rule);
    // I -> 1 in the long-wavelength (unscreened, small-q) limit.
    T dtheta = constants::pi / nrAngleSteps;
    T integral = 0;
    for (SizeType i = 0; i < nrAngleSteps; ++i)
      integral += angularWeight((i + T(0.5)) * dtheta, k);
    integral *= dtheta / constants::pi;
    return md * scatterConst * (1 + 2 * alpha * energy) * integral;
  }

  /// \brief Quasi-elastic, anisotropic 2D scattering: |k'| = |k|, deflection
  /// angle sampled from erfc(q sigma/2)^2 / eps(q)^2 by CDF inversion.
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto valley = this->ptrValley[this->idxValley];
    T k = valley->getNormWaveVec(particle.energy); // elastic: |k'| = |k|

    // cumulative distribution of the deflection magnitude on [0, pi]
    std::array<T, nrAngleSteps + 1> cdf;
    T dtheta = constants::pi / nrAngleSteps;
    cdf[0] = T(0);
    for (SizeType i = 1; i <= nrAngleSteps; ++i)
      cdf[i] = cdf[i - 1] + angularWeight((i - T(0.5)) * dtheta, k);
    T total = cdf[nrAngleSteps];

    T phi = std::atan2(particle.k[1], particle.k[0]);
    T theta;
    if (!(total > T(0))) {
      theta = constants::pi * dist(rng); // degenerate: isotropic magnitude
    } else {
      T target = dist(rng) * total;
      SizeType lo = 1;
      while (lo < nrAngleSteps && cdf[lo] < target)
        ++lo;
      theta = (T(lo) - 1 + (target - cdf[lo - 1]) / (cdf[lo] - cdf[lo - 1])) *
              dtheta;
    }
    if (dist(rng) < T(0.5))
      theta = -theta; // random left / right deflection

    particle.k[0] = k * std::cos(phi + theta);
    particle.k[1] = k * std::sin(phi + theta);
    particle.k[2] = 0;
  }
};

#endif // PIEZOELECTRIC_SINGLE_LAYER_SCATTER_MECHANISM_HPP
