#ifndef EMC_2D_CHARGED_IMPURITY_SCATTER_MECHANISM_HPP
#define EMC_2D_CHARGED_IMPURITY_SCATTER_MECHANISM_HPP

#include <array>
#include <cmath>

#include <ScatterMechanisms/emc2DScreening.hpp>
#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Elastic scattering of 2D carriers off a sheet of charged impurities
 * (Born approximation), for a single-layer semiconductor such as MoS2.
 *
 * A carrier scatters off the screened Coulomb potential of an areal density
 * N_imp of point charges Z e located a distance d from the conducting sheet.
 * The 2D-Fourier-transformed, screened impurity potential is
 *
 *     V(q) = [ Z e^2 / (2 eps0 eps_avg) ] * e^{-q d} / ( q_s + q + r0 q^2 ) ,
 *
 * whose denominator combines the three screening contributions additively:
 *   - the bare 2D Coulomb  1/q                      -> the "q" term,
 *   - the Rytova-Keldysh layer response (1 + r0 q)  -> the "r0 q^2" term,
 *   - the mobile-carrier Thomas-Fermi screening     -> the "q_s" term,
 * where q_s is the 2D free-carrier screening wavevector (emc2DScreening.hpp)
 * evaluated with eps_r = eps_avg. Setting r0 = 0 recovers a purely
 * Thomas-Fermi-screened Coulomb impurity; q_s = 0 recovers the (divergent)
 * Rytova-Keldysh limit -- a finite q_s (finite carrier density) is required for
 * a well-defined total scattering rate, which is also the physical origin of
 * impurity-limited mobility increasing with carrier density (and the
 * metal-insulator transition). References: Ando, Fowler & Stern,
 * Rev. Mod. Phys. 54, 437 (1982); Ma & Jena, PRX 4, 011043 (2014).
 *
 * The scattering is elastic (static potential): |k'| = |k|, q = 2 k sin(theta/2).
 * The total out-scattering rate is
 *
 *     1/tau(E) = N_imp m_c(E) A^2 / (pi hbar^3) * \int_0^pi |V(q)|^2/A^2 dtheta ,
 *     A = Z e^2 / (2 eps0 eps_avg) ,   m_c(E) = conduction effective mass,
 *
 * and the final-state deflection angle is sampled from |V(q(theta))|^2 by CDF
 * inversion.
 *
 * @param impurityDensity areal charged-impurity density N_imp [1 / m^2]
 * @param epsAvg average relative permittivity of the environment (eps_avg)
 * @param screeningWavevector 2D free-carrier screening wavevector q_s [1 / m]
 * @param rytovaKeldyshLength r0 [m] (0 disables the Rytova-Keldysh term)
 * @param remoteDistance d [m], impurity-to-sheet distance (0 = in the interface)
 * @param chargeNumber impurity charge number Z (default 1)
 * @param nameSuffix suffix for filename
 */
template <class T>
class emc2DChargedImpurityScatterMechanism : public emcScatterMechanism<T> {
private:
  T scatterConst; //!< N_imp A^2 / (pi hbar^3)
  T screeningWavevector;
  T rytovaKeldyshLength;
  T remoteDistance;
  mutable std::uniform_real_distribution<T> dist;
  std::string nameSuffix;

  static constexpr SizeType nrAngleSteps = 512;

  /// (unnormalized) |V(q)|^2 / A^2 at scattering angle theta, elastic transfer
  /// q = 2 k sin(theta/2).
  T angularWeight(T theta, T k) const {
    T q = 2 * k * std::sin(theta / 2);
    T denom = screeningWavevector + q + rytovaKeldyshLength * q * q;
    if (denom <= T(0))
      return T(0);
    T v = std::exp(-q * remoteDistance) / denom;
    return v * v;
  }

public:
  emc2DChargedImpurityScatterMechanism() = delete;

  emc2DChargedImpurityScatterMechanism(SizeType inValley, T impurityDensity,
                                       T epsAvg, T inScreeningWavevector,
                                       T inRytovaKeldyshLength = 0,
                                       T inRemoteDistance = 0, T chargeNumber = 1,
                                       std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley),
        screeningWavevector(inScreeningWavevector),
        rytovaKeldyshLength(inRytovaKeldyshLength),
        remoteDistance(inRemoteDistance), dist(0., 1.),
        nameSuffix(inNameSuffix) {
    // prefactor A = Z e^2 / (2 eps0 eps_avg)  [J m]
    T A = chargeNumber * constants::q * constants::q /
          (2 * constants::eps0 * epsAvg);
    scatterConst = impurityDensity * A * A /
                   (constants::pi * std::pow(constants::hbar, 3));
  }

  std::string getName() const { return "ChargedImpurity2D" + nameSuffix; }

  T getScatterRate(T energy, SizeType /*region*/) const {
    auto valley = this->ptrValley[this->idxValley];
    T mc = valley->getEffMassCond(energy);
    T k = valley->getNormWaveVec(energy);
    T dtheta = constants::pi / nrAngleSteps;
    T integral = 0;
    for (SizeType i = 0; i < nrAngleSteps; ++i)
      integral += angularWeight((i + T(0.5)) * dtheta, k);
    integral *= dtheta; // \int_0^pi |V|^2/A^2 dtheta
    return scatterConst * mc * integral;
  }

  /// \brief Elastic, anisotropic 2D scattering: |k'| = |k|, deflection angle
  /// sampled from |V(q)|^2 by CDF inversion.
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto valley = this->ptrValley[this->idxValley];
    T k = valley->getNormWaveVec(particle.energy); // elastic: |k'| = |k|

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

#endif // EMC_2D_CHARGED_IMPURITY_SCATTER_MECHANISM_HPP
