#ifndef EMC_SURFACE_ROUGHNESS_SCATTER_MECHANISM_HPP
#define EMC_SURFACE_ROUGHNESS_SCATTER_MECHANISM_HPP

#include <array>
#include <cmath>

#include <ScatterMechanisms/emc2DScreening.hpp>
#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Interface-roughness-limited scattering of 2D carriers
 * (Prange-Nee / Ando model), for a single-layer semiconductor such as MoS2.
 *
 * Fluctuations of the confining interface couple to the carriers through the
 * normal (out-of-plane) field F they experience; the scattering matrix element
 * is the product of that field and the roughness power spectrum,
 *
 *     |M(q)|^2 = (e F)^2 |Delta(q)|^2 ,
 *     |Delta(q)|^2 = pi Delta^2 Lambda^2 exp( -q^2 Lambda^2 / 4 )   (Gaussian) ,
 *
 * with Delta the RMS roughness amplitude and Lambda the correlation length. For
 * a substrate-supported monolayer the roughness is substrate-conformal and F is
 * the vertical field (from the gate, the carriers, and depletion); a convenient
 * self-field estimate is F = e n_2D / (2 eps0 eps_avg). Because F grows with
 * carrier density, this mechanism strengthens at high density (large gate
 * overdrive) -- the opposite trend to charged-impurity scattering -- and so
 * sets the high-density roll-off of the mobility. The interaction is screened by
 * the mobile carriers via eps(q) (emc2DScreening.hpp). References: Prange & Nee,
 * Phys. Rev. 168, 779 (1968); Ando, Fowler & Stern, RMP 54, 437 (1982).
 *
 * The scattering is elastic (static interface potential): |k'| = |k|,
 * q = 2 k sin(theta/2). The total out-scattering rate is
 *
 *     1/tau(E) = (e F)^2 Delta^2 Lambda^2 m_c(E) / hbar^3
 *                * \int_0^pi exp(-q^2 Lambda^2/4) / eps(q)^2 dtheta ,
 *
 * and the final-state deflection angle is sampled from that distribution.
 *
 * @param effectiveField normal (out-of-plane) field F [V / m]
 * @param roughnessAmplitude RMS roughness Delta [m]
 * @param correlationLength Lambda [m]
 * @param screeningWavevector 2D free-carrier screening wavevector q_s [1 / m]
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcSurfaceRoughnessScatterMechanism : public emcScatterMechanism<T> {
private:
  T scatterConst; //!< (e F)^2 Delta^2 Lambda^2 / hbar^3
  T corrLengthSq; //!< Lambda^2
  T screeningWavevector;
  mutable std::uniform_real_distribution<T> dist;
  std::string nameSuffix;

  static constexpr SizeType nrAngleSteps = 256;

  /// (unnormalized) exp(-q^2 Lambda^2/4)/eps(q)^2 at angle theta, elastic
  /// transfer q = 2 k sin(theta/2).
  T angularWeight(T theta, T k) const {
    T q = 2 * k * std::sin(theta / 2);
    T formFactor = std::exp(-q * q * corrLengthSq / 4);
    return formFactor * twoDScreeningFactor(q, screeningWavevector);
  }

public:
  emcSurfaceRoughnessScatterMechanism() = delete;

  emcSurfaceRoughnessScatterMechanism(SizeType inValley, T effectiveField,
                                      T roughnessAmplitude, T correlationLength,
                                      T inScreeningWavevector,
                                      std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley),
        corrLengthSq(correlationLength * correlationLength),
        screeningWavevector(inScreeningWavevector), dist(0., 1.),
        nameSuffix(inNameSuffix) {
    T eF = constants::q * effectiveField; // [J / m]
    scatterConst = eF * eF * roughnessAmplitude * roughnessAmplitude *
                   corrLengthSq / std::pow(constants::hbar, 3);
  }

  std::string getName() const { return "SurfaceRoughness" + nameSuffix; }

  T getScatterRate(T energy, SizeType /*region*/) const {
    auto valley = this->ptrValley[this->idxValley];
    T mc = valley->getEffMassCond(energy);
    T k = valley->getNormWaveVec(energy);
    T dtheta = constants::pi / nrAngleSteps;
    T integral = 0;
    for (SizeType i = 0; i < nrAngleSteps; ++i)
      integral += angularWeight((i + T(0.5)) * dtheta, k);
    integral *= dtheta; // \int_0^pi exp(-q^2 L^2/4)/eps^2 dtheta
    return scatterConst * mc * integral;
  }

  /// \brief Elastic, anisotropic 2D scattering: |k'| = |k|, deflection angle
  /// sampled from the roughness power spectrum by CDF inversion.
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

#endif // EMC_SURFACE_ROUGHNESS_SCATTER_MECHANISM_HPP
