#ifndef EMC_SCREENED_INTRAVALLEY_OPTICAL_MECHANISM_HPP
#define EMC_SCREENED_INTRAVALLEY_OPTICAL_MECHANISM_HPP

#include <array>
#include <cassert>
#include <cmath>

#include <ScatterMechanisms/emc2DScreening.hpp>
#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Free-carrier-screened intravalley optical-phonon scattering for a
 * single-layer semiconductor (isotropic valley).
 *
 * This is the small-momentum-transfer, intravalley optical (in particular the
 * polar-optical / zone-centre) deformation-potential process treated WITH
 * free-carrier screening, following the approach used in state-of-the-art 2D
 * transport models (e.g. Iglesias et al., 2D Mater. 2023; Ma & Jena,
 * PRX 4, 011043 (2014)): the deformation-potential coupling is kept and the
 * scattering matrix is divided by the 2D static dielectric function eps(q)
 * (emc2DScreening.hpp). This avoids the double counting that would arise from
 * adding an explicit Froehlich term on top of a deformation potential that
 * already folds in the polar mode.
 *
 * Relative to the (isotropic, q-independent) zero-order intervalley mechanism,
 * the total rate acquires the angular screening factor
 *
 *     1/tau^{+-}(E) = m_d(E') scatterConst (1 + 2 alpha E')
 *                     * (1/pi) \int_0^pi dtheta / eps(q)^2 ,
 *
 * with q = sqrt(k^2 + k'^2 - 2 k k' cos theta), k' set by energy conservation
 * (E' = E -+ hbar omega). With eps -> 1 the factor is 1 and the original
 * (unscreened) rate is recovered. The final-state deflection angle is sampled
 * from 1/eps(q)^2 (screening suppresses small-angle scattering). Only intravalley
 * (initial valley == final valley, same subvalley) transitions are handled.
 *
 * @param sigma optical deformation potential [eV / m]
 * @param densityMaterial 2D mass density rho [kg / m^2]
 * @param phononEnergy optical phonon energy [eV]
 * @param isEmission true for phonon emission (E -> E - hbar omega)
 * @param screeningWavevector 2D free-carrier screening wavevector q_s [1 / m]
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcScreenedIntravalleyOpticalMechanism : public emcScatterMechanism<T> {
private:
  T phononEnergy;
  T scatterConst; //!< (sigma q_e/hbar)^2 (N_q + 1/2 -+ 1/2) / (2 rho omega)
  T screeningWavevector;
  bool emission;
  mutable std::uniform_real_distribution<T> dist;
  std::string nameSuffix;

  static constexpr SizeType nrAngleSteps = 128;

  /// (unnormalized) 1/eps(q)^2 at deflection angle theta.
  T angularWeight(T theta, T k, T kPrime) const {
    T q2 = k * k + kPrime * kPrime - 2 * k * kPrime * std::cos(theta);
    T q = std::sqrt(std::max(T(0), q2));
    return twoDScreeningFactor(q, screeningWavevector);
  }

public:
  emcScreenedIntravalleyOpticalMechanism() = delete;

  emcScreenedIntravalleyOpticalMechanism(SizeType inValley, T sigma,
                                         T densityMaterial, T temperature,
                                         T inPhononEnergy, bool inEmission,
                                         T inScreeningWavevector = 0,
                                         std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        screeningWavevector(inScreeningWavevector), emission(inEmission),
        dist(0., 1.), nameSuffix(inNameSuffix) {
    T x = phononEnergy * constants::q / (constants::kB * temperature);
    T omega = phononEnergy * constants::q / constants::hbar;
    T nBose = 1. / (std::exp(x) - 1.);
    T nFactor = emission ? nBose + 1 : nBose;
    scatterConst = std::pow(sigma * constants::q / constants::hbar, 2) *
                   nFactor / (2 * densityMaterial * omega);
  }

  std::string getName() const {
    return std::string("ScreenedIntraOptical") + (emission ? "Em" : "Ab") +
           nameSuffix;
  }

  T getScatterRate(T energy, SizeType /*region*/) const {
    if (emission && energy <= phononEnergy)
      return 0;
    auto valley = this->ptrValley[this->idxValley];
    T finalEnergy = emission ? energy - phononEnergy : energy + phononEnergy;
    T md = valley->getEffMassDOS();
    T alpha = valley->getNonParabolicity();
    T k = valley->getNormWaveVec(energy);
    T kPrime = valley->getNormWaveVec(finalEnergy);
    T dtheta = constants::pi / nrAngleSteps;
    T screenAvg = 0;
    for (SizeType i = 0; i < nrAngleSteps; ++i)
      screenAvg += angularWeight((i + T(0.5)) * dtheta, k, kPrime);
    screenAvg /= nrAngleSteps; // (1/pi) \int_0^pi 1/eps^2 dtheta
    return md * scatterConst * (1 + 2 * alpha * finalEnergy) * screenAvg;
  }

  /// \brief Inelastic intravalley scattering (valley/subvalley unchanged);
  /// deflection angle sampled from 1/eps(q)^2 by CDF inversion.
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto valley = this->ptrValley[this->idxValley];
    T initEnergy = particle.energy;
    particle.energy =
        emission ? initEnergy - phononEnergy : initEnergy + phononEnergy;
    assert(particle.energy > 0);
    T k = valley->getNormWaveVec(initEnergy);
    T kPrime = valley->getNormWaveVec(particle.energy);

    std::array<T, nrAngleSteps + 1> cdf;
    T dtheta = constants::pi / nrAngleSteps;
    cdf[0] = T(0);
    for (SizeType i = 1; i <= nrAngleSteps; ++i)
      cdf[i] = cdf[i - 1] + angularWeight((i - T(0.5)) * dtheta, k, kPrime);
    T total = cdf[nrAngleSteps];

    T phi = std::atan2(particle.k[1], particle.k[0]);
    T psi;
    if (!(total > T(0))) {
      psi = constants::pi * dist(rng);
    } else {
      T target = dist(rng) * total;
      SizeType lo = 1;
      while (lo < nrAngleSteps && cdf[lo] < target)
        ++lo;
      psi = (T(lo) - 1 + (target - cdf[lo - 1]) / (cdf[lo] - cdf[lo - 1])) *
            dtheta;
    }
    if (dist(rng) < T(0.5))
      psi = -psi;

    particle.k[0] = kPrime * std::cos(phi + psi);
    particle.k[1] = kPrime * std::sin(phi + psi);
    particle.k[2] = 0;
  }
};

#endif // EMC_SCREENED_INTRAVALLEY_OPTICAL_MECHANISM_HPP
