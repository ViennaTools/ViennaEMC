#ifndef EMC_REMOTE_SURFACE_OPTICAL_PHONON_MECHANISM_HPP
#define EMC_REMOTE_SURFACE_OPTICAL_PHONON_MECHANISM_HPP

#include <array>
#include <cassert>
#include <cmath>

#include <ScatterMechanisms/emc2DScreening.hpp>
#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Remote surface-optical (SO) phonon scattering of 2D carriers by the
 * polar substrate / gate dielectric, for a single-layer semiconductor.
 *
 * A polar dielectric supports surface-optical phonon modes whose evanescent
 * electric field reaches into the 2D channel and scatters the carriers -- a
 * long-range, Frohlich-like "remote phonon" interaction. Following Ma & Jena,
 * PRX 4, 011043 (2014) (electron-SO Hamiltonian ~ e F_nu e^{-q d}/sqrt(q)), the
 * coupling of mode nu is set by the substrate dielectric step
 *
 *     D_nu = 1/(eps_ox^inf + eps_env) - 1/(eps_ox^0 + eps_env) ,
 *
 * and the total out-scattering rate (this class, per mode and per
 * emission/absorption branch) is
 *
 *     1/tau^{+-}(E) = 2 C m_c(E') (N_q + 1/2 -+ 1/2)
 *                     * \int_0^pi e^{-2 q d} / ( q eps(q)^2 ) dtheta ,
 *     C = e^2 omega_SO D_nu / (4 pi eps0 hbar^2) ,
 *
 * with q = sqrt(k^2 + k'^2 - 2 k k' cos theta) the momentum transfer, k' set by
 * energy conservation (E' = E -+ hbar omega_SO), N_q the Bose occupation, d the
 * carrier-to-surface distance, and eps(q) the optional 2D free-carrier
 * dielectric function (emc2DScreening.hpp). Because the lower-omega_SO polar
 * dielectrics (SiO2, HfO2) have both a larger thermal N_q and a larger D_nu,
 * this mechanism is far stronger for them than for high-omega_SO dielectrics
 * (hBN) -- the physical origin of the mobility gain from hBN encapsulation.
 *
 * NOTE: the finite-thickness (envelope) form factor of the 2D layer is omitted
 * (long-wavelength limit); the dominant q- and temperature-dependence is kept.
 *
 * @param phononEnergy SO phonon energy hbar omega_SO [eV]
 * @param couplingD dielectric step D_nu (dimensionless, >= 0)
 * @param remoteDistance carrier-to-surface distance d [m]
 * @param isEmission true for phonon emission (E -> E - hbar omega_SO)
 * @param screeningWavevector 2D free-carrier screening wavevector q_s [1 / m]
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcRemoteSurfaceOpticalPhononMechanism : public emcScatterMechanism<T> {
private:
  T phononEnergy;
  T remoteDistance;
  T scatterConst; //!< C * (N_q + 1/2 -+ 1/2)
  T screeningWavevector;
  bool emission;
  mutable std::uniform_real_distribution<T> dist;
  std::string nameSuffix;

  static constexpr SizeType nrAngleSteps = 128;

  /// (unnormalized) e^{-2 q d} / (q eps(q)^2) at deflection angle theta.
  T angularWeight(T theta, T k, T kPrime) const {
    T q2 = k * k + kPrime * kPrime - 2 * k * kPrime * std::cos(theta);
    T q = std::sqrt(std::max(T(0), q2));
    if (q <= T(0))
      return T(0);
    return std::exp(-2 * q * remoteDistance) *
           twoDScreeningFactor(q, screeningWavevector) / q;
  }

public:
  emcRemoteSurfaceOpticalPhononMechanism() = delete;

  emcRemoteSurfaceOpticalPhononMechanism(SizeType inValley, T inPhononEnergy,
                                         T couplingD, T inRemoteDistance,
                                         T temperature, bool inEmission,
                                         T inScreeningWavevector = 0,
                                         std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        remoteDistance(inRemoteDistance),
        screeningWavevector(inScreeningWavevector), emission(inEmission),
        dist(0., 1.), nameSuffix(inNameSuffix) {
    T omega = phononEnergy * constants::q / constants::hbar; // [1/s]
    T C = constants::q * constants::q * omega * couplingD /
          (4 * constants::pi * constants::eps0 * constants::hbar *
           constants::hbar);
    T x = phononEnergy * constants::q / (constants::kB * temperature);
    T nBose = 1. / (std::exp(x) - 1.);
    scatterConst = C * (emission ? nBose + 1 : nBose);
  }

  std::string getName() const {
    return std::string("RemoteSO") + (emission ? "Em" : "Ab") + nameSuffix;
  }

  T getScatterRate(T energy, SizeType /*region*/) const {
    if (emission && energy <= phononEnergy)
      return 0;
    auto valley = this->ptrValley[this->idxValley];
    T finalEnergy = emission ? energy - phononEnergy : energy + phononEnergy;
    T k = valley->getNormWaveVec(energy);
    T kPrime = valley->getNormWaveVec(finalEnergy);
    T mc = valley->getEffMassCond(finalEnergy);
    T dtheta = constants::pi / nrAngleSteps;
    T integral = 0;
    for (SizeType i = 0; i < nrAngleSteps; ++i)
      integral += angularWeight((i + T(0.5)) * dtheta, k, kPrime);
    integral *= dtheta;
    return 2 * scatterConst * mc * integral;
  }

  /// \brief Inelastic, anisotropic 2D scattering: |k'| set by E', deflection
  /// angle sampled from e^{-2qd}/(q eps^2) by CDF inversion.
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

#endif // EMC_REMOTE_SURFACE_OPTICAL_PHONON_MECHANISM_HPP
