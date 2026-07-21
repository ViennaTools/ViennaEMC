#ifndef EMC_HUMIDITY_DIELECTRIC_HPP
#define EMC_HUMIDITY_DIELECTRIC_HPP

#include <cmath>

#include <emcConstants.hpp>

/*! \file emcHumidityDielectric.hpp
 * \brief Humidity-dependent effective permittivity of a 2D device surface
 * (ambient channel C3).
 *
 * Adsorbed water raises the effective dielectric constant of the medium above
 * the 2D channel. That permittivity change is the distinctive C3 physics: a
 * higher eps_env
 *   - weakens the (screened) charged-impurity / adsorbate Coulomb scattering
 *     (the "environmental dielectric screening", EDS, effect: |V(q)|^2 falls as
 *     1/eps_env^2), so the mobility RECOVERS with humidity, and
 *   - softens the remote surface-optical-phonon coupling (the dielectric step D
 *     in emcRemoteSurfaceOpticalPhononMechanism uses eps_env),
 * while the water molecules ALSO charge-transfer (C1) and add charged scattering
 * centres (C2). The net humidity response is the competition of these; the model
 * exposes eps_env(RH) so it can be fed into the screening (A2 / B2) and the SO
 * coupling (B1). Because the device capacitance/impedance depends on eps_env,
 * this is also the physical basis for the impedance-spectroscopy / resonant-
 * frequency selectivity idea in the proposal.
 *
 * Water coverage vs relative humidity is taken as a Langmuir isotherm here (a
 * BET / multilayer isotherm, which better captures the steep rise above ~60% RH,
 * is a refinement); the effective permittivity uses a linear (coverage-weighted)
 * effective-medium mix between the dry environment and the wetted value.
 */

/*! \brief Humidity -> effective-permittivity model.
 *
 * @param waterLangmuirConstant K_w so that theta_w = K_w RH / (1 + K_w RH),
 *        with RH the relative humidity as a fraction in [0, 1]
 * @param epsDry effective permittivity of the dry environment (e.g. ~1 for air)
 * @param epsWet effective permittivity at full water coverage (an adsorbed thin
 *        water layer, ~5-30, is far below the bulk-water ~80)
 */
template <class T> class emcHumidityDielectric {
private:
  T waterLangmuirConstant;
  T epsDry;
  T epsWet;

public:
  emcHumidityDielectric() = delete;

  emcHumidityDielectric(T inWaterLangmuirConstant, T inEpsDry, T inEpsWet)
      : waterLangmuirConstant(inWaterLangmuirConstant), epsDry(inEpsDry),
        epsWet(inEpsWet) {}

  /*! \brief Fractional water coverage theta_w(RH), RH in [0, 1]. */
  T waterCoverage(T relativeHumidity) const {
    T krh = waterLangmuirConstant * relativeHumidity;
    return krh / (1 + krh);
  }

  /*! \brief Effective environmental permittivity eps_env(RH) (coverage-weighted
   *  mix between the dry and wetted values). */
  T effectivePermittivity(T relativeHumidity) const {
    return epsDry + waterCoverage(relativeHumidity) * (epsWet - epsDry);
  }
};

#endif // EMC_HUMIDITY_DIELECTRIC_HPP
