#ifndef EMC_ADSORBATE_CHARGE_TRANSFER_HPP
#define EMC_ADSORBATE_CHARGE_TRANSFER_HPP

#include <cmath>

#include <emcConstants.hpp>

/*! \file emcAdsorbateChargeTransfer.hpp
 * \brief Surface-charge-transfer doping of a 2D semiconductor by adsorbed gas
 * molecules (ambient channel C1).
 *
 * The dominant DC effect of the ambient on a 2D conductometric device: gas
 * molecules adsorb on the surface (fractional coverage set by a Langmuir
 * isotherm in the partial pressure / relative humidity) and each adsorbed
 * molecule exchanges charge with the channel. For n-type MoS2, O2 and H2O act
 * as acceptors (they withdraw electrons), so the mobile carrier sheet density
 * n drops as coverage rises, reducing the conductivity sigma = n e mu(n).
 *
 * The three inputs are exactly what the ab-initio workplan supplies:
 *   - adsorption energy E_ads (Langmuir isotherm / desorption barrier) [DFT],
 *   - charge transfer per molecule z (electrons into the channel; z < 0 for an
 *     acceptor) [DFT Bader],
 *   - areal density of adsorption sites N_site (e.g. S-vacancy density) [in-situ].
 *
 * References for the mechanism: Late et al., ACS Nano 6, 5635 (2012); Qiu et
 * al., APL 100, 123104 (2012); "Charge-transfer-based gas sensing", Sci. Rep. 5,
 * 8052 (2015).
 */

/*! \brief Langmuir equilibrium constant K(T) = K0 exp(E_ads / kB T) [1/Pa].
 *  Stronger binding (larger E_ads) -> larger K -> more coverage at a given
 *  pressure. K0 is the pre-exponential (attempt/entropic) factor [1/Pa]. */
template <class T>
T langmuirConstant(T K0, T adsorptionEnergy, T temperature) {
  return K0 * std::exp(adsorptionEnergy * constants::q /
                       (constants::kB * temperature));
}

/*! \brief Surface-charge-transfer doping model for one adsorbate species.
 *
 * @param equilibriumConstant Langmuir constant K [1/Pa] (or dimensionless if the
 *        "pressure" argument is a relative quantity such as relative humidity)
 * @param siteDensity areal density of adsorption sites N_site [1/m^2]
 * @param electronsPerMolecule signed charge transfer per adsorbed molecule z
 *        (electrons added to the channel; z < 0 for an acceptor like O2 / H2O)
 */
template <class T> class emcAdsorbateChargeTransfer {
private:
  T equilibriumConstant;
  T siteDensity;
  T electronsPerMolecule;

public:
  emcAdsorbateChargeTransfer() = delete;

  emcAdsorbateChargeTransfer(T inEquilibriumConstant, T inSiteDensity,
                             T inElectronsPerMolecule)
      : equilibriumConstant(inEquilibriumConstant), siteDensity(inSiteDensity),
        electronsPerMolecule(inElectronsPerMolecule) {}

  /*! \brief Fractional surface coverage theta = K P / (1 + K P) in [0, 1]. */
  T coverage(T pressure) const {
    T kp = equilibriumConstant * pressure;
    return kp / (1 + kp);
  }

  /*! \brief Change in mobile carrier sheet density [1/m^2] induced by the
   *  adsorbate at the given pressure (negative for an acceptor). */
  T carrierDensityChange(T pressure) const {
    return coverage(pressure) * siteDensity * electronsPerMolecule;
  }

  /*! \brief Effective carrier sheet density n = n0 + Delta n [1/m^2], floored at
   *  0 (full depletion). n0 is the pristine (vacuum) carrier density. */
  T effectiveCarrierDensity(T baselineDensity, T pressure) const {
    return std::max(T(0), baselineDensity + carrierDensityChange(pressure));
  }

  /*! \brief Areal density of CHARGED adsorbates [1/m^2] = theta * N_site (each
   *  charge-transferred molecule is a charged scattering centre). Used by the
   *  adsorbate Coulomb-scattering channel (C2). */
  T chargedAdsorbateDensity(T pressure) const {
    return coverage(pressure) * siteDensity;
  }
};

#endif // EMC_ADSORBATE_CHARGE_TRANSFER_HPP
