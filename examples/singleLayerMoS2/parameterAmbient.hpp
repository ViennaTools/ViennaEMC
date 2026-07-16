#include <Ambient/emcAdsorbateChargeTransfer.hpp>
#include <Ambient/emcAdsorbateNoise.hpp>
#include <Ambient/emcHumidityDielectric.hpp>

#include "parameterKaasbjerg.hpp"

/*! \brief Ambient / adsorbate model for single-layer MoS2 (layer "C").
 *
 * This assembles the *ambient-dependent* conductivity model on top of the
 * calibrated intrinsic + supported-device transport base (MoS2Kaasbjerg).
 * It is the direct analogue of parameterKaasbjerg / parameterLi / parameterPilotto:
 * a namespace that turns reusable kernels (include/Ambient/ + the charged-impurity
 * mechanism in include/ScatterMechanisms/) into a coherent material model, with the
 * air/substrate constants gathered in one place.
 *
 * Four ambient channels, all driven by a Langmuir coverage theta(P) / theta(RH):
 *   C1  charge-transfer doping        -> effectiveCarrierDensity(P)
 *   C2  adsorbate Coulomb scattering  -> chargedAdsorbateDensity(P) (fed to B2)
 *   C3  humidity dielectric screening -> effectiveEnvPermittivity(RH)
 *   C4  adsorption/desorption noise    -> makeOnefNoiseModel(...)
 *
 * Default constants reproduce the MoS2/HfO2 FET measured in someData/ (the S13
 * air-vs-N2 experiment): at O2 half-coverage the model gives sigma/sigma0 ~ 0.46,
 * matching the measured I_on(air)/I_on(N2) = 0.463. See MoS2_ambient_note.md.
 *
 * NOTE: this model is MoS2-specific and does not touch the perovskite/silicon
 * models or the other MoS2 parameter sets (Kaasbjerg/Li/Pilotto) — it only *uses*
 * the Kaasbjerg builders.
 */
namespace MoS2Ambient {

using NumType = double;

// --- baseline device (MoS2/HfO2 FET, matches the someData validation) ---------
const NumType baselineDensity = 1e17;        // n0 [m^-2] = 1e13 cm^-2
const NumType substratePermittivity = 20.0;  // HfO2-ish (for the C3 average)
const NumType envScreeningPermittivity = 2.0; // eps for the C2 impurity screening
                                              // (top-surface adsorbate sees ~air)

// The ambient adsorbate layer plays TWO distinct roles. The split + magnitudes are
// calibrated to the literature air/vacuum factor (field-effect mobility 5-10x lower
// in air; water is the PRIMARY Coulomb scatterer, O2 the modest acceptor):
//
//   C2 scattering -- charged Coulomb centres (water + physi-/chemisorbed O2),
//     areal density N_c^max, charge ~1, weakly screened (eps~2). DOMINANT (mu 5-10x).
//   C1 doping     -- O2 charge transfer, a SMALLER acceptor density at charge 0.5.
//
// Both scale with the same Langmuir coverage theta(P) of the ambient layer. The
// calibrated N_c^max ~ 1.5e13 cm^-2 sits in the literature range (1e12-1e13 upper
// end); see MoS2_calibration_note.md / validation/calib_ambient.cpp.

// --- charged Coulomb scatterers (C2, dominant: water + O2) --------------------
const NumType scattererSiteDensity = 1.5e17; // N_c^max [m^-2] = 1.5e13 cm^-2
const NumType scattererCharge = 1.0;         // Z ~ 1 (water-dominated)

// --- O2 acceptor (C1 doping, secondary) --------------------------------------
const NumType o2LangmuirConstant = 1e-4;     // K [1/Pa]
const NumType siteDensity = 5e16;            // N_ct [m^-2] = 5e12 cm^-2 (acceptor)
const NumType electronsPerO2 = -0.5;         // z < 0: acceptor (0.5 e- per molecule)

// --- water adsorbate (drives C3 dielectric screening) -------------------------
const NumType waterLangmuirConstant = 3.0;   // K_w [per RH-fraction]
const NumType permittivityDry = 1.0;         // eps above the flake, dry
const NumType permittivityWet = 30.0;        // eps above the flake, saturated

// factory helpers for the underlying kernels ----------------------------------
inline emcAdsorbateChargeTransfer<NumType> o2Adsorbate() {
  return emcAdsorbateChargeTransfer<NumType>(o2LangmuirConstant, siteDensity,
                                             electronsPerO2);
}
inline emcHumidityDielectric<NumType> waterFilm() {
  return emcHumidityDielectric<NumType>(waterLangmuirConstant, permittivityDry,
                                        permittivityWet);
}

// --- ambient bookkeeping (the "what does the gas do" layer) -------------------
//! C1: O2 coverage at partial pressure P [Pa].
inline NumType coverage(NumType P) { return o2Adsorbate().coverage(P); }

//! C1: carrier density [m^-2] after charge-transfer depletion at pressure P.
inline NumType effectiveCarrierDensity(NumType P) {
  return o2Adsorbate().effectiveCarrierDensity(baselineDensity, P);
}

//! C1: areal density [m^-2] of charged O2 acceptors (the doping count).
inline NumType chargedAdsorbateDensity(NumType P) {
  return o2Adsorbate().chargedAdsorbateDensity(P);
}

//! C2: areal density [m^-2] of charged Coulomb scatterers (water + O2) at
//! coverage theta(P). This is the DOMINANT ambient mobility knob.
inline NumType chargedScattererDensity(NumType P) {
  return coverage(P) * scattererSiteDensity;
}

//! C3: effective environment permittivity seen by the impurity screening at RH,
//! averaged over the flake (adsorbed-water side) and the substrate.
inline NumType effectiveEnvPermittivity(NumType RH) {
  return 0.5 * (waterFilm().effectivePermittivity(RH) + substratePermittivity);
}

// --- the transport model (the "how does that change mobility" layer) ----------
//! Assemble the full scattering stack for a given ambient-derived state:
//!   n       effective carrier density (C1), sets the free-carrier screening,
//!   nImp    areal density of charged impurities (C2); if 0 the term is omitted,
//!   epsAvg  environment permittivity for the impurity screening (C3),
//!   impCharge charge number Z of each impurity (scattering ~ Z^2). Defaults to
//!           the scatterer charge (~1, water-dominated); the C2 density passed in
//!           should then be chargedScattererDensity(P). (For fixed substrate/
//!           interface charges — the humidity demo — pass Z~1 as well.)
//! The intrinsic + polar + piezo terms come straight from the calibrated
//! Kaasbjerg base; only the charged-impurity term carries the ambient state.
template <class DerivedParticleType>
void addAmbientScatterMechanisms(std::unique_ptr<DerivedParticleType> &particleType,
                                 NumType temperature, NumType n, NumType nImp,
                                 NumType epsAvg,
                                 NumType impCharge = scattererCharge) {
  MoS2Kaasbjerg::addValleys(particleType);
  MoS2Kaasbjerg::addAcousticScatterMechanisms(particleType, {0}, temperature);
  MoS2Kaasbjerg::addZeroOrderIntervalleyScatterMechanisms(particleType, {0},
                                                          temperature);
  MoS2Kaasbjerg::addFirstOrderIntervalleyScatterMechanisms(particleType, {0},
                                                           temperature);
  MoS2Kaasbjerg::addFroehlichScatterMechanisms(particleType, {0}, temperature, n,
                                               1.0);
  MoS2Kaasbjerg::addPiezoelectricScatterMechanisms(particleType, {0}, temperature,
                                                   n, 1.0);
  if (nImp > 0) // C2: charged impurities as 2D Coulomb centres at the surface (d=0)
    MoS2Kaasbjerg::addChargedImpurityScatterMechanism(
        particleType, {0}, temperature, nImp, n, epsAvg, 0.0, impCharge);
}

// --- C4: adsorption/desorption conductivity noise -----------------------------
//! Build a McWhorter-style noise model: nSites two-level adsorption sites with
//! switching rates gamma_i log-uniform over [gammaMin, gammaMax] (each at
//! theta = 0.5). The superposition gives 1/f between gammaMin/2pi and gammaMax/2pi
//! and 1/f^2 above; the band width reflects the spread of desorption barriers.
inline emcAdsorbateNoise<NumType>
makeOnefNoiseModel(SizeType nSites, NumType gammaMin, NumType gammaMax,
                   emcRNG &rng) {
  std::vector<NumType> ads(nSites), des(nSites);
  std::uniform_real_distribution<NumType> u(0., 1.);
  const NumType lgMin = std::log(gammaMin), lgMax = std::log(gammaMax);
  for (SizeType i = 0; i < nSites; ++i) {
    NumType gamma = std::exp(lgMin + u(rng) * (lgMax - lgMin));
    ads[i] = des[i] = gamma / 2; // theta_i = 0.5
  }
  return emcAdsorbateNoise<NumType>(ads, des);
}

} // namespace MoS2Ambient
