#ifndef GA2O3_FUNCTIONS_HPP
#define GA2O3_FUNCTIONS_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <ScatterMechanisms/emcAcousticScatterMechanism.hpp>
#include <ScatterMechanisms/emcCoulombScatterMechanism.hpp>
#include <ScatterMechanisms/emcFroehlichInteraction.hpp>
#include <ScatterMechanisms/emcHotPhononFroehlichMechanism.hpp>
#include <ScatterMechanisms/emcScreenedFroehlichInteraction.hpp>
#include <ScatterMechanisms/emcZeroOrderInterValleyScatterMechanism.hpp>
#include <ValleyTypes/emcNonParabolicIsotropValley.hpp>
#include <emcMaterial.hpp>
#include <emcPhononBath.hpp>
#include <emcPlasmonScreening.hpp>
#include <emcUtil.hpp>

/**
 * Material parameters and scatter-mechanism builders for monoclinic β-Ga₂O₃.
 *
 * Band structure
 * --------------
 * The conduction band minimum sits at Γ and is nearly isotropic; the lowest
 * satellite valley lies ≈2.4–2.5 eV above it, so it stays essentially
 * unpopulated for fields below ~450 kV/cm [Ghosh & Singisetti, arXiv:1612.03126].
 * A single non-parabolic isotropic Γ valley is therefore the appropriate
 * analytic band model over the whole field range studied here, and the
 * velocity roll-off above the peak comes from non-parabolicity plus
 * intravalley non-polar scattering rather than from intervalley transfer.
 *
 * Scattering
 * ----------
 *  - polar optical (Fröhlich): dominates the low-field mobility. β-Ga₂O₃ has
 *    12 IR-active modes; following Ma et al. these are lumped into a single
 *    effective mode at 44 meV, which reproduces the measured 300 K mobility.
 *  - intravalley non-polar optical (zero-order deformation potential) and
 *    acoustic: control momentum relaxation at high field. Values come from the
 *    analytic fits Ghosh & Singisetti published to the full-band ab initio
 *    rates expressly "to guide analytical calculations".
 *
 * Parameter sources
 * -----------------
 *  [Ma16]    N. Ma et al., Appl. Phys. Lett. 109, 212101 (2016)
 *            — m*, ε_s, ε_∞, ρ, v_s, D_A, effective POP energy 44 meV.
 *  [Ghosh17] K. Ghosh & U. Singisetti, arXiv:1612.03126 (FBMC)
 *            — analytic DP fits: D₀²/ω₀ = 7.2e18 eV/cm², D_A²/v_s² =
 *              5e-11 eV²s²/cm², satellite valley at ≈2.4 eV, v_peak ≈2e7 cm/s
 *              at 200 kV/cm.
 *  [Guo15]   Z. Guo et al., Appl. Phys. Lett. 106, 111909 (2015)
 *            — anisotropic thermal conductivity 10.9–27 W/mK (motivates the
 *              acoustic reservoir; not a direct input here).
 *
 * Parameters flagged ESTIMATED below are not established in the literature and
 * are intended to be swept, not trusted.
 */
namespace Ga2O3 {

using NumType = double;
using ValleyType = emcNonParabolicIsotropValley<NumType>;
using Acoustic = emcAcousticScatterMechanism<NumType>;
using ZeroOptAb = emcZeroOrderInterValleyAbsorptionScatterMechanism<NumType>;
using ZeroOptEm = emcZeroOrderInterValleyEmissionScatterMechanism<NumType>;
using FroehlichAb = emcFroehlichAbsorption3D<NumType>;
using FroehlichEm = emcFroehlichEmission3D<NumType>;
using HotFroehlichAb = emcHotPhononFroehlichAbsorption3D<NumType>;
using HotFroehlichEm = emcHotPhononFroehlichEmission3D<NumType>;
using ScrFroehlichAb = emcScreenedFroehlichAbsorption3D<NumType>;
using ScrFroehlichEm = emcScreenedFroehlichEmission3D<NumType>;
using ScrHotFroehlichAb = emcScreenedHotPhononFroehlichAbsorption3D<NumType>;
using ScrHotFroehlichEm = emcScreenedHotPhononFroehlichEmission3D<NumType>;
using PhononBath = emcPhononBath<NumType>;
using Screening = emcPlasmonScreening<NumType>;

// ---- material constants -----------------------------------------------------
const NumType epsLo = 10.2;    // static dielectric constant        [Ma16]
const NumType epsHi = 3.573;   // high-frequency dielectric constant[Ma16]
const NumType rho = 5880.;     // mass density [kg / m³]            [Ma16]
const NumType vSound = 6800.;  // sound velocity [m / s]            [Ma16]
const NumType bandGap = 4.85;  // band gap [eV]
// Intrinsic carrier density is ~10^-21 m^-3 at this gap and plays no role in a
// bulk simulation; a nominal non-zero value is used only to keep emcMaterial
// well-formed.
const NumType Ni = 1.;         // [1 / m³], nominal

// ---- Γ valley ---------------------------------------------------------------
const NumType relEffMass = 0.284; // [Ma16]
const SizeType degeneracy = 1;    // single Γ valley
// Kane non-parabolicity α = (1/E_g)(1 - m*/m₀)² = 0.106 1/eV.
const NumType alpha = 0.106;

// ---- polar optical (Fröhlich) ----------------------------------------------
// Effective single-mode POP energy. β-Ga₂O₃'s lowest optical modes span
// 35-48 meV; 44 meV is the value that reproduces the measured 300 K mobility
// [Ma16]. Kept as a variable so a genuine multi-mode set can replace it.
const NumType phononEnergyPOP = 0.044; // [eV]

// ---- multi-mode polar set ---------------------------------------------------
/**
 * Five-mode polar set derived from the ab initio Fröhlich coefficients of
 * Santia et al., AIP Advances 9, 015313 (2019), Table VI, which tabulates
 * C (meV²/a₀²) with |g|² = C/q² for the dominant modes along each of the three
 * reciprocal-lattice directions b₁, b₂, b₃.
 *
 * Reduction procedure:
 *  1. Each mode's coupling weight is taken as C/ħω (the combination that
 *     appears in the Fröhlich scatter constant, which carries one power of ω
 *     and one of the dielectric difference).
 *  2. The 16 tabulated (direction, mode) entries are pooled with weight 1/3
 *     each — a directional average appropriate to the isotropic Γ valley — and
 *     binned into five groups by energy.
 *  3. Weights are renormalized so that Σ_m (1/ε_∞ − 1/ε_s,m) = (1/ε_∞ − 1/ε_s)
 *     from the MEASURED dielectric constants. This enforces the oscillator-
 *     strength sum rule: the total polar coupling is fixed by experiment, and
 *     the ab initio C values only decide how it is DISTRIBUTED over modes.
 *     Each mode's effective static constant follows from
 *     1/ε_s,m = 1/ε_∞ − w_m (1/ε_∞ − 1/ε_s).
 *
 * This matters for hot phonons specifically: the 44 meV single-mode lump
 * concentrates all polar coupling into one bath, whereas the real material
 * spreads it over modes that each accumulate their own, smaller occupation
 * excess. The single-mode reduction is expected to OVERSTATE the hot-phonon
 * effect, and comparing the two is the point of this set.
 *
 * Physics check against the source: the 235 cm⁻¹ (29 meV) mode dominates the
 * phonon-ABSORPTION term and is identified there as the fundamental intrinsic
 * limit on room-temperature mobility, because its energy is near k_BT so its
 * occupation is large. The high-frequency modes (~94 meV) carry the largest C
 * but are barely activated at 300 K — they only matter once carriers heat up,
 * which is exactly the high-field regime studied here.
 */
const std::vector<NumType> modeEnergyMulti = {0.0302, 0.0429, 0.0646, 0.0796,
                                              0.0936}; // [eV]
const std::vector<NumType> modeWeightMulti = {0.1382, 0.0326, 0.2448, 0.1510,
                                              0.4334}; // Σ = 1
//! Per-mode effective static dielectric constant enforcing the sum rule.
inline std::vector<NumType> modeEpsLoMulti() {
  const NumType invHi = 1. / epsHi;
  const NumType total = invHi - 1. / epsLo;
  std::vector<NumType> out;
  for (auto w : modeWeightMulti)
    out.push_back(1. / (invHi - w * total));
  return out;
}

// ---- intravalley non-polar optical (zero-order DP) -------------------------
// [Ghosh17] fit D₀²/ω₀ = 7.2e18 eV/cm², with the two parameters explicitly
// inter-adjustable. Taking ω₀ = 90 meV (their dominant non-polar mode along
// Γ-N/Γ-Z) gives D₀ = sqrt(7.2e18 * 0.090) = 8.05e8 eV/cm = 8.05e10 eV/m.
const NumType phononEnergyNPO = 0.090;  // [eV]
const NumType defPotNPO = 8.05e10;      // [eV / m]

// ---- intravalley acoustic ---------------------------------------------------
// [Ghosh17] zone-centre fit D_A²/v_s² = 5e-11 eV²s²/cm² with v_s = 6.8e5 cm/s
// gives D_A = 4.8 eV; [Ma16] quote 6.9 eV. The two independent estimates agree
// to within 1.4x. The Ghosh value is used for consistency with the high-field
// non-polar fit; sigmaAcAlt is provided for a sensitivity check.
const NumType sigmaAc = 4.8;    // [eV]  [Ghosh17]
const NumType sigmaAcAlt = 6.9; // [eV]  [Ma16]

// ---- phonon bath (hot-phonon bottleneck) ------------------------------------
// ESTIMATED. No measured LO phonon lifetime exists for β-Ga₂O₃. Raman-linewidth
// values for its optical modes are in the ps range and the material is less
// anharmonic than GaN, so ps-scale is the plausible bracket — but this is the
// dominant uncertainty of the whole study and must be swept, not assumed.
const NumType tauLODefault = 5e-12;      // [s]  ESTIMATED — sweep 0.5-20 ps
// Klemens channel LO → 2 LA puts the daughter mode at ħω_LO/2.
const NumType acPhononEnergy = phononEnergyPOP / 2.; // [eV]
// ESTIMATED. β-Ga₂O₃'s thermal conductivity (10.9 W/mK along [100], 27 along
// [010] [Guo15]) is 5-25x below GaN's, so the acoustic reservoir is expected to
// be slow to unload. Sweep.
const NumType tauAcDefault = 20e-12;     // [s]  ESTIMATED — sweep

// q-space discretisation of the phonon bath. Polar scattering is strongly
// small-q peaked, so the emitted phonons pile up near the zone centre where the
// phonon DOS (∝ q²) is small. Resolving q is essential here: a single-bin
// (lumped) bath averages that pile-up away and systematically understates the
// occupation excess.
const SizeType nrPhononBins = 300;
const NumType dqBin = 1e7; // [1 / m] → covers q up to 3e9 1/m

//! Returns the β-Ga₂O₃ material.
template <class T> emcMaterial<T> getGa2O3Material() {
  return emcMaterial<T>{epsLo, rho, Ni, vSound, bandGap};
}

//! Adds the non-parabolic isotropic Γ valley.
template <class DerivedParticleType>
void addGammaValley(std::unique_ptr<DerivedParticleType> &particleType,
                    NumType inAlpha = alpha) {
  particleType->addValley(std::make_unique<ValleyType>(
      relEffMass, particleType->getMass(), degeneracy, inAlpha));
}

//! Adds intravalley acoustic deformation-potential scattering (elastic).
template <class DerivedParticleType, class DeviceType>
void addAcousticScattering(SizeType idxValley,
                           std::unique_ptr<DerivedParticleType> &particleType,
                           const DeviceType &device,
                           const std::vector<int> &idxRegions,
                           NumType sigma = sigmaAc) {
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<Acoustic>(idxValley, sigma, device));
}

//! Adds intravalley non-polar optical (zero-order deformation potential)
//! scattering. The zero-order intervalley mechanism is reused with the final
//! valley equal to the initial one, which is exactly an intravalley optical DP
//! process.
template <class DerivedParticleType, class DeviceType>
void addNonPolarOpticalScattering(
    SizeType idxValley, std::unique_ptr<DerivedParticleType> &particleType,
    const DeviceType &device, const std::vector<int> &idxRegions,
    NumType defPot = defPotNPO, NumType phEnergy = phononEnergyNPO) {
  const std::map<SizeType, std::vector<SizeType>> finalSubValleys = {{0, {0}}};
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<ZeroOptAb>("NPO", idxValley, finalSubValleys,
                                              defPot, phEnergy, device));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<ZeroOptEm>("NPO", idxValley, finalSubValleys,
                                              defPot, phEnergy, device));
}

//! Adds polar optical (Fröhlich) scattering with the *equilibrium*
//! Bose-Einstein phonon occupation — the assumption made by every published
//! β-Ga₂O₃ Monte Carlo study. This is the reference case.
template <class DerivedParticleType>
void addEquilibriumPolarOpticalScattering(
    SizeType idxValley, std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature,
    NumType phEnergy = phononEnergyPOP) {
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<FroehlichAb>(idxValley, phEnergy,
                                                relEffMass, epsHi, epsLo,
                                                temperature, "Ga2O3"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<FroehlichEm>(idxValley, phEnergy,
                                                relEffMass, epsHi, epsLo,
                                                temperature, "Ga2O3"));
}

//! Creates one phonon bath per polar mode.
inline std::vector<std::shared_ptr<PhononBath>>
makePhononBaths(const std::vector<NumType> &modeEnergy, NumType tauLO,
                NumType latticeTemp, NumType Vsim, bool enableAcoustic,
                NumType tauAcoustic, NumType ridleyW = 0.,
                NumType toPhononEnergy = 0., NumType tauTO = 0.) {
  std::vector<std::shared_ptr<PhononBath>> baths;
  for (auto wLO : modeEnergy) {
    baths.push_back(std::make_shared<PhononBath>(
        nrPhononBins, dqBin, tauLO, wLO, latticeTemp, Vsim, enableAcoustic,
        wLO / 2., tauAcoustic, ridleyW, toPhononEnergy, tauTO));
  }
  return baths;
}

//! Adds polar optical (Fröhlich) scattering coupled to a self-consistently
//! evolving, q-resolved phonon population — the hot-phonon case.
template <class DerivedParticleType>
void addHotPolarOpticalScattering(
    SizeType idxValley, std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions,
    const std::vector<NumType> &modeEnergy,
    const std::vector<std::shared_ptr<PhononBath>> &baths) {
  for (SizeType m = 0; m < modeEnergy.size(); ++m) {
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<HotFroehlichAb>(
                        idxValley, modeEnergy[m], relEffMass, epsHi, epsLo,
                        baths[m], "Ga2O3-" + std::to_string(m)));
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<HotFroehlichEm>(
                        idxValley, modeEnergy[m], relEffMass, epsHi, epsLo,
                        baths[m], "Ga2O3-" + std::to_string(m)));
  }
}

//! Adds ionized-impurity (Brooks-Herring) scattering, using the library's
//! existing 3D Coulomb mechanism. Elastic, with screened-Rutherford angular
//! sampling and the non-parabolic DOS correction. The ionized donor density is
//! read from the device doping profile, so it tracks --doping automatically
//! (full ionization assumed — reasonable for Si/Sn donors in n-Ga₂O₃).
//!
//! NOTE: this mechanism screens with the LATTICE temperature (the standard
//! Brooks-Herring assumption of screening by an equilibrium carrier gas),
//! whereas emcPlasmonScreening uses the carrier temperature. Under high field
//! that makes the impurity channel screen slightly more strongly than the polar
//! one. Consistent treatment would need both to use T_e.
template <class DerivedParticleType, class DeviceType>
void addImpurityScattering(SizeType idxValley,
                           std::unique_ptr<DerivedParticleType> &particleType,
                           DeviceType &device,
                           const std::vector<int> &idxRegions) {
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<emcCoulombScatterMechanism<NumType, DeviceType>>(
          idxValley, epsLo, device));
}

//! Creates the shared plasmon-screening state.
//!
//! The background permittivity the carriers screen against is a genuine
//! modelling choice. At the LO frequency the lattice response is the mode
//! being screened, which argues for the optical constant epsHi; the static
//! constant epsLo gives a longer screening length and therefore a weaker
//! effect. Both are reported in the text, so the choice is exposed here
//! rather than fixed.
inline std::shared_ptr<Screening> makeScreening(bool enabled,
                                                NumType epsBackground = epsLo) {
  return std::make_shared<Screening>(epsBackground, enabled);
}

//! Screened polar optical scattering, equilibrium phonon occupation.
//! Pass a disabled screening object to recover the unscreened rates.
template <class DerivedParticleType>
void addScreenedEquilibriumPolarOpticalScattering(
    SizeType idxValley, std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature,
    const std::vector<NumType> &modeEnergy,
    const std::vector<NumType> &modeEpsLo,
    const std::shared_ptr<Screening> &screening) {
  for (SizeType m = 0; m < modeEnergy.size(); ++m) {
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<ScrFroehlichAb>(
                        idxValley, modeEnergy[m], relEffMass, epsHi,
                        modeEpsLo[m], temperature, screening,
                        "Ga2O3-" + std::to_string(m)));
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<ScrFroehlichEm>(
                        idxValley, modeEnergy[m], relEffMass, epsHi,
                        modeEpsLo[m], temperature, screening,
                        "Ga2O3-" + std::to_string(m)));
  }
}

//! Screened polar optical scattering coupled to per-mode hot phonon baths.
template <class DerivedParticleType>
void addScreenedHotPolarOpticalScattering(
    SizeType idxValley, std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions,
    const std::vector<NumType> &modeEnergy,
    const std::vector<NumType> &modeEpsLo,
    const std::vector<std::shared_ptr<PhononBath>> &baths,
    const std::shared_ptr<Screening> &screening, bool qResolved = false,
    bool qResolvedAngle = true) {
  for (SizeType m = 0; m < modeEnergy.size(); ++m) {
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<ScrHotFroehlichAb>(
                        idxValley, modeEnergy[m], relEffMass, epsHi,
                        modeEpsLo[m], baths[m], screening, qResolved,
                        "Ga2O3-" + std::to_string(m), qResolvedAngle));
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<ScrHotFroehlichEm>(
                        idxValley, modeEnergy[m], relEffMass, epsHi,
                        modeEpsLo[m], baths[m], screening, qResolved,
                        "Ga2O3-" + std::to_string(m), qResolvedAngle));
  }
}

} // namespace Ga2O3

#endif // GA2O3_FUNCTIONS_HPP
