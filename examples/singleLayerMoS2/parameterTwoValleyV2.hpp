#ifndef MOS2_PARAMETER_TWO_VALLEY_V2_HPP
#define MOS2_PARAMETER_TWO_VALLEY_V2_HPP

#include <ScatterMechanisms/emc2DChargedImpurityScatterMechanism.hpp>
#include <ScatterMechanisms/emcAcousticSingleLayerScatterMechanism.hpp>
#include <ScatterMechanisms/emcFirstOrderSingleLayerIntervalleyScatterMechanism.hpp>
#include <ScatterMechanisms/emcFroehlichInteractionSingleLayer.hpp>
#include <ScatterMechanisms/emcPiezoelectricSingleLayerScatterMechanism.hpp>
#include <ScatterMechanisms/emcRemoteSurfaceOpticalPhononMechanism.hpp>
#include <ScatterMechanisms/emcSurfaceRoughnessScatterMechanism.hpp>
#include <ScatterMechanisms/emcZeroOrderSingleLayerInterValleyScatterMechanism.hpp>
#include <ValleyTypes/emcNonParabolicAnisotropSingleLayerValley.hpp>
#include <ValleyTypes/emcNonParabolicIsotropSingleLayerValley.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

#include <memory>
#include <string>
#include <vector>

/*! \brief Two-valley (K + Q) first-principles parameter set for monolayer MoS2, VERSION 2
 * (frozen 11 Sep 2026: the kernel of cluster array 196933, kept for the comparison with v3;
 * remaining literature inputs: Kaasbjerg Froehlich constant and width, piezoelectric constant,
 * first-order K-K' couplings, density, sound velocities, phonon energies).
 *
 * Revision of the 2D Materials paper (Sep 2026). Band structure and couplings from
 * the PBE/DFPT chains of this work (norm-conserving PseudoDojo, a = 3.185 A):
 *  - K valley: m* = 0.489 m0 (isotropic, parabolic over the thermal range);
 *  - six Q valleys at 0.50 Gamma-K: m_l = 0.646, m_t = 1.258 m0, E_Q - E_K = 118 meV
 *    (measured value, Fu et al. 2017; our PBE + SOC at the experimental lattice
 *    constant gives 113 meV; the DFPT geometry without SOC gives 190 meV);
 *  - K-valley couplings from the AHC chain at k = K (Xi_LA 6.0, Xi_TA 2.3 eV at the
 *    thermal wavevector; homopolar 3.2e8 eV/cm; K-K' 1.7e8 eV/cm); K -> Q from the
 *    same chain at q = Q - K (five phonon channels);
 *  - Q-valley couplings (intravalley, Q-Q', Q-K') from the AHC chain at k = Q
 *    (chain 1b, 13 wavevectors; review/couplings/README.md).
 * Non-parabolicity from Kane fits E(1+alpha E) = hbar^2 k^2 / 2m* of the PBE bands with the
 * band-edge masses fixed (scratchpad alpha scan, 60 k-points up to 0.5 eV): the K valley is
 * trigonally warped (alpha = -0.1 along K-Gamma, 1.0 along K-M below 0.3 eV), so the
 * isotropic value 0.5 1/eV is the direction average; the Q valley gives 0.2 (transverse)
 * to 0.5 (longitudinal), rounded to 0.3 1/eV. Both are irrelevant at thermal energies
 * (1 + 2 alpha E < 1.03) and matter only for the high-field characteristic.
 * First-order intervalley couplings (K-K'), the 2D Froehlich constant and the
 * piezoelectric constant are the DFPT values of Kaasbjerg et al. (PRB 85, 115317;
 * PRB 87, 235312), unscaled. Extrinsic mechanisms (remote SO phonons, charged
 * impurities, roughness) as in parameterKaasbjerg.hpp, added to both valleys.
 * Function signatures mirror MoS2Kaasbjerg so that drivers can switch by namespace.
 */
namespace MoS2TwoValleyV2 {

using NumType = double;
using SubValleyAssignmentType = std::vector<std::vector<SizeType>>;
using isotropValleyType = emcNonParabolicIsotropSingleLayerValley<NumType>;
using anisotropValleyType = emcNonParabolicAnisotropSingleLayerValley<NumType>;
using acScatterMech = emcAcousticSingleLayerMechanism<NumType>;
using iv0AbScatterMech =
    emcZeroOrderSingleLayerInterValleyAbsorptionScatterMechanism<NumType>;
using iv0EmScatterMech =
    emcZeroOrderSingleLayerInterValleyEmissionScatterMechanism<NumType>;
using iv1AbScatterMech =
    emcFirstOrderSingleLayerInterValleyAbsorptionScatterMechanism<NumType>;
using iv1EmScatterMech =
    emcFirstOrderSingleLayerInterValleyEmissionScatterMechanism<NumType>;
using froehlichEm = emcFroehlichInteractionEmissionSL<NumType>;
using froehlichAb = emcFroehlichInteractionAbsorptionSL<NumType>;
using piezoScatterMech = emcPiezoelectricSingleLayerMechanism<NumType>;
using chargedImpurityMech = emc2DChargedImpurityScatterMechanism<NumType>;
using surfaceRoughnessMech = emcSurfaceRoughnessScatterMechanism<NumType>;
using remoteSOMech = emcRemoteSurfaceOpticalPhononMechanism<NumType>;

// ---------------- band structure ----------------
#ifndef MOS2_TV_EQK
#define MOS2_TV_EQK 0.118 // K-Q separation [eV]
#endif
#ifndef MOS2_TV_ALPHA_K
#define MOS2_TV_ALPHA_K 0.5 // non-parabolicity of the K valley [1/eV] (this work, see below)
#endif
#ifndef MOS2_TV_ALPHA_Q
#define MOS2_TV_ALPHA_Q 0.3 // non-parabolicity of the Q valleys [1/eV] (this work, see below)
#endif
#ifndef MOS2_MEFF_K
#define MOS2_MEFF_K 0.489 // K-valley mass [m0] (this work); override for the sensitivity table
#endif
#ifndef MOS2_HW_LO_GAMMA
#define MOS2_HW_LO_GAMMA 0.048 // LO phonon energy [eV] (reference); 0.0458 = DFPT, this work
#endif
const NumType relEffMassK = MOS2_MEFF_K;
const NumType relEffMassQl = 0.646;
const NumType relEffMassQt = 1.258;
const NumType energyQK = MOS2_TV_EQK;
const NumType dpCalibration = 1.0; // no calibration factor (kept for driver compatibility)
const NumType alphaK = MOS2_TV_ALPHA_K;
const NumType alphaQ = MOS2_TV_ALPHA_Q;

// ---------------- material constants (Kaasbjerg et al.) ----------------
const NumType rho = 3.1e-6;      // kg / m2
const NumType vSoundLA = 6.7e3;  // m / s
const NumType vSoundTA = 4.2e3;  // m / s
const NumType effectiveWidth = 5.41e-10;     // m
const NumType couplingConst = 0.286 * 1e-10; // 2D Froehlich coupling [eV m]
const NumType phononEnergyLOGamma = MOS2_HW_LO_GAMMA; // eV
const NumType piezoConst = 3.0e-11;          // C / m
const NumType rytovaKeldyshLength = 4.0e-9;  // m
const NumType roughnessAmplitude = 3.0e-10;  // m
const NumType roughnessCorrLength = 1.5e-9;  // m
const NumType soRemoteDistance = 5.0e-10;    // m

// ---------------- K-valley couplings (AHC chain, k = K) ----------------
const NumType XiLA_K = 6.0, XiTA_K = 2.3;    // eV, |q| = 0.14 1/A, angle averaged
const NumType D_HO_K = 3.2e10;               // eV / m, homopolar A1' intravalley
const NumType hw_HO = 0.050;                 // eV
const NumType D_KK = 1.7e10;                 // eV / m, K -> K' via the E' mode at K
const NumType hw_KK = 0.041;                 // eV
// first-order K-K' couplings (Kaasbjerg et al., unscaled)
const NumType Iv1DefPotTA = 5.9, Iv1DefPotLA = 3.9, Iv1DefPotTOGamma = 4.0, Iv1DefPotTO = 1.9; // eV
const NumType hw_TA = 0.023, hw_LA = 0.029, hw_TOK = 0.048, hw_TOGamma = 0.048;              // eV

// ---------------- K <-> Q couplings (AHC chain 1, q = Q - K) ----------------
// five phonon channels: {D [eV/m], hbar omega [eV]}; the Q -> K rates use the same
// matrix elements (time reversal), with one final K valley.
struct IvChannel { NumType D; NumType hw; const char *name; };
const std::vector<IvChannel> KQ_near = {
    {1.148e10, 0.0222, "LA"}, {1.436e10, 0.0427, "LO"}, {1.037e10, 0.0493, "HP"},
    {0.567e10, 0.0455, "TA"}, {0.220e10, 0.0165, "TAlow"}};

// ---------------- Q-valley couplings (AHC chain 1b, k = Q1 = 0.515 Gamma-K) ----------------
// Intravalley acoustic couplings at |q| = 0.12 1/A averaged over the longitudinal and
// transverse directions (Xi_L: 2.70 / 3.22, Xi_T: 3.37 / 0.29 eV); they fall to 2.3 / 1.3 eV
// at 0.24 1/A. The flexural (ZA) coupling at Q (2.6 eV-equivalent at 0.12 1/A, transverse q)
// is not part of the deformation-potential kernel. Homopolar A1' at Q: 5.4e8 eV/cm below
// 0.12 1/A (4.3 at 0.24). Intervalley channels {D [eV/m], hbar omega [eV]} from the mode
// decomposition at q = Q1 -> Q2 (60 deg, |q| = 0.68 1/A), Q1 -> Q3 (120 deg, 1.17 1/A),
// Q1 -> Q4 (180 deg, 1.35 1/A) and Q1 -> K' (1.14 1/A); the reverse Q -> K process at
// q = K - Q1 reproduces the K -> Q couplings of chain 1 within 1 %.
const NumType XiL_Q = 3.0, XiT_Q = 2.4;  // eV
const NumType D_HO_Q = 5.4e10;           // eV / m
const std::vector<IvChannel> KQ_far = {   // K -> Q' (three far valleys) and Q -> K'
    {4.95e10, 0.0491, "HP"}, {4.64e10, 0.0285, "LA"}, {2.86e10, 0.0439, "LO"}, {1.56e10, 0.0190, "TA"}};
const std::vector<IvChannel> QQ60 = {     // Q1 -> Q2 / Q6
    {3.87e10, 0.0230, "LA"}, {2.07e10, 0.0492, "HP"}, {1.41e10, 0.0454, "TO"}, {0.78e10, 0.0174, "TA"},
    {0.74e10, 0.0424, "LO"}};
const std::vector<IvChannel> QQ120 = {    // Q1 -> Q3 / Q5
    {3.15e10, 0.0428, "LO"}, {2.97e10, 0.0491, "HP"}, {1.93e10, 0.0285, "LA"}};
const std::vector<IvChannel> QQ180 = {    // Q1 -> Q4
    {5.09e10, 0.0285, "LA"}, {2.92e10, 0.0225, "TA"}, {2.43e10, 0.0474, "HP"}, {1.44e10, 0.0406, "LO"}};

// ---------------- subvalley bookkeeping (six K corners, six Q valleys) ----------------
const SubValleyAssignmentType toSameSub = {{0}, {1}, {2}, {3}, {4}, {5}};
const SubValleyAssignmentType toNextSub = {{1}, {2}, {3}, {4}, {5}, {0}};
const SubValleyAssignmentType toNbrSub = {{1, 5}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}};
const SubValleyAssignmentType toNbrNbrSub = {{2, 4}, {3, 5}, {4, 0}, {5, 1}, {0, 2}, {1, 3}};
const SubValleyAssignmentType toOppositeSub = {{3}, {4}, {5}, {0}, {1}, {2}};
const SubValleyAssignmentType toSlashSub = {{1, 3, 5}, {0, 2, 4}, {1, 3, 5}, {0, 2, 4}, {1, 3, 5}, {0, 2, 4}};
const SubValleyAssignmentType toSameKindSub = {{0, 2, 4}, {1, 3, 5}, {0, 2, 4}, {1, 3, 5}, {0, 2, 4}, {1, 3, 5}};

template <class DerivedParticleType>
void addValleys(std::unique_ptr<DerivedParticleType> &particleType) {
  particleType->addValley(std::make_unique<isotropValleyType>(
      relEffMassK, particleType->getMass(), 6, alphaK));
  const NumType a60 = constants::pi / 3.;
  std::vector<NumType> rotAngles = {0, a60, 2 * a60, 3 * a60, 4 * a60, 5 * a60};
  particleType->addValley(std::make_unique<anisotropValleyType>(
      relEffMassQl, relEffMassQt, particleType->getMass(), 6, alphaQ, rotAngles,
      energyQK));
}

template <class DerivedParticleType, class... TParam>
void helperAddZeroOrderMech(std::unique_ptr<DerivedParticleType> &particleType,
                            const std::vector<int> &idxRegions, TParam... param) {
  particleType->addScatterMechanism(idxRegions, std::make_unique<iv0AbScatterMech>(param...));
  particleType->addScatterMechanism(idxRegions, std::make_unique<iv0EmScatterMech>(param...));
}

template <class DerivedParticleType>
void addAcousticScatterMechanisms(std::unique_ptr<DerivedParticleType> &particleType,
                                  const std::vector<int> &idxRegions, NumType temperature) {
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acScatterMech>(0, XiTA_K, rho, vSoundTA, temperature, "KTA"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acScatterMech>(0, XiLA_K, rho, vSoundLA, temperature, "KLA"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acScatterMech>(1, XiT_Q, rho, vSoundTA, temperature, "QTA"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acScatterMech>(1, XiL_Q, rho, vSoundLA, temperature, "QLA"));
}

/// zero-order optical couplings: intravalley homopolar (K, Q), K-K', K<->Q, Q-Q'
template <class DerivedParticleType>
void addZeroOrderIntervalleyScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {
  // intravalley homopolar A1'
  helperAddZeroOrderMech(particleType, idxRegions, 0, 0, D_HO_K, rho, temperature, hw_HO,
                         toSameSub, "KHomoPolar");
  helperAddZeroOrderMech(particleType, idxRegions, 1, 1, D_HO_Q, rho, temperature, hw_HO,
                         toSameSub, "QHomoPolar");
  // K -> K' (one final valley)
  helperAddZeroOrderMech(particleType, idxRegions, 0, 0, D_KK, rho, temperature, hw_KK,
                         toNextSub, "KtoKprime");
  // K -> Q (three near valleys) and Q -> K (one valley), five channels each
  for (const auto &c : KQ_near) {
    helperAddZeroOrderMech(particleType, idxRegions, 0, 1, c.D, rho, temperature, c.hw,
                           toSameKindSub, std::string("KtoQ") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 1, 0, c.D, rho, temperature, c.hw,
                           toSameSub, std::string("QtoK") + c.name);
  }
  // K -> Q' (three far valleys) and Q -> K' (one valley), four channels each
  for (const auto &c : KQ_far) {
    helperAddZeroOrderMech(particleType, idxRegions, 0, 1, c.D, rho, temperature, c.hw,
                           toSlashSub, std::string("KtoQfar") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 1, 0, c.D, rho, temperature, c.hw,
                           toNextSub, std::string("QtoKprime") + c.name);
  }
  // Q -> Q' : 60 deg (two valleys), 120 deg (two), 180 deg (one)
  for (const auto &c : QQ60)
    helperAddZeroOrderMech(particleType, idxRegions, 1, 1, c.D, rho, temperature, c.hw,
                           toNbrSub, std::string("QtoQ60") + c.name);
  for (const auto &c : QQ120)
    helperAddZeroOrderMech(particleType, idxRegions, 1, 1, c.D, rho, temperature, c.hw,
                           toNbrNbrSub, std::string("QtoQ120") + c.name);
  for (const auto &c : QQ180)
    helperAddZeroOrderMech(particleType, idxRegions, 1, 1, c.D, rho, temperature, c.hw,
                           toOppositeSub, std::string("QtoQ180") + c.name);
}

/// first-order K-K' couplings of Kaasbjerg et al. (valley 0 only)
template <class DerivedParticleType>
void addFirstOrderIntervalleyScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {
  auto add = [&](NumType D1, NumType hw, const char *name) {
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<iv1EmScatterMech>(0, D1, rho, temperature, hw, name));
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<iv1AbScatterMech>(0, D1, rho, temperature, hw, name));
  };
  add(Iv1DefPotTO, hw_TOK, "TO");
  add(Iv1DefPotTOGamma, hw_TOGamma, "TOGamma");
  add(Iv1DefPotTA, hw_TA, "TA");
  add(Iv1DefPotLA, hw_LA, "LA");
}

template <class DerivedParticleType>
void addFroehlichScatterMechanisms(std::unique_ptr<DerivedParticleType> &particleType,
                                   const std::vector<int> &idxRegions, NumType temperature,
                                   NumType carrierDensity = 0, NumType envPermittivity = 1) {
  NumType qs = twoDStaticScreeningWavevector<NumType>(carrierDensity, temperature,
                                                       envPermittivity,
                                                       relEffMassK * constants::me, 4);
  for (SizeType v : {SizeType(0), SizeType(1)}) {
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<froehlichAb>(v, phononEnergyLOGamma, couplingConst,
                                                  effectiveWidth, temperature, "", qs));
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<froehlichEm>(v, phononEnergyLOGamma, couplingConst,
                                                  effectiveWidth, temperature, "", qs));
  }
}

template <class DerivedParticleType>
void addPiezoelectricScatterMechanisms(std::unique_ptr<DerivedParticleType> &particleType,
                                       const std::vector<int> &idxRegions, NumType temperature,
                                       NumType carrierDensity = 0, NumType envPermittivity = 1) {
  NumType qs = twoDStaticScreeningWavevector<NumType>(carrierDensity, temperature,
                                                       envPermittivity,
                                                       relEffMassK * constants::me, 4);
  for (SizeType v : {SizeType(0), SizeType(1)}) {
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<piezoScatterMech>(v, piezoConst, effectiveWidth, rho,
                                                       vSoundTA, temperature, "TA", qs));
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<piezoScatterMech>(v, piezoConst, effectiveWidth, rho,
                                                       vSoundLA, temperature, "LA", qs));
  }
}

template <class DerivedParticleType>
void addChargedImpurityScatterMechanism(std::unique_ptr<DerivedParticleType> &particleType,
                                        const std::vector<int> &idxRegions, NumType temperature,
                                        NumType impurityDensity, NumType carrierDensity,
                                        NumType epsAvg = 1, NumType remoteDistance = 0,
                                        NumType chargeNumber = 1) {
  NumType qs = twoDStaticScreeningWavevector<NumType>(carrierDensity, temperature, epsAvg,
                                                       relEffMassK * constants::me, 4);
  for (SizeType v : {SizeType(0), SizeType(1)})
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<chargedImpurityMech>(v, impurityDensity, epsAvg, qs,
                                                          rytovaKeldyshLength, remoteDistance,
                                                          chargeNumber));
}

template <class DerivedParticleType>
void addSurfaceRoughnessScatterMechanism(std::unique_ptr<DerivedParticleType> &particleType,
                                         const std::vector<int> &idxRegions, NumType temperature,
                                         NumType carrierDensity, NumType epsAvg = 1,
                                         NumType extraField = 0) {
  NumType qs = twoDStaticScreeningWavevector<NumType>(carrierDensity, temperature, epsAvg,
                                                       relEffMassK * constants::me, 4);
  NumType effectiveField =
      constants::q * carrierDensity / (2 * constants::eps0 * epsAvg) + extraField;
  for (SizeType v : {SizeType(0), SizeType(1)})
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<surfaceRoughnessMech>(v, effectiveField, roughnessAmplitude,
                                                           roughnessCorrLength, qs));
}

template <class DerivedParticleType>
void addRemoteSurfaceOpticalPhonon(std::unique_ptr<DerivedParticleType> &particleType,
                                   const std::vector<int> &idxRegions, NumType temperature,
                                   NumType carrierDensity, NumType epsOxInf, NumType epsOx0,
                                   NumType phononEnergy1, NumType phononEnergy2,
                                   NumType epsEnv = 1, NumType remoteDistance = soRemoteDistance) {
  NumType qs = twoDStaticScreeningWavevector<NumType>(carrierDensity, temperature, epsEnv,
                                                       relEffMassK * constants::me, 4);
  NumType dPerMode = 0.5 * (1. / (epsOxInf + epsEnv) - 1. / (epsOx0 + epsEnv));
  for (SizeType v : {SizeType(0), SizeType(1)})
    for (NumType wSO : {phononEnergy1, phononEnergy2}) {
      particleType->addScatterMechanism(
          idxRegions, std::make_unique<remoteSOMech>(v, wSO, dPerMode, remoteDistance,
                                                     temperature, false, qs));
      particleType->addScatterMechanism(
          idxRegions, std::make_unique<remoteSOMech>(v, wSO, dPerMode, remoteDistance,
                                                     temperature, true, qs));
    }
}

} // namespace MoS2TwoValleyV2

#endif // MOS2_PARAMETER_TWO_VALLEY_V2_HPP
