#ifndef MOS2_PARAMETER_TWO_VALLEY_HPP
#define MOS2_PARAMETER_TWO_VALLEY_HPP

#include <ScatterMechanisms/emc2DChargedImpurityScatterMechanism.hpp>
#include <ScatterMechanisms/emcAcousticSingleLayerScatterMechanism.hpp>
#include <ScatterMechanisms/emcAcousticTabulatedSingleLayerScatterMechanism.hpp>
#include <ScatterMechanisms/emcFirstOrderSingleLayerIntervalleyScatterMechanism.hpp>
#include <ScatterMechanisms/emcFroehlichInteractionSingleLayer.hpp>
#include <ScatterMechanisms/emcRemoteSurfaceOpticalPhononMechanism.hpp>
#include <ScatterMechanisms/emcSurfaceRoughnessScatterMechanism.hpp>
#include <ScatterMechanisms/emcZeroOrderSingleLayerInterValleyScatterMechanism.hpp>
#include <ValleyTypes/emcNonParabolicAnisotropSingleLayerValley.hpp>
#include <ValleyTypes/emcNonParabolicIsotropSingleLayerValley.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

/*! \brief Two-valley (K + Q) first-principles parameter set for monolayer MoS2, VERSION 3.
 *
 * Revision of the 2D Materials paper (Sep 2026). Every intrinsic input is a value computed
 * for this work (PBE, norm-conserving PseudoDojo, a = 3.185 A; DFPT chains 1, 1b and 2 of
 * review/couplings/README.md and Supplementary Note S2):
 *  - K valley: isotropic Kane band with m* = 0.446 m0 and alpha = 1.16 1/eV, the least-squares fit of
 *    E(1 + alpha E) = hbar^2 k^2 / 2m* to the direction-averaged (area-equivalent) k(E) of the Wannier
 *    interpolated PBE band over 0 to 0.15 eV (K-M: m 0.458, alpha 1.5; K-Gamma: m 0.44, alpha 0.4 to 0.85;
 *    the fit reproduces the area-equivalent k(E) to 1 % and the band velocities to a few per cent,
 *    whereas the band-edge curvature mass 0.489 with alpha 0.5 reproduces k(E) to 4 % and overestimates
 *    the velocity by up to 17 % at 0.2 eV); the band-edge mass 0.489 (0.477 K-Gamma, 0.530 K-M) is the
 *    quantity reported in the main text;
 *  - Q valleys: m_l = 0.646, m_t = 1.258 m0, alpha = 0.3 1/eV, at E_QK = 118 meV (measured,
 *    Fu et al. 2017; our PBE + SOC at the experimental lattice constant gives 113 meV) and,
 *    with the spin-orbit splitting Delta_Q = 71 meV of Note S1, split into two spin-resolved
 *    sets of three valleys 120 deg apart: Q_low at E_QK and Q_high at E_QK + Delta_Q (for an
 *    electron of given spin the three nearest Q valleys of K and the three farthest ones carry
 *    opposite spin ordering; the EMC has no spin label, so the K -> Q channels are averaged over
 *    spin, D / sqrt(2) into each set, while the Q -> Q and Q -> K channels are spin conserving);
 *  - lattice: mass density 3.03e-6 kg/m^2 from the cell, sound velocities 6.6e3 (LA) and
 *    4.1e3 m/s (TA) from the DFPT acoustic slopes, phonon energies from the DFPT modes;
 *  - K-valley couplings: Xi_LA 6.0, Xi_TA 2.3 eV (AHC chain 1, angle averaged at the thermal
 *    wavevector), K-K' 1.7e8 eV/cm at 40.6 meV, homopolar A1' |g(q)| tabulated from the EPW
 *    interpolation, polar LO |g(q)| tabulated from the EPW interpolation with the 2D long-range
 *    treatment, joined below 0.04 1/A to the strictly 2D limit 0.35 eV / (1 + r0 q) from our Born
 *    charges and dielectric tensor (r0 = 43.5 A); free-carrier screening added as before;
 *  - Q-valley couplings (chain 1b): Xi_L 3.0, Xi_T 2.4 eV, homopolar |g(q)| tabulated;
 *  - intervalley K <-> Q (5 + 4 channels), Q <-> Q' (60, 120, 180 deg) from chains 1 and 1b;
 *  - first-order K-K' couplings from chain 1 with q = K + dq (|dq| = 0.03 to 0.10 1/A, two directions;
 *    linear in |dq| and isotropic): TA 12.0, LA 2.0, TO(K) 4.8 eV and 3.5 eV for the 47.4 meV mode;
 *    intravalley TO(Gamma) 3.6 eV from the EPW TO branch (Kaasbjerg: 5.9 / 3.9 / 1.9 / 4.0).
 * No piezoelectric term (the computed acoustic couplings contain the long-range part; the
 * double count would be below 1 %). Extrinsic mechanisms (remote SO phonons, charged impurities,
 * roughness) as in parameterKaasbjerg.hpp, added to all three valley types.
 * Function signatures mirror MoS2Kaasbjerg so that drivers switch by namespace
 * (parameterKernelSelect.hpp: -DMOS2_TWO_VALLEY; -DMOS2_TWO_VALLEY_V2 selects the frozen v2).
 *  - K-valley acoustic couplings wavevector dependent (xiLA_K_q, xiTA_K_q; MOS2_TV_ACOUSTIC_TAB=0 gives
 *    the constants 6.0 / 2.3 eV); Q valley constants 3.0 / 2.4 eV (two wavevectors only).
 * Sensitivity macros: MOS2_TV_XI_LA_K / MOS2_TV_XI_TA_K (acoustic constants), MOS2_TV_LDA_OPTICAL=1
 * (Kaasbjerg optical set: constant homopolar 4.1e8, K-K' 2.6e8 eV/cm, erfc Froehlich 97 meV),
 * MOS2_TV_FIRST_ORDER=0 (no first-order K-K' couplings), MOS2_TV_EQK, MOS2_TV_DELTAQ (0 = v2 valley model).
 */
namespace MoS2TwoValley {

using NumType = double;
using SubValleyAssignmentType = std::vector<std::vector<SizeType>>;
using isotropValleyType = emcNonParabolicIsotropSingleLayerValley<NumType>;
using anisotropValleyType = emcNonParabolicAnisotropSingleLayerValley<NumType>;
using acScatterMech = emcAcousticSingleLayerMechanism<NumType>;
using acTabScatterMech = emcAcousticTabulatedSingleLayerMechanism<NumType>;
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
using chargedImpurityMech = emc2DChargedImpurityScatterMechanism<NumType>;
using surfaceRoughnessMech = emcSurfaceRoughnessScatterMechanism<NumType>;
using remoteSOMech = emcRemoteSurfaceOpticalPhononMechanism<NumType>;

// ---------------- compile-time switches ----------------
#ifndef MOS2_TV_EQK
#define MOS2_TV_EQK 0.118 // K-Q separation [eV], lower spin subband of Q
#endif
#ifndef MOS2_TV_DELTAQ
#define MOS2_TV_DELTAQ 0.071 // spin-orbit splitting of the Q valleys [eV] (Note S1); 0 = degenerate
#endif
#ifndef MOS2_TV_ALPHA_K
#define MOS2_TV_ALPHA_K 1.16 // non-parabolicity of the K valley [1/eV]: isotropic Kane fit to the Wannier band, 0 to 0.15 eV, both directions (see below)
#endif
#ifndef MOS2_TV_ALPHA_Q
#define MOS2_TV_ALPHA_Q 0.3 // non-parabolicity of the Q valleys [1/eV]
#endif
#ifndef MOS2_MEFF_K
#define MOS2_MEFF_K 0.446 // K-valley transport mass [m0] of the isotropic Kane fit (band-edge curvature mass 0.489); override for the sensitivity table
#endif
#ifndef MOS2_HW_LO_GAMMA
#define MOS2_HW_LO_GAMMA 0.0463 // E' (LO) phonon energy at Gamma [eV], DFPT this work (11.2 THz)
#endif
#ifndef MOS2_TV_XI_LA_K
#define MOS2_TV_XI_LA_K 6.0 // eV
#endif
#ifndef MOS2_TV_XI_TA_K
#define MOS2_TV_XI_TA_K 2.3 // eV
#endif
#ifndef MOS2_TV_ACOUSTIC_TAB
#define MOS2_TV_ACOUSTIC_TAB 1 // 1 = Xi_LA(q), Xi_TA(q) tables for the K valley; 0 = constants XiLA_K, XiTA_K
#endif
#ifndef MOS2_TV_LDA_OPTICAL
#define MOS2_TV_LDA_OPTICAL 0 // 1 = Kaasbjerg optical couplings (sensitivity row)
#endif
#ifndef MOS2_TV_FIRST_ORDER
#define MOS2_TV_FIRST_ORDER 1 // 0 = no first-order K-K' couplings (sensitivity row)
#endif
#ifndef MOS2_TV_IV1_TA
#define MOS2_TV_IV1_TA 12.0 // eV, first-order K-K' TA(K, 22.5 meV): this work (chain 1, q = K + dq, linear and isotropic; Kaasbjerg 5.9)
#endif
#ifndef MOS2_TV_IV1_LA
#define MOS2_TV_IV1_LA 2.0 // eV, first-order K-K' LA(K, 28.5 meV): this work (0.9 to 2.8 eV over directions and |dq|; Kaasbjerg 3.9)
#endif
#ifndef MOS2_TV_IV1_TOK
#define MOS2_TV_IV1_TOK 4.8 // eV, first-order K-K' TO(K, 45.7 meV): this work (4.6 to 5.2 eV; Kaasbjerg 1.9)
#endif
#ifndef MOS2_TV_IV1_A1K
#define MOS2_TV_IV1_A1K 3.5 // eV, first-order K-K' via the 47.4 meV (A1'-like at K) mode: this work (2.9 to 3.9 eV; not in Kaasbjerg's set)
#endif
#ifndef MOS2_TV_IV1_TOG
#define MOS2_TV_IV1_TOG 3.6 // eV, first-order intravalley TO(Gamma): EPW TO branch, this work
#endif

// ---------------- band structure ----------------
const NumType relEffMassK = MOS2_MEFF_K;
const NumType relEffMassQl = 0.646;
const NumType relEffMassQt = 1.258;
const NumType energyQK = MOS2_TV_EQK;
const NumType deltaQ = MOS2_TV_DELTAQ;
const NumType dpCalibration = 1.0; // no calibration factor (kept for driver compatibility)
const NumType alphaK = MOS2_TV_ALPHA_K;
const NumType alphaQ = MOS2_TV_ALPHA_Q;

// ---------------- lattice (this work: a = 3.185 A, DFPT phonons) ----------------
const NumType latticeConstant = 3.185e-10;                                   // m
const NumType cellArea = std::sqrt(3.) / 2 * latticeConstant * latticeConstant; // m^2 (8.79 A^2)
const NumType rho = 3.03e-6;                 // kg / m^2, (M_Mo + 2 M_S) / A_cell
const NumType vSoundLA = 6.6e3;              // m / s, LA slope of the DFPT dispersion (6.58e3 at 0.04 to 0.15 1/A)
const NumType vSoundTA = 4.0e3;              // m / s, TA slope (EPW dispersion 4.02e3 at 0.04 to 0.15 1/A)
const NumType phononEnergyLOGamma = MOS2_HW_LO_GAMMA; // eV, E' at Gamma
const NumType hw_HO = 0.0492;                // eV, A1' (11.90 THz)
const NumType hw_KK = 0.0406;                // eV, E' (LO-like) at K (9.81 THz)
const NumType hw_TA_K = 0.0225, hw_LA_K = 0.0285, hw_TO_K = 0.0457, hw_A1_K = 0.0474; // eV, modes at K (5.44, 6.90, 11.05, 11.47 THz)
const NumType rytovaKeldyshLength = 4.0e-9;  // m (DFPT 2D polarisability, main text Sec. 3)
const NumType roughnessAmplitude = 3.0e-10;  // m   (extrinsic, unchanged)
const NumType roughnessCorrLength = 1.5e-9;  // m   (extrinsic, unchanged)
const NumType soRemoteDistance = 5.0e-10;    // m   (extrinsic, unchanged)

// ---------------- K-valley couplings (AHC chain 1, k = K) ----------------
const NumType XiLA_K = MOS2_TV_XI_LA_K, XiTA_K = MOS2_TV_XI_TA_K; // eV, |q| = 0.14 1/A, angle averaged
// ---------------- wavevector-dependent acoustic couplings of the K valley ----------------
// Xi_LA(q), Xi_TA(q) [eV] over q [1/m] from the EPW interpolation along Gamma-M (2D Coulomb cutoff and
// 2D long-range treatment, Fig. S3a: LA 6.5 eV at q -> 0 falling to 5.3 at 0.3 1/A, TA rising from 0 to
// about 1 eV), times an angular correction from the direct DFPT chain restricted to the directions in
// which that chain is free of its 3D-periodic long-range artefact. The chain without 2D cutoff carries a
// 1/q piezoelectric-like term of about 0.41 eV per (1/A) that appears in the LA channel along Gamma-M and
// in the TA channel along Gamma-K (three-fold pattern; e.g. TA 7.8 / 5.0 / 4.1 / 3.3 eV at 0.05 / 0.08 /
// 0.10 / 0.13 1/A along Gamma-K), so the LA correction uses the Gamma-K LA values (5.6 / 5.4 / 5.3 / 5.0
// / 4.2 eV at 0.05 to 0.26 1/A: rms with the Gamma-M curve, x 0.91 at 0.137 and x 0.90 at 0.30 1/A) and
// the TA correction the artefact-subtracted Gamma-K and 15 deg values (about 1.0 eV at 0.137, 1.8 at 0.26:
// x 1.3 at 0.137, x 1.5 at 0.30); both x 1.0 below 0.05 1/A (isotropic deformation-potential limit) and
// constant beyond 0.8 1/A. The previously used constants 6.0 / 2.3 eV (angle average of the direct chain
// at 0.137 1/A) contained the artefact in the TA channel; MOS2_TV_ACOUSTIC_TAB=0 still selects them.
const std::vector<std::pair<NumType, NumType>> xiLA_K_q = {{0.0000e10, 6.500}, {0.0380e10, 6.495}, {0.0759e10, 6.254}, {0.1139e10, 5.919}, {0.1519e10, 5.645}, {0.1898e10, 5.477}, {0.2278e10, 5.276}, {0.2658e10, 5.050}, {0.3037e10, 4.795}, {0.3417e10, 4.500}, {0.3797e10, 4.129}, {0.4176e10, 3.676}, {0.4556e10, 3.200}, {0.4936e10, 2.830}, {0.5315e10, 2.655}, {0.5695e10, 2.653}, {0.6074e10, 2.755}, {0.6454e10, 2.924}, {0.6834e10, 3.179}, {0.7213e10, 3.613}, {0.7593e10, 4.384}, {0.7973e10, 5.100}};
const std::vector<std::pair<NumType, NumType>> xiTA_K_q = {{0.0000e10, 0.000}, {0.0380e10, 0.122}, {0.0759e10, 0.378}, {0.1139e10, 0.677}, {0.1519e10, 0.962}, {0.1898e10, 1.179}, {0.2278e10, 1.341}, {0.2658e10, 1.441}, {0.3037e10, 1.476}, {0.3417e10, 1.432}, {0.3797e10, 1.410}, {0.4176e10, 1.504}, {0.4556e10, 1.761}, {0.4936e10, 2.093}, {0.5315e10, 2.376}, {0.5695e10, 2.561}, {0.6074e10, 2.651}, {0.6454e10, 2.679}, {0.6834e10, 2.692}, {0.7213e10, 2.757}, {0.7593e10, 2.943}, {0.7973e10, 3.039}};
const NumType D_KK = 1.7e10;                 // eV / m, K -> K' via the E' mode at K (zero order)
// first-order couplings [eV] (matrix element Xi1 sqrt(|k|^2 + |k'|^2), Kaasbjerg's definition)
const NumType Iv1DefPotTA = MOS2_TV_IV1_TA, Iv1DefPotLA = MOS2_TV_IV1_LA;
const NumType Iv1DefPotTO = MOS2_TV_IV1_TOK, Iv1DefPotTOGamma = MOS2_TV_IV1_TOG, Iv1DefPotA1K = MOS2_TV_IV1_A1K;

// ---------------- tabulated intravalley optical couplings |g(q)| ----------------
// Used through the 2D Froehlich mechanism class with a tabulated form factor:
// |g(q)| = C / sqrt(A_cell) * F(q), C = g0 sqrt(A_cell), F(0) = 1; q in 1/m.
// LO (E') polar coupling |g(q)| of the K conduction band: q = 0 to 0.03 1/A from the strictly two-dimensional
// limit g0 / (1 + r0 q) with g0 = 0.349 eV from our Born charges (Z_Mo = -0.99, Z_S = 0.49, E' mode with the
// Mo-S reduced mass, hw 46.2 meV) and r0 = 2 pi alpha_2D = 43.5 A from our dielectric tensor (eps_par 5.31);
// it gives 132 meV at 0.038 1/A against 138 meV from EPW. From 0.038 1/A: EPW with the 2D long-range treatment
// (mode 7 below, mode 6 above the crossing at 0.15 1/A, q along Gamma-M up to 0.8 1/A, held constant beyond).
const NumType g0_LO = 0.3490; // eV, value at q = 0 (form factor 1)
const std::vector<std::pair<NumType, NumType>> ff_LO = {{0.0000e10, 1.0000}, {0.0100e10, 0.6969}, {0.0200e10, 0.5348}, {0.0300e10, 0.4338}, {0.0380e10, 0.3941}, {0.0759e10, 0.2746}, {0.1139e10, 0.2283}, {0.1519e10, 0.2030}, {0.1898e10, 0.1858}, {0.2278e10, 0.1725}, {0.2658e10, 0.1612}, {0.3037e10, 0.1505}, {0.3417e10, 0.1381}, {0.3797e10, 0.1225}, {0.4176e10, 0.1032}, {0.4556e10, 0.0832}, {0.4936e10, 0.0816}, {0.5315e10, 0.0893}, {0.5695e10, 0.0954}, {0.6074e10, 0.1000}, {0.6454e10, 0.1037}, {0.6834e10, 0.1075}, {0.7213e10, 0.1135}, {0.7593e10, 0.1247}, {0.7973e10, 0.1346}};
// homopolar A1' coupling |g(q)| of the K conduction band, EPW mode 8 (q = 0 held at the first point)
const NumType g0_A1K = 0.0774; // eV, value at q = 0 (form factor 1)
const std::vector<std::pair<NumType, NumType>> ff_A1K = {{0.0000e10, 1.0000}, {0.0380e10, 1.0000}, {0.0759e10, 0.8730}, {0.1139e10, 0.7248}, {0.1519e10, 0.5725}, {0.1898e10, 0.4345}, {0.2278e10, 0.3214}, {0.2658e10, 0.2358}, {0.3037e10, 0.1784}, {0.3417e10, 0.1544}, {0.3797e10, 0.1681}, {0.4176e10, 0.2077}, {0.4556e10, 0.2551}, {0.4936e10, 0.2995}, {0.5315e10, 0.3388}, {0.5695e10, 0.3756}, {0.6074e10, 0.4129}, {0.6454e10, 0.4546}, {0.6834e10, 0.5080}, {0.7213e10, 0.0000}, {0.7593e10, 0.0000}, {0.7973e10, 0.0000}};
// homopolar A1' coupling |g(q)| of the Q conduction band, AHC chain 1b (average of longitudinal and transverse q; held constant beyond 0.59 1/A)
const NumType g0_A1Q = 0.0900; // eV, value at q = 0 (form factor 1)
const std::vector<std::pair<NumType, NumType>> ff_A1Q = {{0.0000e10, 1.0000}, {0.0590e10, 1.0000}, {0.1180e10, 0.9833}, {0.2370e10, 0.7767}, {0.3950e10, 0.2844}, {0.5920e10, 0.1244}};
const NumType sqrtCellArea = std::sqrt(cellArea);
// Kaasbjerg reference set (sensitivity row MOS2_TV_LDA_OPTICAL=1 only)
const NumType kaasbjergFroehlichConst = 0.286e-10; // eV m (97 meV)
const NumType kaasbjergEffectiveWidth = 5.41e-10;  // m
const NumType kaasbjergD_HO = 4.1e10, kaasbjergD_KK = 2.6e10; // eV / m

// ---------------- K <-> Q couplings (AHC chain 1, q = Q - K) ----------------
struct IvChannel { NumType D; NumType hw; const char *name; };
const std::vector<IvChannel> KQ_near = {
    {1.148e10, 0.0222, "LA"}, {1.436e10, 0.0427, "LO"}, {1.037e10, 0.0493, "HP"},
    {0.567e10, 0.0455, "TA"}, {0.220e10, 0.0165, "TAlow"}};

// ---------------- Q-valley couplings (AHC chain 1b, k = Q1 = 0.515 Gamma-K) ----------------
const NumType XiL_Q = 3.0, XiT_Q = 2.4;  // eV, |q| = 0.12 1/A, longitudinal/transverse average
const std::vector<IvChannel> KQ_far = {   // K -> Q' (three far valleys) and Q -> K'
    {4.95e10, 0.0491, "HP"}, {4.64e10, 0.0285, "LA"}, {2.86e10, 0.0439, "LO"}, {1.56e10, 0.0190, "TA"}};
const std::vector<IvChannel> QQ60 = {     // 60 deg neighbours
    {3.87e10, 0.0230, "LA"}, {2.07e10, 0.0492, "HP"}, {1.41e10, 0.0454, "TO"}, {0.78e10, 0.0174, "TA"},
    {0.74e10, 0.0424, "LO"}};
const std::vector<IvChannel> QQ120 = {    // 120 deg
    {3.15e10, 0.0428, "LO"}, {2.97e10, 0.0491, "HP"}, {1.93e10, 0.0285, "LA"}};
const std::vector<IvChannel> QQ180 = {    // 180 deg
    {5.09e10, 0.0285, "LA"}, {2.92e10, 0.0225, "TA"}, {2.43e10, 0.0474, "HP"}, {1.44e10, 0.0406, "LO"}};

// ---------------- valley bookkeeping ----------------
// valley 0: K type, six corners (K at even, K' at odd index, 60 deg apart);
// valley 1: Q_low, three valleys at 0, 120, 240 deg (Q_0, Q_2, Q_4 of the six);
// valley 2: Q_high, three valleys at 60, 180, 300 deg (Q_1, Q_3, Q_5).
const SubValleyAssignmentType K_toNextSub = {{1}, {2}, {3}, {4}, {5}, {0}};
const SubValleyAssignmentType K_toAllQ3 = {{0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2}};
const SubValleyAssignmentType Qlow_toK_near = {{0}, {2}, {4}};   // Q_2j -> corner 2j
const SubValleyAssignmentType Qlow_toK_far = {{3}, {5}, {1}};    // Q_2j -> corner 2j + 3
const SubValleyAssignmentType Qhigh_toK_near = {{1}, {3}, {5}};  // Q_2j+1 -> corner 2j + 1
const SubValleyAssignmentType Qhigh_toK_far = {{4}, {0}, {2}};   // Q_2j+1 -> corner 2j + 4
const SubValleyAssignmentType Q3_toOtherTwo = {{1, 2}, {0, 2}, {0, 1}};   // 120 deg, same set
const SubValleyAssignmentType Qlow_toHigh60 = {{0, 2}, {1, 0}, {2, 1}};   // Q_2j -> Q_2j+-1
const SubValleyAssignmentType Qlow_toHigh180 = {{1}, {2}, {0}};           // Q_2j -> Q_2j+3
const SubValleyAssignmentType Qhigh_toLow60 = {{0, 1}, {1, 2}, {2, 0}};   // Q_2j+1 -> Q_2j, Q_2j+2
const SubValleyAssignmentType Qhigh_toLow180 = {{2}, {0}, {1}};           // Q_2j+1 -> Q_2j+4
const SubValleyAssignmentType Q3_toSameSub = {{0}, {1}, {2}};
const SubValleyAssignmentType K_toSameSub = {{0}, {1}, {2}, {3}, {4}, {5}};
const SizeType nrValleyTypes = 3;

template <class DerivedParticleType>
void addValleys(std::unique_ptr<DerivedParticleType> &particleType) {
  particleType->addValley(std::make_unique<isotropValleyType>(
      relEffMassK, particleType->getMass(), 6, alphaK));
  const NumType a60 = constants::pi / 3.;
  std::vector<NumType> anglesLow = {0, 2 * a60, 4 * a60};
  std::vector<NumType> anglesHigh = {a60, 3 * a60, 5 * a60};
  particleType->addValley(std::make_unique<anisotropValleyType>(
      relEffMassQl, relEffMassQt, particleType->getMass(), 3, alphaQ, anglesLow, energyQK));
  particleType->addValley(std::make_unique<anisotropValleyType>(
      relEffMassQl, relEffMassQt, particleType->getMass(), 3, alphaQ, anglesHigh,
      energyQK + deltaQ));
}

template <class DerivedParticleType, class... TParam>
void helperAddZeroOrderMech(std::unique_ptr<DerivedParticleType> &particleType,
                            const std::vector<int> &idxRegions, TParam... param) {
  particleType->addScatterMechanism(idxRegions, std::make_unique<iv0AbScatterMech>(param...));
  particleType->addScatterMechanism(idxRegions, std::make_unique<iv0EmScatterMech>(param...));
}

/// intravalley optical mode with a tabulated |g(q)| (2D Froehlich class), unscreened or
/// screened by the free carriers (screening wavevector qs)
template <class DerivedParticleType>
void helperAddTabulatedOptical(std::unique_ptr<DerivedParticleType> &particleType,
                               const std::vector<int> &idxRegions, SizeType valley, NumType hw,
                               NumType g0, const std::vector<std::pair<NumType, NumType>> &ff,
                               NumType temperature, const std::string &name, NumType qs = 0) {
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<froehlichAb>(valley, hw, g0 * sqrtCellArea, ff, temperature, name, qs));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<froehlichEm>(valley, hw, g0 * sqrtCellArea, ff, temperature, name, qs));
}

template <class DerivedParticleType>
void addAcousticScatterMechanisms(std::unique_ptr<DerivedParticleType> &particleType,
                                  const std::vector<int> &idxRegions, NumType temperature) {
#if MOS2_TV_ACOUSTIC_TAB
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acTabScatterMech>(0, xiTA_K_q, rho, vSoundTA, temperature, "KTA"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acTabScatterMech>(0, xiLA_K_q, rho, vSoundLA, temperature, "KLA"));
#else
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acScatterMech>(0, XiTA_K, rho, vSoundTA, temperature, "KTA"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acScatterMech>(0, XiLA_K, rho, vSoundLA, temperature, "KLA"));
#endif
  for (SizeType v : {SizeType(1), SizeType(2)}) {
    std::string s = v == 1 ? "Qlow" : "Qhigh";
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<acScatterMech>(v, XiT_Q, rho, vSoundTA, temperature, s + "TA"));
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<acScatterMech>(v, XiL_Q, rho, vSoundLA, temperature, s + "LA"));
  }
}

/// zero-order optical couplings: intravalley homopolar (K, Q), K-K', K<->Q, Q-Q'
template <class DerivedParticleType>
void addZeroOrderIntervalleyScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {
  const NumType invSqrt2 = 1. / std::sqrt(2.);
#if MOS2_TV_LDA_OPTICAL
  helperAddZeroOrderMech(particleType, idxRegions, 0, 0, kaasbjergD_HO, rho, temperature, hw_HO,
                         K_toSameSub, "KHomoPolar");
  helperAddZeroOrderMech(particleType, idxRegions, 1, 1, kaasbjergD_HO, rho, temperature, hw_HO,
                         Q3_toSameSub, "QlowHomoPolar");
  helperAddZeroOrderMech(particleType, idxRegions, 2, 2, kaasbjergD_HO, rho, temperature, hw_HO,
                         Q3_toSameSub, "QhighHomoPolar");
  helperAddZeroOrderMech(particleType, idxRegions, 0, 0, kaasbjergD_KK, rho, temperature, hw_KK,
                         K_toNextSub, "KtoKprime");
#else
  // intravalley homopolar A1' with the first-principles q dependence
  helperAddTabulatedOptical(particleType, idxRegions, 0, hw_HO, g0_A1K, ff_A1K, temperature, "KHomoPolar");
  helperAddTabulatedOptical(particleType, idxRegions, 1, hw_HO, g0_A1Q, ff_A1Q, temperature, "QlowHomoPolar");
  helperAddTabulatedOptical(particleType, idxRegions, 2, hw_HO, g0_A1Q, ff_A1Q, temperature, "QhighHomoPolar");
  // K -> K' (one final valley)
  helperAddZeroOrderMech(particleType, idxRegions, 0, 0, D_KK, rho, temperature, hw_KK,
                         K_toNextSub, "KtoKprime");
#endif
  // K -> Q: spin average over the two sets (near and far channels into both, D / sqrt 2)
  for (const auto &c : KQ_near) {
    helperAddZeroOrderMech(particleType, idxRegions, 0, 1, c.D * invSqrt2, rho, temperature, c.hw,
                           K_toAllQ3, std::string("KtoQlow") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 0, 2, c.D * invSqrt2, rho, temperature, c.hw,
                           K_toAllQ3, std::string("KtoQhigh") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 1, 0, c.D, rho, temperature, c.hw,
                           Qlow_toK_near, std::string("QlowtoK") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 2, 0, c.D, rho, temperature, c.hw,
                           Qhigh_toK_near, std::string("QhightoK") + c.name);
  }
  for (const auto &c : KQ_far) {
    helperAddZeroOrderMech(particleType, idxRegions, 0, 1, c.D * invSqrt2, rho, temperature, c.hw,
                           K_toAllQ3, std::string("KtoQlowfar") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 0, 2, c.D * invSqrt2, rho, temperature, c.hw,
                           K_toAllQ3, std::string("KtoQhighfar") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 1, 0, c.D, rho, temperature, c.hw,
                           Qlow_toK_far, std::string("QlowtoKprime") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 2, 0, c.D, rho, temperature, c.hw,
                           Qhigh_toK_far, std::string("QhightoKprime") + c.name);
  }
  // Q -> Q': 120 deg within a spin set, 60 and 180 deg between the sets
  for (const auto &c : QQ120) {
    helperAddZeroOrderMech(particleType, idxRegions, 1, 1, c.D, rho, temperature, c.hw,
                           Q3_toOtherTwo, std::string("QlowtoQlow120") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 2, 2, c.D, rho, temperature, c.hw,
                           Q3_toOtherTwo, std::string("QhightoQhigh120") + c.name);
  }
  for (const auto &c : QQ60) {
    helperAddZeroOrderMech(particleType, idxRegions, 1, 2, c.D, rho, temperature, c.hw,
                           Qlow_toHigh60, std::string("QlowtoQhigh60") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 2, 1, c.D, rho, temperature, c.hw,
                           Qhigh_toLow60, std::string("QhightoQlow60") + c.name);
  }
  for (const auto &c : QQ180) {
    helperAddZeroOrderMech(particleType, idxRegions, 1, 2, c.D, rho, temperature, c.hw,
                           Qlow_toHigh180, std::string("QlowtoQhigh180") + c.name);
    helperAddZeroOrderMech(particleType, idxRegions, 2, 1, c.D, rho, temperature, c.hw,
                           Qhigh_toLow180, std::string("QhightoQlow180") + c.name);
  }
}

/// first-order K-K' (TA, LA, TO at K) and intravalley TO(Gamma) couplings, valley 0 only
template <class DerivedParticleType>
void addFirstOrderIntervalleyScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {
#if MOS2_TV_FIRST_ORDER
  auto add = [&](NumType D1, NumType hw, const char *name) {
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<iv1EmScatterMech>(0, D1, rho, temperature, hw, name));
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<iv1AbScatterMech>(0, D1, rho, temperature, hw, name));
  };
  add(Iv1DefPotTO, hw_TO_K, "TO");
  add(Iv1DefPotTOGamma, phononEnergyLOGamma, "TOGamma");
  add(Iv1DefPotTA, hw_TA_K, "TA");
  add(Iv1DefPotLA, hw_LA_K, "LA");
  add(Iv1DefPotA1K, hw_A1_K, "A1K");
#else
  (void)particleType; (void)idxRegions; (void)temperature;
#endif
}

/// polar-optical (LO) coupling for all valley types, tabulated 2D form, free-carrier screened
template <class DerivedParticleType>
void addFroehlichScatterMechanisms(std::unique_ptr<DerivedParticleType> &particleType,
                                   const std::vector<int> &idxRegions, NumType temperature,
                                   NumType carrierDensity = 0, NumType envPermittivity = 1) {
  NumType qs = twoDStaticScreeningWavevector<NumType>(carrierDensity, temperature,
                                                       envPermittivity,
                                                       relEffMassK * constants::me, 4);
  for (SizeType v = 0; v < nrValleyTypes; ++v) {
#if MOS2_TV_LDA_OPTICAL
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<froehlichAb>(v, phononEnergyLOGamma, kaasbjergFroehlichConst,
                                                  kaasbjergEffectiveWidth, temperature, "", qs));
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<froehlichEm>(v, phononEnergyLOGamma, kaasbjergFroehlichConst,
                                                  kaasbjergEffectiveWidth, temperature, "", qs));
#else
    helperAddTabulatedOptical(particleType, idxRegions, v, phononEnergyLOGamma, g0_LO, ff_LO,
                              temperature, "", qs);
#endif
  }
}

/// No separate piezoelectric term in this kernel: the computed acoustic deformation potentials
/// contain the long-range (piezoelectric) part (Supplementary Note S2). Kept as an empty call so
/// that the drivers can switch kernels by namespace.
template <class DerivedParticleType>
void addPiezoelectricScatterMechanisms(std::unique_ptr<DerivedParticleType> &, const std::vector<int> &,
                                       NumType, NumType = 0, NumType = 1) {}

template <class DerivedParticleType>
void addChargedImpurityScatterMechanism(std::unique_ptr<DerivedParticleType> &particleType,
                                        const std::vector<int> &idxRegions, NumType temperature,
                                        NumType impurityDensity, NumType carrierDensity,
                                        NumType epsAvg = 1, NumType remoteDistance = 0,
                                        NumType chargeNumber = 1) {
  NumType qs = twoDStaticScreeningWavevector<NumType>(carrierDensity, temperature, epsAvg,
                                                       relEffMassK * constants::me, 4);
  for (SizeType v = 0; v < nrValleyTypes; ++v)
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
  for (SizeType v = 0; v < nrValleyTypes; ++v)
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
  for (SizeType v = 0; v < nrValleyTypes; ++v)
    for (NumType wSO : {phononEnergy1, phononEnergy2}) {
      particleType->addScatterMechanism(
          idxRegions, std::make_unique<remoteSOMech>(v, wSO, dPerMode, remoteDistance,
                                                     temperature, false, qs));
      particleType->addScatterMechanism(
          idxRegions, std::make_unique<remoteSOMech>(v, wSO, dPerMode, remoteDistance,
                                                     temperature, true, qs));
    }
}

} // namespace MoS2TwoValley

#endif // MOS2_PARAMETER_TWO_VALLEY_HPP
