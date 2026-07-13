#include <ScatterMechanisms/emc2DChargedImpurityScatterMechanism.hpp>
#include <ScatterMechanisms/emcAcousticSingleLayerScatterMechanism.hpp>
#include <ScatterMechanisms/emcFirstOrderSingleLayerIntervalleyScatterMechanism.hpp>
#include <ScatterMechanisms/emcFroehlichInteractionSingleLayer.hpp>
#include <ScatterMechanisms/emcPiezoelectricSingleLayerScatterMechanism.hpp>
#include <ScatterMechanisms/emcRemoteSurfaceOpticalPhononMechanism.hpp>
#include <ScatterMechanisms/emcSurfaceRoughnessScatterMechanism.hpp>
#include <ScatterMechanisms/emcZeroOrderSingleLayerInterValleyScatterMechanism.hpp>
#include <ValleyTypes/emcParabolicIsotropSingleLayerValley.hpp>

/*! \brief Parameters taken from Paper Kaasbjerg et al.
 *
 * DOI of the paper: https://doi.org/10.1103/PhysRevB.85.115317
 *
 * The band structure of the material is approximated with the
 * multi-valley approach within this paper, the parameters for the
 * valleys of interest (only K-valley are taken from Chapter II A).
 *
 * Phonon scattering and Fröhlich interactions are considered as
 * scatter mechanisms.
 *
 * NOTE: This implementation is still in progress. The adaptation of
 * the state after a scattering event with Fröhlich scattering may has
 * to be adapted.
 */
namespace MoS2Kaasbjerg {

using NumType = double;
using ValleyType = emcParabolicIsotropSingleLayerValley<NumType>;
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

const NumType rho = 3.1e-6; // in kg / m2

// characteristics K-Valley (assume isotrop + parabolic)
const NumType relEffMassK = 0.48;
const SizeType degeneracyK = 6;

// characteristics for Acoustic Scattering
const NumType vSoundLA = 6.7e3; // in m / s
const NumType vSoundTA = 4.2e3; // in m / s
const NumType AcDefPotLA = 2.8; // in eV
const NumType AcDefPotTA = 1.6; // in eV

// characteristics for Zero-order Intervalley Scattering
const NumType Iv0DefPotLO = 2.6e10; // in eV / m
const NumType Iv0DefPotHO = 4.1e10; // in eV / m

// characterisics for First-Order Intervalley Scattering
const NumType Iv1DefPotTA = 5.9;      // in eV
const NumType Iv1DefPotLA = 3.9;      // in eV
const NumType Iv1DefPotTOGamma = 4.0; // in eV
const NumType Iv1DefPotTO = 1.9;      // in eV

// used phonon energies in specific valleys
// PAPER: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.165436
const NumType phononEnergyTA = 0.023;      // in eV
const NumType phononEnergyLA = 0.029;      // in eV
const NumType phononEnergyTOK = 0.048;     // in eV
const NumType phononEnergyTOGamma = 0.048; // in eV
const NumType phononEnergyHP = 0.05;       // in eV
const NumType phononEnergyLOK = 0.041;     // in eV
const NumType phononEnergyLOGamma = 0.048; // in eV

// parameter for Froehlich scattering, found in:
// https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.165436
const NumType couplingConst = 0.286 * 1e-10; // in eV / m
const NumType effectiveWidth = 5.41e-10;     // in m

// piezoelectric constant e11 for monolayer MoS2 (Kaasbjerg et al.,
// PRB 87, 235312 (2013), Table I). effectiveWidth doubles as sigma.
const NumType piezoConst = 3.0e-11; // in C / m

// Rytova-Keldysh 2D polarizability (screening) length for monolayer MoS2
// (~40 A; e.g. Berkelbach et al., PRB 88, 045318 (2013)).
const NumType rytovaKeldyshLength = 4.0e-9; // in m

// representative interface-roughness parameters for MoS2 on an oxide
const NumType roughnessAmplitude = 3.0e-10; // RMS amplitude Delta [m] (~3 A)
const NumType roughnessCorrLength = 1.5e-9; // correlation length Lambda [m]

// carrier-to-substrate distance for remote SO phonons (vdW gap + half-thickness)
const NumType soRemoteDistance = 5.0e-10; // in m

/// assumes that only K-Valley is important for low-field mobility.
/// assumes that K-Valley is parabolic and effective mass is isotrop.
/// also simulates only one subvalley of K!
template <class DerivedParticleType>
void addValleys(std::unique_ptr<DerivedParticleType> &particleType) {
  particleType->addValley(
      std::make_unique<ValleyType>(relEffMassK, particleType->getMass(), 1));
}

template <class DerivedParticleType>
void addAcousticScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acScatterMech>(0, AcDefPotTA, rho, vSoundTA,
                                                  temperature, "TA"));

  particleType->addScatterMechanism(
      idxRegions, std::make_unique<acScatterMech>(0, AcDefPotLA, rho, vSoundLA,
                                                  temperature, "LA"));
}

template <class DerivedParticleType>
void addZeroOrderIntervalleyScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv0EmScatterMech>(
                      0, Iv0DefPotLO, rho, temperature, phononEnergyLOK, "LO"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv0AbScatterMech>(
                      0, Iv0DefPotLO, rho, temperature, phononEnergyLOK, "LO"));
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<iv0EmScatterMech>(0, Iv0DefPotHO, rho, temperature,
                                         phononEnergyHP, "HomoPolar"));
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<iv0AbScatterMech>(0, Iv0DefPotHO, rho, temperature,
                                         phononEnergyHP, "HomoPolar"));
}

template <class DerivedParticleType>
void addFirstOrderIntervalleyScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv1EmScatterMech>(
                      0, Iv1DefPotTO, rho, temperature, phononEnergyTOK, "TO"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv1AbScatterMech>(
                      0, Iv1DefPotTO, rho, temperature, phononEnergyTOK, "TO"));

  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<iv1EmScatterMech>(0, Iv1DefPotTOGamma, rho, temperature,
                                         phononEnergyTOGamma, "TOGamma"));
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<iv1AbScatterMech>(0, Iv1DefPotTOGamma, rho, temperature,
                                         phononEnergyTOGamma, "TOGamma"));

  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv1EmScatterMech>(
                      0, Iv1DefPotTA, rho, temperature, phononEnergyTA, "TA"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv1AbScatterMech>(
                      0, Iv1DefPotTA, rho, temperature, phononEnergyTA, "TA"));

  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv1EmScatterMech>(
                      0, Iv1DefPotLA, rho, temperature, phononEnergyLA, "LA"));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv1AbScatterMech>(
                      0, Iv1DefPotLA, rho, temperature, phononEnergyLA, "LA"));
}

/// Adds the polar-optical (Froehlich) scatter mechanisms. If carrierDensity > 0
/// the interaction is screened by the 2D free-carrier gas (see
/// emc2DScreening.hpp); carrierDensity = 0 (default) leaves it unscreened, so
/// the original free-standing behaviour is reproduced exactly.
///
/// @param carrierDensity mobile electron sheet density n [1 / m^2]
/// @param envPermittivity effective relative permittivity screening the carriers
template <class DerivedParticleType>
void addFroehlichScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature,
    NumType carrierDensity = 0, NumType envPermittivity = 1) {
  // K/K' conduction valleys with spin => degeneracy g = g_spin * g_valley = 4
  NumType screeningWavevector = twoDStaticScreeningWavevector<NumType>(
      carrierDensity, temperature, envPermittivity,
      relEffMassK * constants::me, 4);

  particleType->addScatterMechanism(
      idxRegions, std::make_unique<froehlichAb>(0, phononEnergyLOGamma,
                                                couplingConst, effectiveWidth,
                                                temperature, "",
                                                screeningWavevector));

  particleType->addScatterMechanism(
      idxRegions, std::make_unique<froehlichEm>(0, phononEnergyLOGamma,
                                                couplingConst, effectiveWidth,
                                                temperature, "",
                                                screeningWavevector));
}

/// Adds the piezoelectric acoustic scatter mechanisms (TA + LA branches).
/// If carrierDensity > 0 the (long-range) interaction is screened by the 2D
/// free-carrier gas; carrierDensity = 0 (default) leaves it unscreened.
///
/// @param carrierDensity mobile electron sheet density n [1 / m^2]
/// @param envPermittivity effective relative permittivity screening the carriers
template <class DerivedParticleType>
void addPiezoelectricScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature,
    NumType carrierDensity = 0, NumType envPermittivity = 1) {
  // K/K' conduction valleys with spin => degeneracy g = g_spin * g_valley = 4
  NumType screeningWavevector = twoDStaticScreeningWavevector<NumType>(
      carrierDensity, temperature, envPermittivity,
      relEffMassK * constants::me, 4);

  particleType->addScatterMechanism(
      idxRegions, std::make_unique<piezoScatterMech>(
                      0, piezoConst, effectiveWidth, rho, vSoundTA, temperature,
                      "TA", screeningWavevector));

  particleType->addScatterMechanism(
      idxRegions, std::make_unique<piezoScatterMech>(
                      0, piezoConst, effectiveWidth, rho, vSoundLA, temperature,
                      "LA", screeningWavevector));
}

/// Adds elastic scattering off a sheet of charged impurities. This is an
/// OPTIONAL extrinsic mechanism (not part of the intrinsic setup) enabled when
/// modelling supported/doped/adsorbate-covered films. The free-carrier
/// screening (needed for a finite rate) uses the same 2D Thomas-Fermi
/// wavevector as the polar mechanisms, with eps_r = epsAvg.
///
/// @param impurityDensity areal charged-impurity density N_imp [1 / m^2]
/// @param carrierDensity mobile electron sheet density n [1 / m^2] (sets q_s)
/// @param epsAvg average relative permittivity of the environment
/// @param remoteDistance impurity-to-sheet distance d [m] (0 = interface)
/// @param chargeNumber impurity charge number Z
template <class DerivedParticleType>
void addChargedImpurityScatterMechanism(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature,
    NumType impurityDensity, NumType carrierDensity, NumType epsAvg = 1,
    NumType remoteDistance = 0, NumType chargeNumber = 1) {
  NumType screeningWavevector = twoDStaticScreeningWavevector<NumType>(
      carrierDensity, temperature, epsAvg, relEffMassK * constants::me, 4);

  particleType->addScatterMechanism(
      idxRegions, std::make_unique<chargedImpurityMech>(
                      0, impurityDensity, epsAvg, screeningWavevector,
                      rytovaKeldyshLength, remoteDistance, chargeNumber));
}

/// Adds interface-roughness (Prange-Nee / Ando) scattering. OPTIONAL extrinsic
/// mechanism (not part of the intrinsic setup) for supported/gated films. The
/// normal field is estimated from the carrier self-field F = e n / (2 eps0
/// eps_avg) unless overridden; an extra gate/depletion field can be added.
///
/// @param carrierDensity mobile electron sheet density n [1 / m^2]
/// @param epsAvg average relative permittivity of the environment
/// @param extraField additional normal field (gate/depletion) [V / m]
template <class DerivedParticleType>
void addSurfaceRoughnessScatterMechanism(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature,
    NumType carrierDensity, NumType epsAvg = 1, NumType extraField = 0) {
  NumType screeningWavevector = twoDStaticScreeningWavevector<NumType>(
      carrierDensity, temperature, epsAvg, relEffMassK * constants::me, 4);
  // carrier self-field (half the sheet field) + any external normal field
  NumType effectiveField =
      constants::q * carrierDensity / (2 * constants::eps0 * epsAvg) +
      extraField;

  particleType->addScatterMechanism(
      idxRegions, std::make_unique<surfaceRoughnessMech>(
                      0, effectiveField, roughnessAmplitude, roughnessCorrLength,
                      screeningWavevector));
}

/// Adds remote surface-optical (SO) phonon scattering from a polar substrate /
/// gate dielectric. OPTIONAL extrinsic mechanism (not in the intrinsic setup).
/// Two SO modes are added, each as emission + absorption; the substrate
/// dielectric step is split equally between the modes (a Wang-Mahan
/// decomposition with the intermediate dielectric constant is a refinement).
///
/// Representative substrate parameters (Ma & Jena, PRX 4, 011043 (2014), Tab. I;
/// energies in eV): SiO2  eps_inf=2.5,  eps_0=3.9,  wSO=0.0556, 0.1381;
///                  hBN   eps_inf=4.1,  eps_0=5.09, wSO=0.0931, 0.1791;
///                  HfO2  eps_inf=5.03, eps_0=23.0, wSO=0.0124, 0.0484;
///                  Al2O3 eps_inf=3.2,  eps_0=12.53,wSO=0.0482, 0.0714.
///
/// @param carrierDensity mobile electron sheet density n [1 / m^2] (sets q_s)
/// @param epsOxInf, epsOx0 high-/low-frequency permittivity of the substrate
/// @param phononEnergy1, phononEnergy2 the two SO mode energies [eV]
/// @param epsEnv permittivity of the medium on the other side (e.g. air = 1)
template <class DerivedParticleType>
void addRemoteSurfaceOpticalPhonon(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature,
    NumType carrierDensity, NumType epsOxInf, NumType epsOx0,
    NumType phononEnergy1, NumType phononEnergy2, NumType epsEnv = 1,
    NumType remoteDistance = soRemoteDistance) {
  NumType qs = twoDStaticScreeningWavevector<NumType>(
      carrierDensity, temperature, epsEnv, relEffMassK * constants::me, 4);
  // total dielectric step, split equally between the two SO modes
  NumType dPerMode =
      0.5 * (1. / (epsOxInf + epsEnv) - 1. / (epsOx0 + epsEnv));

  for (NumType wSO : {phononEnergy1, phononEnergy2}) {
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<remoteSOMech>(
                        0, wSO, dPerMode, remoteDistance, temperature,
                        /*emission=*/false, qs)); // absorption
    particleType->addScatterMechanism(
        idxRegions, std::make_unique<remoteSOMech>(
                        0, wSO, dPerMode, remoteDistance, temperature,
                        /*emission=*/true, qs)); // emission
  }
}

} // namespace MoS2Kaasbjerg