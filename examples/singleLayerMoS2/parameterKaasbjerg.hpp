#include <ScatterMechanisms/emcAcousticSingleLayerScatterMechanism.hpp>
#include <ScatterMechanisms/emcFirstOrderSingleLayerIntervalleyScatterMechanism.hpp>
#include <ScatterMechanisms/emcFroehlichInteractionSingleLayer.hpp>
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

template <class DerivedParticleType>
void addFroehlichScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<froehlichAb>(0, phononEnergyLOGamma, couplingConst,
                                    effectiveWidth, temperature));

  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<froehlichEm>(0, phononEnergyLOGamma, couplingConst,
                                    effectiveWidth, temperature));
}

} // namespace MoS2Kaasbjerg