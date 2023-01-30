#include <ScatterMechanisms/emcAcousticSingleLayerScatterMechanism.hpp>
#include <ScatterMechanisms/emcZeroOrderSingleLayerInterValleyScatterMechanism.hpp>
#include <ValleyTypes/emcNonParabolicAnisotropSingleLayerValley.hpp>
#include <ValleyTypes/emcNonParabolicIsotropSingleLayerValley.hpp>
#include <emcConstants.hpp>

/**
 * @brief Parameters taken from Paper Pilotto et al.
 *
 * DOI of the paper: https://doi.org/10.1016/j.sse.2022.108295
 *
 * The band structure of the material is approximated with the
 * multi-valley approach within this paper, the parameter for the
 * valleys of interest (at K- and Q-points) are taken from table 1.
 *
 * As scattering mechanisms, scattering with intrinsic phonons is
 * implemented, the other scatter mechanisms aren't implemented yet.
 * As in the paper the parameter for intrinsic phonon scattering
 * are taken from the paper Li et al.
 */
namespace MoS2Pilotto {
using NumType = double;
using SubValleyAssignmentType = std::vector<std::vector<SizeType>>;

using isotropValleyType = emcNonParabolicIsotropSingleLayerValley<NumType>;
using anisotropValleyType = emcNonParabolicAnisotropSingleLayerValley<NumType>;

using acScatterMech = emcAcousticSingleLayerMechanism<NumType>;
using iv0AbScatterMech =
    emcZeroOrderSingleLayerInterValleyAbsorptionScatterMechanism<NumType>;
using iv0EmScatterMech =
    emcZeroOrderSingleLayerInterValleyEmissionScatterMechanism<NumType>;

const NumType rho = 3.1e-6;    // in kg / m2
const NumType energyQK = 0.16; // in eV
const NumType vSound = 6.6e3;  // in m / s

// characteristics for intra- + intervalley scattering
// average over LA + TA at specific symmetry point
const NumType acPhonEnergyGamma = 0;                 // in eV
const NumType acPhonEnergyK = (23.1 + 29.1) / 2000.; // in eV
const NumType acPhonEnergyM = (19.2 + 29.2) / 2000.; // in eV
const NumType acPhonEnergyQ = (17.9 + 23.6) / 2000.; // in eV
// average over TO + LO + A1 at specific symmetry point
const NumType opPhonEnergyGamma = (48.6 + 48.9 + 50.9) / 3000.; // in eV
const NumType opPhonEnergyK = (46.4 + 42.2 + 51.9) / 3000.;     // in eV
const NumType opPhonEnergyM = (48.2 + 44.3 + 50.1) / 3000.;     // in eV
const NumType opPhonEnergyQ = (48.0 + 44.2 + 52.2) / 3000.;     // in eV

// vectors that assign to each initial subvalley (first idx) the possible
// final subvalleys
const SubValleyAssignmentType toSameSub = {{0}, {1}, {2}, {3}, {4}, {5}};
const SubValleyAssignmentType toNextSub = {{1}, {2}, {3}, {4}, {5}, {0}};
const SubValleyAssignmentType toNbrSub = {{1, 5}, {0, 2}, {1, 3},
                                          {2, 4}, {3, 5}, {4, 0}};
const SubValleyAssignmentType toNbrNbrSub = {{2, 4}, {3, 5}, {4, 0},
                                             {5, 1}, {0, 2}, {1, 3}};
const SubValleyAssignmentType toOppositeSub = {{3}, {4}, {5}, {0}, {1}, {2}};
const SubValleyAssignmentType toSlashSub = {{1, 3, 5}, {0, 2, 4}, {1, 3, 5},
                                            {0, 2, 4}, {1, 3, 5}, {0, 2, 4}};
const SubValleyAssignmentType toSameKindSub = {{0, 2, 4}, {1, 3, 5}, {0, 2, 4},
                                               {1, 3, 5}, {0, 2, 4}, {0, 2, 4}};

template <class DerivedParticleType>
void addValleys(std::unique_ptr<DerivedParticleType> &particleType) {
  // K - Valleys (non parabolic, isotrop valleys)
  particleType->addValley(std::make_unique<isotropValleyType>(
      0.47, particleType->getMass(), 6, 0.94));

  // Q - Valleys (non parabolic, anisotrop valleys)
  NumType angle60dg = constants::pi / 3.;
  std::vector<NumType> rotAngles = {
      0, angle60dg, 2 * angle60dg, constants::pi, 4 * angle60dg, 5 * angle60dg};
  particleType->addValley(std::make_unique<anisotropValleyType>(
      0.54, 1.14, particleType->getMass(), 6, 1.16, rotAngles, energyQK));
}

template <class DerivedParticleType>
void addAcousticScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {

  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<acScatterMech>(0, 4.5, rho, vSound, temperature, "K"));
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<acScatterMech>(1, 2.8, rho, vSound, temperature, "Q"));
}

template <class DerivedParticleType, class... TParam>
void helperAddZeroOrderMech(std::unique_ptr<DerivedParticleType> &particleType,
                            const std::vector<int> &idxRegions,
                            TParam... param) {
  // Absorption
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv0AbScatterMech>(param...));
  // Emission
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<iv0EmScatterMech>(param...));
}

template <class DerivedParticleType>
void addZeroOrderIntervalleyScatterMechanisms(
    std::unique_ptr<DerivedParticleType> &particleType,
    const std::vector<int> &idxRegions, NumType temperature) {

  // K -> K (Gamma - phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 0, 0, 5.8e10, rho,
                         temperature, opPhonEnergyGamma, toSameSub,
                         "KOpGammaToK");

  // K -> K' (K-phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 0, 0, 1.4e10, rho,
                         temperature, acPhonEnergyK, toNextSub, "KAcKToKSlash");
  helperAddZeroOrderMech(particleType, idxRegions, 0, 0, 2.0e10, rho,
                         temperature, opPhonEnergyK, toNextSub, "KOpKToKSlash");

  // K -> Q (Q-phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 0, 1, 0.93e9, rho,
                         temperature, acPhonEnergyQ, toSameKindSub, "KAcQToQ");
  helperAddZeroOrderMech(particleType, idxRegions, 0, 1, 1.9e10, rho,
                         temperature, opPhonEnergyQ, toSameKindSub, "KOpQToQ");

  // K -> Q' (M-phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 0, 1, 4.4e10, rho,
                         temperature, acPhonEnergyM, toSlashSub, "KAcMToQ");
  helperAddZeroOrderMech(particleType, idxRegions, 0, 1, 5.6e10, rho,
                         temperature, opPhonEnergyM, toSlashSub, "KOpMToQ");

  // Q1 -> Q1 (Gamma - phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 1, 1, 7.1e10, rho,
                         temperature, opPhonEnergyGamma, toSameSub, "QOpQToQ1");

  // Q1 -> Q2 / Q6 (Q - phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 1, 1, 2.1e10, rho,
                         temperature, acPhonEnergyQ, toNbrSub, "QAcQToQ2");
  helperAddZeroOrderMech(particleType, idxRegions, 1, 1, 4.8e10, rho,
                         temperature, opPhonEnergyQ, toNbrSub, "QOpQToQ2");

  // Q -> Q3 / Q5 (M - phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 1, 1, 2.0e10, rho,
                         temperature, acPhonEnergyM, toNbrNbrSub, "QAcMToQ3");
  helperAddZeroOrderMech(particleType, idxRegions, 1, 1, 4.0e10, rho,
                         temperature, opPhonEnergyM, toNbrNbrSub, "QOpMToQ3");

  // Q -> Q4 (K - phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 1, 1, 4.8e10, rho,
                         temperature, acPhonEnergyK, toOppositeSub, "QAcKToQ4");
  helperAddZeroOrderMech(particleType, idxRegions, 1, 1, 6.5e10, rho,
                         temperature, opPhonEnergyK, toOppositeSub, "QOpKToQ4");

  // Q -> K (Q - phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 1, 0, 1.5e10, rho,
                         temperature, acPhonEnergyQ, toSameSub, "QAcQToK");
  helperAddZeroOrderMech(particleType, idxRegions, 1, 0, 2.4e10, rho,
                         temperature, opPhonEnergyQ, toSameSub, "QOpQToK");

  // Q -> K' (M - phonons)
  helperAddZeroOrderMech(particleType, idxRegions, 1, 0, 4.4e10, rho,
                         temperature, acPhonEnergyM, toNextSub, "QAcMToKSlash");
  helperAddZeroOrderMech(particleType, idxRegions, 1, 0, 6.6e10, rho,
                         temperature, opPhonEnergyM, toNextSub, "QOpMToKSlash");
}

} // namespace MoS2Pilotto