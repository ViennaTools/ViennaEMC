#include <ScatterMechanisms/emcAcousticScatterMechanism.hpp>
#include <ScatterMechanisms/emcCoulombScatterMechanism.hpp>
#include <ScatterMechanisms/emcFirstOrderInterValleyScatterMechanism.hpp>
#include <ScatterMechanisms/emcZeroOrderInterValleyScatterMechanism.hpp>
#include <ValleyTypes/emcNonParabolicAnistropValley.hpp>
#include <emcMaterial.hpp>
#include <emcUtil.hpp>

#include <map>

using NumType = double;
using ValleyType = emcNonParabolicAnisotropValley<NumType>;
using Acoustic = emcAcousticScatterMechanism<NumType>;
using ZeroIvAb = emcZeroOrderInterValleyAbsorptionScatterMechanism<NumType>;
using ZeroIvEm = emcZeroOrderInterValleyEmissionScatterMechanism<NumType>;
using FirstIvAb = emcFirstOrderInterValleyAbsorptionScatterMechanism<NumType>;
using FirstIvEm = emcFirstOrderInterValleyEmissionScatterMechanism<NumType>;

namespace Silicon {
//! material characteristics
static const NumType epsR = 11.8;
const NumType rho = 2329.;    // in kg / m3
const NumType Ni = 1.45e16;   // in 1. / m3
const NumType vSound = 9040;  // in m / s
const NumType bandgap = 1.15; // in eV

//! \brief X-valley characteristics.
//! X - Valley consists of 6 conduction
//! bands that are grouped in pairs.
const NumType relEffMassLong = 0.916;  // Vasileska et al.
const NumType relEffMassTrans = 0.196; // Vasileska et al.
const SizeType degeneracy = 3;
const NumType alpha = 0.5; // in 1 / eV

//! characteristics for scattering (Vasileska et al.)
const NumType sigmaAc = 9;            // in eV
const NumType defPot0f = 5.23e10;     // in eV / m
const NumType defPot0g = 5.23e10;     // in eV / m
const NumType defPot1f = 2.5;         // in eV / m
const NumType defPot1g = 4.;          // in eV / m
const NumType phononEnergy0f = 0.06;  // in eV
const NumType phononEnergy0g = 0.06;  // in eV
const NumType phononEnergy1f = 0.023; // in eV
const NumType phononEnergy1g = 0.018; // in eV

//! maps that relate the initial to the final valleys for
//! intervalley scattering.
const std::map<SizeType, std::vector<SizeType>> finalGValleys = {
    {0, {0}}, {1, {1}}, {2, {2}}};
const std::map<SizeType, std::vector<SizeType>> finalFValleys = {
    {0, {1, 1, 2, 2}}, {1, {0, 0, 2, 2}}, {2, {0, 0, 1, 1}}};

//! Returns Material Silicon
template <class T> emcMaterial<T> getSiliconMaterial() {
  return emcMaterial<T>{epsR, rho, Ni, vSound, bandgap};
}

//! Adds X-Valleys of silicon to the given particleType.
//! As described above the 6 conduction band valleys are
//! paired.
template <class DerivedParticleType>
void addXValley(std::unique_ptr<DerivedParticleType> &particleType) {
  std::unique_ptr<ValleyType> siliconValley(
      new ValleyType({relEffMassLong, relEffMassTrans, relEffMassTrans},
                     particleType->getMass(), degeneracy, alpha));
  siliconValley->setSubValleyEllipseCoordSystem(0, {1, 0, 0}, {0, 1, 0},
                                                {0, 0, 1});
  siliconValley->setSubValleyEllipseCoordSystem(1, {0, 1, 0}, {1, 0, 0},
                                                {0, 0, 1});
  siliconValley->setSubValleyEllipseCoordSystem(2, {0, 0, 1}, {0, 1, 0},
                                                {1, 0, 0});
  particleType->addValley(std::move(siliconValley));
}

//! Adds Acoustic Scattering using the silicon characteristics to the
//! given particleType.
template <class DerivedParticleType, class DeviceType>
void addAcousticScattering(SizeType idxValley,
                           std::unique_ptr<DerivedParticleType> &particleType,
                           const DeviceType &device,
                           const std::vector<int> &idxRegions) {
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<Acoustic>(idxValley, sigmaAc, device));
}

//! Adds Coulomb Scattering using the silicon characteristics to the
//! given particleType.
template <class DerivedParticleType, class DeviceType>
void addCoulombScattering(SizeType idxValley,
                          std::unique_ptr<DerivedParticleType> &particleType,
                          DeviceType &device,
                          const std::vector<int> &idxRegions) {
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<emcCoulombScatterMechanism<NumType, DeviceType>>(
          idxValley, epsR, device));
}

//! Adds Zero Order (Optical) Intervalley Scattering using the silicon
//! characteristics to the given particleType. Here we distinguish
//! between f- and g-scattering and between the emission and the
//! absorption of a phonon.
template <class DerivedParticleType, class DeviceType>
void addZeroOrderInterValleyScattering(
    SizeType idxValley, std::unique_ptr<DerivedParticleType> &particleType,
    const DeviceType &device, const std::vector<int> &idxRegions) {
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<ZeroIvAb>("F", idxValley, finalFValleys,
                                             defPot0f, phononEnergy0f, device));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<ZeroIvEm>("F", idxValley, finalFValleys,
                                             defPot0f, phononEnergy0f, device));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<ZeroIvAb>("G", idxValley, finalGValleys,
                                             defPot0g, phononEnergy0g, device));
  particleType->addScatterMechanism(
      idxRegions, std::make_unique<ZeroIvEm>("G", idxValley, finalGValleys,
                                             defPot0g, phononEnergy0g, device));
}

//! Adds First Order (Optical) Intervalley Scattering using the silicon
//! characteristics to the given particleType. Here we distinguish
//! between f- and g-scattering and between the emission and the
//! absorption of a phonon.
template <class DerivedParticleType, class DeviceType>
void addFirstOrderInterValleyScattering(
    SizeType idxValley, std::unique_ptr<DerivedParticleType> &particleType,
    const DeviceType &device, const std::vector<int> &idxRegions) {
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<FirstIvAb>("F", idxValley, finalFValleys, defPot1f,
                                  phononEnergy1f, device));
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<FirstIvEm>("F", idxValley, finalFValleys, defPot1f,
                                  phononEnergy1f, device));
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<FirstIvAb>("G", idxValley, finalGValleys, defPot1g,
                                  phononEnergy1g, device));
  particleType->addScatterMechanism(
      idxRegions,
      std::make_unique<FirstIvEm>("G", idxValley, finalGValleys, defPot1g,
                                  phononEnergy1g, device));
}
} // namespace Silicon
