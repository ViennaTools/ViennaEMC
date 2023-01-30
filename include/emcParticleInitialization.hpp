#ifndef EMC_PARTICLE_INITIALIZATION_HPP
#define EMC_PARTICLE_INITIALIZATION_HPP

#include <array>
#include <random>

#include <emcConstants.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

/// initialize particle positions.
/// random position near coordinate is assigned.
template <class T, SizeType Dim>
std::array<T, Dim> initParticlePos(const std::array<SizeType, Dim> &coord,
                                   const std::array<SizeType, Dim> &extent,
                                   const std::array<T, Dim> &spacing,
                                   emcRNG &rng) {
  std::array<T, Dim> posPart;
  std::uniform_real_distribution<T> distr(0., 1.);
  for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
    if (coord[idxDim] == extent[idxDim] - 1)
      posPart[idxDim] = (coord[idxDim] - distr(rng) * 0.5) * spacing[idxDim];
    else if (coord[idxDim] == 0)
      posPart[idxDim] = distr(rng) * 0.5 * spacing[idxDim];
    else
      posPart[idxDim] = (coord[idxDim] + distr(rng) - 0.5) * spacing[idxDim];
  }
  return posPart;
}

/// initialize energy and wave-vector of particle.
/// energy is assigned using maxwellian distribution.
/// direction of wave-vector is assigned randomly.
template <class T, SizeType Dim, template <class, SizeType> class DeviceType,
          class ValleyType>
void initParticleKSpaceMaxwellian(emcParticle<T> &part,
                                  const std::array<SizeType, Dim> &coord,
                                  const DeviceType<T, Dim> &device,
                                  const ValleyType *valley, emcRNG &rng) {
  std::uniform_real_distribution<T> dist(0., 1.);
  std::uniform_real_distribution<T> distForLog(1e-6, 1.);
  part.energy = -1.5 * device.getThermalVoltage() * std::log(distForLog(rng));
  part.k = initRandomDirection(valley->getNormWaveVec(part.energy), dist(rng),
                               dist(rng));
  auto extent = device.getGridExtent();
  for (SizeType idxDir = 0; idxDir < Dim; idxDir++) {
    if (((coord[idxDir] == 0) && (part.k[idxDir] < 0)) ||
        ((coord[idxDir] == extent[idxDir] - 1) && (part.k[idxDir] > 0)))
      part.k[idxDir] *= -1;
  }
}

/// used for injection of particles at ohmic contacts.
/// initialize energy and wave-vector of particle.
/// energy is assigned using maxwellian distribution transversal to a near
/// boundary and vel weighted maxwellian normal to boundary.
template <class T, SizeType Dim, template <class, SizeType> class DeviceType,
          class ValleyType>
void initParticleKSpaceVelWeightedMaxwellian(
    emcParticle<T> &part, const std::array<SizeType, Dim> &coord,
    const DeviceType<T, Dim> &device, const ValleyType *valley, emcRNG &rng) {
  initParticleKSpaceMaxwellian(part, coord, device, valley, rng);
  // adapt wave-vec if particle at boundary in this direction
  std::uniform_real_distribution<T> distForLog(1e-6, 1.);
  auto extent = device.getGridExtent();
  for (SizeType idxDir = 0; idxDir < Dim; idxDir++) {
    if (coord[idxDir] == 0 || coord[idxDir] == extent[idxDir] - 1) {
      T factor = (coord[idxDir] == 0) ? 1 : -1;
      part.k[idxDir] =
          factor / constants::hbar *
          std::sqrt(-2 * valley->getEffMassCond() * device.getThermalVoltage() *
                    constants::q * std::log(distForLog(rng)));
    }
    part.energy = valley->getEnergy(part.k);
  }
}

#endif // EMC_PARTICLE_INITIALIZATION_HPP