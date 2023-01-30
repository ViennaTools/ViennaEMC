#ifndef PARTICLE_DRIFT_HPP
#define PARTICLE_DRIFT_HPP

#include <emcConstants.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

/// function lets particle drift for time dt under the influence of the
/// given force.
/// NOTE: doesn't adapt position if particle leaves domain!
template <class T, class ValleyType, SizeType Dim>
void drift(T dt, emcParticle<T> &part, const ValleyType *valley,
           std::array<T, Dim> &pos, const std::array<T, 3> &force) {
  auto vogtFactor = valley->getVogtTransformationFactor();
  auto kOldECS = valley->transformToEllipseCoord(part.subValley, part.k);
  auto kNewECS = kOldECS;
  auto forceECS = valley->transformToEllipseCoord(part.subValley, force);
  // calculate dK = - force * dt / hbar * vogtFactor in ECS
  for (SizeType idxDim = 0; idxDim < 3; idxDim++) {
    kNewECS[idxDim] +=
        forceECS[idxDim] * dt * vogtFactor[idxDim] / constants::hbar;
  }
  part.k = valley->transformToDeviceCoord(part.subValley, kNewECS);
  part.energy = valley->getEnergy(part.k);
  // calculate new position (leapfrog scheme) in ECS
  std::array<T, 3> dPos;
  for (SizeType idxDim = 0; idxDim < 3; idxDim++) {
    T avgK = (kNewECS[idxDim] + kOldECS[idxDim]) / 2;
    dPos[idxDim] = constants::hbar * vogtFactor[idxDim] * avgK * dt /
                   valley->getEffMassCond(part.energy);
  }
  dPos = valley->transformToDeviceCoord(part.subValley, dPos);
  std::transform(pos.begin(), pos.end(), dPos.begin(), pos.begin(),
                 std::plus<>{});
}

/// checks position of particle, if it is out of bounds particle is
/// either marked to be removed (if it left through ohmic contact) or
/// it is reflected (if it left through gate / artificial boundary)
template <class T, class DeviceType, SizeType Dim, class ParticleType>
bool handleParticleAtBoundary(emcParticle<T> &part, std::array<T, Dim> &pos,
                              ParticleType *particleType, DeviceType &device,
                              emcRNG &rng) {
  bool remove = false;
  if (device.isOutOfBounds(pos)) {
    // calculate nearest boundary position
    std::array<T, Dim> posBoundary;
    auto maxPos = device.getMaxPos();
    std::transform(pos.begin(), pos.end(), maxPos.begin(), posBoundary.begin(),
                   [](const T &pos, const T &max) {
                     return std::max(0., std::min(pos, max));
                   });
    // get type of nearest boundary, handle particle at that kind of boundary
    std::array<SizeType, Dim> coordBoundary = device.posToCoord(posBoundary);
    auto &surface = device.getSurface();
    if (surface.isOhmicContact(coordBoundary)) {
      remove = true;
      pos = posBoundary;
    } else { // scatter particle
      particleType->scatterParticleAtBoundary(
          surface.getBoundaryPos(coordBoundary), part, device, pos, rng);
    }
  }
  return remove;
}

#endif // PARTICLE_DRIFT_HPP