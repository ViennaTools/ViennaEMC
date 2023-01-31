#ifndef EMC_CONSTANT_SURFACE_SCATTER_MECHANISM_HPP
#define EMC_CONSTANT_SURFACE_SCATTER_MECHANISM_HPP

#include <SurfaceScatterMechanisms/emcSurfaceScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

/**
 * @brief Assumes constant probability of diffuse and specular scattering.
 *
 * @param specularityParam constant parameter between 0 and 1 representing
 * the probability of a particle being scattered specularly (if param higher,
 * the roughness of the surface is higher)
 * 1 ... perfectly smooth surface (only specular scattering)
 * 0 ... perfectly rough surface (only diffuse scattering)
 */
template <class T, class DeviceType, SizeType Dim = DeviceType::Dimension>
class emcConstantSurfaceScatterMechanism
    : public emcSurfaceScatterMechanism<T, DeviceType> {

  T specularityParam;

public:
  emcConstantSurfaceScatterMechanism() = delete;

  emcConstantSurfaceScatterMechanism(T inSpecularityParam,
                                     std::array<T, Dim> inMaxPos)
      : emcSurfaceScatterMechanism<T, DeviceType>(inMaxPos),
        specularityParam(inSpecularityParam) {}

private:
  T getDiffScatterProb(emcParticle<T> & /*particle*/) const {
    return 1 - specularityParam;
  }

  void scatterParticleDiffusely(emcParticle<T> &particle,
                                std::array<T, DeviceType::Dimension> &pos,
                                emcRNG &rng) const {
    T theta = std::asin(std::sqrt(this->dist(rng)));
    T phi = 2 * constants::pi * this->dist(rng);
    this->calculateAndAssignKAndPos(particle, pos, theta, phi);
  }
};

#endif