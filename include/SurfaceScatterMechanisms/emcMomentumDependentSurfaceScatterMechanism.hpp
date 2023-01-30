#ifndef EMC_MOMENTUM_DEPENDENT_SURFACE_SCATTER_MECHANISM_HPP
#define EMC_MOMENTUM_DEPENDENT_SURFACE_SCATTER_MECHANISM_HPP

#include <SurfaceScatterMechanisms/emcSurfaceScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

/// @param roughnessHeight RMS roughness height of the surface
template <class T, class DeviceType, SizeType Dim = DeviceType::Dimension>
class emcMomentumDependentSurfaceScatterMechanism
    : public emcSurfaceScatterMechanism<T, DeviceType> {

  T roughnessHeight;

public:
  emcMomentumDependentSurfaceScatterMechanism(T inRoughnessHeight,
                                              std::array<T, Dim> inMaxPos)
      : emcSurfaceScatterMechanism<T, DeviceType>(inMaxPos),
        roughnessHeight(inRoughnessHeight) {}

  T getDiffScatterProb(emcParticle<T> &particle) const {
    return 1 -
           (std::exp(-std::pow(2 * roughnessHeight *
                                   particle.k[this->getIndexFromBoundaryPos()],
                               2)));
  }

private:
  void scatterParticleDiffusely(emcParticle<T> &particle,
                                std::array<T, DeviceType::Dimension> &pos,
                                emcRNG &rng) const {
    T speed = norm(particle.k);
    T theta = solveForTheta(this->dist(rng), roughnessHeight, speed);
    T phi = 2 * constants::pi * this->dist(rng);
    this->calculateAndAssignKAndPos(particle, pos, theta, phi);
  }

  T solveForTheta(T r, T height, T speed) const {
    T c = std::pow(2 * height * speed, 2);
    T e = std::exp(-c);
    T startEstimation = std::sqrt(r * (1 / (1 - e) - 1 / c));
    T error = 1;

    T tol = 1e-10;
    int it = 0;
    int max_it = 10000;

    T x1;
    T x = startEstimation;

    while (error > tol && it < max_it) {
      T cos = std::cos(x);
      T sin = std::sin(x);
      T ecos = std::exp(-c * std::pow(cos, 2));
      T esin = std::exp(-c * std::pow(sin, 2));

      x1 = x - (((e - ecos) / c + std::pow(sin, 2) - r * (1 - (1 - e) / c)) *
                (2 * esin * c * sin * cos * (ecos - 1)) / (e * (c - 1) - 1));

      error = std::fabs(x1 - x);
      x = x1;

      it++;
    }

    return x;
  }
};

#endif