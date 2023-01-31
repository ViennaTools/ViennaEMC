#ifndef EMC_PARTICLE_HPP
#define EMC_PARTICLE_HPP

#include <array>

#include <emcUtil.hpp>

/// class that holds the information that is needed for a
/// simulated particle
template <class T> struct emcParticle {
  std::array<T, 3> k = {0, 0, 0}; /**< wave-vector */
  T energy = 0;                   /**< energy [eV] */
  T tau = 1;                      /**< rem. free flight time [s] */
  T grainTau = 1; /**< free flight time between grain scattering [s] */
  SizeType valley = 0;
  SizeType subValley = 0;
  SizeType region = 0;
};

#endif // EMC_PARTICLE_HPP