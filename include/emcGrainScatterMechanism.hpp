#ifndef EMC_GRAIN_SCATTER_MECHANISM_HPP
#define EMC_GRAIN_SCATTER_MECHANISM_HPP

#include <emcParticle.hpp>

/*! \brief This class handles the scattering of a particle at grains.
 *
 * It is assumed that the mean rate at which particles scatter at grains
 * is given. Additionally, the transmission probability of a particle through
 * the grain is given.
 *
 * @param transmissionProb probability which determines whether the particle
 * transmits through the grain or is reflected
 * @param scatterRate rate at which particles scatter at grains (in average)
 */
template <class T> struct emcGrainScatterMechanism {
protected:
  mutable std::uniform_real_distribution<T> dist;
  T transmissionProb;
  T scatterRate;

public:
  /// @brief default scatter mechanism, with predefined mean scatter
  /// rate and a transmission probability of 50%
  emcGrainScatterMechanism() : emcGrainScatterMechanism(0.5, 1e14) {}

  emcGrainScatterMechanism(T inTransmissionProb, T inScatterRate)
      : dist(0., 1.), transmissionProb(inTransmissionProb),
        scatterRate(inScatterRate) {}

  /// \brief returns the rate at which the particle scatters at the given
  /// grains.
  T getScatterRate() const { return scatterRate; }

  /*! \brief function which scatters the given particle. It determines whether
   * the particle is reflected or transmitted.
   *
   * @param particle characteristics of the particle which scatters at the grain
   */
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    if (dist(rng) > transmissionProb) {
      scatterReflectingParticle(particle, rng);
    } else {
      scatterTransmittingParticle(particle, rng);
    }
  }

  /*! \brief function which reflects the given particle of the grain.
   *
   * It is assumed that the reflection is elastic and diffusive (randomized
   * direction in opposite hemisphere, theta bigger than 90 degree).
   *
   * @param particle characteristics of the corresponding particle
   */
  void scatterReflectingParticle(emcParticle<T> &particle, emcRNG &rng) const {
    T rand = dist(rng);
    if (rand < 0.5)
      rand += 0.5;
    particle.k = initRandomDirectionWithRespectToCurrentK(
        particle.k, 1 - 2 * rand, dist(rng));
  };

  /*! \brief function which transmits the given particle through the grain.
   *
   * It is assumed that the transmission is elastic and diffusive (randomized
   * direction within the same hemisphere, theta smaller than 90 degree).
   *
   * @param particle characteristics of the corresponding particle
   */
  void scatterTransmittingParticle(emcParticle<T> &particle,
                                   emcRNG &rng) const {
    T rand = dist(rng);
    if (rand > 0.5)
      rand -= 0.5;
    particle.k = initRandomDirectionWithRespectToCurrentK(
        particle.k, 1 - 2 * rand, dist(rng));
  }
};

#endif // EMC_GRAIN_SCATTER_MECHANISM_HPP