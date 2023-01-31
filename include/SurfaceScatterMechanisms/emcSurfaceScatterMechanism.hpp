#ifndef EMC_SURFACE_SCATTER_MECHANISM_HPP
#define EMC_SURFACE_SCATTER_MECHANISM_HPP

#include <emcBoundaryPos.hpp>
#include <emcParticle.hpp>

/*! \brief Generic Surface Scatter Mechanism. This class handles Scattering of a
 * particle at the boundary of the device.
 *
 * The derived classes determine at which rate and at which angles the particle
 * scatter
 *
 * @param boundaryPos Boundary of the device for which the scattering is handled
 * @param maxPos position of the maximal extent of the device (needed for
 * determination of position of MAX-boundPos)
 */
template <class T, class DeviceType, SizeType Dim = DeviceType::Dimension>
struct emcSurfaceScatterMechanism {
protected:
  emcBoundaryPos boundaryPos;
  std::array<T, Dim> maxPos;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcSurfaceScatterMechanism() = delete;

  emcSurfaceScatterMechanism(std::array<T, Dim> inMaxPos)
      : dist(0., 1.), maxPos(inMaxPos) {}

  virtual ~emcSurfaceScatterMechanism() = default;

  /*! \brief function which scatters the given particle. It determines whether
   * the scattering is defusive or specularly.
   *
   * @param particle for which the scatter-rate should be calculated
   * @param pos position-vector to assign the new position of the particle to
   */
  void scatterParticle(emcParticle<T> &particle, std::array<T, Dim> &pos,
                       emcRNG &rng) const {
    if (dist(rng) < getDiffScatterProb(particle))
      scatterParticleDiffusely(particle, pos, rng);
    else
      scatterParticleSpecularly(particle, pos);
  }

  // sets the boundary of the device for which the scattering is handled
  void setBoundaryPosition(emcBoundaryPos inBoundaryPos) {
    boundaryPos = inBoundaryPos;
  }

protected:
  /*! \brief function which returns the probability of a given particle
   * scattering diffusely. It needs to be overridden in the concrete classes.
   *
   * @param particle for which the scatter-rate should be calculated
   */
  virtual T getDiffScatterProb(emcParticle<T> &particle) const = 0;

  /*! \brief function which scatters the given particle diffusely. It needs to
   * be overridden in the concrete classes.
   *
   * @param particle which should be scattered
   * @param pos position-vector to assign the new position of the particle to
   */
  virtual void scatterParticleDiffusely(emcParticle<T> &particle,
                                        std::array<T, Dim> &pos,
                                        emcRNG &rng) const = 0;

  /*! \brief function which scatters the given particle specularly.
   *
   * @param particle which should be scattered
   * @param pos position-vector to assign the new position of the particle to
   */
  void scatterParticleSpecularly(emcParticle<T> &particle,
                                 std::array<T, Dim> &pos) const {
    for (SizeType idx = 0; idx < Dim; idx++)
      if (pos[idx] < 0) {
        pos[idx] *= -1;
        if (particle.k[idx] < 0)
          particle.k[idx] *= -1;
      } else if (pos[idx] > maxPos[idx]) {
        pos[idx] = 2 * maxPos[idx] - pos[idx];
        if (particle.k[idx] > 0)
          particle.k[idx] *= -1;
      }
  };

  /// Helper function, that determines the absolut position of the
  /// max-boundary of the given axis-index in case the position passed is
  /// outside of the device (MAX-direction) or 0 incase it didn't
  T getAbsolutPositionOfBoundary(SizeType index,
                                 std::array<T, Dim> &pos) const {
    if (pos[index] > maxPos[index]) {
      return maxPos[index];
    } else {
      return 0;
    }
  }

  /// Helper function that checks whether the position is out of bounds of the
  /// device for the given axis-index
  T isOutOfBounds(SizeType index, std::array<T, Dim> &pos) const {
    return (pos[index] > maxPos[index] || pos[index] < 0);
  }

  /// Helper function, that obtains the index of a position or k vector which is
  /// perpendicular to the assigned boundary position
  int getIndexFromBoundaryPos(int delta = 0) const {
    auto index = toUnderlying(boundaryPos);
    return (index / 2 + delta) % 3;
  }

  /// Returns whether the assigned boundary is at a max position of the device
  bool isScatteringAtMaxPos() const {
    auto index = toUnderlying(boundaryPos);
    return index % 2 == 1;
  }

  /*! \brief function which reassigns the values of the position and k vector
   * after a scattering event with the given polar and azimuthal angle and
   * adapts the position (s.t. it is inside the device) and the wave vector
   * (s.t. it points into the device).
   *
   * @param particle particle characterics from scattered particle
   * @param pos position-vector of the particle
   * @param theta polar angle for k after the scattering
   * @param phi azimuthal angle for k after the scattering
   */
  void calculateAndAssignKAndPos(emcParticle<T> &particle,
                                 std::array<T, Dim> &pos, T theta,
                                 T phi) const {
    // set random direction into the device
    T speed = norm(particle.k);
    auto perpendicularIndex = getIndexFromBoundaryPos();
    T signFactor = 1;
    if (isScatteringAtMaxPos())
      signFactor = -1;

    particle.k[perpendicularIndex] = speed * std::cos(theta) * signFactor;
    particle.k[getIndexFromBoundaryPos(1)] =
        speed * std::sin(theta) * std::cos(phi);
    particle.k[getIndexFromBoundaryPos(2)] =
        speed * std::sin(theta) * std::sin(phi);

    // adapt particle position + adapt wave vector direction (if required)
    scatterParticleSpecularly(particle, pos);
  }
};

#endif // EMC_SURFACE_SCATTER_MECHANISM_HPP