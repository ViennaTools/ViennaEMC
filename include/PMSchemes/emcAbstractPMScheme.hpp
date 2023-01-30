#ifndef EMC_ABSTRACT_PM_SCHEME_HPP
#define EMC_ABSTRACT_PM_SCHEME_HPP

#include <emcGrid.hpp>
#include <emcUtil.hpp>

#include <array>
#include <vector>

/*! \brief Abstract Base Class for Particle-Mesh-Scheme (PMScheme).
 *
 * A PMScheme relates the discrete representation of the device and
 * the continuous particle positions. The assignment of particles
 * to a given mesh (discrete representation of a device), the
 * interpolation of discrete characteristics to continuous ones
 * and the discrete derivation of the potential to obtain the
 * electric field are handled by objects of this class.
 *
 */
template <class T, class DeviceType> class emcAbstractPMScheme {
public:
  static const SizeType Dim = DeviceType::Dimension;

  /*! \brief assigns one simulated particle (representing nrCarriers carriers)
   * to the given mesh (single particle assignment)
   *
   * @param pos position of the particle
   * @param nrCarriers how many carriers are represented by that particle (> 1
   * if it is superparticle)
   * @param spacing spacing of the discrete representation
   * @param gridNrParticles grid, to which the particle is assigned (will be
   * adapted)
   */
  virtual void assignToMesh(const std::array<T, Dim> &pos, SizeType nrCarriers,
                            const std::array<T, Dim> &spacing,
                            emcGrid<T, Dim> &gridNrParticles) const = 0;

  /*! \brief assigns all particles from container to the given mesh (each
   * represents nrCarriers particles) (multiple particles assignment)
   *
   * @param position container that holds the continuous position of all
   * particles
   * @param nrCarriers how many carriers are represented by one particle (> 1 if
   * superparticles are used)
   * @param spacing spacing of the discrete representation
   * @param gridNrParticles grid, to which the particles are assigned (will be
   * adpated)
   */
  virtual void assignToMesh(const std::vector<std::array<T, Dim>> &position,
                            SizeType nrCarriers,
                            const std::array<T, Dim> &spacing,
                            emcGrid<T, Dim> &gridNrParticles) const = 0;

  /*! \brief interpolates the given electric field to the particle position
   * and calculates the force on the particle
   *
   * @param eField el. field at discrete grid points
   * @param pos continuous position of the particle
   * @param spacing spacing of the discrete representation (grid)
   * @param charge charge of the given particle
   * @return force that is seen by particle at current position
   */
  virtual std::array<T, 3>
  interpolateForce(const std::vector<emcGrid<T, Dim>> &eField,
                   const std::array<T, Dim> &pos,
                   const std::array<T, Dim> &spacing, T charge) const = 0;

  /*! \brief calculates the electric field from the given potential
   *
   * @param eField el. field at discrete grid points (will be overwritten)
   * @param potential given potential
   * @param device given device
   */
  virtual void calcEField(std::vector<emcGrid<T, Dim>> &eField,
                          const emcGrid<T, Dim> &potential,
                          const DeviceType &device) const = 0;
};

#endif // EMC_ABSTRACT_PM_SCHEME_HPP