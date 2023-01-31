#ifndef EMC_ABSTRACT_SOLVER_HPP
#define EMC_ABSTRACT_SOLVER_HPP

#include <emcGrid.hpp>
#include <emcUtil.hpp>

///\brief Base Class for PoissonSolver.
template <class T, class DeviceType, class ParticleHandler>
struct emcAbstractSolver {
  typedef emcGrid<T, DeviceType::Dimension> GridType;

  /*! \brief Calculates the equilibrium potential.
   *
   * Assumes that particle concentrations are in equilibrium and
   * discards applied voltage at contacts.
   *
   * @param pot grid that holds initial guess for potential, is overwritten
   * with newly calculated potential
   * @param device used device for the potential calculation
   * @param resetBC boolean telling if the boundary conditions have to be reset
   */
  virtual void calcEquilibriumPotential(GridType &pot, const DeviceType &device,
                                        bool resetBC = true) = 0;

  /*! \brief Calculates the non-equilibrium potential.
   *
   * Assumes that electron concentration is given at all grid points,
   * hole concentration is in equilibrium ans uses the applied voltage
   * at the contacts.
   *
   * @param pot grid that holds initial guess for potential, is overwritten
   * with newly calculated potential
   * @param device used device for the potential calculation
   * @param eConc electron concentration at each grid point
   * @param resetBC boolean telling if the boundary conditions have to be reset
   *
   * // TODO loose assumption that only electron concentration is given.
   */
  virtual void calcNonEquilibriumPotential(GridType &pot,
                                           const DeviceType &device,
                                           const GridType &eConc,
                                           bool resetBC = true) = 0;

  /*! \brief Calculates the bakcground potential.
   *
   * Assumes that particle concentrations are 0 (because they are handled by
   * FMM) adapts BC so that carrier effects are discarded.
   *
   * @param pot grid that holds initial guess for potential, is overwritten
   * with newly calculated potential
   * @param device used device for the potential calculation
   * @param particleHandler holds information on all the particles.
   * @param resetBC boolean telling if the boundary conditions have to be reset
   */
  virtual void calcBackgroundPotential(GridType &pot, const DeviceType &device,
                                       ParticleHandler &handler,
                                       bool resetBC = true) = 0;
};

#endif // EMC_ABSTRACT_SOLVER_HPP