#ifndef EMC_ABSTRACT_VALLEY_HPP
#define EMC_ABSTRACT_VALLEY_HPP

#include <array>

#include <emcMessage.hpp>
#include <emcUtil.hpp>

/*! \brief Abstract Class used for the implementation of Valley in the
 * Bandstructure of a Semiconductor.
 *
 * Valleys are simplifications of the bandstructure, using the effective mass
 * approach (analytical approximation for bandstructure). Note: One Valley in
 * this implementation can consist of multiple subvalleys that all have the same
 * effective masses and bottom energies, but can have different orientations of
 * their ellipse coordinate systems.
 *
 * For more information see README-file in the ValleyType (current) folder.
 */
template <class T> class emcAbstractValley {
public:
  virtual ~emcAbstractValley() = default;

  /// returns the density of states effective mass
  /// @param energy energy of particle [eV]
  virtual T getEffMassDOS(T energy = 0) const = 0;

  /// returns the conductive effective mass
  /// @param energy energy of particle [ev]
  virtual T getEffMassCond(T energy = 0) const = 0;

  /// returns alpha (non-Parabolicity Factor) [1 / eV]
  /// return 0 if valley is parabolic
  virtual T getNonParabolicity() const = 0;

  /// return the energy at the bottom of the valley [eV]
  virtual T getBottomEnergy() const = 0;

  /// returns degeneracy factor of valley
  /// = nr of equivallent subvalleys
  virtual SizeType getDegeneracyFactor() const = 0;

  /// returns the norm of the wave-vector with given energy
  /// @param energy energy of particle [eV]
  virtual T getNormWaveVec(T energy) const = 0;

  /// returns energy (in eV)
  /// @param k wave-vector of particle
  virtual T getEnergy(const std::array<T, 3> &k) const = 0;

  /// returns Gamma = h^2 * k^2 / 2 m (in eV)
  /// @param energy energy of particle [eV]
  virtual T getGamma(T energy) const = 0;

  /// returns velocity of a particle (in m / s)
  /// @param k wave-vector of particle
  /// @param energy energy of particle [eV]
  /// @param idxSubValley index of subValley of particle
  virtual std::array<T, 3> getVelocity(const std::array<T, 3> &k, T energy,
                                       SizeType idxSubValley) const = 0;

  /// returns an array with the factors of the Herring-Vogt Transformation
  /// matrix (diagonal elements in effective mass tensor in ellipse
  /// coordinate system), set to (1,1,1) if effective mass is isotropical
  virtual const std::array<T, 3> &getVogtTransformationFactor() const = 0;

  /// helper that transforms vec from device to ellipse coordinates
  /// of specific subvalley
  virtual std::array<T, 3>
  transformToDeviceCoord(SizeType idxSubValley,
                         const std::array<T, 3> &vec) const = 0;

  /// helper that transforms vec from ellipse to device coordinates
  virtual std::array<T, 3>
  transformToEllipseCoord(SizeType idxSubValley,
                          const std::array<T, 3> &vec) const = 0;

  /// helper that checks if a valley-object uses realistic parameter
  void check() const {
    if (getDegeneracyFactor() < 1) {
      emcMessage::getInstance()
          .addError("Degeneracy Factor of a valley has to be at least 1.")
          .print();
    }
    if (getNonParabolicity() < 0) {
      emcMessage::getInstance()
          .addError("NonParabolicity of a valley has to be positive.")
          .print();
    }
  }
};

#endif // EMC_ABSTRACT_VALLEY