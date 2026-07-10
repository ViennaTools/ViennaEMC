#ifndef EMC_ACCEPTOR_HPP
#define EMC_ACCEPTOR_HPP

#include <ParticleType/emcParticleType.hpp>
#include <emcConstants.hpp>

/**
 * @brief Non-moving ionised acceptor ion.
 *
 * Symmetric counterpart to emcDonor: charge = -q, moved = false.
 * Assumes every acceptor atom has accepted an additional electron and is
 * therefore negatively charged. Placed at the same positions and density
 * as the hole population it charge-compensates.
 *
 * Used in conjunction with emcHole to model p-type regions or to provide
 * the fixed background charge for photo-generated hole populations in
 * device simulations with FMM particle-particle interactions.
 *
 * @tparam T          numeric type (float or double)
 * @tparam DeviceType device type
 */
template <class T, class DeviceType>
struct emcAcceptor : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::ValueVec ValueVec;
  typedef typename DeviceType::SizeVec SizeVec;

  static const SizeType Dim = DeviceType::Dimension;

  emcAcceptor() : emcParticleType<T, DeviceType>(1, 1) {}

  std::string getName() const { return "Acceptor"; }

  T getCharge() const { return -constants::q; }

  bool isMoved() const { return false; }

  bool isInjected() const { return false; }

  T getInitialNrParticles(const SizeVec &coord, const DeviceType &device,
                          const emcGrid<T, Dim> & /*potential*/) {
    T hDensity = std::fabs(device.getDopingProfile().getDoping(coord));
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 ||
          coord[idxDim] == device.getGridExtent()[idxDim] - 1)
        hDensity *= T(0.5);
    }
    return hDensity * device.getCellVolume();
  }
};

#endif // EMC_ACCEPTOR_HPP
