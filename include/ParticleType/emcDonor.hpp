#ifndef EMC_DONOR_HPP
#define EMC_DONOR_HPP

#include <ParticleType/emcParticleType.hpp>

/** @brief Non-moving donor ion.
 *
 * Assumes that every donor atom has given additional electron away.
 */
template <class T, class DeviceType>
struct emcDonor : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::ValueVec ValueVec;
  typedef typename DeviceType::SizeVec SizeVec;

  static const SizeType Dim = DeviceType::Dimension;

  emcDonor() : emcParticleType<T, DeviceType>(1, 1) {}

  std::string getName() const { return "Donor"; }

  T getCharge() const { return constants::q; }

  bool isMoved() const { return false; }

  bool isInjected() const { return false; }

  // same nr. as electrons at each coordinate
  T getInitialNrParticles(const SizeVec &coord, const DeviceType &device,
                          const emcGrid<T, Dim> & /*potential*/) {
    T eDensity = device.getDopingProfile().getDoping(coord);
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 ||
          coord[idxDim] == device.getGridExtent()[idxDim] - 1)
        eDensity *= 0.5;
    }
    return eDensity * device.getCellVolume();
  }
};

#endif