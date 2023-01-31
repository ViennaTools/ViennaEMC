#ifndef EMC_DOPING_PROFILE_HPP
#define EMC_DOPING_PROFILE_HPP

#include <map>
#include <memory>

#include <emcGrid.hpp>
#include <emcMessage.hpp>
#include <emcUtil.hpp>

/*! \brief Represents doping profile of a device.
 *
 * Stores the doping characteristics of the device on a grid.
 * The profile is grouped in different indexed regions; currently
 * each region possesses a constant doping region. Additionally,
 * the doping concentration can be normalized by a given normalization
 * parameter.
 *
 * @param nrRegions stores the number of doping regions
 * @param idxDopingRegions stores the index of the doping region at each coord
 * @param doping stores the amount of doping at each coord
 * @param idxToDopingMap map which relates the index of the region and the
 * associated doping
 * @param normalizationParam normalization parameter for the doping
 * concentration
 */
template <class T, SizeType Dim> class emcDopingProfile {
  /// array containing values of type SizeType, with dimension
  /// of device; used to represent coordinates in the grid.
  typedef std::array<SizeType, Dim> SizeVec;

  SizeType nrRegions;
  emcGrid<int, Dim> idxDopingRegions;
  emcGrid<T, Dim> doping;
  std::map<int, T> idxToDopingMap;
  T normalizationParam;

public:
  emcDopingProfile() = delete;

  /*! \brief Constructor which initializes the extent of the grid
   * and the intrinsic doping region. The normalization parameter is
   * set to 1.
   *
   * @param gridExtent extent of the grid that represents the device
   * @param Ni intrinsic doping [1 / m ^ 3]
   */
  emcDopingProfile(const SizeVec &gridExtent, T Ni)
      : emcDopingProfile(gridExtent, Ni, 1.) {
    idxToDopingMap[-1] = Ni;
  }

  /*! \brief Constructor.
   *
   * @param gridExtent extent of the grid that represents the device
   * @param Ni intrinsic doping [1 / m ^ 3]
   * @param inNormParam normalization parameter for the doping
   */
  emcDopingProfile(const SizeVec &gridExtent, T Ni, T inNormParam)
      : nrRegions(0), idxDopingRegions(gridExtent, -1), doping(gridExtent, Ni),
        idxToDopingMap(), normalizationParam(inNormParam) {
    idxToDopingMap[-1] = Ni;
  }

  /*! \brief Adds a rectangular doping region with a constant amount of doping.
   *
   * Doping region is defined by inMinCoord, which is the minimum corner of the
   * rectangular region and inMaxCoord, which is the maximum corner.
   *
   * @param inMinCoord minimal corner of the doping region
   * @param inMaxCorner maximal corner of the doping region
   * @param inDoping amount of doping
   */
  void addConstantDopingRegion(SizeVec inMinCoord, SizeVec inMaxCoord,
                               T inDoping) {
    checkAndAdaptRegionCoords(inMinCoord, inMaxCoord);
    idxDopingRegions.fill(nrRegions, inMinCoord, inMaxCoord);
    doping.fill(inDoping, inMinCoord, inMaxCoord);
    idxToDopingMap[nrRegions] = inDoping;
    nrRegions++;
  }

  /*! \brief returns a grid that represents the device and
   * that stores the index of the doping region at each
   * grid point.
   */
  const emcGrid<int, Dim> &getDopingRegionIdx() const {
    return idxDopingRegions;
  }

  /*! \brief returns the region index at the given coordinate.
   * @param coord coordinate of the grid
   */
  int getDopingRegionIdx(const SizeVec &coord) const {
    return idxDopingRegions[coord];
  }

  /*! \brief returns a grid that represents the device and that
   * stores the doping at each grid.
   * @param normalized determines if the amount of doping should be
   * normalized by normalization parameter.
   */
  emcGrid<T, Dim> getDoping(bool normalized = false) const {
    if (!normalized)
      return doping;
    else {
      auto normedDoping = doping;
      for (auto &amountDoping : normedDoping)
        amountDoping /= normalizationParam;
      return normedDoping;
    }
  }

  /*! \brief returns the doping at a given coordinate.
   * @param coord coordinate of the device
   * @param normalized determines if the amount of doping should be
   * normalized by normalization parameter.
   */
  T getDoping(const SizeVec &coord, bool normalized = false) const {
    if (!normalized)
      return doping[coord];
    else
      return doping[coord] / normalizationParam;
  }

  /*! \brief returns the doping in the given doping region.
   * @param idxRegion index of the doping region of interest
   * @param normalized determines if the amount of doping should be
   * normalized by the normalization parameter.
   */
  T getDoping(int idxRegion, bool normalized = false) const {
    checkIdxRegion(idxRegion);
    T doping = idxToDopingMap.find(idxRegion)->second;
    if (normalized)
      doping /= normalizationParam;
    return doping;
  }

  /*! \brief returns the doping in the given doping region.
   * @param idxRegion index of the doping region of interest
   * @param normalized determines if the amount of doping should be
   * normalized.
   */
  T getDoping(SizeType idxRegion, bool normalized = false) const {
    return getDoping(static_cast<int>(idxRegion), normalized);
  }

  SizeType getNrDopingRegions() const { return nrRegions; }

private:
  void checkIdxRegion(int idxRegion) const {
    if (idxRegion < -1 || idxRegion >= (int)getNrDopingRegions()) {
      emcMessage::getInstance()
          .addError("Idx for Doping Region is out of bounds.")
          .print();
    }
  }

  /// checks if added DopingRegion is in Geometry
  /// adapts parts outside of geometry + throws error whole region is outside
  void checkAndAdaptRegionCoords(SizeVec &minCoord, SizeVec &maxCoord) {
    bool posOutOfBounds = false;
    bool maxSmallerThanMin = false;
    auto deviceExtent = idxDopingRegions.getExtent();

    for (int idxDim = 0; idxDim < Dim; idxDim++) {
      // check if doping region is completely out of bounds
      if (std::min(minCoord[idxDim], maxCoord[idxDim]) > deviceExtent[idxDim] ||
          std::max(minCoord[idxDim], maxCoord[idxDim]) < 0) {
        emcMessage::getInstance()
            .addError("Added Region is completely out of bounds. \n Please "
                      "adapt MaxPos and MinPos.")
            .print();
      }
      // check if part of region is out of bounds, ignore that part that is out
      // of bounds
      if (std::min(minCoord[idxDim], maxCoord[idxDim]) < 0 ||
          std::max(maxCoord[idxDim], minCoord[idxDim]) > deviceExtent[idxDim]) {
        posOutOfBounds = true;
        SizeType minExtent = 0;
        minCoord[idxDim] = std::max(minExtent, minCoord[idxDim]);
        maxCoord[idxDim] = std::max(minExtent, maxCoord[idxDim]);
        minCoord[idxDim] = std::min(deviceExtent[idxDim], minCoord[idxDim]);
        maxCoord[idxDim] = std::min(deviceExtent[idxDim], maxCoord[idxDim]);
      }
      // check if maxCoord is bigger than minCoord in one direction
      if (maxCoord[idxDim] < minCoord[idxDim]) {
        maxSmallerThanMin = true;
        std::swap(minCoord[idxDim], maxCoord[idxDim]);
      }
    }

    if (posOutOfBounds) {
      emcMessage::getInstance()
          .addWarning("Dopingregion was partially out of bounds. "
                      "Region-boundary that was outside is set to value of "
                      "device-boundary.")
          .print();
    }
    if (maxSmallerThanMin) {
      emcMessage::getInstance()
          .addWarning("Minimal position of region in one direction was "
                      "bigger than Maximal one. Values for MinPos and MaxPos "
                      "were swapped.")
          .print();
    }
  }

  template <class, SizeType> friend class emcDevice;
};

#endif // EMC_DOPING_PROFILE_HPP