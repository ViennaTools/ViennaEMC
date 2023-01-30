#ifndef EMC_SURFACE_HPP
#define EMC_SURFACE_HPP

#include <algorithm>
#include <map>
#include <memory>

#include <emcBoundaryPos.hpp>
#include <emcContact.hpp>
#include <emcGrid.hpp>
#include <emcMessage.hpp>
#include <emcUtil.hpp>

/*! \brief Class that handles the surface of the device, with all
 * contacts and their characteristics.
 *
 * @tparam T Numerical Type
 * @tparam DimSurface dimension of the surface (= dimension of device - 1)
 * @param contacts vector holding pointers to all the contacts on the surface
 * @param idxContactGrid vector that contains a grid of each boundary (facet).
 * Each grid contains the index to the contact that is placed at the specific
 * coordinate. If no contact is placed at a coordinate -1 is stored in the grid.
 * @param allBoundaryPositions vector that holds all possible values for
 * emcBoundaryPos for that device
 * @param voltageNormParam normalization parameter for the voltage
 */
template <class T, SizeType DimSurface> class emcSurface {
  static_assert(DimSurface == 1 || DimSurface == 2,
                "Wrong Dimension for Surface, possible dimensions are: {1,2}.");

public:
  static const SizeType DimDevice = DimSurface + 1;
  typedef std::array<SizeType, DimSurface> SizeVecSurface;
  typedef std::array<T, DimSurface> ValueVecSurface;

private:
  std::vector<std::unique_ptr<emcContact<T>>> contacts;
  std::vector<emcGrid<int, DimSurface>> idxContactGrid;
  std::vector<emcBoundaryPos> allBoundaryPositions;
  T voltageNormParam;

public:
  emcSurface() = delete;

  /*! \brief Constructor.
   *
   * Prepares idxContactGrid and allBoundaryPositions based on the
   * given extent and dimension.
   *
   * @param deviceGridExtent extent of the grid of the device.
   * @param inVoltageNormalizationParam parameter which is used for
   * the normalization of the voltage
   */
  emcSurface(const std::array<SizeType, DimDevice> &deviceGridExtent,
             T inVoltageNormalizationParam = 1.)
      : contacts(), idxContactGrid(), allBoundaryPositions(),
        voltageNormParam(inVoltageNormalizationParam) {
    initAllBoundaryPositions();
    for (auto &boundPos : allBoundaryPositions) {
      auto extentBound = getCoordBoundary(deviceGridExtent, boundPos);
      emcGrid<int, DimSurface> tmpGrid{extentBound, -1};
      idxContactGrid.push_back(tmpGrid);
    }
  }

  SizeType getNrContacts() const { return contacts.size(); }

  /**
   * @brief returns the index of the contact at the given
   * boundary and coordinate on the boundary.
   *
   * @param coord coordinate on the given boundary
   * @param boundPos boundary position (enum emcBoundaryPos)
   */
  int getContactIdx(const SizeVecSurface &coord,
                    emcBoundaryPos boundPos) const {
    return idxContactGrid[toUnderlying(boundPos)][coord];
  }

  /**
   * @brief returns the index of the contact at this coordinate.
   * If coordinate is not on boundary, also -1 is returned.
   *
   * NOTE: this index is not unique for grid points that lie on
   * multiply boundary positions (at corners, edges).
   */
  int getContactIdx(const std::array<SizeType, DimDevice> &deviceCoord) const {
    emcBoundaryPos boundPos = getBoundaryPos(deviceCoord);
    if (boundPos != emcBoundaryPos::INVALID) {
      return idxContactGrid[toUnderlying(boundPos)]
                           [getCoordBoundary(deviceCoord, boundPos)];
    } else
      return -1;
  }

  /**
   * @brief returns the voltage of a specific contact.
   *
   * @param idxContact index of the contact
   * @param normalized determines whether the voltage should be normalized
   * by normalization parameter.
   */
  T getContactVoltage(SizeType idxContact, bool normalized = false) const {
    checkContactIdx(idxContact);
    if (!normalized)
      return contacts[idxContact]->getAppliedVoltage();
    else
      return contacts[idxContact]->getAppliedVoltage() / voltageNormParam;
  }

  /**
   * @brief returns the voltage of the specific contact.
   *
   * @param coord coordinate on the given boundary
   * @param boundPos boundary position (enum emcBoundaryPos)
   * @param normalized determines whether the voltage should be normalized
   * by normalization parameter
   */
  T getContactVoltage(const SizeVecSurface &coord, emcBoundaryPos boundPos,
                      bool normalized = false) const {
    return getContactVoltage(getContactIdx(coord, boundPos), normalized);
  }

  /**
   * @brief returns the type of the given contact.
   *
   * @param idxContact index of the contact
   */
  emcContactType getContactType(SizeType idxContact) const {
    checkContactIdx(idxContact);
    return contacts[idxContact]->getType();
  }

  /**
   * @brief returns type of the given contact.
   *
   * @param coord coordinate on the given boundary
   * @param boundPos boundary position (enum emcBoundaryPos)
   */
  emcContactType getContactType(const SizeVecSurface &coord,
                                emcBoundaryPos boundPos) const {
    return getContactType(getContactIdx(coord, boundPos));
  }

  /**
   * @brief returns additional parameter of a given contact.
   *
   * @param idxContact index of the contact
   * @param idxInformation index of the additional parameter
   */
  T getContactFurtherParameter(SizeType idxContact,
                               SizeType idxInformation) const {
    checkContactIdx(idxContact);
    return contacts[idxContact]->getFurtherParameter(idxInformation);
  }

  /**
   * @brief returns additional parameter of a given contact.
   *
   * @param coord coordinate on the given boundary
   * @param boundPos boundary position (enum emcBoundaryPos)
   * @param idxInformation index of the additional parameter
   */
  T getContactFurtherParameter(const SizeVecSurface &coord,
                               emcBoundaryPos boundPos,
                               const SizeType idxInformation) const {
    return getContactFurtherParameter(getContactIdx(coord, boundPos),
                                      idxInformation);
  }

  /**
   * @brief returns the index of the ohmic contact.
   *
   * NOTE: this index is only unique for ohmic contacts as only for this
   * type of contact all occurences of the grid point on all boundary
   * positions are updated!
   */
  int getOhmicContactIdx(const std::array<SizeType, DimDevice> &coord) const {
    return getContactIdx(coord);
  }

  const std::vector<emcBoundaryPos> &getAllBoundaryPos() const {
    return allBoundaryPositions;
  }

  /*! \brief adds a contact to a specific boundary from minPos to maxPos
   * and stores the contact in the contacts-vector.
   *
   * @param boundaryPos boundary position of the contact (enum emcBoundPos)
   * @param type type of the contact (enum emcContactType)
   * @param minCoord coordinate of the minimal corner of the contact
   * @param maxCoord coordinate of the maximal corner of the contact
   * @param epsOxide relative dielectric constant of the oxide (needed for gate
   * contact)
   * @param thickness thickness of the oxide (needed for gate contact)
   * @param barrierHeight barrier height at the contact (needed for gate
   * contact)
   */
  void addContact(emcBoundaryPos boundaryPos, emcContactType type,
                  T appliedVoltage, const SizeVecSurface &minCoord,
                  const SizeVecSurface &maxCoord, T epsOxide = 0,
                  T thickness = 0, T barrierHeight = 0) {
    checkContactPos(boundaryPos, minCoord, maxCoord);
    idxContactGrid[toUnderlying(boundaryPos)].fill(contacts.size(), minCoord,
                                                   maxCoord);

    switch (type) {
    case emcContactType::GATE:
      contacts.push_back(std::make_unique<emcGateContact<T>>(
          appliedVoltage, epsOxide, thickness, barrierHeight));
      break;
    case emcContactType::OHMIC:
      updateAllOccurences(boundaryPos, minCoord, maxCoord, contacts.size());
      contacts.push_back(std::make_unique<emcOhmicContact<T>>(appliedVoltage));
      break;
    case emcContactType::SCHOTTKY:
      updateAllOccurences(boundaryPos, minCoord, maxCoord, contacts.size());
      contacts.push_back(
          std::make_unique<emcSchottkyContact<T>>(appliedVoltage));
      break;
    }
  }

  /// @brief helper that checks if artificial boundary is present at the
  /// boundary coordinate coord
  bool isArtificialBoundary(const SizeVecSurface &coord,
                            emcBoundaryPos boundPos) const {
    return getContactIdx(coord, boundPos) == -1;
  }

  /// @brief helper that checks if specific contact-type is present at the
  /// boundary coordinate coord
  bool isContactType(const SizeVecSurface &coord, emcBoundaryPos boundPos,
                     emcContactType type) const {
    T idxContact = getContactIdx(coord, boundPos);
    return !(idxContact == -1) && getContactType(idxContact) == type;
  }

  /// @brief helper that checks if a given coord is at an ohmic contact
  bool isOhmicContact(const std::array<SizeType, DimDevice> &coord) const {
    auto boundPos = getBoundaryPos(coord);
    if (boundPos != emcBoundaryPos::INVALID)
      return isContactType(getCoordBoundary(coord, boundPos), boundPos,
                           emcContactType::OHMIC);
    return false;
  }

  /// @brief prints idxContactGrid (for debugging)
  void print(std::ostream &out = std::cout) const {
    for (auto &boundary : idxContactGrid) {
      boundary.print(out);
      out << "\n";
    }
  }

  /// @brief prints the indices of contacts of a specific boundary given by
  /// boundPos. (for debugging)
  void print(emcBoundaryPos boundPos, std::ostream &out = std::cout) const {
    if (toUnderlying(boundPos) < 2 * DimDevice)
      idxContactGrid.at(toUnderlying(boundPos)).print(out);
  }

  /// @brief helper that checks if deviceCoord is on Boundary
  bool isOnBoundary(const std::array<SizeType, DimDevice> &deviceCoord) const {
    bool onBoundary = false;
    for (auto &boundaryPos : allBoundaryPositions) {
      if (isOnBoundary(deviceCoord, boundaryPos)) {
        onBoundary = true;
        break;
      }
    }
    return onBoundary;
  }

  /// @brief helper that checks if deviceCoord is on a specific boundary
  bool isOnBoundary(const std::array<SizeType, DimDevice> &deviceCoord,
                    emcBoundaryPos boundPos) const {
    SizeType idxDim = toUnderlying(boundPos) / 2;
    if (toUnderlying(boundPos) % 2 == 1)
      return deviceCoord[idxDim] == (getDeviceExtent(idxDim) - 1);
    else
      return deviceCoord[idxDim] == 0;
  }

  /// @brief calculates the coordinate on the surface from the coordinate on
  /// the contact and the boundary the coordinate is on
  SizeVecSurface
  getCoordBoundary(const std::array<SizeType, DimDevice> &deviceCoord,
                   emcBoundaryPos idxBoundary) const {
    checkBoundaryIdx(idxBoundary);
    // TODO onBoundary(coord,idxBoundary)
    SizeVecSurface coordBound;
    SizeType idxFixed = toUnderlying(idxBoundary) / 2;
    std::copy(deviceCoord.begin(), deviceCoord.begin() + idxFixed,
              coordBound.begin());
    std::copy(deviceCoord.begin() + idxFixed + 1, deviceCoord.end(),
              coordBound.begin() + idxFixed);
    return coordBound;
  }

  /// @brief calculates the device-coordinate from a given boundary-coordinate
  std::array<SizeType, DimDevice>
  getCoordDevice(const SizeVecSurface &boundaryCoord,
                 emcBoundaryPos idxBoundary) const {
    checkBoundaryIdx(idxBoundary);
    std::array<SizeType, DimDevice> deviceCoord;
    // fill value from fixed dimension on surface
    SizeType idxFixed = toUnderlying(idxBoundary) / 2;
    if (toUnderlying(idxBoundary) % 2 == 0)
      deviceCoord[idxFixed] = 0;
    else
      deviceCoord[idxFixed] = getDeviceExtent(idxFixed) - 1;
    // copy rest of the coordinate
    std::copy(boundaryCoord.begin(), boundaryCoord.begin() + idxFixed,
              deviceCoord.begin());
    std::copy(boundaryCoord.begin() + idxFixed, boundaryCoord.end(),
              deviceCoord.begin() + idxFixed + 1);

    return deviceCoord;
  }

  /// returns a boundaryPos on which the given device-coordinate can be found
  /// NOTE: not unique for coordinates that are on multiple boundaries!
  emcBoundaryPos
  getBoundaryPos(const std::array<SizeType, DimDevice> &deviceCoord) const {
    emcBoundaryPos boundPos;
    for (const auto &boundPos : allBoundaryPositions) {
      if (isOnBoundary(deviceCoord, boundPos)) {
        return boundPos;
      }
    }
    return emcBoundaryPos::INVALID;
  }

  /// @brief initializes the coordinate at the boundary.
  /// Used for the iteration over the surface.
  void initCoord(emcBoundaryPos &boundPos, SizeVecSurface &coord) const {
    boundPos = emcBoundaryPos::XMIN;
    coord.fill(0);
  }

  /// @brief advances coordinates at surface.
  /// Used for the iteration over the surface.
  void advanceCoord(emcBoundaryPos &boundPos, SizeVecSurface &coord) const {
    idxContactGrid[toUnderlying(boundPos)].advanceCoord(coord);
    if (idxContactGrid[toUnderlying(boundPos)].isEndCoord(coord) &&
        boundPos != allBoundaryPositions.back()) {
      coord.fill(0);
      SizeType idxBoundary = toUnderlying(boundPos);
      boundPos = static_cast<emcBoundaryPos>(idxBoundary + 1);
    }
  }

  /// @brief checks if coordinate is the last coordinate in grid.
  /// Used for the iteration over the surface.
  bool isEndCoord(emcBoundaryPos &boundPos, SizeVecSurface &coord) const {
    return (boundPos == allBoundaryPositions.back() &&
            idxContactGrid.back().isEndCoord(coord));
  }

private:
  /// @brief tests if the added contact is on multiple boundaries, if it is all
  /// occurences in idxContactGrid are updated.
  void updateAllOccurences(emcBoundaryPos boundPos,
                           const SizeVecSurface &minCoord,
                           const SizeVecSurface &maxCoord, int idxContact) {
    for (SizeType idxDim = 0; idxDim < DimSurface; idxDim++) {
      auto startMin = minCoord;
      auto startMax = maxCoord;
      while (startMin[idxDim] <= maxCoord[idxDim]) {
        auto deviceCoordFromMin = getCoordDevice(startMin, boundPos);
        for (auto &boundaryPos : allBoundaryPositions) {
          if (isOnBoundary(deviceCoordFromMin, boundaryPos)) {
            auto boundaryCoord =
                getCoordBoundary(deviceCoordFromMin, boundaryPos);
            idxContactGrid[toUnderlying(boundaryPos)].fill(
                idxContact, boundaryCoord, boundaryCoord);
          }
          if (DimSurface == 2) {
            auto deviceCoordFromMax = getCoordDevice(startMax, boundPos);
            if (isOnBoundary(deviceCoordFromMax, boundaryPos)) {
              auto boundaryCoord =
                  getCoordBoundary(deviceCoordFromMax, boundaryPos);
              idxContactGrid[toUnderlying(boundaryPos)].fill(
                  idxContact, boundaryCoord, boundaryCoord);
            }
          }
        }
        startMin[idxDim]++;
        startMax[idxDim]--;
      }
    }
  }

  /// @brief helper that returns the extent of the grid of the device
  /// in a specific dimension
  SizeType getDeviceExtent(SizeType idxDim) const {
    if (idxDim != 0)
      return idxContactGrid[0].getSize(idxDim - 1);
    else
      return idxContactGrid[2].getSize(0);
  }

  /// @brief returns extent of a specific boundary-grid
  SizeVecSurface getExtentBoundary(emcBoundaryPos idxBoundary) const {
    checkBoundaryIdx(idxBoundary);
    return idxContactGrid[toUnderlying(idxBoundary)].getExtent();
  }

  /// @brief appends all existing boundaryPositions to this vector
  void initAllBoundaryPositions() {
    allBoundaryPositions = {emcBoundaryPos::XMIN, emcBoundaryPos::XMAX,
                            emcBoundaryPos::YMIN, emcBoundaryPos::YMAX};
    if (DimDevice == 3) {
      allBoundaryPositions.insert(allBoundaryPositions.end(),
                                  {emcBoundaryPos::ZMIN, emcBoundaryPos::ZMAX});
    }
  }

  /// @brief helper that checks if index for contact vector id valid
  void checkContactIdx(const SizeType idxContact) const {
    if (idxContact >= getNrContacts())
      emcMessage::getInstance()
          .addError("Index for Contact is out of bounds.")
          .print();
  }

  /// @brief helper that checks if idx for boundary is valid
  void checkBoundaryIdx(const emcBoundaryPos idxBoundary) const {
    if (toUnderlying(idxBoundary) >= 2 * DimDevice)
      emcMessage::getInstance()
          .addError("Index for Boundary is out of bounds.")
          .print();
  }

  /// @brief helper that checks if contact position is valid
  void checkContactPos(emcBoundaryPos boundaryPos,
                       const SizeVecSurface &minCoord,
                       const SizeVecSurface &maxCoord) const {
    // check if coordinates are valid
    if (!idxContactGrid[toUnderlying(boundaryPos)].isValid(minCoord)) {
      emcMessage::getInstance()
          .addError("MinPos for Contact is out of bounds.")
          .print();
    }
    if (!idxContactGrid[toUnderlying(boundaryPos)].isValid(maxCoord)) {
      emcMessage::getInstance()
          .addError("MaxCoord for Contact is out of bounds.")
          .print();
    }
    // check if minPos[idxDim] <= maxPos[idxDim]
    for (int idxDim = 0; idxDim < DimSurface; idxDim++) {
      if (maxCoord[idxDim] < minCoord[idxDim]) {
        emcMessage::getInstance()
            .addError("Component of MaxPos is bigger than component of MinPos "
                      "in at least one direction for added contact.")
            .print();
      }
    }
  }
};

#endif // EMC_SURFACE_HPP