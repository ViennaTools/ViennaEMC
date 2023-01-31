#ifndef EMC_DEVICE_HPP
#define EMC_DEVICE_HPP

#include <math.h>
#include <numeric>

#include <emcConstants.hpp>
#include <emcDopingProfile.hpp>
#include <emcGrid.hpp>
#include <emcMaterial.hpp>
#include <emcMessage.hpp>
#include <emcSurface.hpp>
#include <emcUtil.hpp>

/**
 * @brief Class which defines the semiconductor device.
 *
 * Assumes:
 *  - device is a box defined by minimal and maximal position
 *  - minimal position is at origin
 *  - fixed spacing of grid in each dimension
 *
 * @tparam T NumericType
 * @tparam Dim Dimension of the device
 * @param material semiconductor material
 * @param maxPos position of maximal corner, describes the size of the device
 * [m]
 * @param spacing spacing of the grid in each direction [m]
 * @param temperature temperature of the device [K]
 * @param dopingProfile describes the net amount of doping at each grid point,
 * summarized in regions
 * @param surface object which describes the surface of the device, with the
 * position, types and characteristics of all contacts
 * @param deviceWidth for 2D devices, describes the ectent in the 3rd dimension,
 * not needed for 3D dimensions
 */
template <class T, SizeType Dim> class emcDevice {
  static_assert(Dim == 2 || Dim == 3,
                "Wrong Dimension for Device, possible dimensions are: {2,3}.");

public:
  static const SizeType Dimension = Dim;

  typedef T ValueType;
  typedef emcMaterial<T> MaterialType;
  typedef emcSurface<T, Dim - 1> SurfaceType;
  typedef emcDopingProfile<T, Dim> DopingProfileType;
  /// vector containing Numeric Type with Dimension of device
  typedef std::array<T, Dim> ValueVec;
  /// vector containing SizeType with Dimension of device
  typedef typename emcDopingProfile<T, Dim>::SizeVec SizeVec;
  /// vector containing Numeric Type with Dim of surface (deviceDim-1)
  typedef typename emcSurface<T, Dim - 1>::SizeVecSurface SizeVecSurface;
  /// vector containing SizeType with Dim of surface (deviceDim-1)
  typedef typename emcSurface<T, Dim - 1>::ValueVecSurface ValueVecSurface;

private:
  MaterialType material;
  ValueVec maxPos;
  ValueVec spacing;
  T temperature;
  T deviceWidth{1e-6};

  /// pre-calculated characteristics of the device.
  T cellVolume, normedCellVolume;
  T thermalVoltage, debyeLength;

  DopingProfileType dopingProfile;
  SurfaceType surface;

public:
  emcDevice() = delete;

  /**
   * @brief Construct a new emcDevice object.
   *
   * @param inMaterial material of the device
   * @param inMaxPos position of maximal corner of the device [m]
   * @param inSpacing spacing of the grid in each direction [m]
   * @param inTemperature temperature of the device [K]
   */
  emcDevice(MaterialType inMaterial, ValueVec inMaxPos, ValueVec inSpacing,
            T inTemperature = 300)
      : material(inMaterial), maxPos(inMaxPos), spacing(inSpacing),
        temperature(inTemperature),
        thermalVoltage(constants::kB / constants::q * temperature),
        debyeLength(std::sqrt(constants::eps0 * material.getEpsR() *
                              thermalVoltage / constants::q /
                              material.getNi())),
        dopingProfile(maxPosToExtent(inMaxPos, inSpacing), material.getNi(),
                      material.getNi()),
        surface(maxPosToExtent(inMaxPos, inSpacing), thermalVoltage) {
    calcCellVolume();
  }

  /**
   * @brief Sets the device dimension in the 3rd dimension.
   */
  template <SizeType Dimension = Dim>
  typename std::enable_if<(Dimension == 2)>::type
  setDeviceWidth(T inDeviceWidth) {
    deviceWidth = inDeviceWidth;
    calcCellVolume();
  }

  const MaterialType &getMaterial() const { return material; }

  const DopingProfileType &getDopingProfile() const { return dopingProfile; }

  const SurfaceType &getSurface() const { return surface; }

  /// returns spacing (optional normalization by debye-length)
  ValueVec getSpacing(bool normalized = false) const {
    if (!normalized)
      return spacing;
    else
      return normalizeLength(spacing);
  }

  T getTemperature() const { return temperature; }

  /// returns thermal voltage Vt = kb * T / q
  T getThermalVoltage() const { return thermalVoltage; }

  /// returns intrinsic debyelength = sqrt(Vt * eps / q / Ni)
  T getDebyeLength() const { return debyeLength; }

  /// returns maxPos (optional normalization by debye-length)
  ValueVec getMaxPos(bool normalized = false) const {
    if (!normalized)
      return maxPos;
    else
      return normalizeLength(maxPos);
  }

  /// get the cellVolume (optional normalization)
  T getCellVolume(bool normalized = false) const {
    if (normalized)
      return normedCellVolume;
    else
      return cellVolume;
  }

  /// \brief returns a vector that holds the extent of the grid in each
  /// direction
  SizeVec getGridExtent() const { return dopingProfile.doping.getExtent(); }

  /**
   * \brief adds new constant doping region (box defined by minimal + maximal
   * corner)
   *
   * @param inMinPos position of minimal corner [m]
   * @param inMaxPos position of maximal corner [m]
   * @param inDoping net doping in this region [1 / m^3], with n-doping being
   * defined as positive net doping
   */
  void addConstantDopingRegion(ValueVec inMinPos, ValueVec inMaxPos,
                               T inDoping) {
    dopingProfile.addConstantDopingRegion(posToCoord(inMinPos),
                                          posToCoord(inMaxPos),
                                          inDoping /* + material.getNi() */);
  }

  /**
   * @brief add Ohmic contact.
   *
   * @param boundaryPos position of the boundary (enum emcBoundaryPos)
   * @param voltage applied Voltage at contact [V]
   * @param minPos position of minimal corner of the contact [m]
   * @param maxPos position of maximal corner of the contact [m]
   */
  void addOhmicContact(emcBoundaryPos boundaryPos, T voltage,
                       ValueVecSurface minPos, ValueVecSurface maxPos) {
    surface.addContact(boundaryPos, emcContactType::OHMIC, voltage,
                       posToCoord(boundaryPos, minPos),
                       posToCoord(boundaryPos, maxPos));
  }

  /**
   * @brief add Gate Contact.
   *
   * @param boundaryPos position of the boundary (enum emcBoundaryPos)
   * @param voltage applied Voltage at contact [V]
   * @param minPos position of minimal corner of the contact [m]
   * @param maxPos position of maximal corner of the contact [m]
   * @param epsRoxide relative dielectric constant of the oxide
   * @param thickness thickness of the oxide [m]
   * @param barrierHeight barrier height at the contact [V]
   */
  void addGateContact(emcBoundaryPos boundaryPos, T voltage,
                      ValueVecSurface minPos, ValueVecSurface maxPos,
                      T epsRoxide, T thickness, T barrierHeight) {
    surface.addContact(boundaryPos, emcContactType::GATE, voltage,
                       posToCoord(boundaryPos, minPos),
                       posToCoord(boundaryPos, maxPos), epsRoxide, thickness,
                       barrierHeight);
  }

  // TODO implement addSchottkyContact

  /// @brief helper function that tests if a position is out of bounds
  /// @param position position of interest
  bool isOutOfBounds(ValueVec &position) const {
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (position[idxDim] < 0 || position[idxDim] > maxPos[idxDim])
        return true;
    }
    return false;
  }

  T normalizeLength(T length) const { return length / debyeLength; }

  ValueVec normalizeLength(const ValueVec &length) const {
    ValueVec normedLength = length;
    for (auto &val : normedLength)
      val /= debyeLength;
    return normedLength;
  }

  T normalizeDoping(T doping) const { return doping / material.getNi(); }

  T normalizeVoltage(T voltage) const { return voltage / thermalVoltage; }

  T undoNormalizeLength(T normLength) const { return normLength * debyeLength; }

  T undoNormalizeDoping(T normDoping) const {
    return normDoping * material.getNi();
  }

  T undoNormalizeVoltage(T normVoltage) const {
    return normVoltage * thermalVoltage;
  }

  /**
   * @brief helper function that advances coordinate.
   * Used to iterate through device.
   *
   * @param coord current coordinate
   */
  void advanceCoord(SizeVec &coord) const {
    dopingProfile.doping.advanceCoord(coord);
  }

  /**
   * @brief helper function that checks if coord is the last coordinate in the
   * device. Used for the iteration through the device, to determine if the
   * iteration should be stopped.
   *
   * @param coord current coordinate
   */
  bool isEndCoord(const SizeVec &coord) const {
    return dopingProfile.doping.isEndCoord(coord);
  }

  /// @brief helper which calculates the position in device from the coordinate
  ValueVec coordToPos(const SizeVec &coord) const {
    ValueVec pos;
    std::transform(
        coord.begin(), coord.end(), spacing.begin(), pos.begin(),
        [](const SizeType &idx, const T &delta) { return idx * delta; });
    return pos;
  }

  /// @brief helper which calculates coordinates of grid from position in device
  SizeVec posToCoord(const ValueVec &position) const {
    SizeVec coord;
    std::transform(position.begin(), position.end(), spacing.begin(),
                   coord.begin(), [](const T &pos, const T &delta) -> SizeType {
                     return std::round(pos / delta);
                   });
    return coord;
  }

  /// @brief helper which calculates coordinate on surface from position on
  /// surface
  template <SizeType Dimension = Dim>
  typename std::enable_if<(Dimension == 2), SizeVecSurface>::type
  posToCoord(emcBoundaryPos boundPos, const ValueVecSurface &position) const {
    SizeVecSurface coord;
    if (toUnderlying(boundPos) / 2 == 0)
      coord[0] = std::round(position[0] / spacing[1]);
    else
      coord[0] = std::round(position[0] / spacing[0]);
    return coord;
  }

  /// @brief helper that calculates coordinate on surface from position on
  /// surface
  template <SizeType Dimension = Dim>
  typename std::enable_if<(Dimension == 3), SizeVecSurface>::type
  posToCoord(emcBoundaryPos boundPos, const ValueVecSurface &position) const {
    SizeVecSurface coord;
    if (toUnderlying(boundPos) / 2 == 0) {
      coord[0] = std::round(position[0] / spacing[1]);
      coord[1] = std::round(position[1] / spacing[2]);
    } else if (toUnderlying(boundPos) / 2 == 1) {
      coord[0] = std::round(position[0] / spacing[0]);
      coord[1] = std::round(position[1] / spacing[2]);
    } else {
      coord[0] = std::round(position[0] / spacing[0]);
      coord[1] = std::round(position[1] / spacing[1]);
    }
    return coord;
  }

private:
  /// @brief helper function that calculates the volume of the cell in the grid
  void calcCellVolume() {
    cellVolume = std::accumulate(spacing.begin(), spacing.end(), 1.,
                                 std::multiplies<T>());
    ValueVec normedSpacing = getSpacing(true);
    normedCellVolume = std::accumulate(
        normedSpacing.begin(), normedSpacing.end(), 1., std::multiplies<T>());
    if (Dim == 2) {
      cellVolume *= deviceWidth;
      normedCellVolume *= normalizeLength(deviceWidth);
    }
  }
};

#endif // EMC_DEVICE_HPP