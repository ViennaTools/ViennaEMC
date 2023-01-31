#ifndef EMC_GRID_HPP
#define EMC_GRID_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <emcMessage.hpp>
#include <emcUtil.hpp>

/*! \brief A class that stores a grid of a specific dimension as 1D-vector.
 *
 * @tparam T type contained in the grid
 * @tparam Dim Dimension of the grid
 */
template <typename T, SizeType Dim> class emcGrid {
  static_assert(Dim > 0 && Dim < 4,
                "Wrong Dimension of Grid, possible dimensions are: {1,2,3}.");

public:
  /// array-type that holds coordinate of grid point
  typedef std::array<SizeType, Dim> CoordVec;

private:
  std::vector<T> data; // container, that stores the data in 1D
  CoordVec extent;     // nr of grid points in each dimension

public:
  // --- constructors --------------------------------------------------------
  emcGrid() = default;

  emcGrid(const emcGrid<T, Dim> &rhs) : data(), extent(rhs.getExtent()) {
    resize(extent);
    std::copy(rhs.begin(), rhs.end(), begin());
  }

  emcGrid(emcGrid<T, Dim> &&rhs)
      : data(std::move(rhs.data)), extent(std::move(rhs.extent)) {}

  /*! Constructs an empty grid with the given extent.
   * @param inExtent array that tells extension of grid in each dimension
   */
  emcGrid(const CoordVec &inExtent) : data(), extent(inExtent) {
    resize(extent);
  }

  /*! Constructs a grid with the given extent and fills it with the given
   * value.
   * @param inExtent array that tells extension of grid in each dimension
   * @param defaultValue value with which the grid is filled
   */
  emcGrid(const CoordVec &inExtent, T defaultValue) : data(), extent(inExtent) {
    resize(extent);
    fill(defaultValue);
  }

  // --- get-functions --------------------------------------------------------
  /*! Returns the size of a specific dimension
   * @param idxDim index of dimension
   */
  SizeType getSize(SizeType idxDim) const {
    checkDimIdx(idxDim);
    return extent[idxDim];
  }

  /// returns total size ( = nr. of grid points)
  SizeType getSize() const { return data.size(); }

  CoordVec getExtent() const { return extent; }

  /// returns next value in the direction of idxDim (if it exists)
  const T &getNextValue(const CoordVec &currCoord, int idxDim) const {
    checkNbrExistence(currCoord, idxDim, 1);
    CoordVec nextCoord{currCoord};
    nextCoord[idxDim]++;
    return data[coordToIdx(nextCoord)];
  }

  /// returns prev value in the direction of idxDim (if it exists)
  const T &getPrevValue(const CoordVec &currCoord, int idxDim) const {
    checkNbrExistence(currCoord, idxDim, -1);
    CoordVec prevCoord(currCoord);
    prevCoord[idxDim]--;
    return data[coordToIdx(prevCoord)];
  }

  // --- further functions ----------------------------------------------
  /*! Fills whole container with same value.
   * @param value value that is inserted in container
   */
  void fill(const T value) { std::fill(data.begin(), data.end(), value); }

  /*! Fills box from minCoord to maxCoord  with same value.
   * @param value value that is inserted in container
   * @param minCoord minimum corner of the box
   * @param maxCoord maximum corner of the box
   */
  void fill(const T value, const CoordVec &minCoord, const CoordVec &maxCoord) {
    CoordVec currCoord;
    currCoord.fill(0);
    currCoord[0] = -1;
    std::replace_if(
        begin(), end(),
        [this, &currCoord, &minCoord, &maxCoord](T &entry) {
          advanceCoord(currCoord);
          for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
            if (minCoord[idxDim] > currCoord[idxDim] ||
                maxCoord[idxDim] < currCoord[idxDim])
              return false;
          }
          return true;
        },
        value);
  }

  /// helper that advances coordinate to go through 1D
  /// data structure the same way that iterator does
  void advanceCoord(CoordVec &currCoord) const {
    currCoord[0]++;
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (currCoord[idxDim] >= extent[idxDim] && idxDim != Dim - 1) {
        currCoord[idxDim + 1]++;
        for (SizeType idxResetDim = 0; idxResetDim <= idxDim; idxResetDim++)
          currCoord[idxResetDim] = 0;
      }
    }
  }

  /// check if coord is the coord after the last coord (like iterator end())
  /// enables iterating through grid combined with advanceCoord()
  bool isEndCoord(const CoordVec &coord) const {
    return coord[Dim - 1] == extent[Dim - 1];
  }

  // fill data structure with sequence of values (for debugging)
  void iota(const T value) { std::iota(data.begin(), data.end(), value); }

  /// print grid (for debugging)
  void print(std::ostream &out = std::cout) const {
    CoordVec currCoord;
    currCoord.fill(0);
    out << extent << "\n";
    std::for_each(begin(), end(), [this, &out, &currCoord](const T &entry) {
      for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
        // if idx out of bounds adapt indices
        if (currCoord[idxDim] >= extent[idxDim]) {
          currCoord[idxDim + 1]++;
          for (SizeType idxResetDim = 0; idxResetDim <= idxDim; idxResetDim++)
            currCoord[idxResetDim] = 0;
          out << "\n";
          if (idxDim == 1)
            out << "\n";
        }
      }
      out << entry;
      if (currCoord[0] != extent[0] - 1)
        out << " ";
      currCoord[0]++;
    });
    out << "\n";
  }

  // --- iterator ----------------------------------------------------------
  typename std::vector<T>::const_iterator begin() const { return data.begin(); }
  typename std::vector<T>::const_iterator end() const { return data.end(); }
  typename std::vector<T>::iterator begin() { return data.begin(); }
  typename std::vector<T>::iterator end() { return data.end(); }

  // --- operator ----------------------------------------------------------
  /// accesses data structure at specific coordinate
  T &operator[](const CoordVec &coord) {
    checkCoord(coord);
    return data[coordToIdx(coord)];
  }

  /// accesses data structure at specific coordinate
  const T &operator[](const CoordVec &coord) const {
    checkCoord(coord);
    return data[coordToIdx(coord)];
  }

  // elementwise addition and substraction of data
  emcGrid<T, Dim> operator+=(const emcGrid<T, Dim> &rhs) {
    checkExtent(rhs);
    std::transform(rhs.begin(), rhs.end(), begin(), begin(), std::plus<T>());
    return *this;
  }

  emcGrid<T, Dim> operator-=(emcGrid<T, Dim> &rhs) {
    checkExtent(rhs);
    std::transform(rhs.begin(), rhs.end(), begin(), begin(), std::minus<T>());
    return *this;
  }

  // elementwise assignment of data
  emcGrid<T, Dim> &operator=(const emcGrid<T, Dim> &other) {
    checkExtent(other);
    std::copy(other.begin(), other.end(), begin());
    return *this;
  }

  /// helper that tells if coordinate is on boundary
  bool onBoundary(const CoordVec &coord) const {
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] == 0 || coord[idxDim] == extent[idxDim] - 1)
        return true;
    }
    return false;
  }

  /// helper that checks if both grids have same extent
  bool hasSameExtent(const emcGrid<T, Dim> &other) const {
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (other.getSize(idxDim) != extent[idxDim])
        return false;
    }
    return true;
  }

  /// helper that checks if a coordinate is on the grid (valid)
  bool isValid(const CoordVec &coord) const {
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      if (coord[idxDim] >= extent[idxDim])
        return false;
    }
    return true;
  }

private:
  /// resizes the container to a given extension
  void resize(const CoordVec &newExtent) {
    const auto newsize = std::accumulate(std::begin(extent), std::end(extent),
                                         SizeType{1}, std::multiplies<>{});
    data.resize(newsize);
  }

  /// helper that calculates 1D-index from a given coordinate
  SizeType coordToIdx(const CoordVec &coord) const {
    SizeType idx1D = 0;
    SizeType tmp;
    for (int idxDim = Dim - 1; idxDim >= 0; idxDim--) {
      tmp = coord[idxDim];
      for (int idxOtherDim = idxDim - 1; idxOtherDim >= 0; idxOtherDim--) {
        tmp *= extent[idxOtherDim];
      }
      idx1D += tmp;
    }
    assert(idx1D < getSize() || "Index of Grid out of bounds.");
    return idx1D;
  }

  /// helper that checks if index for dimension is valid
  void checkDimIdx(const SizeType idxDim) const {
    if (idxDim >= Dim) {
      emcMessage::getInstance()
          .addError("Index for Dimension is out of bounds!")
          .print();
    }
  }

  /// helper that checks if coordinate is valid
  void checkCoord(const CoordVec &coord) const {
    if (!isValid(coord)) {
      emcMessage::getInstance()
          .addError("Coordinate is out of bounds!")
          .print();
    }
  }

  /// helper that checks if two grids have same extent
  void checkExtent(const emcGrid<T, Dim> &other) const {
    if (!hasSameExtent(other)) {
      emcMessage::getInstance()
          .addError("Extents of the Grids don't match, can't perform operation "
                    "on them!")
          .print();
    }
  }

  /// helper that checks if requested neighbour exists
  void checkNbrExistence(const CoordVec &currCoord, SizeType idxDim,
                         int Dir) const {
    bool nbrExists = true;
    if (Dir == -1 && currCoord[idxDim] == 0)
      nbrExists = false;
    else if (Dir == 1 && currCoord[idxDim] == extent[idxDim] - 1)
      nbrExists = false;
    if (!nbrExists) {
      emcMessage::getInstance()
          .addError("Neighbour of coordinate is out of bounds!")
          .print();
    }
  }
}; // class emcGrid

// binary operators (+, -) componentwise
template <typename T, SizeType Dim>
emcGrid<T, Dim> operator+(const emcGrid<T, Dim> &lhs,
                          const emcGrid<T, Dim> &rhs) {
  if (!lhs.hasSameExtent(rhs)) {
    emcMessage::getInstance()
        .addError("Extents of the Grids don't match, can't perform operation "
                  "on them!")
        .print();
  }
  emcGrid<T, Dim> res = lhs;
  std::transform(rhs.begin(), rhs.end(), lhs.begin(), res.begin(),
                 std::plus<T>());
  return res;
}

template <typename T, SizeType Dim>
emcGrid<T, Dim> operator-(const emcGrid<T, Dim> &lhs,
                          const emcGrid<T, Dim> &rhs) {
  if (!lhs.hasSameExtent(rhs)) {
    emcMessage::getInstance()
        .addError("Extents of the Grids don't match, can't perform operation "
                  "on them!")
        .print();
  }
  emcGrid<T, Dim> res = lhs;
  std::transform(rhs.begin(), rhs.end(), lhs.begin(), res.begin(),
                 std::minus<T>());
  return res;
}

#endif // EMC_GRID_HPP