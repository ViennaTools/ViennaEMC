#ifndef EFIELD_CALCULATION_HPP
#define EFIELD_CALCULATION_HPP

#include <vector>

#include <emcGrid.hpp>
#include <emcSurface.hpp>
#include <emcUtil.hpp>

/// calculates negative derivative of potential (el. field)
/// at grid points
template <template <typename, SizeType> class DeviceType, class T, SizeType Dim>
void calcEFieldAtGridPts(std::vector<emcGrid<T, Dim>> &eField,
                         const emcGrid<T, Dim> &potential,
                         const DeviceType<T, Dim> &device) {
  auto h = device.getSpacing();
  for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
    std::array<SizeType, Dim> coord;
    for (coord.fill(0); !device.isEndCoord(coord); device.advanceCoord(coord)) {
      if (coord[idxDim] != 0 &&
          coord[idxDim] != eField[idxDim].getSize(idxDim) - 1) {
        eField[idxDim][coord] =
            device.undoNormalizeVoltage(potential.getPrevValue(coord, idxDim) -
                                        potential.getNextValue(coord, idxDim)) /
            (2 * h[idxDim]);
      }
    }
    setEFieldBoundaryValues(eField, idxDim, device);
  }
}

/// calculates negative derivative of potential (el. field)
/// at edge midpoints (needed for NEC)
template <template <typename, SizeType> class DeviceType, class T, SizeType Dim>
void calcEFieldAtEdgeMidPts(std::vector<emcGrid<T, Dim>> &eField,
                            const emcGrid<T, Dim> &potential,
                            const DeviceType<T, Dim> &device) {
  auto h = device.getSpacing();
  for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
    std::array<SizeType, Dim> coord;
    for (coord.fill(0); !device.isEndCoord(coord); device.advanceCoord(coord)) {
      if (coord[idxDim] != 0 &&
          coord[idxDim] != eField[idxDim].getSize(idxDim) - 1) {
        eField[idxDim][coord] =
            device.undoNormalizeVoltage(potential[coord] -
                                        potential.getNextValue(coord, idxDim)) /
            h[idxDim];
      }
    }
    setEFieldBoundaryValues(eField, idxDim, device);
  }
}

/// sets the boundary values for the electric field.
/// normal to artificial boundary -> E = 0.
/// normal to contact -> E = E(nearest grid point).
template <template <typename, SizeType> class DeviceType, class T, SizeType Dim>
void setEFieldBoundaryValues(std::vector<emcGrid<T, Dim>> &eField,
                             SizeType idxDim,
                             const DeviceType<T, Dim> &device) {
  for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
    std::array<SizeType, Dim - 1> coordSurf;
    emcBoundaryPos boundPos;
    auto &surface = device.getSurface();
    for (surface.initCoord(boundPos, coordSurf);
         !surface.isEndCoord(boundPos, coordSurf);
         surface.advanceCoord(boundPos, coordSurf)) {
      auto coordDev = surface.getCoordDevice(coordSurf, boundPos);
      auto idxDim = toUnderlying(boundPos) / 2;
      if (surface.isArtificialBoundary(coordSurf, boundPos)) {
        eField[idxDim][coordDev] = 0;
      } else if (toUnderlying(boundPos) % 2 == 0) {
        eField[idxDim][coordDev] =
            eField[idxDim].getNextValue(coordDev, idxDim);
      } else {
        eField[idxDim][coordDev] =
            eField[idxDim].getPrevValue(coordDev, idxDim);
      }
    }
  }
}

#endif // EFIELD_CALCULATION_HPP