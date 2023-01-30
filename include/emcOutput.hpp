#ifndef EMC_OUTPUT_HPP
#define EMC_OUTPUT_HPP

#include <fstream>
#include <string>
#include <vector>

#include <emcGrid.hpp>
#include <emcUtil.hpp>

/// writes the given grid to a file in build-folder
/// NOTE: fileName should be given without extension
template <class T, SizeType Dim>
void writeToFile(const emcGrid<T, Dim> &grid, std::string fileName) {
  std::ofstream resultFile;
  resultFile.open(fileName + ".txt");
  grid.print(resultFile);
  resultFile.close();
}

/// writes the transformed grid to a file in build-folder
/// NOTE: fileName should be given without extension
template <class T, class DeviceType, SizeType Dim>
void writeToFile(const emcGrid<T, Dim> &grid, std::string fileName,
                 T (*func)(const T &, const DeviceType &device),
                 const DeviceType &device) {
  std::ofstream out;
  out.open(fileName + ".txt");
  std::array<T, Dim> currCoord;
  auto extent = grid.getExtent();
  currCoord.fill(0);
  out << extent << "\n";
  for (const auto &entry : grid) {
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
    out << func(entry, device);
    if (currCoord[0] != extent[0] - 1)
      out << " ";
    currCoord[0]++;
  }
  out << "\n";
  out.close();
}

/// writes input to a file in build-folder
/// NOTE: fileName should be given without extension
template <class T>
void writeToFile(const std::vector<std::vector<std::vector<int>>> &nettoNrPart,
                 const std::vector<std::vector<std::vector<T>>> &current,
                 SizeType idxType, T stepTime, T startTime,
                 std::string fileName) {
  std::ofstream out;
  out.open(fileName + ".txt");
  for (SizeType nrStep = 1; nrStep <= current.size(); nrStep++) {
    out << startTime + nrStep * stepTime;
    for (auto &nrPart : nettoNrPart[nrStep - 1][idxType])
      out << " " << nrPart;
    for (auto &curr : current[nrStep - 1][idxType])
      out << " " << curr;
    out << "\n";
  }
  out.close();
}

/// standard grid transform functions
template <class T, class DeviceType>
T undoNormalizationPotential(const T &val, const DeviceType &device) {
  return device.undoNormalizeVoltage(val);
}

template <class T, class DeviceType>
T undoNormalizationConcentration(const T &val, const DeviceType &device) {
  return device.undoNormalizeDoping(val);
}

#endif // EMC_OUTPUT_HPP