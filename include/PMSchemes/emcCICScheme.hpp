#ifndef EMC_CIC_SCHEME_HPP
#define EMC_CIC_SCHEME_HPP

#include <PMSchemes/emcAbstractPMScheme.hpp>
#include <PMSchemes/emcEFieldCalculation.hpp>
#include <emcGrid.hpp>
#include <emcMessage.hpp>
#include <emcUtil.hpp>

#include <array>
#include <vector>

/*! \brief Cloud - In - Cell Scheme
 *
 * Assigns parts of the particle to nearest points in each direction (4 nearest
 * points for 2D and 8 for 3D), weighted by distance to the points. Force is
 * also calculated as weighted average of the nearest points. Electric field is
 * calculated at grid points.
 *
 * Note: Only implemented for 2D and 3D!
 *
 */
template <class T, class DeviceType>
class emcCICScheme : public emcAbstractPMScheme<T, DeviceType> {
public:
  static const SizeType Dim = DeviceType::Dimension;

  static_assert(Dim == 2 || Dim == 3,
                "PMScheme CIC is only implemented for 2D and 3D.");

  // 2D single particle assignment to 4 nearest points weighted by distance
  void assignToMesh(const std::array<T, 2> &pos, SizeType nrCarriers,
                    const std::array<T, 2> &spacing,
                    emcGrid<T, 2> &gridNrParticles) const {
    SizeType minX = std::floor(pos[0] / spacing[0]);
    SizeType minY = std::floor(pos[1] / spacing[1]);
    T wX = pos[0] / spacing[0] - minX;
    T wY = pos[1] / spacing[1] - minY;
    gridNrParticles[{minX, minY}] += wX * wY * nrCarriers;
    gridNrParticles[{minX + 1, minY}] += (1 - wX) * wY * nrCarriers;
    gridNrParticles[{minX, minY + 1}] += wX * (1 - wY) * nrCarriers;
    gridNrParticles[{minX + 1, minY + 1}] += (1 - wX) * (1 - wY) * nrCarriers;
  }

  // 3D single particle assignment to 8 nearest points weighted by distance
  void assignToMesh(const std::array<T, 3> &pos, SizeType nrCarriers,
                    const std::array<T, 3> &spacing,
                    emcGrid<T, 3> &gridNrParticles) const {
    SizeType minX = std::floor(pos[0] / spacing[0]);
    SizeType minY = std::floor(pos[1] / spacing[1]);
    SizeType minZ = std::floor(pos[2] / spacing[2]);
    T wX = pos[0] / spacing[0] - minX;
    T wY = pos[1] / spacing[1] - minY;
    T wZ = pos[2] / spacing[2] - minZ;
    gridNrParticles[{minX, minY, minZ}] += wX * wY * wZ * nrCarriers;
    gridNrParticles[{minX + 1, minY, minZ}] += (1 - wX) * wY * wZ * nrCarriers;
    gridNrParticles[{minX, minY + 1, minZ}] += wX * (1 - wY) * wZ * nrCarriers;
    gridNrParticles[{minX + 1, minY + 1, minZ}] +=
        (1 - wX) * (1 - wY) * wZ * nrCarriers;
    gridNrParticles[{minX, minY, minZ + 1}] += wX * wY * (1 - wZ) * nrCarriers;
    gridNrParticles[{minX + 1, minY, minZ + 1}] +=
        (1 - wX) * wY * (1 - wZ) * nrCarriers;
    gridNrParticles[{minX, minY + 1, minZ + 1}] +=
        wX * (1 - wY) * (1 - wZ) * nrCarriers;
    gridNrParticles[{minX + 1, minY + 1, minZ + 1}] +=
        (1 - wX) * (1 - wY) * (1 - wZ) * nrCarriers;
  }

  // 2D multiple particle assignment to 4 nearest grid points weighted by
  // distance
  void assignToMesh(const std::vector<std::array<T, 2>> &position,
                    SizeType nrCarriers, const std::array<T, 2> &spacing,
                    emcGrid<T, 2> &gridNrParticles) const {
    SizeType minX, minY;
    T wX, wY;
    for (const auto &pos : position) {
      minX = std::floor(pos[0] / spacing[0]);
      minY = std::floor(pos[1] / spacing[1]);
      wX = pos[0] / spacing[0] - minX;
      wY = pos[1] / spacing[1] - minY;
      gridNrParticles[{minX, minY}] += wX * wY * nrCarriers;
      gridNrParticles[{minX + 1, minY}] += (1 - wX) * wY * nrCarriers;
      gridNrParticles[{minX, minY + 1}] += wX * (1 - wY) * nrCarriers;
      gridNrParticles[{minX + 1, minY + 1}] += (1 - wX) * (1 - wY) * nrCarriers;
    }
  }

  // 3D multiple particle assignment to 8 nearest grid points weighted by
  // distance
  void assignToMesh(const std::vector<std::array<T, 3>> &position,
                    SizeType nrCarriers, const std::array<T, 3> &spacing,
                    emcGrid<T, 3> &gridNrParticles) const {
    SizeType minX, minY, minZ;
    T wX, wY, wZ;
    for (const auto &pos : position) {
      minX = std::floor(pos[0] / spacing[0]);
      minY = std::floor(pos[1] / spacing[1]);
      minZ = std::floor(pos[2] / spacing[2]);
      wX = pos[0] / spacing[0] - minX;
      wY = pos[1] / spacing[1] - minY;
      wZ = pos[2] / spacing[2] - minZ;
      gridNrParticles[{minX, minY, minZ}] += wX * wY * wZ * nrCarriers;
      gridNrParticles[{minX + 1, minY, minZ}] +=
          (1 - wX) * wY * wZ * nrCarriers;
      gridNrParticles[{minX, minY + 1, minZ}] +=
          wX * (1 - wY) * wZ * nrCarriers;
      gridNrParticles[{minX + 1, minY + 1, minZ}] +=
          (1 - wX) * (1 - wY) * wZ * nrCarriers;
      gridNrParticles[{minX, minY, minZ + 1}] +=
          wX * wY * (1 - wZ) * nrCarriers;
      gridNrParticles[{minX + 1, minY, minZ + 1}] +=
          (1 - wX) * wY * (1 - wZ) * nrCarriers;
      gridNrParticles[{minX, minY + 1, minZ + 1}] +=
          wX * (1 - wY) * (1 - wZ) * nrCarriers;
      gridNrParticles[{minX + 1, minY + 1, minZ + 1}] +=
          (1 - wX) * (1 - wY) * (1 - wZ) * nrCarriers;
    }
  }

  /// 2D force interpolation using eField-value of 4 nearest grid points
  /// weighted by distance
  std::array<T, 3> interpolateForce(const std::vector<emcGrid<T, 2>> &eField,
                                    const std::array<T, 2> &pos,
                                    const std::array<T, 2> &spacing,
                                    T charge) const {
    SizeType minX = std::floor(pos[0] / spacing[0]);
    SizeType minY = std::floor(pos[1] / spacing[1]);
    T wX = pos[0] / spacing[0] - minX;
    T wY = pos[1] / spacing[1] - minY;

    std::array<T, 3> force;
    for (SizeType idxDim = 0; idxDim < 2; idxDim++) {
      force[idxDim] =
          charge * (eField[idxDim][{minX, minY}] * wX * wY +
                    eField[idxDim][{minX + 1, minY}] * (1 - wX) * wY +
                    eField[idxDim][{minX, minY + 1}] * wX * (1 - wY) +
                    eField[idxDim][{minX + 1, minY + 1}] * (1 - wY) * (1 - wY));
    }
    force[2] = 0;
    return force;
  }

  /// 3D force interpolation using eField-value of 8 nearest grid points
  /// weighted by distance
  std::array<T, 3> interpolateForce(const std::vector<emcGrid<T, 3>> &eField,
                                    const std::array<T, 3> &pos,
                                    const std::array<T, 3> &spacing,
                                    T charge) const {
    SizeType minX = std::floor(pos[0] / spacing[0]);
    SizeType minY = std::floor(pos[1] / spacing[1]);
    SizeType minZ = std::floor(pos[2] / spacing[2]);
    T wX = pos[0] / spacing[0] - minX;
    T wY = pos[1] / spacing[1] - minY;
    T wZ = pos[2] / spacing[2] - minZ;
    std::array<T, 3> force;
    for (SizeType idxDim = 0; idxDim < 3; idxDim++) {
      force[idxDim] =
          (eField[idxDim][{minX, minY, minZ}] * wX * wY * wZ +
           eField[idxDim][{minX + 1, minY, minZ}] * (1 - wX) * wY * wZ +
           eField[idxDim][{minX, minY + 1, minZ}] * wX * (1 - wY) * wZ +
           eField[idxDim][{minX + 1, minY + 1, minZ}] * (1 - wY) * (1 - wY) *
               wZ +
           eField[idxDim][{minX, minY, minZ + 1}] * wX * wY * (1 - wZ) +
           eField[idxDim][{minX + 1, minY, minZ + 1}] * (1 - wX) * wY *
               (1 - wZ) +
           eField[idxDim][{minX, minY + 1, minZ + 1}] * wX * (1 - wY) *
               (1 - wZ) +
           eField[idxDim][{minX + 1, minY + 1, minZ + 1}] * (1 - wY) *
               (1 - wY) * (1 - wZ)) *
          charge;
    }
    return force;
  }

  /// calculation of EField at grid points (from potential)
  void calcEField(std::vector<emcGrid<T, Dim>> &eField,
                  const emcGrid<T, Dim> &potential,
                  const DeviceType &device) const {
    calcEFieldAtGridPts(eField, potential, device);
  }
};

#endif // EMC_CIC_SCHEME_HPP