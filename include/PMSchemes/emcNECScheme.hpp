#ifndef EMC_NEC_SCHEME_HPP
#define EMC_NEC_SCHEME_HPP

#include <PMSchemes/emcAbstractPMScheme.hpp>
#include <PMSchemes/emcEFieldCalculation.hpp>
#include <emcGrid.hpp>
#include <emcMessage.hpp>
#include <emcUtil.hpp>

#include <array>
#include <vector>

/*! \brief Nearest - Element - Center Scheme
 *
 * Assigns parts of the particle equally to points of the element
 * with the nearest center. Force is interpolated by averaging the
 * force of the points of the element with the neares center.
 * Electric field has to be calculated at edge midpoints.
 *
 * Note: Only implemented for 2D!
 *
 */
template <class T, class DeviceType>
class emcNECScheme : public emcAbstractPMScheme<T, DeviceType> {
public:
  static const SizeType Dim = DeviceType::Dimension;

  static_assert(Dim == 2, "PMScheme NEC is only implemented for 2D.");

  // 2D single particle assignment equally to 4 nearest grid points
  void assignToMesh(const std::array<T, 2> &pos, SizeType nrCarriers,
                    const std::array<T, 2> &spacing,
                    emcGrid<T, 2> &gridNrParticles) const {
    SizeType minX = std::floor(pos[0] / spacing[0]);
    SizeType minY = std::floor(pos[1] / spacing[1]);
    T weight = 0.25 * nrCarriers;
    gridNrParticles[{minX, minY}] += weight;
    gridNrParticles[{minX + 1, minY}] += weight;
    gridNrParticles[{minX, minY + 1}] += weight;
    gridNrParticles[{minX + 1, minY + 1}] += weight;
  }

  // 3D single particle assignment equally to 8 nearest grid points
  void assignToMesh(const std::array<T, 3> &pos, SizeType nrCarriers,
                    const std::array<T, 3> &spacing,
                    emcGrid<T, 3> &gridNrParticles) const {
    SizeType minX = std::floor(pos[0] / spacing[0]);
    SizeType minY = std::floor(pos[1] / spacing[1]);
    SizeType minZ = std::floor(pos[2] / spacing[2]);
    T weight = 0.125 * nrCarriers;
    gridNrParticles[{minX, minY, minZ}] += weight;
    gridNrParticles[{minX + 1, minY, minZ}] += weight;
    gridNrParticles[{minX, minY + 1, minZ}] += weight;
    gridNrParticles[{minX + 1, minY + 1, minZ}] += weight;
    gridNrParticles[{minX, minY, minZ + 1}] += weight;
    gridNrParticles[{minX + 1, minY, minZ + 1}] += weight;
    gridNrParticles[{minX, minY + 1, minZ + 1}] += weight;
    gridNrParticles[{minX + 1, minY + 1, minZ + 1}] += weight;
  }

  // 2D multiple particle assignment equally to 4 nearest grid points
  void assignToMesh(const std::vector<std::array<T, 2>> &position,
                    SizeType nrCarriers, const std::array<T, 2> &spacing,
                    emcGrid<T, 2> &gridNrParticles) const {
    SizeType minX, minY;
    T weight = 0.25 * nrCarriers;
    for (const auto &pos : position) {
      minX = std::floor(pos[0] / spacing[0]);
      minY = std::floor(pos[1] / spacing[1]);
      gridNrParticles[{minX, minY}] += weight;
      gridNrParticles[{minX + 1, minY}] += weight;
      gridNrParticles[{minX, minY + 1}] += weight;
      gridNrParticles[{minX + 1, minY + 1}] += weight;
    }
  }

  // 3D multiple particle assignment to 8 nearest grid points
  void assignToMesh(const std::vector<std::array<T, 3>> &position,
                    SizeType nrCarriers, const std::array<T, 3> &spacing,
                    emcGrid<T, 3> &gridNrParticles) const {
    SizeType minX, minY, minZ;
    T weight = 0.125 * nrCarriers;
    for (const auto &pos : position) {
      minX = std::floor(pos[0] / spacing[0]);
      minY = std::floor(pos[1] / spacing[1]);
      minZ = std::floor(pos[2] / spacing[2]);
      gridNrParticles[{minX, minY, minZ}] += weight;
      gridNrParticles[{minX + 1, minY, minZ}] += weight;
      gridNrParticles[{minX, minY + 1, minZ}] += weight;
      gridNrParticles[{minX + 1, minY + 1, minZ}] += weight;
      gridNrParticles[{minX, minY, minZ + 1}] += weight;
      gridNrParticles[{minX + 1, minY, minZ + 1}] += weight;
      gridNrParticles[{minX, minY + 1, minZ + 1}] += weight;
      gridNrParticles[{minX + 1, minY + 1, minZ + 1}] += weight;
    }
  }

  // 2D force interpolation
  std::array<T, 3> interpolateForce(const std::vector<emcGrid<T, 2>> &eField,
                                    const std::array<T, 2> &pos,
                                    const std::array<T, 2> &spacing,
                                    T charge) const {
    SizeType minX = std::floor(pos[0] / spacing[0]);
    SizeType minY = std::floor(pos[1] / spacing[1]);

    std::array<T, 3> force;
    force[0] =
        charge * (eField[0][{minX, minY}] + eField[0][{minX, minY + 1}]) / 2;
    force[1] =
        charge * (eField[1][{minX, minY}] + eField[1][{minX + 1, minY}]) / 2;
    force[2] = 0;
    return force;
  }

  /// calculation of EField at edge midpoints (from potential)
  void calcEField(std::vector<emcGrid<T, Dim>> &eField,
                  const emcGrid<T, Dim> &potential,
                  const DeviceType &device) const {
    calcEFieldAtEdgeMidPts(eField, potential, device);
  }
};

#endif // EMC_NEC_SCHEME_HPP