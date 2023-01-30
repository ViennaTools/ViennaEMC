#include <PMSchemes/emcAbstractPMScheme.hpp>
#include <emcGrid.hpp>
#include <emcMessage.hpp>
#include <emcUtil.hpp>

#include <array>
#include <vector>

/*! \brief Nearest - Element - Center Scheme as done in ViennaWD.
 *
 * Adaptations to normal emcNECScheme.hpp:
 *  - boundary conditions of electric field different
 *  (- rounding of x-position of force)
 *
 * Note: Only implemented for 2D!
 *
 */
template <class T, class DeviceType>
class emcNECSchemeVWD : public emcAbstractPMScheme<T, DeviceType> {
public:
  static const SizeType Dim = DeviceType::Dimension;

  static_assert(Dim == 2, "PMScheme NEC-VWD is only implemented for 2D.");

  /// @brief 2D single particle assignment
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

  /// @brief 2D multiple particle assignment
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

  /// @brief force interpolation as done in VWD.
  /// Diff to emc: round x-position, instead of floor
  /// (assigned to next element center in x-direction)
  std::array<T, 3> interpolateForce(const std::vector<emcGrid<T, 2>> &eField,
                                    const std::array<T, 2> &pos,
                                    const std::array<T, 2> &spacing,
                                    T charge) const {
    SizeType minX = std::round(pos[0] / spacing[0]);
    // SizeType minX = std::floor(pos[0] / spacing[0]); // EMC-version
    SizeType minY = std::floor(pos[1] / spacing[1]);

    std::array<T, 3> force;
    force[0] =
        charge * (eField[0][{minX, minY}] + eField[0][{minX, minY + 1}]) / 2;
    if (minX != eField[1].getSize(0) - 1) {
      force[1] =
          charge * (eField[1][{minX, minY}] + eField[1][{minX + 1, minY}]) / 2;
    } else {
      force[1] = charge * eField[1][{minX, minY}];
    }
    force[2] = 0;
    return force;
  }

  /// @brief calculation of electric field as performed in VWD.
  /// Diff to emc: boundary conditions for E-Field
  /// eField at minBoundaries calculated normally.
  /// eField at maxBoundaries set to 0.
  void calcEField(std::vector<emcGrid<T, Dim>> &eField,
                  const emcGrid<T, Dim> &potential,
                  const DeviceType &device) const {
    auto h = device.getSpacing();
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
      std::array<SizeType, Dim> coord;
      for (coord.fill(0); !device.isEndCoord(coord);
           device.advanceCoord(coord)) {
        if (coord[idxDim] != eField[idxDim].getSize(idxDim) - 1) {
          eField[idxDim][coord] =
              device.undoNormalizeVoltage(
                  potential[coord] - potential.getNextValue(coord, idxDim)) /
              h[idxDim];
        } else
          eField[idxDim][coord] = 0;
      }
    }
  }
};