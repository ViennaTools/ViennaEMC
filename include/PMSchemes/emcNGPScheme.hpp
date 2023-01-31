#ifndef EMC_NGP_SCHEME_HPP
#define EMC_NGP_SCHEME_HPP

#include <PMSchemes/emcAbstractPMScheme.hpp>
#include <PMSchemes/emcEFieldCalculation.hpp>
#include <emcGrid.hpp>
#include <emcUtil.hpp>

#include <array>
#include <vector>

/*! \brief Nearest - Grid - Point Scheme
 *
 * Assigns each particle to nearest grid point, uses
 * force of nearest grid point anc calculates el field
 * at grid points.
 */
template <class T, class DeviceType>
class emcNGPScheme : public emcAbstractPMScheme<T, DeviceType> {
public:
  static const SizeType Dim = DeviceType::Dimension;

  /// single particle assignment to nearest grid point
  void assignToMesh(const std::array<T, Dim> &pos, SizeType nrCarriers,
                    const std::array<T, Dim> &spacing,
                    emcGrid<T, Dim> &gridNrParticles) const {
    std::array<SizeType, Dim> coord;
    std::transform(pos.begin(), pos.end(), spacing.begin(), coord.begin(),
                   [](const T &val, const T &h) -> SizeType {
                     return std::round(val / h);
                   });
    gridNrParticles[coord] += nrCarriers;
  }

  /// multiple particle assignment to nearest grid point
  void assignToMesh(const std::vector<std::array<T, Dim>> &position,
                    SizeType nrCarriers, const std::array<T, Dim> &spacing,
                    emcGrid<T, Dim> &gridNrParticles) const {
    std::array<SizeType, Dim> coord;
    for (const auto &pos : position) {
      std::transform(pos.begin(), pos.end(), spacing.begin(), coord.begin(),
                     [](const T &val, const T &h) -> SizeType {
                       return std::round(val / h);
                     });
      gridNrParticles[coord] += nrCarriers;
    }
  }

  /// force interpolation at particle position using eField-value of nearest
  /// grid point
  std::array<T, 3> interpolateForce(const std::vector<emcGrid<T, Dim>> &eField,
                                    const std::array<T, Dim> &pos,
                                    const std::array<T, Dim> &spacing,
                                    T charge) const {
    std::array<T, 3> force;
    if (Dim == 2)
      force[2] = 0;
    std::array<SizeType, Dim> coord;
    std::transform(pos.begin(), pos.end(), spacing.begin(), coord.begin(),
                   [spacing](const T &val, const T &h) -> SizeType {
                     return std::round(val / h);
                   });
    for (SizeType idxDim = 0; idxDim < Dim; idxDim++)
      force[idxDim] = charge * eField[idxDim][coord];
    return force;
  }

  /// calculation of EField at grid points (from potential)
  void calcEField(std::vector<emcGrid<T, Dim>> &eField,
                  const emcGrid<T, Dim> &potential,
                  const DeviceType &device) const {
    calcEFieldAtGridPts(eField, potential, device);
  }
};

#endif // EMC_NGP_SCHEME_HPP