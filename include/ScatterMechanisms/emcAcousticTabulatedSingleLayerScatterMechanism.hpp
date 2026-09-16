#ifndef ACOUSTIC_TABULATED_SINGLE_LAYER_SCATTER_MECHANISM_HPP
#define ACOUSTIC_TABULATED_SINGLE_LAYER_SCATTER_MECHANISM_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \brief Intravalley acoustic phonon scattering of a single layer with a
 * wavevector-dependent deformation potential Xi(q).
 *
 * Same physics as emcAcousticSingleLayerMechanism (elastic, equipartition,
 * 2D, isotropic parabolic or non-parabolic valley), but the deformation
 * potential is a table Xi(q) [eV] over the phonon wavevector q [1/m] (linear
 * interpolation, constant beyond both ends), as obtained from first-principles
 * electron-phonon matrix elements (e.g. the Wannier interpolation of the LA and
 * TA couplings of monolayer MoS2, which vanish or saturate with q). For elastic
 * scattering on a circle of radius k the momentum transfer is
 * q = 2 k sin(theta / 2), so the rate is the angular average of Xi(q)^2 and
 * the deflection angle theta is sampled from the distribution ~ Xi(q(theta))^2
 * by numerical inversion of the cumulative distribution. With a constant table
 * this class reduces exactly to emcAcousticSingleLayerMechanism.
 *
 * @param scatterConst pre-computed constant (rate for Xi = 1 eV)
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcAcousticTabulatedSingleLayerMechanism : public emcScatterMechanism<T> {
private:
  T scatterConst;
  std::vector<T> tabQ, tabXi;
  mutable std::uniform_real_distribution<T> dist;
  std::string nameSuffix;

  /// Xi(q) [eV] from the table
  T deformationPotential(T q) const {
    if (q <= tabQ.front())
      return tabXi.front();
    if (q >= tabQ.back())
      return tabXi.back();
    auto it = std::upper_bound(tabQ.begin(), tabQ.end(), q);
    SizeType i = it - tabQ.begin();
    T w = (q - tabQ[i - 1]) / (tabQ[i] - tabQ[i - 1]);
    return tabXi[i - 1] + w * (tabXi[i] - tabXi[i - 1]);
  }

public:
  emcAcousticTabulatedSingleLayerMechanism() = delete;

  /*! \brief Constructor.
   * @param inValley index of the valley the mechanism is assigned to
   * @param table (q [1/m], Xi [eV]) pairs in ascending q
   * @param densityMaterial mass density [kg / m^2]
   * @param velSound sound velocity of the branch [m / s]
   * @param temperature lattice temperature [K]
   * @param inNameSuffix suffix for filename
   */
  emcAcousticTabulatedSingleLayerMechanism(
      SizeType inValley, const std::vector<std::pair<T, T>> &table,
      T densityMaterial, T velSound, T temperature,
      std::string inNameSuffix = "")
      : dist(0., 1.), nameSuffix(inNameSuffix), emcScatterMechanism<T>(inValley) {
    for (const auto &qx : table) {
      tabQ.push_back(qx.first);
      tabXi.push_back(qx.second);
    }
    T cL = densityMaterial * std::pow(velSound, 2);
    scatterConst = std::pow(constants::q, 2) * constants::kB * temperature /
                   (cL * std::pow(constants::hbar, 3));
  }

  std::string getName() const { return "AcousticTabSL" + nameSuffix; }

  /// rate = (angular average of Xi(q)^2) x rate for Xi = 1 eV
  T getScatterRate(T energy, SizeType /*region*/) const {
    auto valley = this->ptrValley[this->idxValley];
    T md = valley->getEffMassDOS();
    T alpha = valley->getNonParabolicity();
    T k = valley->getNormWaveVec(energy);
    constexpr SizeType N = 64;
    T sum = 0;
    for (SizeType i = 0; i < N; ++i) {
      T theta = (i + T(0.5)) * constants::pi / N; // symmetric in theta -> [0, pi]
      T xi = deformationPotential(2 * k * std::sin(theta / 2));
      sum += xi * xi;
    }
    return md * scatterConst * (1 + 2 * alpha * energy) * sum / N;
  }

  /// \brief Elastic intravalley scattering; the deflection angle follows
  /// Xi(q(theta))^2, isotropic for a constant table.
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto currValley = this->ptrValley[particle.valley];
    T k = currValley->getNormWaveVec(particle.energy);
    constexpr SizeType N = 64;
    std::array<T, N + 1> cdf;
    cdf[0] = 0;
    for (SizeType i = 1; i <= N; ++i) {
      T theta = (i - T(0.5)) * constants::pi / N;
      T xi = deformationPotential(2 * k * std::sin(theta / 2));
      cdf[i] = cdf[i - 1] + xi * xi;
    }
    T psi;
    if (!(cdf[N] > 0)) {
      psi = 2 * constants::pi * dist(rng);
    } else {
      T target = dist(rng) * cdf[N];
      SizeType lo = 1;
      while (lo < N && cdf[lo] < target)
        ++lo;
      T frac = (target - cdf[lo - 1]) / (cdf[lo] - cdf[lo - 1]);
      T psiMag = (T(lo) - 1 + frac) * constants::pi / N;
      psi = (dist(rng) < T(0.5)) ? psiMag : (2 * constants::pi - psiMag);
    }
    // rotate the in-plane wavevector by psi (isotropic valley frame)
    T phi = std::atan2(particle.k[1], particle.k[0]) + psi;
    auto vogtFactor = currValley->getVogtTransformationFactor();
    std::array<T, 3> kk{std::cos(phi) / vogtFactor[0], std::sin(phi) / vogtFactor[1], 0};
    T factor = 1. / std::sqrt(kk[0] * kk[0] + kk[1] * kk[1]);
    kk[0] *= factor * k;
    kk[1] *= factor * k;
    particle.k = kk;
  }
};

#endif // ACOUSTIC_TABULATED_SINGLE_LAYER_SCATTER_MECHANISM_HPP
