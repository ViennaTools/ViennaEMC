#ifndef EMC_MULTIVALLEY_PAULI_EXCLUSION_HPP
#define EMC_MULTIVALLEY_PAULI_EXCLUSION_HPP

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include <emcConstants.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

/**
 * Per-valley / per-subvalley Lugli-Ferry k-space occupancy grid for
 * Pauli-exclusion rejection in a MULTIVALLEY, single-layer (2D) EMC.
 *
 * Why a dedicated class: in the analytic multivalley band model every valley
 * measures its wave vector from its OWN minimum (all valleys are "centred at
 * k = 0" in particle.k). A single global occupancy grid (emcPauliExclusion)
 * therefore conflates all valleys/subvalleys into the same k-slots and blocks
 * scattering ACROSS valleys, which is unphysical. This class keeps a SEPARATE
 * 2D (kx, ky) occupancy grid for each (valley, subvalley), so the Pauli
 * blocking (Burstein-Moss / band filling) is applied within each valley's own
 * Fermi sea. Carriers have k_z = 0, so the grid is 2D.
 *
 * Capacity (max SIMULATED particles, each = one real carrier, per subvalley per
 * k-slot of area dk^2), including spin:
 *
 *   N_c = floor( 2 * A_sim * dk^2 / (2 pi)^2 )
 *
 * N_c must be >= ~5 for meaningful blocking statistics, and the simulated
 * carrier count per subvalley must be consistent with the intended density so
 * that the Fermi sea fills the slots to ~N_c (otherwise the grid over- or
 * under-blocks).
 *
 * Reference: Lugli & Ferry, IEEE TED 32 (1985) 2431.
 *
 * @tparam T floating-point type
 */
template <class T> class emcMultiValleyPauliExclusion {
public:
  T dk;       // k-bin side length [m^-1]
  int nBins;  // grid dimension per in-plane axis
  int offset; // index offset so bin(k) = floor(k/dk) + offset >= 0
  int N_c;    // max simulated particles per (subvalley, k-slot)

  SizeType nRejected = 0;
  SizeType nScattered = 0;

  emcMultiValleyPauliExclusion() = delete;

  /**
   * @param inDk k-bin width [m^-1]
   * @param kMax maximum |kx|,|ky| to represent [m^-1]
   * @param Asim simulation AREA [m^2]
   * @param valleyDegeneracies number of subvalleys in each valley (e.g. {6, 6})
   */
  emcMultiValleyPauliExclusion(T inDk, T kMax, T Asim,
                               const std::vector<SizeType> &valleyDegeneracies)
      : dk(inDk) {
    int halfBins = static_cast<int>(std::ceil(kMax / dk)) + 1;
    nBins = 2 * halfBins + 1;
    offset = halfBins;

    // per-subvalley 2D capacity (spin factor 2)
    T cap = T(2) * Asim * dk * dk / (T(4) * constants::pi * constants::pi);
    N_c = std::max(1, static_cast<int>(std::floor(cap)));

    // linear (valley, subvalley) key offsets
    keyOffset.assign(valleyDegeneracies.size(), 0);
    SizeType total = 0;
    for (SizeType v = 0; v < valleyDegeneracies.size(); ++v) {
      keyOffset[v] = total;
      total += valleyDegeneracies[v];
    }
    occ.assign(total, std::vector<int>(
                           static_cast<std::size_t>(nBins) * nBins, 0));

    std::cout << "MultiValleyPauliExclusion: dk = " << dk * 1e-7 << "e7 m^-1"
              << "  nBins = " << nBins << "^2"
              << "  N_c = " << N_c << " per (subvalley, slot)"
              << "  nKeys = " << total << "\n";
    if (N_c < 5)
      std::cerr << "  WARNING: N_c = " << N_c
                << " < 5 — box too small for reliable Pauli blocking.\n";
  }

  /** Rebuild all occupancy grids from the current particles. */
  template <class PartVec> void buildGrid(const PartVec &particles) {
    for (auto &g : occ)
      std::fill(g.begin(), g.end(), 0);
    for (const auto &p : particles)
      occ[key(p.valley, p.subValley)][flatIdx(p.k)] += 1;
  }

  /** True if the (valley, subvalley) k-slot of kNew is at/above capacity. */
  bool isBlocked(SizeType valley, SizeType subValley,
                 const std::array<T, 3> &kNew) const {
    return occ[key(valley, subValley)][flatIdx(kNew)] >= N_c;
  }

  /** Move one particle between (valley, subvalley, k) slots after an accepted
   *  scatter. */
  void update(SizeType valleyOld, SizeType subOld,
              const std::array<T, 3> &kOld, SizeType valleyNew, SizeType subNew,
              const std::array<T, 3> &kNew) {
    int &o = occ[key(valleyOld, subOld)][flatIdx(kOld)];
    o = std::max(0, o - 1);
    occ[key(valleyNew, subNew)][flatIdx(kNew)] += 1;
  }

  T getRejectionRate() const {
    return nScattered > 0 ? T(nRejected) / T(nScattered) : T(0);
  }
  void resetCounters() {
    nRejected = 0;
    nScattered = 0;
  }

private:
  std::vector<SizeType> keyOffset;   // start key index for each valley
  std::vector<std::vector<int>> occ; // occ[key][ix*nBins + iy]

  SizeType key(SizeType valley, SizeType subValley) const {
    return keyOffset[valley] + subValley;
  }

  int axisIdx(T kcomp) const {
    int i = static_cast<int>(std::floor(kcomp / dk)) + offset;
    return std::max(0, std::min(nBins - 1, i));
  }

  int flatIdx(const std::array<T, 3> &k) const {
    return axisIdx(k[0]) * nBins + axisIdx(k[1]);
  }
};

#endif // EMC_MULTIVALLEY_PAULI_EXCLUSION_HPP
