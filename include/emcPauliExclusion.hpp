#ifndef EMC_PAULI_EXCLUSION_HPP
#define EMC_PAULI_EXCLUSION_HPP

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include <emcConstants.hpp>
#include <emcUtil.hpp>

/**
 * Lugli-Ferry k-space occupancy grid for Pauli-exclusion rejection in EMC.
 *
 * After each phonon-scatter event the destination k-slot is checked: if it
 * already holds N_c simulated particles (the Pauli capacity of that slot) the
 * scatter is rejected and treated as self-scatter, implementing the
 * Burstein-Moss / band-filling effect.
 *
 * Algorithm: Lugli & Ferry, IEEE TED 32 (1985) 2431; applied to MHPs in
 * Faber, van de Ven, Filipovic, Loi, Koster (2025, in preparation).
 *
 * Grid sizing
 * ──────────
 * The k-space is divided into cubic bins of side dk [m^-1].  The maximum
 * number of SIMULATED particles (each representing one real carrier) that
 * can occupy a bin is:
 *
 *   N_c = floor( 2 · V_sim · dk³ / (8π³) )    [eq. 4 of Faber et al. 2025]
 *
 * where the factor 2 accounts for spin degeneracy.  N_c must be ≥ 5 for
 * meaningful blocking statistics; this requires V_sim ≥ 20π³/dk³.
 *
 * For dk = 4×10^7 m^-1 (= 4×10^5 cm^-1, Paper 3 value) the minimum box is
 * ≈ (270 nm)³.  The recommended setup for band-filling runs is V = (464 nm)³
 * with N_sim ≈ 10^5 particles (same as Paper 3).
 *
 * WARNING: do NOT use the default (100 nm)³ box of hotCarrierMHP.cpp with
 * USE_BAND_FILLING = true — set USE_BAND_FILLING which switches to the
 * larger box automatically.
 *
 * @tparam T  floating-point type (float or double)
 */
template <class T>
class emcPauliExclusion {
public:
  T dk;      // k-bin side length [m^-1]
  int nBins; // grid dimension per axis (covers -kMax … +kMax)
  int offset; // index offset so bin(k) = floor(k/dk) + offset ≥ 0
  int N_c;   // max simulated particles per k-slot

  std::vector<int> occ; // flattened 3D occupancy [ix*nBins² + iy*nBins + iz]

  SizeType nRejected  = 0; // Pauli-rejected scatter events since last reset
  SizeType nScattered = 0; // total scatter events since last reset

  emcPauliExclusion() = delete;

  /**
   * @param inDk      k-bin width [m^-1]  (Paper 3: 4×10^7 m^-1)
   * @param kMax      maximum |k| component to represent [m^-1];
   *                  set to sqrt(2·m*·E_max·q)/ħ for the heaviest carrier
   * @param Vsim      simulation volume [m^3]
   */
  emcPauliExclusion(T inDk, T kMax, T Vsim)
      : dk(inDk) {
    // Number of bins per axis: enough to cover ±kMax with one extra on each
    // side to absorb occasional out-of-range k-vectors.
    int halfBins = static_cast<int>(std::ceil(kMax / dk)) + 1;
    nBins  = 2 * halfBins + 1;
    offset = halfBins;

    T cap = T(2) * Vsim * dk * dk * dk
            / (T(8) * constants::pi * constants::pi * constants::pi);
    N_c = std::max(1, static_cast<int>(std::floor(cap)));

    occ.assign(static_cast<std::size_t>(nBins) * nBins * nBins, 0);

    std::cout << "PauliExclusion: dk = " << dk * 1e-7 << "e7 m^-1"
              << "  nBins = " << nBins << "^3"
              << "  N_c = " << N_c << " per slot\n";

    if (N_c < 5)
      std::cerr << "  WARNING: N_c = " << N_c
                << " < 5 — simulation box is too small for reliable "
                   "Pauli blocking.  Use V >= (270 nm)^3 with dk = 4e7 m^-1.\n";
  }

  // ── grid management ─────────────────────────────────────────────────────

  /** Rebuild occupancy from scratch.  Call once at the start of each timestep
   *  (after recombination and ESC extraction, before moveParticleType). */
  template <class PartVec>
  void buildGrid(const PartVec &particles) {
    std::fill(occ.begin(), occ.end(), 0);
    for (const auto &p : particles)
      occ[flatIdx(p.k)] += 1;
  }

  /** Return true if the k-slot for kNew is already at or above capacity. */
  bool isBlocked(const std::array<T, 3> &kNew) const {
    return occ[flatIdx(kNew)] >= N_c;
  }

  /** Move one particle from kOld slot to kNew slot.
   *  Call ONLY after a scatter has been accepted (isBlocked returned false). */
  void update(const std::array<T, 3> &kOld, const std::array<T, 3> &kNew) {
    occ[flatIdx(kOld)] = std::max(0, occ[flatIdx(kOld)] - 1);
    occ[flatIdx(kNew)] += 1;
  }

  /** Remove one particle from the grid (on recombination / ESC extraction).
   *  Called externally only if the Pauli grid is kept live between timesteps;
   *  otherwise use buildGrid() at the start of each step. */
  void remove(const std::array<T, 3> &k) {
    occ[flatIdx(k)] = std::max(0, occ[flatIdx(k)] - 1);
  }

  /** Add one particle (inverse of remove). */
  void add(const std::array<T, 3> &k) { occ[flatIdx(k)] += 1; }

  // ── diagnostics ─────────────────────────────────────────────────────────

  T getRejectionRate() const {
    return nScattered > 0 ? T(nRejected) / T(nScattered) : T(0);
  }

  void resetCounters() {
    nRejected  = 0;
    nScattered = 0;
  }

private:
  // Convert k-vector component to grid index, clamped to [0, nBins-1].
  int axisIdx(T kcomp) const {
    int i = static_cast<int>(std::floor(kcomp / dk)) + offset;
    return std::max(0, std::min(nBins - 1, i));
  }

  int flatIdx(const std::array<T, 3> &k) const {
    int ix = axisIdx(k[0]);
    int iy = axisIdx(k[1]);
    int iz = axisIdx(k[2]);
    return ix * nBins * nBins + iy * nBins + iz;
  }
};

#endif // EMC_PAULI_EXCLUSION_HPP
