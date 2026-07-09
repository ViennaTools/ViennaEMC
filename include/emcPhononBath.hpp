#ifndef EMC_PHONON_BATH_HPP
#define EMC_PHONON_BATH_HPP

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include <emcConstants.hpp>
#include <emcUtil.hpp>

/**
 * Tracks the non-equilibrium LO phonon population N_q for each phonon
 * wavevector bin following the Lugli (1987) algorithm, as used in
 * Faber, Filipovic, Koster, J. Phys. Chem. Lett. 2024, 15, 12601.
 *
 * Evolution equation per timestep Δt (eq. 3-5 of Faber et al.):
 *   N_q^{t+Δt} = N_q^t + g(q)*Δt - (Δt/τ_LO)*(N_q^t - N_0)
 * where
 *   g(q)  = (1/D_ph(q)) * [n_em(q) - n_abs(q)] / Δt
 *   D_ph  = q^2 * Δq / (2π^2 * V_sim)            [phonon DOS]
 *
 * Usage in the simulation loop:
 *   1. Carrier-phonon scatter events call recordEmission(|q|) or
 *      recordAbsorption(|q|) with the phonon wavevector magnitude.
 *   2. After all particles are moved, call update(dt) to evolve N_q
 *      and reset the event counters.
 *   3. Call scatterHandler.reinitScatterTables() to rebuild the scatter
 *      rate look-up tables with the updated phonon occupation.
 *
 * @tparam T  floating-point type (float or double)
 *
 * @param nrBins   number of phonon wavevector bins
 * @param dq       bin width [1/m]
 * @param tauLO    LO phonon lifetime (Klemens / Ridley decay) [s]
 * @param N0       equilibrium Bose-Einstein occupation
 * @param Nq       current phonon occupation per bin (size nrBins)
 * @param nEm      phonon emission counts per bin since last update()
 * @param nAbs     phonon absorption counts per bin since last update()
 * @param Vsim     simulation volume [m^3]
 */
template <class T>
class emcPhononBath {
public:
  SizeType nrBins;
  T dq;    // [1/m]
  T tauLO; // [s]
  T N0;    // equilibrium Bose-Einstein occupation (same for all bins,
           // single-mode approximation)
  T Vsim;  // [m^3]

  std::vector<T> Nq;   // current occupation per bin
  std::vector<T> nEm;  // emission event count per bin (reset each update)
  std::vector<T> nAbs; // absorption event count per bin (reset each update)

  emcPhononBath() = delete;

  /**
   * @param inNrBins     number of q-bins (try 100-300)
   * @param inDq         bin width [1/m] (typical: 1e6 m^{-1})
   * @param inTauLO      LO phonon lifetime [s] (typical: 0.1-1 ps for MHPs)
   * @param phononEnergy LO phonon energy [eV]
   * @param latticeTemp  lattice temperature [K]
   * @param inVsim       simulation volume [m^3]
   */
  emcPhononBath(SizeType inNrBins, T inDq, T inTauLO, T phononEnergy,
                T latticeTemp, T inVsim)
      : nrBins(inNrBins), dq(inDq), tauLO(inTauLO), Vsim(inVsim),
        Nq(inNrBins), nEm(inNrBins, T(0)), nAbs(inNrBins, T(0)) {
    T exponent = constants::q * phononEnergy / (constants::kB * latticeTemp);
    N0 = T(1) / (std::exp(exponent) - T(1));
    std::fill(Nq.begin(), Nq.end(), N0);
  }

  /// Centre wavevector of bin i [1/m]
  T qCentre(SizeType i) const { return (T(i) + T(0.5)) * dq; }

  /// Bin index for phonon wavevector magnitude q (clamped)
  SizeType binOf(T q) const {
    SizeType idx = static_cast<SizeType>(q / dq);
    return (idx >= nrBins) ? nrBins - 1 : idx;
  }

  /// Current N_q for a given phonon wavevector magnitude q [1/m]
  T getNq(T q) const { return Nq[binOf(q)]; }

  /// Record a phonon emission event (carrier emits phonon with |q|)
  void recordEmission(T q) { nEm[binOf(q)] += T(1); }

  /// Record a phonon absorption event (carrier absorbs phonon with |q|)
  void recordAbsorption(T q) { nAbs[binOf(q)] += T(1); }

  /**
   * Advance phonon populations by one timestep Δt and reset event counters.
   * Implements eq. (3)–(5) of Faber, Filipovic, Koster, JPCL 2024.
   *
   * @param dt  simulation timestep [s]
   */
  void update(T dt) {
    for (SizeType i = 0; i < nrBins; i++) {
      T q = qCentre(i);
      // number of phonon modes in this q-bin = q² dq Vsim / (2π²)
      T Dph = q * q * dq * Vsim /
              (T(2) * constants::pi * constants::pi);
      T g = (Dph > T(0)) ? (nEm[i] - nAbs[i]) / (Dph * dt) : T(0);
      // evolution equation (eq. 3)
      Nq[i] += g * dt - (dt / tauLO) * (Nq[i] - N0);
      if (Nq[i] < T(0))
        Nq[i] = T(0); // occupation cannot go negative
      nEm[i] = T(0);
      nAbs[i] = T(0);
    }
  }

  /// Weighted-average phonon occupation (weight = phonon DOS ∝ q²)
  T getMeanNq() const {
    T sumW = T(0), sumWN = T(0);
    for (SizeType i = 0; i < nrBins; i++) {
      T q = qCentre(i);
      T w = q * q;
      sumW += w;
      sumWN += w * Nq[i];
    }
    return (sumW > T(0)) ? sumWN / sumW : N0;
  }

  T getN0() const { return N0; }
  T getTauLO() const { return tauLO; }
};

#endif // EMC_PHONON_BATH_HPP
