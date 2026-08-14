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
 * Base evolution equation per timestep Δt (eq. 3-5 of Faber et al.):
 *   N_q^{t+Δt} = N_q^t + g(q)*Δt - (Δt/τ_LO)*(N_q^t - N_0)
 * where
 *   g(q)  = (1/D_ph(q)) * [n_em(q) - n_abs(q)] / Δt
 *   D_ph  = q^2 * Δq / (2π^2 * V_sim)            [phonon DOS]
 *
 * Optional acoustic phonon bath (TO/LA/TA branches + Klemens cascade)
 * ─────────────────────────────────────────────────────────────────────
 * The Klemens channel LO → LA + TA deposits energy into acoustic modes.
 * If those modes heat up, the reverse process suppresses the effective LO
 * decay — the acoustic phonon bottleneck. This is modelled by replacing the
 * fixed equilibrium target N_0 with the LO occupation at the current
 * acoustic temperature T_ac:
 *
 *   N_q decay target:  N_LO*(T_ac) = 1 / (exp(ħω_LO / kB T_ac) - 1)
 *
 * A single effective acoustic mode (combining LA + TA at the zone-boundary
 * Klemens frequency ħω_ac ≈ ħω_LO / 2) is tracked:
 *
 *   dN_ac/dt = <(N_q - N_LO*(T_ac))>_DOS / τ_LO  - (N_ac - N_ac_eq) / τ_ac
 *
 * where <·>_DOS is the phonon-DOS-weighted average over q-bins, τ_ac is the
 * acoustic phonon–lattice coupling time, and N_ac_eq = 1/(exp(ħω_ac/kB T_lat)-1).
 * T_ac is back-computed from N_ac via the Planck inversion each step.
 *
 * When N_ac → N_ac_eq the target N_LO*(T_ac) → N_0 and the original
 * single-bath dynamics are fully recovered (backward-compatible).
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
 */
template <class T>
class emcPhononBath {
public:
  SizeType nrBins;
  T dq;    // [1/m]
  T tauLO; // [s]
  // Optional per-bin LO lifetime profile [s]. Empty means the uniform tauLO
  // above, and the update loop then follows the original code path exactly.
  // Set via setTauLOProfile(); used to test sensitivity to a wavevector-
  // dependent lifetime, for which no measurement exists in these materials.
  std::vector<T> tauLOProfile;
  T N0;    // equilibrium Bose-Einstein occupation of LO mode at T_lattice
  T Vsim;  // [m^3]

  std::vector<T> Nq;   // current occupation per bin
  std::vector<T> nEm;  // emission event count per bin (reset each update)
  std::vector<T> nAbs; // absorption event count per bin (reset each update)

  // Prefix sums for getNqInWindow(): cumW[i] = Σ_{j<i} 1/q_j,
  // cumWN[i] = Σ_{j<i} N_j/q_j. Both have nrBins+1 entries and are refreshed
  // by rebuildWindowSums() whenever Nq changes.
  std::vector<T> cumW;
  std::vector<T> cumWN;

  // Screening wavevector squared [1/m^2] used in the prefix-sum weight. Zero
  // reproduces the unscreened weight 1/q exactly, so this defaults to the
  // unscreened case and only changes behavior once a caller sets it.
  T qs2Screen = T(0);

private:
  T loPhononEnergyEV; // ħω_LO [eV]
  T latTempK;         // lattice temperature [K]

  // Acoustic bath (optional)
  bool hasAcousticBath = false;
  T acPhononEnergyEV = T(0); // ħω_ac [eV]  (≈ ħω_LO / 2 for Klemens channel)
  T tauAcoustic      = T(0); // acoustic–lattice coupling time [s]
  T N_ac             = T(0); // current acoustic phonon occupation
  T N_ac_eq          = T(0); // equilibrium acoustic occupation at T_lattice

  // Ridley channel (optional): LO → LA + TO, competing with Klemens LO → 2 LA.
  // A fraction wRidley of the anharmonic decay feeds a transverse-optical
  // reservoir instead of the acoustic one. Both reservoirs back-heat the LO
  // population; the LO decay target is their occupation-weighted blend, so
  // wRidley = 0 recovers the Klemens-only model exactly.
  bool hasRidley     = false;
  T wRidley          = T(0); // branching ratio into the TO reservoir [0,1]
  T toPhononEnergyEV = T(0); // ħω_TO [eV]
  T tauTO            = T(0); // TO–lattice coupling time [s]
  T N_TO             = T(0); // current TO phonon occupation
  T N_TO_eq          = T(0); // equilibrium TO occupation at T_lattice

  /// Refresh the prefix sums backing getNqInWindow(). Called after every change
  /// to Nq, so the windowed query stays O(1).
  void rebuildWindowSums() {
    cumW.assign(nrBins + 1, T(0));
    cumWN.assign(nrBins + 1, T(0));
    for (SizeType i = 0; i < nrBins; i++) {
      const T q = (T(i) + T(0.5)) * dq;
      // Transition weight for the (screened) Froehlich vertex:
      //   |M|^2 x phase space  ~  [1/(q^2 + q_s^2)] x q dq
      // which reduces to the unscreened 1/q when q_s = 0. The SAME weight sets
      // both the rate (via the window-averaged occupation) and the final-state
      // angle, so the two stay consistent by construction.
      const T denom = q * q + qs2Screen;
      const T w = (denom > T(0)) ? q / denom : T(0);
      cumW[i + 1] = cumW[i] + w;
      cumWN[i + 1] = cumWN[i] + w * Nq[i];
    }
  }

  // Planck occupation for given energy [eV] and temperature [K]
  T planck(T energyEV, T tempK) const {
    T x = constants::q * energyEV / (constants::kB * tempK);
    return T(1) / (std::exp(x) - T(1));
  }

  // Temperature [K] from Planck occupation (inversion of N = 1/(exp(x)-1))
  T planckTemp(T energyEV, T N) const {
    if (N <= T(0))
      return T(0);
    return constants::q * energyEV / (constants::kB * std::log(T(1) + T(1) / N));
  }

public:
  emcPhononBath() = delete;

  /**
   * @param inNrBins            number of q-bins (try 100-300)
   * @param inDq                bin width [1/m] (typical: 1e6 m^{-1})
   * @param inTauLO             LO phonon lifetime [s]
   * @param phononEnergy        LO phonon energy [eV]
   * @param latticeTemp         lattice temperature [K]
   * @param inVsim              simulation volume [m^3]
   * @param inEnableAcoustic    enable TO/LA/TA acoustic bath (default: false)
   * @param inAcPhononEnergy    acoustic phonon energy ħω_ac [eV]
   *                            (≈ phononEnergy / 2 for Klemens LO → LA+TA)
   * @param inTauAcoustic       acoustic–lattice thermalization time [s]
   *                            (typical: 10–100 ps for MHPs)
   * @param inWRidley           Ridley branching ratio [0,1]: fraction of the
   *                            anharmonic decay following LO → LA + TO instead
   *                            of the Klemens channel LO → 2 LA. 0 disables it.
   * @param inToPhononEnergy    ħω_TO [eV] of the Ridley daughter TO mode
   * @param inTauTO             TO–lattice thermalization time [s]
   */
  emcPhononBath(SizeType inNrBins, T inDq, T inTauLO, T phononEnergy,
                T latticeTemp, T inVsim,
                bool inEnableAcoustic  = false,
                T inAcPhononEnergy     = T(0),
                T inTauAcoustic        = T(0),
                T inWRidley            = T(0),
                T inToPhononEnergy     = T(0),
                T inTauTO              = T(0))
      : nrBins(inNrBins), dq(inDq), tauLO(inTauLO), Vsim(inVsim),
        Nq(inNrBins), nEm(inNrBins, T(0)), nAbs(inNrBins, T(0)),
        loPhononEnergyEV(phononEnergy), latTempK(latticeTemp) {
    N0 = planck(phononEnergy, latticeTemp);
    std::fill(Nq.begin(), Nq.end(), N0);

    if (inEnableAcoustic) {
      hasAcousticBath    = true;
      acPhononEnergyEV   = inAcPhononEnergy;
      tauAcoustic        = inTauAcoustic;
      N_ac_eq            = planck(inAcPhononEnergy, latticeTemp);
      N_ac               = N_ac_eq;
    }

    if (inEnableAcoustic && inWRidley > T(0) && inToPhononEnergy > T(0)) {
      hasRidley        = true;
      wRidley          = (inWRidley > T(1)) ? T(1) : inWRidley;
      toPhononEnergyEV = inToPhononEnergy;
      tauTO            = inTauTO;
      N_TO_eq          = planck(inToPhononEnergy, latticeTemp);
      N_TO             = N_TO_eq;
    }

    rebuildWindowSums();
  }

  /**
   * Set the screening wavevector entering the transition weight.
   *
   * Rebuilds the prefix sums only when the value actually changes, so calling
   * this every timestep is cheap. MUST be called from serial code (it mutates
   * the prefix sums that the parallel particle loop reads).
   *
   * @param inQs2 squared inverse screening length [1/m^2]; 0 = unscreened
   */
  void setScreeningQ2(T inQs2) {
    if (inQs2 == qs2Screen)
      return;
    qs2Screen = inQs2;
    rebuildWindowSums();
  }

  T getScreeningQ2() const { return qs2Screen; }

  /// Per-bin LO lifetime profile [s]; size must equal nrBins. An empty vector
  /// restores the uniform lifetime.
  void setTauLOProfile(std::vector<T> inProfile) {
    tauLOProfile = std::move(inProfile);
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

  /// Record a phonon emission event (carrier emits phonon with |q|).
  /// Called from inside the OpenMP-parallel particle-move loop, so the
  /// counter increment must be atomic: without it, concurrent threads lose
  /// read-modify-write updates on the shared bin, and the number of lost
  /// events depends on thread contention. That makes N_q -- and hence the
  /// whole hot-phonon bottleneck -- depend on machine load rather than on
  /// physics (observed: 6% drift in <N_LO> between an idle and a loaded host).
  void recordEmission(T q) {
    const SizeType i = binOf(q);
#ifdef _OPENMP
#pragma omp atomic
#endif
    nEm[i] += T(1);
  }

  /// Record a phonon absorption event (carrier absorbs phonon with |q|).
  /// Atomic for the same reason as recordEmission().
  void recordAbsorption(T q) {
    const SizeType i = binOf(q);
#ifdef _OPENMP
#pragma omp atomic
#endif
    nAbs[i] += T(1);
  }

  /**
   * Advance phonon populations by one timestep Δt and reset event counters.
   *
   * Without acoustic bath: implements eq. (3)-(5) of Faber et al. JPCL 2024.
   * With acoustic bath: replaces the fixed N_0 target with N_LO*(T_ac) and
   * evolves N_ac to track the energy cascade into acoustic modes.
   *
   * @param dt  simulation timestep [s]
   */
  void update(T dt) {
    // Determine LO decay target (= N0 unless a daughter reservoir is hot).
    // With the Ridley channel active the target is the branching-weighted
    // blend of the two daughter reservoirs' back-heating, so that wRidley = 0
    // reproduces the Klemens-only target exactly.
    T N_target = N0;
    if (hasAcousticBath) {
      T N_klemens = N0;
      if (N_ac > N_ac_eq) {
        T T_ac    = planckTemp(acPhononEnergyEV, N_ac);
        N_klemens = planck(loPhononEnergyEV, T_ac);
      }
      if (hasRidley) {
        T N_ridley = N0;
        if (N_TO > N_TO_eq) {
          T T_TO   = planckTemp(toPhononEnergyEV, N_TO);
          N_ridley = planck(loPhononEnergyEV, T_TO);
        }
        N_target = (T(1) - wRidley) * N_klemens + wRidley * N_ridley;
      } else {
        N_target = N_klemens;
      }
    }

    // Update each LO q-bin and accumulate DOS-weighted decay for acoustic.
    // With a lifetime profile the acoustic source becomes <(Nq-N*)/tau(q)>_DOS,
    // reducing to <Nq-N*>_DOS / tauLO for a uniform profile.
    const bool hasProfile = !tauLOProfile.empty();
    T sumW      = T(0);
    T sumWdN    = T(0); // Σ w * (Nq_old - N_target)            (uniform path)
    T sumWdNtau = T(0); // Σ w * (Nq_old - N_target) / tau(q)   (profile path)

    for (SizeType i = 0; i < nrBins; i++) {
      T q   = qCentre(i);
      T Dph = q * q * dq * Vsim / (T(2) * constants::pi * constants::pi);
      T g   = (Dph > T(0)) ? (nEm[i] - nAbs[i]) / (Dph * dt) : T(0);
      const T tau_i = hasProfile ? tauLOProfile[i] : tauLO;

      // Accumulate before updating (use values at start of step)
      if (hasAcousticBath) {
        T w   = q * q; // DOS weight ∝ q²
        sumW  += w;
        if (hasProfile)
          sumWdNtau += w * (Nq[i] - N_target) / tau_i;
        else
          sumWdN += w * (Nq[i] - N_target);
      }

      // LO evolution equation (eq. 3 of Faber et al., generalised target)
      Nq[i] += g * dt - (dt / tau_i) * (Nq[i] - N_target);
      if (Nq[i] < T(0))
        Nq[i] = T(0);

      nEm[i]  = T(0);
      nAbs[i] = T(0);
    }

    // Evolve daughter reservoirs. The decayed LO population is split between
    // the Klemens (acoustic) and Ridley (TO) channels by the branching ratio:
    //   dN_ac/dt = (1-w_R) <ΔN>_DOS / τ_LO - (N_ac - N_ac_eq) / τ_ac
    //   dN_TO/dt =    w_R  <ΔN>_DOS / τ_LO - (N_TO - N_TO_eq) / τ_TO
    if (hasAcousticBath) {
      T fKlemens = hasRidley ? (T(1) - wRidley) : T(1);
      if (!hasProfile) {
        // Uniform lifetime: original expressions, kept verbatim so the
        // default path is unchanged operation for operation.
        T dN_LO_mean = (sumW > T(0)) ? sumWdN / sumW : T(0);
        N_ac += dt * fKlemens * dN_LO_mean / tauLO
              - dt * (N_ac - N_ac_eq) / tauAcoustic;
        if (N_ac < N_ac_eq)
          N_ac = N_ac_eq; // reservoir cannot cool below the lattice
        if (hasRidley) {
          N_TO += dt * wRidley * dN_LO_mean / tauLO
                - dt * (N_TO - N_TO_eq) / tauTO;
          if (N_TO < N_TO_eq)
            N_TO = N_TO_eq;
        }
      } else {
        // Lifetime profile: the reservoir source is <(Nq - N*)/tau(q)>_DOS.
        T decayFlux = (sumW > T(0)) ? sumWdNtau / sumW : T(0);
        N_ac += dt * fKlemens * decayFlux
              - dt * (N_ac - N_ac_eq) / tauAcoustic;
        if (N_ac < N_ac_eq)
          N_ac = N_ac_eq;
        if (hasRidley) {
          N_TO += dt * wRidley * decayFlux
                - dt * (N_TO - N_TO_eq) / tauTO;
          if (N_TO < N_TO_eq)
            N_TO = N_TO_eq;
        }
      }
    }

    rebuildWindowSums();
  }

  /**
   * Coupling-weighted mean occupation over the wavevector window a transition
   * can actually reach: [qMin, qMax] = [|k_i − k_f|, k_i + k_f].
   *
   * For the Fröhlich interaction the probability of a transition at wavevector
   * q goes as |M(q)|² × (phase space) ∝ (1/q²) × q dq = dq/q, so the population
   * a carrier samples is the 1/q-weighted mean over the allowed window — NOT
   * the DOS-weighted (q²) mean returned by getMeanNq(). The distinction is
   * large whenever the non-equilibrium excess is concentrated at small q, which
   * is the generic polar case: emission piles up near q_min, exactly where the
   * 1/q² coupling is strongest and the q² phonon DOS is smallest.
   *
   * O(1) via prefix sums maintained by update(), so this is cheap enough to
   * call once per energy level on every scatter-table rebuild.
   *
   * Falls back to getMeanNq() when the window is empty or degenerate (e.g. a
   * single-bin bath), so lumped-bath users are unaffected.
   */
  T getNqInWindow(T qMin, T qMax) const {
    if (nrBins < 2 || qMax <= qMin)
      return getMeanNq();
    // Bin range covering the window (inclusive), clamped to the grid.
    long lo = static_cast<long>(std::floor(qMin / dq));
    long hi = static_cast<long>(std::ceil(qMax / dq));
    if (lo < 0)
      lo = 0;
    if (hi > static_cast<long>(nrBins))
      hi = static_cast<long>(nrBins);
    if (hi - lo < 1)
      return getMeanNq();
    const T wSum = cumW[hi] - cumW[lo];
    if (wSum <= T(0))
      return getMeanNq();
    return (cumWN[hi] - cumWN[lo]) / wSum;
  }

  /**
   * Sample a phonon wavevector from the q-resolved transition distribution.
   *
   * When the occupation varies across the allowed window, the FINAL-STATE ANGLE
   * must be sampled from the same distribution that sets the rate, or the two
   * are inconsistent: the rate would count small-q transitions as more likely
   * while the angular sampler still treated the window as uniformly occupied.
   *
   * For the unscreened Fröhlich vertex, with q² = k_i² + k_f² − 2k_i k_f cosθ,
   *
   *     P(cosθ) dcosθ ∝ [f(q)/(q²+q_s²)] dcosθ  →  P(q) dq ∝ f(q) q dq/(q²+q_s²)
   *
   * where f(q) = N_q (absorption) or N_q + 1 (emission). Inverting that CDF
   * needs only the prefix sums already maintained for getNqInWindow():
   * Σ N_j/q_j for the absorption part and Σ 1/q_j for the spontaneous part.
   * Cost is O(log nrBins) per event.
   *
   * The caller converts back with cosθ = (k_i² + k_f² − q²) / (2 k_i k_f).
   *
   * Screening is handled through setScreeningQ2(): the weight becomes
   * q/(q² + q_s²), which reduces to 1/q when q_s = 0. Both the rate and this
   * sampler read the same prefix sums, so they cannot drift out of step.
   *
   * @param qMin,qMax allowed window [1/m]
   * @param emission  true for phonon emission (f = N_q + 1), false absorption
   * @param r         uniform random number in [0,1)
   */
  T sampleQ(T qMin, T qMax, bool emission, T r) const {
    if (nrBins < 2 || qMax <= qMin)
      return qMin;
    long lo = static_cast<long>(std::floor(qMin / dq));
    long hi = static_cast<long>(std::ceil(qMax / dq));
    if (lo < 0)
      lo = 0;
    if (hi > static_cast<long>(nrBins))
      hi = static_cast<long>(nrBins);
    if (hi - lo < 1)
      return qMin;

    auto S = [&](long i) {
      return emission ? (cumWN[i] + cumW[i]) : cumWN[i];
    };
    const T sLo = S(lo), sHi = S(hi);
    const T span = sHi - sLo;
    if (!(span > T(0)))
      return T(0.5) * (qMin + qMax); // degenerate: fall back to window centre

    const T target = sLo + r * span;
    // binary search for the bin whose cumulative sum brackets the target
    long a = lo, b = hi;
    while (b - a > 1) {
      const long mid = (a + b) / 2;
      if (S(mid) <= target)
        a = mid;
      else
        b = mid;
    }
    // linear interpolation inside bin a
    const T sA = S(a), sB = S(a + 1);
    const T frac = (sB > sA) ? (target - sA) / (sB - sA) : T(0.5);
    const T q = (T(a) + frac) * dq;
    return std::max(qMin, std::min(qMax, q));
  }

  /// Weighted-average LO phonon occupation (weight = phonon DOS ∝ q²)
  T getMeanNq() const {
    T sumW = T(0), sumWN = T(0);
    for (SizeType i = 0; i < nrBins; i++) {
      T q = qCentre(i);
      T w = q * q;
      sumW  += w;
      sumWN += w * Nq[i];
    }
    return (sumW > T(0)) ? sumWN / sumW : N0;
  }

  /// Current acoustic phonon occupation (0 if acoustic bath disabled)
  T getMeanNac() const { return N_ac; }

  /// Equilibrium acoustic occupation at lattice temperature
  T getAcousticN0() const { return N_ac_eq; }

  /// Acoustic phonon temperature [K] back-computed from N_ac; returns
  /// lattice temperature when acoustic bath is disabled or at equilibrium.
  T getAcousticTemp() const {
    if (!hasAcousticBath || N_ac <= N_ac_eq)
      return latTempK;
    return planckTemp(acPhononEnergyEV, N_ac);
  }

  /// Current TO (Ridley daughter) phonon occupation
  T getMeanNTO() const { return N_TO; }

  /// TO phonon temperature [K]; lattice temperature when Ridley is disabled.
  T getTOTemp() const {
    if (!hasRidley || N_TO <= N_TO_eq)
      return latTempK;
    return planckTemp(toPhononEnergyEV, N_TO);
  }

  bool acousticBathEnabled() const { return hasAcousticBath; }
  bool ridleyEnabled() const { return hasRidley; }
  T getRidleyBranching() const { return wRidley; }
  T getN0()    const { return N0; }
  T getTauLO() const { return tauLO; }
  T getLOEnergy() const { return loPhononEnergyEV; }
};

#endif // EMC_PHONON_BATH_HPP
