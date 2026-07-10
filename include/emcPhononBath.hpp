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
  T N0;    // equilibrium Bose-Einstein occupation of LO mode at T_lattice
  T Vsim;  // [m^3]

  std::vector<T> Nq;   // current occupation per bin
  std::vector<T> nEm;  // emission event count per bin (reset each update)
  std::vector<T> nAbs; // absorption event count per bin (reset each update)

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

    // Update each LO q-bin and accumulate DOS-weighted decay for acoustic
    T sumW   = T(0);
    T sumWdN = T(0); // Σ w * (Nq_old - N_target)

    for (SizeType i = 0; i < nrBins; i++) {
      T q   = qCentre(i);
      T Dph = q * q * dq * Vsim / (T(2) * constants::pi * constants::pi);
      T g   = (Dph > T(0)) ? (nEm[i] - nAbs[i]) / (Dph * dt) : T(0);

      // Accumulate before updating (use values at start of step)
      if (hasAcousticBath) {
        T w   = q * q; // DOS weight ∝ q²
        sumW  += w;
        sumWdN += w * (Nq[i] - N_target);
      }

      // LO evolution equation (eq. 3 of Faber et al., generalised target)
      Nq[i] += g * dt - (dt / tauLO) * (Nq[i] - N_target);
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
      T dN_LO_mean = (sumW > T(0)) ? sumWdN / sumW : T(0);
      T fKlemens   = hasRidley ? (T(1) - wRidley) : T(1);

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
    }
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
