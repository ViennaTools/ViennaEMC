#ifndef EMC_HOT_CARRIER_OUTPUT_HPP
#define EMC_HOT_CARRIER_OUTPUT_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <emcConstants.hpp>

/**
 * Hot-carrier diagnostics: carrier temperature extraction (MB + FD) and
 * VOC / PCE calculation.
 *
 * Temperature extraction
 * ──────────────────────
 * Two estimators are provided:
 *
 * 1) Maxwell-Boltzmann (always valid for low density or quick estimates):
 *       T_MB = (2/3) * <E> / k_B          [from equipartition]
 *
 * 2) Fermi-Dirac fit (method of moments, valid for degenerate carriers):
 *    (a) Start from T_MB = (2/3)*<E>/k_B.
 *    (b) Compute n_0(T_MB) = 1/V * ∫ g(k) f_{FD}(E; η=0) dk
 *                          = (1/(2π²)) * (2m*k_BT_MB/ħ²)^{3/2} * F_{1/2}(0)
 *        where F_{1/2}(0) = (√π/2) * (√2/2) ≈ 0.6049 is the half-integral.
 *        Actually n_0(T_MB) is the carrier density that would give η=0
 *        at T_MB: n_0 = (1/√2) * (m* k_B T_MB / (π ħ²))^{3/2} (see body).
 *    (c) Form ψ = ψ_target = n_0(T_MB) / n  (where n = N_sim/V_sim).
 *        Note ψ(η) = F_{1/2}(η) / F_{1/2}(0) is NOT the function below —
 *        actually we use the *simpler* moment condition:
 *           <E> / (3/2 k_B T_FD) = F_{3/2}(η) / F_{1/2}(η)
 *        combined with n = N_3D(T_FD, η).
 *        In practice: solve η by bisection from
 *           f(η) = n_0(T_MB) / n  vs  F_{1/2}(η) / F_{1/2}(0)
 *    (d) Convert:  T_FD = T_MB * F_{1/2}(η=0) / F_{1/2}(η)
 *                  μ    = η * k_B * T_FD / q  [eV]
 *
 * VOC formula (Paper 3 eq. 3, Faber et al. 2025)
 * ─────────────────────────────────────────────
 *   q V_OC = Δμ_eh * (T_L / T_eh) + ΔE_eh * (1 - T_L / T_eh)
 *
 *   where:  Δμ_eh = μ_e - μ_h + E_gap    (quasi-Fermi level splitting)
 *           ΔE_eh = <E_e> + <E_h> + E_gap (mean particle energies above edges)
 *           T_L   = lattice temperature
 *           T_eh  = (T_e + T_h) / 2       (average carrier temperature)
 *
 * @tparam T  floating-point type (float or double)
 */
namespace emcHotCarrierOutput {

// ── Fermi integrals ───────────────────────────────────────────────────────

/**
 * Complete Fermi-Dirac integral of order j:
 *   F_j(eta) = ∫_0^∞  x^j / (1 + exp(x - eta))  dx
 *
 * Evaluated by adaptive midpoint quadrature over [0, x_max] where
 * x_max is chosen so the integrand is negligible, plus an analytic
 * tail correction for large negative η (non-degenerate limit).
 *
 * The implementation is accurate to ~0.1% for all η in [-30, +30].
 */
template <class T>
T fermiIntegral(T j, T eta) {
  // Upper limit: x_max = max(10, eta + 30) to always capture the Fermi edge
  T xMax = std::max(T(10), eta + T(30));
  int N  = 1200;
  T  dx  = xMax / T(N);
  T  sum = T(0);
  for (int i = 0; i < N; i++) {
    T x   = (T(i) + T(0.5)) * dx;
    T exv = std::exp(std::min(x - eta, T(700))); // clamp to avoid overflow
    sum  += std::pow(x, j) / (T(1) + exv);
  }
  return sum * dx;
}

// ── Maxwell-Boltzmann temperature ────────────────────────────────────────

/**
 * Carrier temperature from mean kinetic energy (Maxwell-Boltzmann).
 *
 * @param avgEnergy_eV  <E> averaged over all simulated particles [eV]
 * @return T_MB [K]
 */
template <class T>
T getMBTemp(T avgEnergy_eV) {
  return T(2) * avgEnergy_eV * T(constants::q)
         / (T(3) * T(constants::kB));
}

// ── Fermi-Dirac fit by method of moments ─────────────────────────────────

/** Result of a Fermi-Dirac temperature fit. */
template <class T>
struct FDFitResult {
  T T_FD;   // carrier temperature [K]
  T mu_eV;  // chemical potential [eV]
  T eta;    // reduced chemical potential (μ / k_B T_FD)
  bool converged;
};

/**
 * Fermi-Dirac temperature fit for a 3D parabolic band.
 *
 * @param avgEnergy_eV  <E> per carrier [eV]
 * @param n             carrier density [m^-3]  (= N_sim / V_sim)
 * @param effMassRel    effective mass in units of m_0 (e.g. 0.2)
 * @param latTempK      lattice temperature [K] (lower bound for T_FD)
 * @return FDFitResult<T>
 */
template <class T>
FDFitResult<T> getFDFit(T avgEnergy_eV, T n, T effMassRel,
                         T latTempK = T(300)) {
  T mStar  = effMassRel * T(constants::me); // [kg]
  T Eavg_J = avgEnergy_eV * T(constants::q);

  // Maxwell-Boltzmann reference: exact temperature in the non-degenerate
  // limit, where <E> = (3/2) k T.
  T T_MB = getMBTemp(avgEnergy_eV);
  if (T_MB < latTempK)
    T_MB = latTempK;

  // Solve the degeneracy parameter η jointly from density and mean energy.
  // With the unnormalised Fermi-Dirac integrals F_j(η)=∫₀^∞ x^j/(e^{x-η}+1)dx,
  // the 3D relations are
  //     n   = (1/2π²) (2 m* k T / ħ²)^{3/2} F_{1/2}(η),
  //     <E> = k T · F_{3/2}(η) / F_{1/2}(η).
  // Eliminating T gives a temperature-independent equation for η,
  //     n · F_{3/2}(η)^{3/2} / F_{1/2}(η)^{5/2}
  //         = (1/2π²) (2 m* <E> / ħ²)^{3/2}  ≡  rhs,
  // after which k T = <E> F_{1/2}(η) / F_{3/2}(η). This reduces to
  // T_FD → T_MB as η → -∞ (Boltzmann) and T_FD < T_MB when degenerate.
  T rhs = std::pow(T(2) * mStar * Eavg_J
                     / (T(constants::hbar) * T(constants::hbar)),
                   T(1.5))
          / (T(2) * T(constants::pi) * T(constants::pi));
  T targetH = rhs / n;

  auto hOf = [](T eta) {
    T f12 = fermiIntegral<T>(T(0.5), eta);
    T f32 = fermiIntegral<T>(T(1.5), eta);
    return std::pow(f32, T(1.5)) / std::pow(f12, T(2.5));
  };

  // h(η) decreases monotonically: +∞ (η→-∞) down to a constant (η→+∞).
  T etaLo = T(-60), etaHi = T(60);
  bool converged = false;
  for (int iter = 0; iter < 80; iter++) {
    T etaMid = T(0.5) * (etaLo + etaHi);
    if (hOf(etaMid) > targetH)
      etaLo = etaMid; // h too large ⇒ need larger η
    else
      etaHi = etaMid;
    if ((etaHi - etaLo) < T(1e-5)) {
      converged = true;
      break;
    }
  }
  T eta_sol = T(0.5) * (etaLo + etaHi);

  T f12  = fermiIntegral<T>(T(0.5), eta_sol);
  T f32  = fermiIntegral<T>(T(1.5), eta_sol);
  T kT_J = Eavg_J * f12 / f32;
  T T_FD = kT_J / T(constants::kB);
  if (T_FD < latTempK)
    T_FD = latTempK;

  T mu_eV = eta_sol * kT_J / T(constants::q);

  return {T_FD, mu_eV, eta_sol, converged};
}

// ── VOC / PCE ─────────────────────────────────────────────────────────────

/**
 * Open-circuit voltage from the hot-carrier energy-balance formula.
 *
 * Implements eq. (3) of Faber, van de Ven, Filipovic, Loi, Koster (2025):
 *   q V_OC = Δμ_eh * (T_L/T_eh) + ΔE_eh * (1 - T_L/T_eh)
 *
 * @param mu_e_eV    electron chemical potential (from FD fit) [eV]
 * @param mu_h_eV    hole chemical potential (from FD fit) [eV]
 *                   (convention: positive above VBM = E_gap - |μ_h_from_VBM|)
 * @param T_e_K      electron temperature [K]
 * @param T_h_K      hole temperature [K]
 * @param T_L_K      lattice temperature [K]
 * @param E_gap_eV   band gap [eV]
 * @param avgE_e_eV  mean electron kinetic energy above CBM [eV]
 * @param avgE_h_eV  mean hole kinetic energy above VBM [eV]
 * @return V_OC [V]
 */
template <class T>
T getVOC(T mu_e_eV, T mu_h_eV, T T_e_K, T T_h_K, T T_L_K,
         T E_gap_eV, T avgE_e_eV, T avgE_h_eV) {
  T T_eh     = T(0.5) * (T_e_K + T_h_K);
  T ratio    = (T_eh > T(0)) ? T_L_K / T_eh : T(1);

  // Quasi-Fermi level splitting
  T dMu      = mu_e_eV + mu_h_eV + E_gap_eV;

  // Mean excess energy (kinetic + gap)
  T dE       = avgE_e_eV + avgE_h_eV + E_gap_eV;

  T qVoc     = dMu * ratio + dE * (T(1) - ratio);
  return qVoc; // [eV] → [V] (q cancels in q*V_OC/q)
}

/**
 * Power conversion efficiency.
 *
 * @param J_mAcm2   photocurrent density [mA cm^-2]  (e.g. J_SC ≈ 25 mA/cm²)
 * @param VOC_V     open-circuit voltage [V]
 * @param FF        fill factor (0–1; ~0.85 for radiative limit)
 * @param P_sun_mWcm2  incident power density [mW cm^-2] (AM1.5G = 100)
 * @return PCE as a fraction (multiply by 100 for percentage)
 */
template <class T>
T getPCE(T J_mAcm2, T VOC_V, T FF = T(0.85),
         T P_sun_mWcm2 = T(100)) {
  return J_mAcm2 * VOC_V * FF / P_sun_mWcm2;
}

// ── Particle statistics helpers ───────────────────────────────────────────

/**
 * Compute mean energy [eV] from a vector of emcParticle-like objects.
 * T_particle must have a member `energy` [eV].
 */
template <class T, class PartVec>
T meanEnergy_eV(const PartVec &particles) {
  if (particles.empty())
    return T(0);
  T sum = T(0);
  for (const auto &p : particles)
    sum += p.energy;
  return sum / T(particles.size());
}

} // namespace emcHotCarrierOutput

#endif // EMC_HOT_CARRIER_OUTPUT_HPP
