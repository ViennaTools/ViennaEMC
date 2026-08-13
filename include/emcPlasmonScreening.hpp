#ifndef EMC_PLASMON_SCREENING_HPP
#define EMC_PLASMON_SCREENING_HPP

#include <cmath>

#include <emcConstants.hpp>
#include <emcUtil.hpp>

/**
 * Free-carrier (plasmon) screening state shared between scatter mechanisms.
 *
 * Holds the inverse screening length q_s, which the simulation loop refreshes
 * from the current carrier density and carrier temperature. Polar-optical
 * mechanisms that hold a pointer to this object then screen their Coulomb-like
 * 1/q² vertex as
 *
 *     1/q²  →  1/(q² + q_s²)
 *
 * which turns the Fröhlich angular integral from ln(q₊/q₋) into
 * ½ ln((q₊² + q_s²)/(q₋² + q_s²)) — see emcScreenedFroehlichInteraction.hpp.
 *
 * Debye (non-degenerate) form:
 *
 *     q_s² = n e² / (ε_s ε₀ k_B T_e)
 *
 * evaluated at the *carrier* temperature T_e, not the lattice temperature,
 * because under field the carriers are hot and screen more weakly.
 *
 * Two caveats that matter when interpreting results:
 *
 *  - STATIC screening over-screens the polar interaction. The plasma cannot
 *    respond on a timescale shorter than 1/ω_p, so at ω ≈ ω_LO the true
 *    (dynamic) screening is weaker than the static limit. A proper treatment
 *    diagonalises the coupled LO-plasmon modes. Using the static limit
 *    therefore yields a CONSERVATIVE bound: any effect that survives it is
 *    robust against screening.
 *  - The Debye form assumes a non-degenerate gas. Above ~4e18 cm^-3 in a
 *    m* = 0.28 band at 300 K the gas is degenerate and the Thomas-Fermi form
 *    applies instead; Debye then overestimates screening at low carrier
 *    temperature. Under high field T_e is large, which pushes the gas back
 *    toward non-degenerate, but the high-density low-field corner should be
 *    read with this in mind.
 *
 * @tparam T floating-point type
 */
template <class T> class emcPlasmonScreening {
private:
  T qs2 = T(0);   // squared inverse screening length [1/m²]
  T epsStatic;    // static relative dielectric constant
  bool enabled;

public:
  emcPlasmonScreening() = delete;

  /**
   * @param inEpsStatic static relative dielectric constant
   * @param inEnabled   when false, getQs2() stays 0 and every screened
   *                    mechanism reduces exactly to its unscreened form
   */
  emcPlasmonScreening(T inEpsStatic, bool inEnabled = true)
      : epsStatic(inEpsStatic), enabled(inEnabled) {}

  /**
   * Refresh the screening wavevector.
   *
   * @param density     carrier density [1 / m³]
   * @param carrierTemp carrier temperature [K]
   */
  void update(T density, T carrierTemp) {
    if (!enabled || density <= T(0) || carrierTemp <= T(0)) {
      qs2 = T(0);
      return;
    }
    qs2 = density * constants::q * constants::q /
          (epsStatic * constants::eps0 * constants::kB * carrierTemp);
  }

  /// Set the screening wavevector directly [1/m²] (for tests / overrides).
  void setQs2(T inQs2) { qs2 = enabled ? inQs2 : T(0); }

  /// Squared inverse screening length [1/m²]; 0 disables screening.
  T getQs2() const { return qs2; }

  /// Inverse screening length [1/m]
  T getQs() const { return std::sqrt(qs2); }

  /// Screening length [m]; infinity when disabled.
  T getScreeningLength() const {
    return (qs2 > T(0)) ? T(1) / std::sqrt(qs2)
                        : std::numeric_limits<T>::infinity();
  }

  bool isEnabled() const { return enabled; }
};

#endif // EMC_PLASMON_SCREENING_HPP
