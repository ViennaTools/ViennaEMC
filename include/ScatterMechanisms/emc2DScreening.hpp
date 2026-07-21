#ifndef EMC_2D_SCREENING_HPP
#define EMC_2D_SCREENING_HPP

#include <cmath>

#include <emcConstants.hpp>

/*! \file emc2DScreening.hpp
 * \brief Free-carrier (static) screening for 2D semiconductors.
 *
 * Long-wavelength 2D Thomas-Fermi / static-RPA screening of the long-range
 * (Froehlich polar-optical, charged-impurity) interactions by the mobile
 * carrier gas. The screened matrix element is
 *
 *     |g_scr(q)|^2 = |g_bare(q)|^2 / eps(q)^2 ,   eps(q) = 1 + q_s / q ,
 *
 * with the 2D screening wavevector
 *
 *     q_s = e^2 (dn/dmu) / (2 eps0 eps_r) ,
 *     dn/dmu = D0 * (1 - exp(-n / (D0 kB T))) ,   D0 = g m* / (2 pi hbar^2) .
 *
 * dn/dmu (the thermally-averaged density of states) interpolates between the
 * degenerate limit (dn/dmu -> D0, giving the standard 2D Thomas-Fermi
 * wavevector q_s = g/a_B*) and the non-degenerate/Debye limit
 * (dn/dmu -> n/kB T), so q_s is valid across the relevant density/temperature
 * range. This is the mechanism behind the experimentally observed increase of
 * MoS2 mobility with carrier density (screening of the polar-optical and
 * charged-impurity scattering).
 *
 * References: Stern, PRL 18, 546 (1967); Kaasbjerg, Thygesen & Jacobsen,
 * PRB 85, 115317 (2012); Ma & Jena, PRX 4, 011043 (2014).
 *
 * NOTE: eps(q) = 1 + q_s/q uses the long-wavelength form for all q; the
 * q > 2 k_F reduction of the static Lindhard function is neglected (it slightly
 * over-screens large-angle / large-q scattering).
 */

/*! \brief 2D static (Thomas-Fermi / long-wavelength RPA) screening wavevector.
 *
 * @param carrierDensity mobile 2D carrier sheet density n [1 / m^2]
 * @param temperature    carrier temperature [K]
 * @param envPermittivity effective relative permittivity of the dielectric
 *        environment screening the carriers (eps_r)
 * @param dosEffMass     density-of-states effective mass m* [kg]
 * @param degeneracy     total band degeneracy g = g_spin * g_valley
 *        (e.g. 4 for the MoS2 K/K' conduction valleys with spin)
 * @return screening wavevector q_s [1 / m] (0 if there are no carriers)
 */
template <class T>
T twoDStaticScreeningWavevector(T carrierDensity, T temperature,
                                T envPermittivity, T dosEffMass,
                                T degeneracy = T(4)) {
  if (carrierDensity <= T(0) || temperature <= T(0) || dosEffMass <= T(0))
    return T(0);
  T D0 = degeneracy * dosEffMass /
         (2 * constants::pi * constants::hbar * constants::hbar); // [1/(J m^2)]
  T kBT = constants::kB * temperature;                            // [J]
  T dndmu = D0 * (1 - std::exp(-carrierDensity / (D0 * kBT)));    // [1/(J m^2)]
  return constants::q * constants::q * dndmu /
         (2 * constants::eps0 * envPermittivity); // [1/m]
}

/*! \brief 2D static dielectric function eps(q) = 1 + q_s / q (>= 1).
 *  Returns 1 (unscreened) for a non-positive momentum transfer guard. */
template <class T> T twoDStaticDielectric(T q, T screeningWavevector) {
  if (q <= T(0) || screeningWavevector <= T(0))
    return T(1);
  return T(1) + screeningWavevector / q;
}

/*! \brief Screening reduction factor 1 / eps(q)^2 applied to a bare |g(q)|^2. */
template <class T> T twoDScreeningFactor(T q, T screeningWavevector) {
  T eps = twoDStaticDielectric(q, screeningWavevector);
  return T(1) / (eps * eps);
}

#endif // EMC_2D_SCREENING_HPP
