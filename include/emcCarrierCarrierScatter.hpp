#ifndef EMC_CARRIER_CARRIER_SCATTER_HPP
#define EMC_CARRIER_CARRIER_SCATTER_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

#include <emcConstants.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

/**
 * @brief Binary carrier-carrier scattering via the random-pairing method.
 *
 * Models stochastic two-body Coulomb collisions between mobile carriers using
 * Debye-screened differential cross-sections (Lugli & Ferry, Solid-State
 * Electron. 28, 1223, 1985). Each call to scatter() randomly pairs the
 * carrier ensemble, evaluates a collision probability per pair from the
 * screened Coulomb cross-section, and — when a collision fires — rotates the
 * relative wave vector in the centre-of-mass (CM) frame while exactly
 * conserving the total momentum and total kinetic energy of the pair.
 *
 * Collision probability per pair per time step:
 *
 *   P = C * n * |k_rel| * dt / [ q_D^2 * (q_D^2 + 4|k_rel|^2) ]
 *
 * where
 *   C    = q^4 * m_eff / (pi * eps^2 * hbar^3)   [SI, s^-1]
 *   q_D^2 = n * q^2 / (eps0 * eps_r * k_B * T)   [m^-2]
 *   k_rel = k_i - k_j                             [m^-1]
 *
 * Deflection angle in the CM frame is sampled from the screened Rutherford
 * distribution:
 *
 *   cos(theta) = [ 2r(1+beta) - beta ] / (2r + beta),  beta = q_D^2 / (4|k_rel|^2)
 *
 * Energy is recalculated from the updated wave vectors using the parabolic
 * dispersion E = hbar^2 |k|^2 / (2 m_eff q), where q converts J to eV.
 * This is exact for parabolic bands and introduces only an O(alpha) correction
 * for weakly non-parabolic bands.
 *
 * @tparam T  floating-point precision (float or double)
 */
template <class T>
class emcCarrierCarrierScatter {
private:
  T effMass;    // effective mass [kg]
  T Vsim;       // simulation volume [m^3]
  T latTempK;   // lattice temperature [K] (used for q_D)
  T probConst;  // q^4 * m_eff / (pi * eps^2 * hbar^3)  [s^-1]
  T qD2perNT;   // q^2 / (eps0 * eps_r * k_B)  so that q_D^2 = n * qD2perNT / T

  mutable std::uniform_real_distribution<T> dist{T(0), T(1)};

  // Parabolic kinetic energy in eV from wave vector [m^-1]
  T kineticEnergy(const std::array<T, 3> &k) const {
    T k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
    return constants::hbar * constants::hbar * k2 /
           (T(2) * effMass * constants::q);
  }

public:
  /**
   * @param inEpsR      relative permittivity for Debye screening (use
   *                    optical / high-frequency value eps_inf for hot carriers)
   * @param relEffMass  relative effective mass (m_star / m_e)
   * @param inVsim      simulation box volume [m^3]
   * @param inTempK     lattice temperature [K]
   */
  emcCarrierCarrierScatter(T inEpsR, T relEffMass, T inVsim, T inTempK)
      : effMass(relEffMass * constants::me), Vsim(inVsim), latTempK(inTempK) {
    T eps = constants::eps0 * inEpsR;
    probConst = constants::q * constants::q * constants::q * constants::q *
                effMass /
                (constants::pi * eps * eps *
                 constants::hbar * constants::hbar * constants::hbar);
    qD2perNT = constants::q * constants::q / (eps * constants::kB);
  }

  /**
   * @brief Perform one round of carrier-carrier binary collisions.
   *
   * @param particles  particle vector for one carrier type (modified in place)
   * @param dt         time step [s]
   * @param rng        random number generator
   */
  void scatter(std::vector<emcParticle<T>> &particles, T dt,
               emcRNG &rng) const {
    SizeType N = particles.size();
    if (N < 2)
      return;

    T n   = T(N) / Vsim;
    T qD2 = n * qD2perNT / latTempK; // Debye wave vector squared [m^-2]

    // Randomly pair all carriers
    std::vector<SizeType> idx(N);
    std::iota(idx.begin(), idx.end(), SizeType(0));
    std::shuffle(idx.begin(), idx.end(), rng);

    SizeType nPairs = N / 2;
    for (SizeType p = 0; p < nPairs; ++p) {
      auto &pi = particles[idx[2 * p]];
      auto &pj = particles[idx[2 * p + 1]];

      // Centre-of-mass and relative wave vectors
      std::array<T, 3> kCM, kRel;
      for (int d = 0; d < 3; ++d) {
        kCM[d]  = T(0.5) * (pi.k[d] + pj.k[d]);
        kRel[d] = pi.k[d] - pj.k[d];
      }
      T kRelSq   = kRel[0]*kRel[0] + kRel[1]*kRel[1] + kRel[2]*kRel[2];
      T kRelNorm = std::sqrt(kRelSq);

      // Skip pairs with nearly identical momenta (avoids 0/0 in angle formula)
      if (kRelNorm < T(1e4))
        continue;

      // Collision probability P = probConst * n * |k_rel| * dt
      //                           / (q_D^2 * (q_D^2 + 4|k_rel|^2))
      T denom = qD2 * (qD2 + T(4) * kRelSq);
      T P = probConst * n * kRelNorm * dt / denom;

      if (dist(rng) >= std::min(P, T(1)))
        continue;

      // Sample polar deflection angle from screened Rutherford distribution:
      //   cos(theta) = (2r(1+beta) - beta) / (2r + beta)
      //   beta = q_D^2 / (4|k_rel|^2)
      T beta = qD2 / (T(4) * kRelSq);
      T r    = dist(rng);
      T cosTheta = (T(2) * r * (T(1) + beta) - beta) / (T(2) * r + beta);
      cosTheta = std::clamp(cosTheta, T(-1), T(1));

      // Rotate k_rel by (theta, phi_random) — preserves |k_rel|
      std::array<T, 3> kRelNew =
          initRandomDirectionWithRespectToCurrentK(kRel, cosTheta, dist(rng));

      // Update wave vectors: k_i' = k_CM + k_rel'/2, k_j' = k_CM - k_rel'/2
      for (int d = 0; d < 3; ++d) {
        pi.k[d] = kCM[d] + T(0.5) * kRelNew[d];
        pj.k[d] = kCM[d] - T(0.5) * kRelNew[d];
      }

      // Recompute kinetic energies (parabolic dispersion, result in eV)
      pi.energy = kineticEnergy(pi.k);
      pj.energy = kineticEnergy(pj.k);
    }
  }
};

#endif // EMC_CARRIER_CARRIER_SCATTER_HPP
