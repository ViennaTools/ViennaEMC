#ifndef EMC_INTER_CARRIER_SCATTER_HPP
#define EMC_INTER_CARRIER_SCATTER_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>

#include <emcConstants.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

/**
 * @brief Inter-species (electron-hole) binary carrier-carrier scattering.
 *
 * Uses the same Lugli-Ferry random-pairing algorithm as emcCarrierCarrierScatter
 * but handles two different effective masses, so the CM frame is asymmetric.
 *
 * The key quantities that generalise the same-mass formulas are:
 *
 *   fe = m_e / (m_e + m_h)          hole-to-total mass fraction
 *   fh = m_h / (m_e + m_h)
 *   m_r = m_e * m_h / (m_e + m_h)  reduced mass
 *   m_eff = 2 * m_r                 effective mass for the probability formula
 *
 * Effective relative wave vector (reduces to k_i - k_j for equal masses):
 *
 *   k_rel_eff = 2 * k_rel_CM = 2 * (fh * k_e - fe * k_h)
 *
 * With these definitions the collision probability formula is identical in
 * structure to the same-mass case:
 *
 *   P = C * n_eff * |k_rel_eff| * dt / (q_D^2 * (q_D^2 + 4|k_rel_eff|^2))
 *
 * where C = q^4 * m_eff / (pi * eps^2 * hbar^3) and
 *   n_eff = max(n, p) [m^-3]   (so that min(N_e,N_h) pairs reproduce the
 *                                physical N_e * N_h * sigma * v_rel rate)
 *
 * Debye screening uses the total carrier density n + p so both species
 * contribute to screening.
 *
 * Wave-vector update (exactly conserves crystal momentum and kinetic energy):
 *
 *   k_e' = fe * (k_e + k_h) + k_rel_eff_new / 2
 *   k_h' = fh * (k_e + k_h) - k_rel_eff_new / 2
 *
 * For equal masses (fe = fh = 0.5) these reduce to the familiar symmetric
 * k_CM +/- k_rel/2 update used in emcCarrierCarrierScatter.
 *
 * @tparam T  floating-point precision (float or double)
 */
template <class T>
class emcInterCarrierScatter {
private:
  T effMassE;  // electron effective mass [kg]
  T effMassH;  // hole effective mass [kg]
  T fe;        // m_e / (m_e + m_h)
  T fh;        // m_h / (m_e + m_h)
  T Vsim;      // simulation volume [m^3]
  T latTempK;  // lattice temperature [K]
  T probConst; // q^4 * 2*m_r / (pi * eps^2 * hbar^3)
  T qD2perNT;  // q^2 / (eps0 * eps_r * kB)

  mutable std::uniform_real_distribution<T> dist{T(0), T(1)};

  T kineticEnergyE(const std::array<T, 3> &k) const {
    T k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    return constants::hbar * constants::hbar * k2 /
           (T(2) * effMassE * constants::q);
  }
  T kineticEnergyH(const std::array<T, 3> &k) const {
    T k2 = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    return constants::hbar * constants::hbar * k2 /
           (T(2) * effMassH * constants::q);
  }

public:
  /**
   * @param inEpsR       relative permittivity for Debye screening
   *                     (use optical / high-frequency eps_inf for hot carriers)
   * @param relEffMassE  electron relative effective mass (m_star / m_e)
   * @param relEffMassH  hole relative effective mass (m_star / m_e)
   * @param inVsim       simulation box volume [m^3]
   * @param inTempK      lattice temperature [K]
   */
  emcInterCarrierScatter(T inEpsR, T relEffMassE, T relEffMassH,
                         T inVsim, T inTempK)
      : effMassE(relEffMassE * constants::me),
        effMassH(relEffMassH * constants::me),
        fe(relEffMassE / (relEffMassE + relEffMassH)),
        fh(relEffMassH / (relEffMassE + relEffMassH)),
        Vsim(inVsim), latTempK(inTempK) {
    T eps  = constants::eps0 * inEpsR;
    T mRed = relEffMassE * relEffMassH /
             (relEffMassE + relEffMassH) * constants::me; // reduced mass [kg]
    // m_eff = 2 * m_r so the formula matches the intra-species convention
    T mEff = T(2) * mRed;
    probConst = constants::q * constants::q * constants::q * constants::q *
                mEff /
                (constants::pi * eps * eps *
                 constants::hbar * constants::hbar * constants::hbar);
    qD2perNT = constants::q * constants::q / (eps * constants::kB);
  }

  /**
   * @brief Perform one round of electron-hole binary collisions.
   *
   * Randomly pairs min(N_e, N_h) electron-hole pairs. For each pair the
   * collision probability is evaluated; accepted collisions rotate the
   * effective relative wave vector while exactly conserving the total crystal
   * momentum and total kinetic energy.
   *
   * @param electrons  electron particle vector (modified in place)
   * @param holes      hole particle vector (modified in place)
   * @param dt         time step [s]
   * @param rng        random number generator
   */
  void scatter(std::vector<emcParticle<T>> &electrons,
               std::vector<emcParticle<T>> &holes,
               T dt, emcRNG &rng) const {
    SizeType Ne = electrons.size();
    SizeType Nh = holes.size();
    if (Ne == 0 || Nh == 0)
      return;

    T n = T(Ne) / Vsim;
    T p = T(Nh) / Vsim;

    // Both species contribute to Debye screening
    T qD2 = (n + p) * qD2perNT / latTempK;

    // n_eff = max(n,p) so min(N_e,N_h) paired events reproduce N_e*N_h*sigma*v
    T nEff = std::max(n, p);

    // Random pairing: shuffle both index arrays independently
    std::vector<SizeType> idxE(Ne), idxH(Nh);
    std::iota(idxE.begin(), idxE.end(), SizeType(0));
    std::iota(idxH.begin(), idxH.end(), SizeType(0));
    std::shuffle(idxE.begin(), idxE.end(), rng);
    std::shuffle(idxH.begin(), idxH.end(), rng);

    SizeType nPairs = std::min(Ne, Nh);
    for (SizeType pair = 0; pair < nPairs; ++pair) {
      auto &pe = electrons[idxE[pair]];
      auto &ph = holes[idxH[pair]];

      // k_rel_eff = 2 * k_rel_CM = 2 * (fh*k_e - fe*k_h)
      // For equal masses this reduces to k_e - k_h (the intra-species k_rel).
      std::array<T, 3> kRelEff;
      for (int d = 0; d < 3; ++d)
        kRelEff[d] = T(2) * (fh * pe.k[d] - fe * ph.k[d]);

      T kRelEffSq   = kRelEff[0]*kRelEff[0] + kRelEff[1]*kRelEff[1]
                    + kRelEff[2]*kRelEff[2];
      T kRelEffNorm = std::sqrt(kRelEffSq);

      if (kRelEffNorm < T(1e4))
        continue;

      T denom = qD2 * (qD2 + T(4) * kRelEffSq);
      T P     = probConst * nEff * kRelEffNorm * dt / denom;

      if (dist(rng) >= std::min(P, T(1)))
        continue;

      // Screened Rutherford deflection angle
      T beta     = qD2 / (T(4) * kRelEffSq);
      T r        = dist(rng);
      T cosTheta = (T(2) * r * (T(1) + beta) - beta) / (T(2) * r + beta);
      cosTheta   = std::clamp(cosTheta, T(-1), T(1));

      std::array<T, 3> kRelEffNew =
          initRandomDirectionWithRespectToCurrentK(kRelEff, cosTheta, dist(rng));

      // Update wave vectors (conserves crystal momentum and kinetic energy):
      //   k_e' = fe*(k_e+k_h) + k_rel_eff_new/2
      //   k_h' = fh*(k_e+k_h) - k_rel_eff_new/2
      for (int d = 0; d < 3; ++d) {
        T kTot    = pe.k[d] + ph.k[d];
        pe.k[d] = fe * kTot + T(0.5) * kRelEffNew[d];
        ph.k[d] = fh * kTot - T(0.5) * kRelEffNew[d];
      }

      pe.energy = kineticEnergyE(pe.k);
      ph.energy = kineticEnergyH(ph.k);
    }
  }
};

#endif // EMC_INTER_CARRIER_SCATTER_HPP
