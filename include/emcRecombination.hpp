#ifndef EMC_RECOMBINATION_HPP
#define EMC_RECOMBINATION_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include <emcConstants.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

/**
 * @brief Radiative and Auger band-to-band recombination for bipolar EMC.
 *
 * Each call to recombine() draws the number of events from a Poisson
 * distribution and removes the corresponding electron-hole pairs from the
 * particle vectors. For Auger events the released recombination energy is
 * deposited into a randomly chosen third carrier.
 *
 * Rates (per unit volume per unit time):
 *
 *   Radiative:  R_rad  = B * n * p
 *   Auger CCCH: R_ccch = C_n * n^2 * p   (energy to third electron)
 *   Auger CHHS: R_chhs = C_p * n * p^2   (energy to third hole)
 *
 * Particles are removed using swap-with-last O(1) deletion so the vectors
 * stay compact. Index arithmetic accounts for the index shift that occurs
 * when the third Auger carrier happens to be the last element.
 *
 * Typical MAPbI3 parameters (SI units):
 *   B   ~ 1e-16 m^3/s  (= 1e-10 cm^3/s)
 *   C_n = C_p ~ 1e-39 m^6/s  (= 1e-27 cm^6/s)
 *   E_gap ~ 1.6 eV
 *
 * At n = p = 1e24 m^-3 these give tau_rad ~ tau_Auger ~ 1-10 ns, well
 * beyond the 10 ps hot-carrier cooling window. The Poisson sampler handles
 * the resulting sub-integer expected event count correctly.
 *
 * @tparam T  floating-point precision (float or double)
 */
template <class T>
class emcRecombination {
private:
  T B;        // radiative coefficient [m^3/s]
  T C_n;      // Auger CCCH coefficient [m^6/s]
  T C_p;      // Auger CHHS coefficient [m^6/s]
  T E_gap;    // band gap [eV]
  T effMassE; // electron effective mass [kg]
  T effMassH; // hole effective mass [kg]
  T Vsim;     // simulation volume [m^3]

  SizeType nRadiative = 0;
  SizeType nAuger = 0;

  mutable std::uniform_real_distribution<T> dist{T(0), T(1)};

  // Knuth's algorithm for Poisson(lambda); normal approximation for large lambda
  SizeType samplePoisson(T lambda, emcRNG &rng) const {
    if (lambda <= T(0))
      return 0;
    if (lambda > T(30)) {
      std::normal_distribution<T> nd(lambda, std::sqrt(lambda));
      T s = nd(rng);
      return s > T(0) ? SizeType(std::round(s)) : SizeType(0);
    }
    T L = std::exp(-lambda);
    SizeType k = 0;
    T p = T(1);
    do {
      p *= dist(rng);
      k++;
    } while (p > L);
    return k - 1;
  }

  // O(1) particle removal: swap with last element, then pop.
  template <class PartVec, class PosVec>
  static void removeParticle(PartVec &parts, PosVec &pos, SizeType idx) {
    if (idx != parts.size() - 1) {
      std::swap(parts[idx], parts.back());
      std::swap(pos[idx], pos.back());
    }
    parts.pop_back();
    pos.pop_back();
  }

  // Set carrier to new energy E_new [eV] with a uniformly random k direction.
  void energizeCarrier(emcParticle<T> &part, T effMass, T E_new,
                       emcRNG &rng) const {
    E_new = std::max(E_new, T(0));
    T kNew = std::sqrt(T(2) * effMass * E_new * constants::q) / constants::hbar;
    part.k = initRandomDirection(kNew, dist(rng), dist(rng));
    part.energy = E_new;
  }

public:
  /**
   * @param inB          radiative recombination coefficient [m^3/s]
   * @param inC_n        Auger coefficient for CCCH process [m^6/s]
   * @param inC_p        Auger coefficient for CHHS process [m^6/s]
   * @param inE_gap      band gap [eV]
   * @param relEffMassE  electron relative effective mass (m_star / m_e)
   * @param relEffMassH  hole relative effective mass (m_star / m_e)
   * @param inVsim       simulation volume [m^3]
   */
  emcRecombination(T inB, T inC_n, T inC_p, T inE_gap, T relEffMassE,
                   T relEffMassH, T inVsim)
      : B(inB), C_n(inC_n), C_p(inC_p), E_gap(inE_gap),
        effMassE(relEffMassE * constants::me),
        effMassH(relEffMassH * constants::me), Vsim(inVsim) {}

  /**
   * @brief Apply one time step of radiative and Auger recombination.
   *
   * @param electrons  electron particle vector (modified in place)
   * @param holes      hole particle vector (modified in place)
   * @param posE       electron position vector (kept in sync with electrons)
   * @param posH       hole position vector (kept in sync with holes)
   * @param dt         time step [s]
   * @param rng        random number generator
   */
  template <class PosVec>
  void recombine(std::vector<emcParticle<T>> &electrons,
                 std::vector<emcParticle<T>> &holes, PosVec &posE,
                 PosVec &posH, T dt, emcRNG &rng) {
    if (electrons.empty() || holes.empty())
      return;

    T n = T(electrons.size()) / Vsim;
    T p = T(holes.size()) / Vsim;

    // ── Radiative: e + h → photon ────────────────────────────────────
    {
      SizeType N = samplePoisson(B * n * p * Vsim * dt, rng);
      for (SizeType ev = 0;
           ev < N && !electrons.empty() && !holes.empty(); ++ev) {
        SizeType ie = SizeType(dist(rng) * electrons.size()) % electrons.size();
        SizeType ih = SizeType(dist(rng) * holes.size()) % holes.size();
        removeParticle(electrons, posE, ie);
        // holes unchanged by the electron removal
        removeParticle(holes, posH, ih);
        ++nRadiative;
      }
    }

    // ── Auger CCCH: e1 + h → e3 (third electron gains energy) ───────
    {
      SizeType N = samplePoisson(C_n * n * n * p * Vsim * dt, rng);
      for (SizeType ev = 0;
           ev < N && electrons.size() >= 2 && !holes.empty(); ++ev) {
        SizeType Ne = electrons.size();
        SizeType Nh = holes.size();

        SizeType ie1 = SizeType(dist(rng) * Ne) % Ne;
        SizeType ie3 = SizeType(dist(rng) * Ne) % Ne;
        if (ie3 == ie1)
          ie3 = (ie1 + 1) % Ne; // ensure distinct
        SizeType ih = SizeType(dist(rng) * Nh) % Nh;

        T E_recomb = E_gap + electrons[ie1].energy + holes[ih].energy;
        T E3_new   = electrons[ie3].energy + E_recomb;

        // If ie3 is the last element, removeParticle(ie1) will move it to ie1
        SizeType ie3_after = ie3;
        if (ie1 != Ne - 1 && ie3 == Ne - 1)
          ie3_after = ie1;

        removeParticle(electrons, posE, ie1);
        removeParticle(holes, posH, ih);
        energizeCarrier(electrons[ie3_after], effMassE, E3_new, rng);
        ++nAuger;
      }
    }

    // ── Auger CHHS: e + h1 → h3 (third hole gains energy) ───────────
    {
      SizeType N = samplePoisson(C_p * n * p * p * Vsim * dt, rng);
      for (SizeType ev = 0;
           ev < N && !electrons.empty() && holes.size() >= 2; ++ev) {
        SizeType Ne = electrons.size();
        SizeType Nh = holes.size();

        SizeType ie  = SizeType(dist(rng) * Ne) % Ne;
        SizeType ih1 = SizeType(dist(rng) * Nh) % Nh;
        SizeType ih3 = SizeType(dist(rng) * Nh) % Nh;
        if (ih3 == ih1)
          ih3 = (ih1 + 1) % Nh;

        T E_recomb = E_gap + electrons[ie].energy + holes[ih1].energy;
        T E3_new   = holes[ih3].energy + E_recomb;

        SizeType ih3_after = ih3;
        if (ih1 != Nh - 1 && ih3 == Nh - 1)
          ih3_after = ih1;

        removeParticle(electrons, posE, ie);
        // electrons removal does not affect hole indices
        removeParticle(holes, posH, ih1);
        energizeCarrier(holes[ih3_after], effMassH, E3_new, rng);
        ++nAuger;
      }
    }
  }

  SizeType getNrRadiative() const { return nRadiative; }
  SizeType getNrAuger()     const { return nAuger; }
};

#endif // EMC_RECOMBINATION_HPP
