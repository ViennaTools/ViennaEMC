#ifndef EMC_ENERGY_SELECTIVE_CONTACT_HPP
#define EMC_ENERGY_SELECTIVE_CONTACT_HPP

#include <algorithm>
#include <cassert>
#include <vector>

#include <emcConstants.hpp>
#include <emcParticle.hpp>

/**
 * @brief Energy-selective contact (ESC) for hot carrier solar cells.
 *
 * Models a contact that transmits only carriers whose kinetic energy (above
 * the band edge) falls in a narrow window [E_ex - deltaE/2, E_ex + deltaE/2].
 * Carriers outside the window are unaffected. The transmission rate is
 * controlled by a single time constant tau_ex:
 *
 *   P_extract = min(dt / tau_ex, 1)  per eligible carrier per step
 *
 * A second ESC object of the same type is used for holes, with independent
 * extraction energy and window width so the e/h quasi-Fermi levels can be
 * set independently (Würfel's HCSC principle).
 *
 * Implementation details
 * ──────────────────────
 * Particles are removed via swap-with-last (O(1)). The loop runs backward
 * so that after a removal at index i the element moved from the back is at
 * index i, which will not be revisited — correct because that element has a
 * higher original index and was already decided upon in a previous iteration.
 *
 * Diagnostics
 * ───────────
 * getNrExtracted()          cumulative extraction count
 * getMeanExtractedEnergy()  mean kinetic energy of extracted carriers [eV]
 * getCurrentDensity()       time-averaged current density [A/m^2]
 * resetCounters()           resets counters (e.g. for per-interval rates)
 *
 * Typical MAPbI3 HCSC starting point:
 *   E_ex   ~ 0.4 eV    (near mean carrier energy at T_c ~ 3000 K)
 *   deltaE ~ 0.05 eV   (narrow window maximises extracted chemical potential)
 *   tau_ex ~ 1 ps      (comparable to LO phonon scattering time)
 *
 * @tparam T  floating-point precision (float or double)
 */
template <class T>
class emcEnergySelectiveContact {
private:
  T E_ex;       // centre of extraction window [eV]
  T halfDeltaE; // half-width of extraction window [eV]
  T tauEx;      // extraction time constant [s]

  SizeType nExtracted = 0;
  T sumExtractedEnergy = T(0);

  mutable std::uniform_real_distribution<T> dist{T(0), T(1)};

  template <class PartVec, class PosVec>
  static void removeParticle(PartVec &parts, PosVec &pos, SizeType idx) {
    if (idx != parts.size() - 1) {
      std::swap(parts[idx], parts.back());
      std::swap(pos[idx], pos.back());
    }
    parts.pop_back();
    pos.pop_back();
  }

public:
  /**
   * @param inE_ex     centre of the extraction energy window [eV]
   * @param inDeltaE   full width of the extraction window [eV]
   * @param inTauEx    extraction time constant [s]; 1/tau_ex is the rate at
   *                   which an eligible carrier is removed per unit time
   */
  emcEnergySelectiveContact(T inE_ex, T inDeltaE, T inTauEx)
      : E_ex(inE_ex), halfDeltaE(T(0.5) * inDeltaE), tauEx(inTauEx) {
    assert(inDeltaE > T(0) && "ESC window width must be positive");
    assert(inTauEx  > T(0) && "ESC extraction time constant must be positive");
  }

  /**
   * @brief Stochastically extract carriers in the energy window for one step.
   *
   * Iterates backward so swap-with-last removal is safe without index fixup.
   *
   * @param particles  carrier vector (modified in place)
   * @param pos        position vector (kept in sync with particles)
   * @param dt         time step [s]
   * @param rng        random number generator
   */
  template <class PosVec>
  void extract(std::vector<emcParticle<T>> &particles, PosVec &pos, T dt,
               emcRNG &rng) {
    if (particles.empty())
      return;

    T P = std::min(dt / tauEx, T(1));

    // Backward iteration: when element at i is removed (swap with back + pop),
    // the old back element lands at i and is never revisited. Since we process
    // indices from high to low, that element was already evaluated and kept —
    // so skipping it again is correct.
    for (int i = static_cast<int>(particles.size()) - 1; i >= 0; --i) {
      T E = particles[static_cast<SizeType>(i)].energy;
      if (std::fabs(E - E_ex) <= halfDeltaE && dist(rng) < P) {
        sumExtractedEnergy += E;
        ++nExtracted;
        removeParticle(particles, pos, static_cast<SizeType>(i));
      }
    }
  }

  //! Cumulative number of carriers extracted since construction (or last reset).
  SizeType getNrExtracted() const { return nExtracted; }

  //! Mean kinetic energy [eV] of all extracted carriers; 0 if none extracted.
  T getMeanExtractedEnergy() const {
    return nExtracted > 0 ? sumExtractedEnergy / static_cast<T>(nExtracted)
                          : T(0);
  }

  /**
   * @brief Time-averaged extracted current density [A/m^2].
   *
   * @param contactArea  physical contact area [m^2]
   * @param elapsedTime  total simulated time over which extraction occurred [s]
   */
  T getCurrentDensity(T contactArea, T elapsedTime) const {
    if (nExtracted == 0 || elapsedTime <= T(0))
      return T(0);
    return static_cast<T>(nExtracted) * constants::q /
           (contactArea * elapsedTime);
  }

  //! Reset cumulative counters (e.g. between output intervals for rate tracking).
  void resetCounters() {
    nExtracted = 0;
    sumExtractedEnergy = T(0);
  }

  T getExtractionEnergy()  const { return E_ex; }
  T getWindowWidth()       const { return T(2) * halfDeltaE; }
  T getExtractionTimeCst() const { return tauEx; }
};

#endif // EMC_ENERGY_SELECTIVE_CONTACT_HPP
