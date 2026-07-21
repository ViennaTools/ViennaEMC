#ifndef EMC_ADSORBATE_NOISE_HPP
#define EMC_ADSORBATE_NOISE_HPP

#include <cmath>
#include <random>
#include <vector>

#include <emcConstants.hpp>
#include <emcUtil.hpp>

/*! \file emcAdsorbateNoise.hpp
 * \brief Kinetic Monte Carlo of dynamic adsorption/desorption for
 * low-frequency conductivity noise (ambient channel C4).
 *
 * Each surface adsorption site switches stochastically between empty and
 * occupied: an empty site adsorbs at rate r_ads = k_ads * P (proportional to
 * the partial pressure), an occupied site desorbs at rate
 * r_des = nu * exp(-E_des / kB T) (Arrhenius, E_des the desorption barrier from
 * DFT / TPD). The occupied count N_ads(t) therefore fluctuates in time; through
 * charge transfer (C1) and adsorbate Coulomb scattering (C2) this modulates the
 * carrier density and mobility, so the conductivity fluctuates:
 * S_sigma(f) / sigma^2 = (d ln sigma / d N_ads)^2 * S_{N_ads}(f).
 *
 * Spectral shapes:
 *   - a single site (one time constant) -> a random-telegraph signal with a
 *     Lorentzian power spectrum, corner f_c = (r_ads + r_des) / 2 pi;
 *   - many sites with a DISTRIBUTION of desorption barriers (hence a
 *     log-uniform spread of time constants) -> a superposition of Lorentzians
 *     that gives 1/f noise over the corresponding frequency band
 *     (McWhorter mechanism).
 *
 * The equilibrium coverage of site i is theta_i = r_ads / (r_ads + r_des,i),
 * consistent with the Langmuir isotherm of emcAdsorbateChargeTransfer (C1).
 *
 * References: McWhorter (1957); low-frequency noise in MoS2 FETs — Nano Lett.
 * (2013), APL Mater. 2, 092515 (2014); IEEE T-ED (2018).
 */

/*! \brief Arrhenius desorption rate r_des = nu exp(-E_des / kB T) [1/s]. */
template <class T>
T arrheniusDesorptionRate(T attemptFrequency, T desorptionEnergy,
                          T temperature) {
  return attemptFrequency * std::exp(-desorptionEnergy * constants::q /
                                     (constants::kB * temperature));
}

/*! \brief KMC of a set of two-state adsorption sites (per-site rates).
 *
 * @param adsorptionRates per-site adsorption rate r_ads [1/s] (physically
 *        r_ads = k_ads P, common to all sites; a per-site vector also allows the
 *        theta = 0.5 McWhorter idealization)
 * @param desorptionRates per-site desorption rate r_des [1/s] (size sets the
 *        number of sites; one value for RTN, a distribution for 1/f)
 */
template <class T> class emcAdsorbateNoise {
private:
  std::vector<T> adsorptionRates;
  std::vector<T> desorptionRates;
  std::vector<char> occupied;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcAdsorbateNoise() = delete;

  emcAdsorbateNoise(const std::vector<T> &inAdsorptionRates,
                    const std::vector<T> &inDesorptionRates)
      : adsorptionRates(inAdsorptionRates), desorptionRates(inDesorptionRates),
        occupied(inDesorptionRates.size(), 0), dist(0., 1.) {}

  /*! \brief Convenience: common adsorption rate r_ads for all sites. */
  emcAdsorbateNoise(T commonAdsorptionRate,
                    const std::vector<T> &inDesorptionRates)
      : emcAdsorbateNoise(
            std::vector<T>(inDesorptionRates.size(), commonAdsorptionRate),
            inDesorptionRates) {}

  SizeType nSites() const { return occupied.size(); }

  /*! \brief Initialize each site to its equilibrium occupation probability
   *  theta_i = r_ads / (r_ads + r_des), so the run starts in steady state. */
  void seedEquilibrium(emcRNG &rng) {
    for (SizeType i = 0; i < occupied.size(); ++i) {
      T theta = adsorptionRates[i] / (adsorptionRates[i] + desorptionRates[i]);
      occupied[i] = (dist(rng) < theta) ? 1 : 0;
    }
  }

  /*! \brief Advance all sites by dt (must be << 1 / max rate) and return the
   *  number of occupied sites N_ads. */
  SizeType step(T dt, emcRNG &rng) {
    SizeType n = 0;
    for (SizeType i = 0; i < occupied.size(); ++i) {
      if (occupied[i]) {
        if (dist(rng) < desorptionRates[i] * dt)
          occupied[i] = 0;
      } else {
        if (dist(rng) < adsorptionRates[i] * dt)
          occupied[i] = 1;
      }
      n += occupied[i];
    }
    return n;
  }
};

#endif // EMC_ADSORBATE_NOISE_HPP
