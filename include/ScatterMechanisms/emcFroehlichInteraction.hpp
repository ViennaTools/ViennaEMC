#ifndef EMC_FROEHLICH_INTERACTION_HPP
#define EMC_FROEHLICH_INTERACTION_HPP

#include <cassert>
#include <cmath>
#include <string>

#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcUtil.hpp>

/**
 * 3D Fröhlich (polar optical phonon) scatter rate for bulk semiconductors.
 *
 * Scatter rate for Kane non-parabolic bands (Jacoboni & Reggiani 1983,
 * Rev. Mod. Phys. 55, 645; eq. 5.122 of Jacoboni & Lugli 2011):
 *
 *   Γ(E) = C * N_q / k_i * ln((k_i + k_f) / |k_i - k_f|)
 *
 * where
 *   C     = e^2 w0 m_star / (4pi hbar^2) * (1/eps_inf - 1/eps_s) / eps0_SI   [m^{-1} s^{-1}]
 *   k_i   = √(2 m* γ(E) q) / ħ        [m^{-1}]
 *   k_f   = √(2 m* γ(E') q) / ħ       [m^{-1}]  (E' = E ± ħω₀)
 *   γ(E)  = E(1 + α E)                [eV, Kane gamma function]
 *   N_q   = equilibrium Bose-Einstein occupation = 1 / (exp(ħω₀/kT) - 1)
 *
 * For absorption N_q appears directly; for emission the factor is N_q + 1.
 *
 * The scatterParticle() method samples the 3D Fröhlich angular distribution
 * (Fröhlich 1939; see also Ferry "Semiconductor Transport" p. 220):
 *   cos θ = (1 + f - (1+2f)^r) / f
 * where f = 2√(E_i E_f) / (√E_i - √E_f)² and r ∈ [0,1) is a random number.
 *
 * @param phononEnergy   LO phonon energy [eV]
 * @param effMass        band-edge effective mass [kg]
 * @param eps_hi         optical (high-frequency) relative dielectric constant
 * @param eps_lo         static (low-frequency) relative dielectric constant
 * @param scatterConst   precomputed rate prefactor [m^{-1} s^{-1}]
 * @param N_bose         equilibrium Bose-Einstein phonon occupation
 */

// ---- helper: precompute common Fröhlich constant ----------------------------
// C = e² ω₀ m* / (4π ħ²) * (1/ε_hi - 1/ε_lo) / ε₀   [m^{-1} s^{-1}]
template <class T>
T froehlichScatterConst3D(T phononEnergy, T effMass, T eps_hi, T eps_lo) {
  T omega0 = phononEnergy * constants::q / constants::hbar;
  return constants::q * constants::q * omega0 * effMass /
         (4. * constants::pi * constants::hbar * constants::hbar) *
         (T(1) / eps_hi - T(1) / eps_lo) / constants::eps0;
}

// ---- helper: equilibrium Bose-Einstein occupation ---------------------------
template <class T>
T boseEinstein(T phononEnergy, T temperature) {
  T x = constants::q * phononEnergy / (constants::kB * temperature);
  return T(1) / (std::exp(x) - T(1));
}

// =============================================================================
/** 3D Fröhlich absorption: carrier absorbs an LO phonon (E → E + ħω₀). */
template <class T>
class emcFroehlichAbsorption3D : public emcScatterMechanism<T> {
private:
  T phononEnergy; // [eV]
  T effMass;      // [kg]
  T scatterConst; // C [m^{-1} s^{-1}]
  T N_bose;       // equilibrium Bose-Einstein occupation
  std::string nameSuffix;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcFroehlichAbsorption3D() = delete;

  /**
   * @param inValley       index of the valley (same for initial and final)
   * @param inPhononEnergy LO phonon energy [eV]
   * @param relEffMass     relative effective mass (m_star / m_e)
   * @param eps_hi         optical (high-frequency) relative dielectric constant
   * @param eps_lo         static (low-frequency) relative dielectric constant
   * @param temperature    lattice temperature [K]
   * @param inNameSuffix   suffix for output file name
   */
  emcFroehlichAbsorption3D(SizeType inValley, T inPhononEnergy, T relEffMass,
                            T eps_hi, T eps_lo, T temperature,
                            std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        effMass(relEffMass * constants::me), nameSuffix(inNameSuffix),
        dist(0., 1.) {
    scatterConst = froehlichScatterConst3D(phononEnergy, effMass, eps_hi, eps_lo);
    N_bose = boseEinstein(phononEnergy, temperature);
  }

  std::string getName() const { return "FroehlichAbsorption3D" + nameSuffix; }

  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    auto valley = this->ptrValley[this->idxValley];
    T gammaI = valley->getGamma(energy);
    T energyFinal = energy + phononEnergy;
    T gammaF = valley->getGamma(energyFinal);
    if (gammaI <= T(0) || gammaF <= T(0))
      return T(0);
    T kI = std::sqrt(T(2) * effMass * gammaI * constants::q) / constants::hbar;
    T kF = std::sqrt(T(2) * effMass * gammaF * constants::q) / constants::hbar;
    // kF > kI for absorption; log argument is always > 1
    T lnFactor = std::log((kI + kF) / (kF - kI));
    return scatterConst * N_bose / kI * lnFactor;
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    T initEnergy = particle.energy;
    particle.energy += phononEnergy;
    T finalEnergy = particle.energy;

    // sample scattering angle from 3D Fröhlich distribution
    T f = T(2) * std::sqrt(initEnergy * finalEnergy) /
          (std::sqrt(initEnergy) - std::sqrt(finalEnergy)) /
          (std::sqrt(initEnergy) - std::sqrt(finalEnergy));
    T cosTheta = (T(1) + f - std::pow(T(1) + T(2) * f, dist(rng))) / f;
    cosTheta = std::max(T(-1), std::min(T(1), cosTheta));

    T kNew = this->ptrValley[this->idxValley]->getNormWaveVec(particle.energy);
    particle.k = initRandomDirectionWithRespectToCurrentK(particle.k, cosTheta,
                                                          dist(rng));
    // rescale to correct magnitude
    T kCurr = norm(particle.k);
    if (kCurr > T(0))
      particle.k = scale(particle.k, kNew / kCurr);
  }
};

// =============================================================================
/** 3D Fröhlich emission: carrier emits an LO phonon (E → E − ħω₀). */
template <class T>
class emcFroehlichEmission3D : public emcScatterMechanism<T> {
private:
  T phononEnergy;
  T effMass;
  T scatterConst;
  T N_bose;
  std::string nameSuffix;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcFroehlichEmission3D() = delete;

  /**
   * @param inValley       index of the valley
   * @param inPhononEnergy LO phonon energy [eV]
   * @param relEffMass     relative effective mass (m_star / m_e)
   * @param eps_hi         optical (high-frequency) relative dielectric constant
   * @param eps_lo         static (low-frequency) relative dielectric constant
   * @param temperature    lattice temperature [K]
   * @param inNameSuffix   suffix for output file name
   */
  emcFroehlichEmission3D(SizeType inValley, T inPhononEnergy, T relEffMass,
                          T eps_hi, T eps_lo, T temperature,
                          std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        effMass(relEffMass * constants::me), nameSuffix(inNameSuffix),
        dist(0., 1.) {
    scatterConst = froehlichScatterConst3D(phononEnergy, effMass, eps_hi, eps_lo);
    N_bose = boseEinstein(phononEnergy, temperature);
  }

  std::string getName() const { return "FroehlichEmission3D" + nameSuffix; }

  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    if (energy <= phononEnergy)
      return T(0);
    auto valley = this->ptrValley[this->idxValley];
    T gammaI = valley->getGamma(energy);
    T energyFinal = energy - phononEnergy;
    T gammaF = valley->getGamma(energyFinal);
    if (gammaI <= T(0) || gammaF <= T(0))
      return T(0);
    T kI = std::sqrt(T(2) * effMass * gammaI * constants::q) / constants::hbar;
    T kF = std::sqrt(T(2) * effMass * gammaF * constants::q) / constants::hbar;
    // kI > kF for emission; log argument is always > 1
    T lnFactor = std::log((kI + kF) / (kI - kF));
    return scatterConst * (N_bose + T(1)) / kI * lnFactor;
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    T initEnergy = particle.energy;
    particle.energy -= phononEnergy;
    assert(particle.energy > T(0));
    T finalEnergy = particle.energy;

    T f = T(2) * std::sqrt(initEnergy * finalEnergy) /
          (std::sqrt(initEnergy) - std::sqrt(finalEnergy)) /
          (std::sqrt(initEnergy) - std::sqrt(finalEnergy));
    T cosTheta = (T(1) + f - std::pow(T(1) + T(2) * f, dist(rng))) / f;
    cosTheta = std::max(T(-1), std::min(T(1), cosTheta));

    T kNew = this->ptrValley[this->idxValley]->getNormWaveVec(particle.energy);
    particle.k = initRandomDirectionWithRespectToCurrentK(particle.k, cosTheta,
                                                          dist(rng));
    T kCurr = norm(particle.k);
    if (kCurr > T(0))
      particle.k = scale(particle.k, kNew / kCurr);
  }
};

#endif // EMC_FROEHLICH_INTERACTION_HPP
