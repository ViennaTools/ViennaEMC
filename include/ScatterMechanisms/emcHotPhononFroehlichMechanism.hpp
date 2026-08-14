#ifndef EMC_HOT_PHONON_FROEHLICH_MECHANISM_HPP
#define EMC_HOT_PHONON_FROEHLICH_MECHANISM_HPP

#include <cassert>
#include <cmath>
#include <memory>
#include <string>

#include <ScatterMechanisms/emcFroehlichInteraction.hpp>
#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcPhononBath.hpp>
#include <emcUtil.hpp>

/**
 * 3D Fröhlich (polar optical) absorption with hot phonon bottleneck (HPB).
 *
 * The scatter rate uses the non-equilibrium phonon occupation N_q from a
 * shared emcPhononBath rather than the equilibrium Bose-Einstein value.
 * After each scatter event the phonon wavevector magnitude |q| = |k_i - k_f|
 * is recorded in the bath for the population-evolution bookkeeping.
 *
 * Algorithm: Lugli (Solid-State Electron. 1988, 31, 667), as implemented in
 * Faber, Filipovic, Koster, J. Phys. Chem. Lett. 2024, 15, 12601.
 *
 * Simulation loop pseudo-code:
 * @code
 *   for each timestep dt:
 *     handler.moveParticles(dt);           // scatter events recorded in bath
 *     phononBath->update(dt);              // evolve N_q, reset counters
 *     scatterHandler.reinitScatterTables(); // rebuild tables with new N_q
 * @endcode
 *
 * @param phononEnergy   LO phonon energy [eV]
 * @param effMass        band-edge effective mass [kg]
 * @param scatterConst   rate prefactor [m^{-1} s^{-1}]
 * @param phononBath     shared phonon bath (mutable: records scatter events)
 */
template <class T>
class emcHotPhononFroehlichAbsorption3D : public emcScatterMechanism<T> {
private:
  T phononEnergy;
  T effMass;
  T scatterConst;
  std::shared_ptr<emcPhononBath<T>> phononBath;
  std::string nameSuffix;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcHotPhononFroehlichAbsorption3D() = delete;

  /**
   * @param inValley       valley index
   * @param inPhononEnergy LO phonon energy [eV]
   * @param relEffMass     relative effective mass (m_star / m_e)
   * @param eps_hi         optical (high-frequency) relative dielectric constant
   * @param eps_lo         static (low-frequency) relative dielectric constant
   * @param inPhononBath   shared phonon bath (owns N_q state)
   * @param inNameSuffix   output file name suffix
   */
  emcHotPhononFroehlichAbsorption3D(
      SizeType inValley, T inPhononEnergy, T relEffMass, T eps_hi, T eps_lo,
      std::shared_ptr<emcPhononBath<T>> inPhononBath,
      std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        effMass(relEffMass * constants::me), phononBath(std::move(inPhononBath)),
        nameSuffix(inNameSuffix), dist(0., 1.) {
    scatterConst =
        froehlichScatterConst3D(phononEnergy, effMass, eps_hi, eps_lo);
  }

  std::string getName() const {
    return "HotPhononFroehlichAbsorption3D" + nameSuffix;
  }

  /// Scatter rate using current non-equilibrium N_q from the phonon bath.
  /// Uses the DOS-weighted average N_q since the rate table is energy-resolved
  /// but not wavevector-resolved (single-mode approximation).
  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    auto valley = this->ptrValley[this->idxValley];
    T gammaI = valley->getGamma(energy);
    T energyFinal = energy + phononEnergy;
    T gammaF = valley->getGamma(energyFinal);
    if (gammaI <= T(0) || gammaF <= T(0))
      return T(0);
    T kI = std::sqrt(T(2) * effMass * gammaI * constants::q) / constants::hbar;
    T kF = std::sqrt(T(2) * effMass * gammaF * constants::q) / constants::hbar;
    T lnFactor = std::log((kI + kF) / (kF - kI));
    T Nq = phononBath->getMeanNq();
    return scatterConst * Nq / kI * lnFactor;
  }

  /// Scatter particle and record the phonon wavevector in the bath.
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    std::array<T, 3> kOld = particle.k;
    T initEnergy = particle.energy;
    particle.energy += phononEnergy;
    T finalEnergy = particle.energy;

    T f = T(2) * std::sqrt(initEnergy * finalEnergy) /
          (std::sqrt(initEnergy) - std::sqrt(finalEnergy)) /
          (std::sqrt(initEnergy) - std::sqrt(finalEnergy));
    T cosTheta = (T(1) + f - std::pow(T(1) + T(2) * f, dist(rng))) / f;
    cosTheta = std::max(T(-1), std::min(T(1), cosTheta));

    T kNew = this->ptrValley[this->idxValley]->getNormWaveVec(particle.energy);
    particle.k =
        initRandomDirectionWithRespectToCurrentK(particle.k, cosTheta, dist(rng));
    T kCurr = norm(particle.k);
    if (kCurr > T(0))
      particle.k = scale(particle.k, kNew / kCurr);

    // record which phonon mode was involved: |q| = |k_new - k_old|
    std::array<T, 3> qVec = subtract(particle.k, kOld);
    T qMag = norm(qVec);
    phononBath->recordAbsorption(qMag); // carrier absorbs → phonon population --
  }
};

// =============================================================================
/** 3D Fröhlich emission with non-equilibrium phonon bath (HPB). */
template <class T>
class emcHotPhononFroehlichEmission3D : public emcScatterMechanism<T> {
private:
  T phononEnergy;
  T effMass;
  T scatterConst;
  std::shared_ptr<emcPhononBath<T>> phononBath;
  std::string nameSuffix;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcHotPhononFroehlichEmission3D() = delete;

  /**
   * @param inValley       valley index
   * @param inPhononEnergy LO phonon energy [eV]
   * @param relEffMass     relative effective mass (m_star / m_e)
   * @param eps_hi         optical (high-frequency) relative dielectric constant
   * @param eps_lo         static (low-frequency) relative dielectric constant
   * @param inPhononBath   shared phonon bath
   * @param inNameSuffix   output file name suffix
   */
  emcHotPhononFroehlichEmission3D(
      SizeType inValley, T inPhononEnergy, T relEffMass, T eps_hi, T eps_lo,
      std::shared_ptr<emcPhononBath<T>> inPhononBath,
      std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        effMass(relEffMass * constants::me), phononBath(std::move(inPhononBath)),
        nameSuffix(inNameSuffix), dist(0., 1.) {
    scatterConst =
        froehlichScatterConst3D(phononEnergy, effMass, eps_hi, eps_lo);
  }

  std::string getName() const {
    return "HotPhononFroehlichEmission3D" + nameSuffix;
  }

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
    T lnFactor = std::log((kI + kF) / (kI - kF));
    T Nq = phononBath->getMeanNq();
    return scatterConst * (Nq + T(1)) / kI * lnFactor;
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    std::array<T, 3> kOld = particle.k;
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
    particle.k =
        initRandomDirectionWithRespectToCurrentK(particle.k, cosTheta, dist(rng));
    T kCurr = norm(particle.k);
    if (kCurr > T(0))
      particle.k = scale(particle.k, kNew / kCurr);

    // record emitted phonon wavevector
    std::array<T, 3> qVec = subtract(particle.k, kOld);
    T qMag = norm(qVec);
    phononBath->recordEmission(qMag); // carrier emits → phonon population ++
  }
};

#endif // EMC_HOT_PHONON_FROEHLICH_MECHANISM_HPP
