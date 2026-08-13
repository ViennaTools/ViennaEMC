#ifndef EMC_SCREENED_FROEHLICH_INTERACTION_HPP
#define EMC_SCREENED_FROEHLICH_INTERACTION_HPP

#include <cassert>
#include <cmath>
#include <memory>
#include <string>

#include <ScatterMechanisms/emcFroehlichInteraction.hpp>
#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>
#include <emcPhononBath.hpp>
#include <emcPlasmonScreening.hpp>
#include <emcUtil.hpp>

/**
 * 3D Fröhlich (polar optical) scattering with free-carrier (plasmon) screening,
 * in both equilibrium-phonon and non-equilibrium (hot-phonon) variants.
 *
 * Screening replaces the bare Coulomb-like vertex 1/q² by 1/(q² + q_s²), with
 * q_s supplied by a shared emcPlasmonScreening object that the simulation loop
 * refreshes each step. Two things change, and BOTH matter for transport:
 *
 * 1. THE RATE. The angular integral
 *        ∫ q dq / q²  = ln(q₊/q₋)
 *    becomes
 *        ∫ q dq / (q² + q_s²) = ½ ln((q₊² + q_s²)/(q₋² + q_s²))
 *    with q₊ = k_i + k_f and q₋ = |k_i − k_f|.
 *
 * 2. THE ANGLE. Screening suppresses small-q (small-angle) events
 *    preferentially, so the scattering becomes more isotropic and each event
 *    relaxes more momentum. Screening the rate alone would get the transport
 *    effect qualitatively wrong, so the final-state sampling is screened too.
 *
 *    With A' = k_i² + k_f² + q_s² and B = 2 k_i k_f, inverting the cumulative
 *    distribution of P(cosθ) ∝ 1/(A' − B cosθ) gives
 *
 *        cosθ = [A' − (A' + B) ((A' − B)/(A' + B))^r] / B,   r ∈ [0,1)
 *
 *    which reduces exactly to the standard unscreened Fröhlich sampling
 *    cosθ = [1 + f − (1 + 2f)^(1−r)]/f,  f = B/(A − B),  when q_s → 0.
 *
 * Setting q_s = 0 (or disabling the screening object) reproduces the unscreened
 * mechanisms in emcFroehlichInteraction.hpp / emcHotPhononFroehlichMechanism.hpp
 * bit-for-bit in the rate and statistically in the angle.
 *
 * See emcPlasmonScreening.hpp for why the static screening limit used here is a
 * CONSERVATIVE (over-screening) bound on the polar interaction.
 */

// ---- helper: screened angular log factor ------------------------------------
// ½ ln((q₊² + q_s²)/(q₋² + q_s²)); reduces to ln(q₊/q₋) when q_s = 0.
template <class T>
T screenedFroehlichLogFactor(T kI, T kF, T qs2) {
  const T qPlus2 = (kI + kF) * (kI + kF);
  const T qMinus2 = (kI - kF) * (kI - kF);
  if (qs2 <= T(0)) {
    if (qMinus2 <= T(0))
      return T(0);
    return T(0.5) * std::log(qPlus2 / qMinus2);
  }
  return T(0.5) * std::log((qPlus2 + qs2) / (qMinus2 + qs2));
}

// ---- helper: screened Fröhlich final-state polar angle -----------------------
template <class T>
T sampleScreenedFroehlichCosTheta(T kI, T kF, T qs2, T r) {
  const T B = T(2) * kI * kF;
  if (B <= T(0))
    return T(1) - T(2) * r; // isotropic fallback
  const T Ap = kI * kI + kF * kF + qs2;
  const T num = Ap - B;
  const T den = Ap + B;
  if (num <= T(0) || den <= T(0))
    return T(1) - T(2) * r;
  const T cosTheta = (Ap - den * std::pow(num / den, r)) / B;
  return std::max(T(-1), std::min(T(1), cosTheta));
}


// ---- helper: q-resolved final-state angle -----------------------------------
// When the phonon occupation varies across the allowed window, the angle must
// be drawn from the SAME distribution that sets the rate. Samples q from
// P(q) ~ f(q)/q via the bath's prefix sums, then converts to cos(theta).
// Screening enters through the bath's weight (setScreeningQ2), so this is
// valid screened or not, provided the caller keeps the bath in sync.
template <class T>
T sampleQResolvedCosTheta(const emcPhononBath<T> &bath, T kI, T kF,
                          bool emission, T r) {
  const T B = T(2) * kI * kF;
  if (B <= T(0))
    return T(1) - T(2) * r;
  const T q =
      bath.sampleQ(std::fabs(kI - kF), kI + kF, emission, r);
  const T cosTheta = (kI * kI + kF * kF - q * q) / B;
  return std::max(T(-1), std::min(T(1), cosTheta));
}

// =============================================================================
/** Screened 3D Fröhlich absorption with equilibrium phonon occupation. */
template <class T>
class emcScreenedFroehlichAbsorption3D : public emcScatterMechanism<T> {
private:
  T phononEnergy, effMass, scatterConst, N_bose;
  std::shared_ptr<emcPlasmonScreening<T>> screening;
  std::string nameSuffix;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcScreenedFroehlichAbsorption3D() = delete;

  emcScreenedFroehlichAbsorption3D(
      SizeType inValley, T inPhononEnergy, T relEffMass, T eps_hi, T eps_lo,
      T temperature, std::shared_ptr<emcPlasmonScreening<T>> inScreening,
      std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        effMass(relEffMass * constants::me), screening(std::move(inScreening)),
        nameSuffix(inNameSuffix), dist(0., 1.) {
    scatterConst =
        froehlichScatterConst3D(phononEnergy, effMass, eps_hi, eps_lo);
    N_bose = boseEinstein(phononEnergy, temperature);
  }

  std::string getName() const {
    return "ScreenedFroehlichAbsorption3D" + nameSuffix;
  }

  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    auto valley = this->ptrValley[this->idxValley];
    T gammaI = valley->getGamma(energy);
    T gammaF = valley->getGamma(energy + phononEnergy);
    if (gammaI <= T(0) || gammaF <= T(0))
      return T(0);
    T kI = std::sqrt(T(2) * effMass * gammaI * constants::q) / constants::hbar;
    T kF = std::sqrt(T(2) * effMass * gammaF * constants::q) / constants::hbar;
    return scatterConst * N_bose / kI *
           screenedFroehlichLogFactor(kI, kF, screening->getQs2());
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto valley = this->ptrValley[this->idxValley];
    T kI = valley->getNormWaveVec(particle.energy);
    particle.energy += phononEnergy;
    T kF = valley->getNormWaveVec(particle.energy);

    T cosTheta =
        sampleScreenedFroehlichCosTheta(kI, kF, screening->getQs2(), dist(rng));
    particle.k =
        initRandomDirectionWithRespectToCurrentK(particle.k, cosTheta, dist(rng));
    T kCurr = norm(particle.k);
    if (kCurr > T(0))
      particle.k = scale(particle.k, kF / kCurr);
  }
};

// =============================================================================
/** Screened 3D Fröhlich emission with equilibrium phonon occupation. */
template <class T>
class emcScreenedFroehlichEmission3D : public emcScatterMechanism<T> {
private:
  T phononEnergy, effMass, scatterConst, N_bose;
  std::shared_ptr<emcPlasmonScreening<T>> screening;
  std::string nameSuffix;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcScreenedFroehlichEmission3D() = delete;

  emcScreenedFroehlichEmission3D(
      SizeType inValley, T inPhononEnergy, T relEffMass, T eps_hi, T eps_lo,
      T temperature, std::shared_ptr<emcPlasmonScreening<T>> inScreening,
      std::string inNameSuffix = "")
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        effMass(relEffMass * constants::me), screening(std::move(inScreening)),
        nameSuffix(inNameSuffix), dist(0., 1.) {
    scatterConst =
        froehlichScatterConst3D(phononEnergy, effMass, eps_hi, eps_lo);
    N_bose = boseEinstein(phononEnergy, temperature);
  }

  std::string getName() const {
    return "ScreenedFroehlichEmission3D" + nameSuffix;
  }

  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    if (energy <= phononEnergy)
      return T(0);
    auto valley = this->ptrValley[this->idxValley];
    T gammaI = valley->getGamma(energy);
    T gammaF = valley->getGamma(energy - phononEnergy);
    if (gammaI <= T(0) || gammaF <= T(0))
      return T(0);
    T kI = std::sqrt(T(2) * effMass * gammaI * constants::q) / constants::hbar;
    T kF = std::sqrt(T(2) * effMass * gammaF * constants::q) / constants::hbar;
    return scatterConst * (N_bose + T(1)) / kI *
           screenedFroehlichLogFactor(kI, kF, screening->getQs2());
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto valley = this->ptrValley[this->idxValley];
    T kI = valley->getNormWaveVec(particle.energy);
    particle.energy -= phononEnergy;
    assert(particle.energy > T(0));
    T kF = valley->getNormWaveVec(particle.energy);

    T cosTheta =
        sampleScreenedFroehlichCosTheta(kI, kF, screening->getQs2(), dist(rng));
    particle.k =
        initRandomDirectionWithRespectToCurrentK(particle.k, cosTheta, dist(rng));
    T kCurr = norm(particle.k);
    if (kCurr > T(0))
      particle.k = scale(particle.k, kF / kCurr);
  }
};

// =============================================================================
/** Screened 3D Fröhlich absorption with non-equilibrium (hot) phonon bath. */
template <class T>
class emcScreenedHotPhononFroehlichAbsorption3D
    : public emcScatterMechanism<T> {
private:
  T phononEnergy, effMass, scatterConst;
  std::shared_ptr<emcPhononBath<T>> phononBath;
  std::shared_ptr<emcPlasmonScreening<T>> screening;
  bool qResolved;
  // Separate switch for the final-state angle. Setting this false while
  // qResolved is true reproduces the INCONSISTENT treatment discussed in the
  // text: the rate follows the resolved distribution while the angle is drawn
  // as if the occupation were uniform. Provided so that combination can be
  // measured rather than argued about; it is not physically correct.
  bool qResolvedAngle;
  std::string nameSuffix;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcScreenedHotPhononFroehlichAbsorption3D() = delete;

  emcScreenedHotPhononFroehlichAbsorption3D(
      SizeType inValley, T inPhononEnergy, T relEffMass, T eps_hi, T eps_lo,
      std::shared_ptr<emcPhononBath<T>> inPhononBath,
      std::shared_ptr<emcPlasmonScreening<T>> inScreening,
      bool inQResolved = false, std::string inNameSuffix = "",
      bool inQResolvedAngle = true)
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        effMass(relEffMass * constants::me), phononBath(std::move(inPhononBath)),
        screening(std::move(inScreening)), qResolved(inQResolved),
        qResolvedAngle(inQResolvedAngle), nameSuffix(inNameSuffix),
        dist(0., 1.) {
    scatterConst =
        froehlichScatterConst3D(phononEnergy, effMass, eps_hi, eps_lo);
  }

  std::string getName() const {
    return "ScreenedHotPhononFroehlichAbsorption3D" + nameSuffix;
  }

  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    auto valley = this->ptrValley[this->idxValley];
    T gammaI = valley->getGamma(energy);
    T gammaF = valley->getGamma(energy + phononEnergy);
    if (gammaI <= T(0) || gammaF <= T(0))
      return T(0);
    T kI = std::sqrt(T(2) * effMass * gammaI * constants::q) / constants::hbar;
    T kF = std::sqrt(T(2) * effMass * gammaF * constants::q) / constants::hbar;
    const T Nq = qResolved
                     ? phononBath->getNqInWindow(std::fabs(kI - kF), kI + kF)
                     : phononBath->getMeanNq();
    return scatterConst * Nq / kI *
           screenedFroehlichLogFactor(kI, kF, screening->getQs2());
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto valley = this->ptrValley[this->idxValley];
    std::array<T, 3> kOld = particle.k;
    T kI = valley->getNormWaveVec(particle.energy);
    particle.energy += phononEnergy;
    T kF = valley->getNormWaveVec(particle.energy);

    // Angle must come from the same distribution as the rate. The bath's
    // transition weight already carries q_s (setScreeningQ2), so the q-resolved
    // sampler is correct screened or unscreened.
    const T qs2 = screening->getQs2();
    T cosTheta =
        (qResolved && qResolvedAngle)
            ? sampleQResolvedCosTheta(*phononBath, kI, kF, false, dist(rng))
            : sampleScreenedFroehlichCosTheta(kI, kF, qs2, dist(rng));
    particle.k =
        initRandomDirectionWithRespectToCurrentK(particle.k, cosTheta, dist(rng));
    T kCurr = norm(particle.k);
    if (kCurr > T(0))
      particle.k = scale(particle.k, kF / kCurr);

    phononBath->recordAbsorption(norm(subtract(particle.k, kOld)));
  }
};

// =============================================================================
/** Screened 3D Fröhlich emission with non-equilibrium (hot) phonon bath. */
template <class T>
class emcScreenedHotPhononFroehlichEmission3D : public emcScatterMechanism<T> {
private:
  T phononEnergy, effMass, scatterConst;
  std::shared_ptr<emcPhononBath<T>> phononBath;
  std::shared_ptr<emcPlasmonScreening<T>> screening;
  bool qResolved;
  // Separate switch for the final-state angle. Setting this false while
  // qResolved is true reproduces the INCONSISTENT treatment discussed in the
  // text: the rate follows the resolved distribution while the angle is drawn
  // as if the occupation were uniform. Provided so that combination can be
  // measured rather than argued about; it is not physically correct.
  bool qResolvedAngle;
  std::string nameSuffix;
  mutable std::uniform_real_distribution<T> dist;

public:
  emcScreenedHotPhononFroehlichEmission3D() = delete;

  emcScreenedHotPhononFroehlichEmission3D(
      SizeType inValley, T inPhononEnergy, T relEffMass, T eps_hi, T eps_lo,
      std::shared_ptr<emcPhononBath<T>> inPhononBath,
      std::shared_ptr<emcPlasmonScreening<T>> inScreening,
      bool inQResolved = false, std::string inNameSuffix = "",
      bool inQResolvedAngle = true)
      : emcScatterMechanism<T>(inValley), phononEnergy(inPhononEnergy),
        effMass(relEffMass * constants::me), phononBath(std::move(inPhononBath)),
        screening(std::move(inScreening)), qResolved(inQResolved),
        qResolvedAngle(inQResolvedAngle), nameSuffix(inNameSuffix),
        dist(0., 1.) {
    scatterConst =
        froehlichScatterConst3D(phononEnergy, effMass, eps_hi, eps_lo);
  }

  std::string getName() const {
    return "ScreenedHotPhononFroehlichEmission3D" + nameSuffix;
  }

  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    if (energy <= phononEnergy)
      return T(0);
    auto valley = this->ptrValley[this->idxValley];
    T gammaI = valley->getGamma(energy);
    T gammaF = valley->getGamma(energy - phononEnergy);
    if (gammaI <= T(0) || gammaF <= T(0))
      return T(0);
    T kI = std::sqrt(T(2) * effMass * gammaI * constants::q) / constants::hbar;
    T kF = std::sqrt(T(2) * effMass * gammaF * constants::q) / constants::hbar;
    const T Nq = qResolved
                     ? phononBath->getNqInWindow(std::fabs(kI - kF), kI + kF)
                     : phononBath->getMeanNq();
    return scatterConst * (Nq + T(1)) / kI *
           screenedFroehlichLogFactor(kI, kF, screening->getQs2());
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    auto valley = this->ptrValley[this->idxValley];
    std::array<T, 3> kOld = particle.k;
    T kI = valley->getNormWaveVec(particle.energy);
    particle.energy -= phononEnergy;
    assert(particle.energy > T(0));
    T kF = valley->getNormWaveVec(particle.energy);

    const T qs2 = screening->getQs2();
    T cosTheta =
        (qResolved && qResolvedAngle)
            ? sampleQResolvedCosTheta(*phononBath, kI, kF, true, dist(rng))
            : sampleScreenedFroehlichCosTheta(kI, kF, qs2, dist(rng));
    particle.k =
        initRandomDirectionWithRespectToCurrentK(particle.k, cosTheta, dist(rng));
    T kCurr = norm(particle.k);
    if (kCurr > T(0))
      particle.k = scale(particle.k, kF / kCurr);

    phononBath->recordEmission(norm(subtract(particle.k, kOld)));
  }
};

#endif // EMC_SCREENED_FROEHLICH_INTERACTION_HPP
