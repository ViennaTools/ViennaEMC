#ifndef EMC_ZERO_ORDER_SINGLE_LAYER_INTERVALLEY_SCATTER_MECHANISM_HPP
#define EMC_ZERO_ORDER_SINGLE_LAYER_INTERVALLEY_SCATTER_MECHANISM_HPP

#include <map>

#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>

/*! \brief Zero Order Intervalley Absorption Scatter Mechanism for Single Layer
 * of Material.
 *
 * For single layer it is assumed that there is no movement in z-direction.
 * Formula for scattering from Paper of Kaasbjerg et al. with given link:
 * (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.115317)
 *
 * @param scatterConst pre-computed constant needed for calculation of
 * scatter rate
 * @param phononEnergy energy of the involved phonon [eV]
 * @param idxFinalValley index of the final valley
 * @param nrFinalValleys nr of possible final subvalleys
 * @param finalSubValleys vector that stores all possible final subvalleys for
 * all initial valleys
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcZeroOrderSingleLayerInterValleyAbsorptionScatterMechanism
    : public emcScatterMechanism<T> {
private:
  T scatterConst;
  T phononEnergy;
  SizeType idxFinalValley;
  SizeType nrFinalValleys;
  std::string nameSuffix;
  std::vector<std::vector<SizeType>> finalSubValleys;

  mutable std::uniform_real_distribution<T> dist;

public:
  emcZeroOrderSingleLayerInterValleyAbsorptionScatterMechanism() = delete;

  /// \brief Constructor for Kaasbjerg et al. comparison.
  /// Assumes that only one valley is used with one subvalley.
  emcZeroOrderSingleLayerInterValleyAbsorptionScatterMechanism(
      SizeType inValley, T sigma, T densityMaterial, T temperature,
      T inPhononEnergy, std::string inNameSuffix = "")
      : emcZeroOrderSingleLayerInterValleyAbsorptionScatterMechanism(
            inValley, inValley, sigma, densityMaterial, temperature,
            inPhononEnergy, std::vector<std::vector<SizeType>>(),
            inNameSuffix) {}

  /*! \brief Constructor.
   *
   * @param inValley index of of initial valley (valley to which scatter
   * mechanism is assigned)
   * @param inFinalValley index of final valley
   * @param sigma deformation potential [eV / m]
   * @param densityMaterial density of material [1 / m^2]
   * @param temperature temperature of material (device)
   * @param inPhononEnergy energy of the involved phonon [eV]
   * @param inFinalSubValleys assignment of scattering from initial subvalley to
   * final ones.
   * @param inNameSuffix suffix for filename
   */
  emcZeroOrderSingleLayerInterValleyAbsorptionScatterMechanism(
      SizeType inValley, SizeType inFinalValley, T sigma, T densityMaterial,
      T temperature, T inPhononEnergy,
      std::vector<std::vector<SizeType>> inFinalSubValleys,
      std::string inNameSuffix = "")
      : dist(0., 1.), idxFinalValley(inFinalValley), nameSuffix(inNameSuffix),
        phononEnergy(inPhononEnergy), emcScatterMechanism<T>(inValley),
        finalSubValleys(inFinalSubValleys) {

    if (finalSubValleys.empty())
      nrFinalValleys = 1; // for Kaasbjerg comparison
    else {
      nrFinalValleys = finalSubValleys[0].size();
    }

    T exponent = phononEnergy * constants::q / (constants::kB * temperature);
    T omega = phononEnergy * constants::q / constants::hbar;

    scatterConst = nrFinalValleys *
                   std::pow(sigma * constants::q / constants::hbar, 2) /
                   (2 * densityMaterial * omega * (std::exp(exponent) - 1));
  }

  std::string getName() const {
    return "ZeroInterValleyAbsorptionSL" + nameSuffix;
  }

  /// \brief scatter-rate calculation can be used for parabolic /
  /// nonparabolic and istrop / anistrop valleys.
  T getScatterRate(T energy, SizeType /*region*/) const {
    auto initValley = this->ptrValley[this->idxValley];
    auto finalValley = this->ptrValley[idxFinalValley];
    T deltaValley =
        finalValley->getBottomEnergy() - initValley->getBottomEnergy();
    T eFinal = energy + phononEnergy - deltaValley;
    if (eFinal > 0) {
      T md = finalValley->getEffMassDOS();
      T alpha = finalValley->getNonParabolicity();
      return md * scatterConst * (1 + 2 * alpha * eFinal);
    }
    return 0;
  }

  /// \brief Inelastic, Isotrop Intervalley or Intravalley Scattering.
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    // adapt valley and optionally choose random final subvalley
    particle.valley = idxFinalValley;
    if (!finalSubValleys.empty()) {
      particle.subValley =
          finalSubValleys[particle.subValley]
                         [std::floor(dist(rng) * nrFinalValleys)];
    }

    // adapt energy of particle
    auto initValley = this->ptrValley[this->idxValley];
    auto finalValley = this->ptrValley[idxFinalValley];
    T deltaValley =
        finalValley->getBottomEnergy() - initValley->getBottomEnergy();
    particle.energy += phononEnergy - deltaValley;
    assert(particle.energy > 0);

    T angle = 2 * constants::pi * dist(rng);
    auto vogtFactor = finalValley->getVogtTransformationFactor();

    // random direction weighted by vogtFactor (see Esseni p. 276)
    std::array<T, 3> k{std::cos(angle) / vogtFactor[0],
                       std::sin(angle) / vogtFactor[1], 0};
    T factor = 1. / (std::sqrt(k[0] * k[0] + k[1] * k[1]));
    k[0] *= factor;
    k[1] *= factor;

    particle.k =
        finalValley->transformToDeviceCoord(particle.subValley, particle.k);
    T normK = finalValley->getNormWaveVec(particle.energy);
    k[0] *= normK;
    k[1] *= normK;
    particle.k = k;
  }

  /// @brief checks if final valley and all subvalleys are valid.
  void check() final {
    std::string msg;
    if (idxFinalValley >= this->ptrValley.size()) {
      msg = getName() + ": idxFinalValley " + std::to_string(idxFinalValley) +
            " is not valid.";
      emcMessage::getInstance().addError(msg).print();
    }

    auto initValley = this->ptrValley[this->idxValley];
    auto finalValley = this->ptrValley[idxFinalValley];
    SizeType initDeg = initValley->getDegeneracyFactor();
    SizeType finalDeg = finalValley->getDegeneracyFactor();

    if (!finalSubValleys.empty()) {
      for (SizeType idxInitSub = 0; idxInitSub < initDeg; idxInitSub++) {
        if (finalSubValleys.at(idxInitSub).size() != nrFinalValleys) {
          msg = getName() + ": Nr. of final subvalleys not consistent.";
          emcMessage::getInstance().addWarning(msg).print();
        }
        for (const auto &idxFinalSub : finalSubValleys.at(idxInitSub)) {
          if (idxFinalSub >= finalDeg) {
            msg = getName() + ": Used idx " + std::to_string(idxFinalSub) +
                  " for valley of degeneracy " + std::to_string(finalDeg) +
                  " is not valid.";
            emcMessage::getInstance().addError(msg).print();
          }
        }
      }
    }
  }
};

/*! \brief Zero Order Intervalley Emission Scatter Mechanism for Single Layer
 * of Material.
 *
 * For single layer it is assumed that there is no movement in z-direction.
 * Formula for scattering from Paper of Kaasbjerg et al. with given link:
 * (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.115317)
 *
 * @param scatterConst pre-computed constant needed for calculation of
 * scatter rate
 * @param phononEnergy energy of the involved phonon [eV]
 * @param idxFinalValley index of the final valley
 * @param nrFinalValleys nr of possible final subvalleys
 * @param finalSubValleys vector that stores all possible final subvalleys for
 * all initial valleys
 * @param nameSuffix suffix for filename
 */
template <class T>
class emcZeroOrderSingleLayerInterValleyEmissionScatterMechanism
    : public emcScatterMechanism<T> {
private:
  T scatterConst;
  mutable std::uniform_real_distribution<T> dist;
  T phononEnergy;
  SizeType idxFinalValley;
  SizeType nrFinalValleys;
  std::string nameSuffix;
  std::vector<std::vector<SizeType>> finalSubValleys;

public:
  emcZeroOrderSingleLayerInterValleyEmissionScatterMechanism() = delete;

  /// \brief Constructor for Kaasbjerg et al. comparison.
  /// Assumes that only one valley is used with one subvalley.
  emcZeroOrderSingleLayerInterValleyEmissionScatterMechanism(
      SizeType inValley, T sigma, T densityMaterial, T temperature,
      T inPhononEnergy, std::string inNameSuffix = "")
      : emcZeroOrderSingleLayerInterValleyEmissionScatterMechanism(
            inValley, inValley, sigma, densityMaterial, temperature,
            inPhononEnergy, std::vector<std::vector<SizeType>>(),
            inNameSuffix) {}

  /*! \brief Constructor.
   *
   * @param inValley index of of initial valley (valley to which scatter
   * mechanism is assigned)
   * @param inFinalValley index of final valley
   * @param sigma deformation potential [eV / m]
   * @param densityMaterial density of material [1 / m^2]
   * @param temperature temperature of material (device)
   * @param inPhononEnergy energy of the involved phonon [eV]
   * @param inFinalSubValleys assignment of scattering from initial subvalley to
   * final ones.
   * @param inNameSuffix suffix for filename
   */
  emcZeroOrderSingleLayerInterValleyEmissionScatterMechanism(
      SizeType inValley, SizeType inFinalValley, T sigma, T densityMaterial,
      T temperature, T inPhononEnergy,
      std::vector<std::vector<SizeType>> inFinalSubValleys,
      std::string inNameSuffix = "")
      : dist(0., 1.), idxFinalValley(inFinalValley), nameSuffix(inNameSuffix),
        phononEnergy(inPhononEnergy), emcScatterMechanism<T>(inValley),
        finalSubValleys(inFinalSubValleys) {
    if (finalSubValleys.empty())
      nrFinalValleys = 1; // for Kaasbjerg comparison
    else {
      nrFinalValleys = finalSubValleys[0].size();
    }

    T exponent = (phononEnergy * constants::q) / (constants::kB * temperature);
    T omega = (phononEnergy * constants::q) / constants::hbar;

    scatterConst = nrFinalValleys *
                   std::pow(sigma * constants::q / constants::hbar, 2) *
                   std::exp(exponent) /
                   (2 * densityMaterial * omega * (std::exp(exponent) - 1));
  }

  std::string getName() const {
    return "ZeroInterValleyEmissionSL" + nameSuffix;
  }

  /// \brief scatter-rate calculation can be used for parabolic /
  /// nonparabolic and istrop / anistrop valleys.
  T getScatterRate(T energy, SizeType /*region*/) const {
    auto initValley = this->ptrValley[this->idxValley];
    auto finalValley = this->ptrValley[idxFinalValley];
    T deltaValley =
        finalValley->getBottomEnergy() - initValley->getBottomEnergy();
    T eFinal = energy - phononEnergy - deltaValley;
    if (eFinal > 0) {
      T md = finalValley->getEffMassDOS();
      T alpha = finalValley->getNonParabolicity();
      return md * scatterConst * (1 + 2 * alpha * eFinal);
    }
    return 0;
  }

  /// \brief Inelastic, Isotrop Intervalley or Intravalley Scattering.
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    particle.valley = idxFinalValley;
    if (!finalSubValleys.empty()) {
      particle.subValley =
          finalSubValleys[particle.subValley]
                         [std::floor(dist(rng) * nrFinalValleys)];
    }

    // adapt energy of particle
    auto initValley = this->ptrValley[this->idxValley];
    auto finalValley = this->ptrValley[idxFinalValley];
    T deltaValley =
        finalValley->getBottomEnergy() - initValley->getBottomEnergy();
    particle.energy -= (deltaValley + phononEnergy);
    assert(particle.energy >= 0);

    T angle = 2 * constants::pi * dist(rng);
    auto vogtFactor = finalValley->getVogtTransformationFactor();

    // random direction weighted by vogtFactor (see Esseni p. 276)
    std::array<T, 3> k{std::cos(angle) / vogtFactor[0],
                       std::sin(angle) / vogtFactor[1], 0};
    T factor = 1. / (std::sqrt(k[0] * k[0] + k[1] * k[1]));
    k[0] *= factor;
    k[1] *= factor;

    particle.k =
        finalValley->transformToDeviceCoord(particle.subValley, particle.k);
    T normK = finalValley->getNormWaveVec(particle.energy);
    k[0] *= normK;
    k[1] *= normK;
    particle.k = k;
  }

  /// @brief checks if final valley and all subvalleys are valid.
  void check() final {
    std::string msg;
    if (idxFinalValley >= this->ptrValley.size()) {
      msg = getName() + ": idxFinalValley " + std::to_string(idxFinalValley) +
            " is not valid.";
      emcMessage::getInstance().addError(msg).print();
    }

    auto initValley = this->ptrValley[this->idxValley];
    auto finalValley = this->ptrValley[idxFinalValley];
    SizeType initDeg = initValley->getDegeneracyFactor();
    SizeType finalDeg = finalValley->getDegeneracyFactor();

    if (!finalSubValleys.empty()) {
      for (SizeType idxInitSub = 0; idxInitSub < initDeg; idxInitSub++) {
        if (finalSubValleys.at(idxInitSub).size() != nrFinalValleys) {
          msg = getName() + ": Nr. of final subvalleys not consistent.";
          emcMessage::getInstance().addWarning(msg).print();
        }
        for (const auto &idxFinalSub : finalSubValleys.at(idxInitSub)) {
          if (idxFinalSub >= finalDeg) {
            msg = getName() + ": Used idx " + std::to_string(idxFinalSub) +
                  " for valley of degeneracy " + std::to_string(finalDeg) +
                  " is not valid.";
            emcMessage::getInstance().addError(msg).print();
          }
        }
      }
    }
  }
};

#endif // EMC_ZERO_ORDER_SINGLE_LAYER_INTERVALLEY_SCATTER_MECHANISM_HPP