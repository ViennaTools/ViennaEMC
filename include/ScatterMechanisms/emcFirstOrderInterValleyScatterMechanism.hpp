#ifndef EMC_FIRST_ORDER_INTERVALLEY_SCATTER_MECHANISM_HPP
#define EMC_FIRST_ORDER_INTERVALLEY_SCATTER_MECHANISM_HPP

#include <map>

#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <emcConstants.hpp>

// helper function that calculates the scatter constant needed for
// first order intervalley scattering
template <class T>
T getFirstOrderScatterConst(T defPot, T phEnergy, T rho, T temp,
                            SizeType nrFValleys, bool isAbsorption = true) {
  T result = nrFValleys * std::sqrt(2) * pow(constants::q, 5. / 2.) *
             pow(defPot, 2) / (M_PI * rho * pow(constants::hbar, 4) * phEnergy);
  T nrPh =
      1. / (std::exp(constants::q * phEnergy / (constants::kB * temp)) - 1.);
  result = isAbsorption ? result * nrPh : result * (nrPh + 1);
  return result;
}

/*! \brief First Order Optical Intervalley Absorption Scatter Mechanism.
 *
 * Assumes that mechanism is inelastic and isotropical.
 *
 * @param nameSuffix suffix for filename
 * @param phononEnergy energy of the involved phonon [eV]
 * @param scatterConst pre-computed constant needed for calculation of
 * scatter rate
 * @param nrFValleys nr of possible final subvalleys
 * @param idxFinalValley index of the final valley
 * @param finalSubValleys vector that stores all possible final subvalleys for
 * all initial valleys
 * @param deltaValley energy difference between initial and final valley
 */
template <class T>
class emcFirstOrderInterValleyAbsorptionScatterMechanism
    : public emcScatterMechanism<T> {
protected:
  std::string nameSuffix;
  T phononEnergy;
  T scatterConst;
  SizeType nrFValleys;
  SizeType idxFinalValley;
  std::map<SizeType, std::vector<SizeType>> finalSubValleys;
  mutable std::uniform_real_distribution<T> dist;
  mutable T deltaValley;

public:
  emcFirstOrderInterValleyAbsorptionScatterMechanism() = delete;

  /*! \brief Constructor in case final and initial valley are the same.
   *
   * @param inNameSuffix suffix for filename
   * @param inIdxValley index of of initial valley (valley to which scatter
   * mechanism is assigned)
   * @param inFinalSubValleys assignment of scattering from initial subvalley to
   * final ones.
   * @param defPotential deformation potential [eV / m]
   * @param inPhononEnergy energy of the involved phonon [eV]
   * @param device used device
   */
  template <class DeviceType>
  emcFirstOrderInterValleyAbsorptionScatterMechanism(
      std::string inNameSuffix, SizeType inIdxValley,
      std::map<SizeType, std::vector<SizeType>> inFinalSubValleys,
      T defPotential, T inPhononEnergy, DeviceType &device)
      : emcFirstOrderInterValleyAbsorptionScatterMechanism(
            inNameSuffix, inIdxValley, inIdxValley, inFinalSubValleys,
            defPotential, inPhononEnergy, device) {}

  /*! \brief Constructor in case final and initial valley are different.
   *
   * @param inNameSuffix suffix for filename
   * @param inIdxValley index of of initial valley (valley to which scatter
   * mechanism is assigned)
   * @param inIdxFinalValley index of final valley
   * @param inFinalSubValleys assignment of scattering from initial subvalley to
   * final ones.
   * @param defPotential deformation potential [eV / m]
   * @param inPhononEnergy energy of the involved phonon [eV]
   * @param device used device
   */
  template <class DeviceType>
  emcFirstOrderInterValleyAbsorptionScatterMechanism(
      std::string inNameSuffix, SizeType inIdxValley, SizeType inIdxFinalValley,
      std::map<SizeType, std::vector<SizeType>> inFinalSubValleys,
      T defPotential, T inPhononEnergy, const DeviceType &device)
      : emcScatterMechanism<T>(inIdxValley), idxFinalValley(inIdxFinalValley),
        dist(0., 1.), nameSuffix(inNameSuffix), phononEnergy(inPhononEnergy),
        finalSubValleys(inFinalSubValleys), deltaValley(0),
        nrFValleys(inFinalSubValleys.at(0).size()) {
    scatterConst = getFirstOrderScatterConst(
        defPotential, phononEnergy, device.getMaterial().getRho(),
        device.getTemperature(), nrFValleys, true);
  }

  std::string getName() const {
    return "FirstInterValleyAbsorption" + nameSuffix;
  }

  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    auto initValley = this->ptrValley[this->idxValley];
    auto finalValley = this->ptrValley[idxFinalValley];
    deltaValley =
        finalValley->getBottomEnergy() - initValley->getBottomEnergy();
    T energyFinal = energy + phononEnergy - deltaValley;
    if (energyFinal > 0) {
      T md = finalValley->getEffMassDOS();
      T alpha = finalValley->getNonParabolicity();
      T gamma = finalValley->getGamma(energy);
      T gammaF = finalValley->getGamma(energyFinal);
      return scatterConst * pow(md, 5. / 2.) * std::sqrt(gammaF) *
             (2 * alpha * energyFinal + 1.0) * (gamma + gammaF);
    } else
      return 0;
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    // adapt valley + subvalley
    particle.valley = idxFinalValley;
    particle.subValley =
        finalSubValleys.at(particle.subValley)[rng() % nrFValleys];

    // adapt energy + wave-vector
    particle.energy += (phononEnergy - deltaValley);
    T knew = this->ptrValley[idxFinalValley]->getNormWaveVec(particle.energy);
    particle.k = initRandomDirection(knew, dist(rng), dist(rng));
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

    for (SizeType idxInitSub = 0; idxInitSub < initDeg; idxInitSub++) {
      if (finalSubValleys.find(idxInitSub) == finalSubValleys.end()) {
        msg = getName() +
              ": No final subvalleys given for initial subvalley index " +
              std::to_string(idxInitSub);
        emcMessage::getInstance().addError(msg).print();
      }
      if (finalSubValleys.at(idxInitSub).size() != nrFValleys) {
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
};

/*! \brief First Order Optical Intervalley Emission Scatter Mechanism.
 *
 * Assumes that mechanism is inelastic and isotropical.
 *
 * @param nameSuffix suffix for filename
 * @param phononEnergy energy of the involved phonon [eV]
 * @param scatterConst pre-computed constant needed for calculation of
 * scatter rate
 * @param nrFValleys nr of possible final subvalleys
 * @param idxFinalValley index of the final valley
 * @param finalSubValleys vector that stores all possible final subvalleys for
 * all initial valleys
 * @param deltaValley energy difference between initial and final valley
 */
template <class T>
class emcFirstOrderInterValleyEmissionScatterMechanism
    : public emcScatterMechanism<T> {
protected:
  std::string nameSuffix;
  T phononEnergy;
  T scatterConst;
  SizeType nrFValleys;
  SizeType idxFinalValley;
  std::map<SizeType, std::vector<SizeType>> finalSubValleys;
  mutable std::uniform_real_distribution<T> dist;
  mutable T deltaValley;

public:
  emcFirstOrderInterValleyEmissionScatterMechanism() = delete;

  /*! \brief Constructor in case final and initial valley are the same.
   *
   * @param inNameSuffix suffix for filename
   * @param inIdxValley index of of initial valley (valley to which scatter
   * mechanism is assigned)
   * @param inFinalSubValleys assignment of scattering from initial subvalley to
   * final ones.
   * @param defPotential deformation potential [eV / m]
   * @param inPhononEnergy energy of the involved phonon [eV]
   * @param device used device
   */
  template <class DeviceType>
  emcFirstOrderInterValleyEmissionScatterMechanism(
      std::string inNameSuffix, SizeType inIdxValley,
      std::map<SizeType, std::vector<SizeType>> inFinalSubValleys,
      T defPotential, T inPhononEnergy, const DeviceType &device)
      : emcFirstOrderInterValleyEmissionScatterMechanism(
            inNameSuffix, inIdxValley, inIdxValley, inFinalSubValleys,
            defPotential, inPhononEnergy, device) {}

  /*! \brief Constructor in case final and initial valley are different.
   *
   * @param inNameSuffix suffix for filename
   * @param inIdxValley index of of initial valley (valley to which scatter
   * mechanism is assigned)
   * @param inIdxFinalValley index of final valley
   * @param inFinalSubValleys assignment of scattering from initial subvalley to
   * final ones.
   * @param defPotential deformation potential [eV / m]
   * @param inPhononEnergy energy of the involved phonon [eV]
   * @param device used device
   */
  template <class DeviceType>
  emcFirstOrderInterValleyEmissionScatterMechanism(
      std::string inNameSuffix, SizeType inIdxValley, SizeType inIdxFinalValley,
      std::map<SizeType, std::vector<SizeType>> inFinalSubValleys,
      T defPotential, T inPhononEnergy, const DeviceType &device)
      : emcScatterMechanism<T>(inIdxValley), idxFinalValley(inIdxFinalValley),
        dist(0., 1.), nameSuffix(inNameSuffix), phononEnergy(inPhononEnergy),
        finalSubValleys(inFinalSubValleys), deltaValley(0),
        nrFValleys(inFinalSubValleys.at(0).size()) {
    scatterConst = getFirstOrderScatterConst(
        defPotential, phononEnergy, device.getMaterial().getRho(),
        device.getTemperature(), nrFValleys, false);
  }

  std::string getName() const {
    return "FirstInterValleyEmission" + nameSuffix;
  }

  T getScatterRate(T energy, SizeType /*idxRegion*/) const {
    auto initValley = this->ptrValley[this->idxValley];
    auto finalValley = this->ptrValley[idxFinalValley];
    deltaValley =
        finalValley->getBottomEnergy() - initValley->getBottomEnergy();
    T energyFinal = energy - phononEnergy - deltaValley;
    if (energyFinal > 0) {
      T md = finalValley->getEffMassDOS();
      T alpha = finalValley->getNonParabolicity();
      T gamma = finalValley->getGamma(energy);
      T gammaF = finalValley->getGamma(energyFinal);
      return scatterConst * pow(md, 5. / 2.) * std::sqrt(gammaF) *
             (2 * alpha * energyFinal + 1.0) * (gamma + gammaF);
    } else
      return 0;
  }

  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    // adapt valley + subvalley
    particle.valley = idxFinalValley;
    particle.subValley =
        finalSubValleys.at(particle.subValley)[rng() % nrFValleys];

    // adapt energy + wave-vector
    particle.energy -= (phononEnergy + deltaValley);
    T knew = this->ptrValley[idxFinalValley]->getNormWaveVec(particle.energy);
    particle.k = initRandomDirection(knew, dist(rng), dist(rng));
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

    for (SizeType idxInitSub = 0; idxInitSub < initDeg; idxInitSub++) {
      if (finalSubValleys.find(idxInitSub) == finalSubValleys.end()) {
        msg = getName() +
              ": No final subvalleys given for initial subvalley index " +
              std::to_string(idxInitSub);
        emcMessage::getInstance().addError(msg).print();
      }
      if (finalSubValleys.at(idxInitSub).size() != nrFValleys) {
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
};

#endif // EMC_FIRST_ORDER_INTERVALLEY_SCATTER_MECHANISM_HPP