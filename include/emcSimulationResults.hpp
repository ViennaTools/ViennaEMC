#ifndef EMC_SIMULATION_RESULTS_HPP
#define EMC_SIMULATION_RESULTS_HPP

#include <vector>

#include <emcConstants.hpp>
#include <emcGrid.hpp>
#include <emcOutput.hpp>
#include <emcSimulationParameter.hpp>
#include <emcUtil.hpp>

/**
 * @brief Class that holds the simulation results: current
 * at each contact, average potential and particle concentration.
 *
 * Current is calculated via the number of particles that leave and enter
 * the contacts at each step.
 *
 * @tparam T Numeric Type
 * @tparam DeviceType Type of the used device
 * @param param the parameter used in the current simulation
 * @param currPot current potential profile
 * @param avgPot average potential profile (tracked in the last steps)
 * @param currConc current particle concentration of each particle type
 * @param avgConc average particle concentration of each particle type
 * @param nettoPart counts the netto number of particles that leave the
 * device at each non-transient step (accessable through idxStep, idxType,
 * idxContact), (removed particles - injected particles)
 * @param current mean current that was measured at each non-transient step
 * (accessible through idxStep, idxType, idxContact)
 * @param nrPart number of particles that are currently assigned to each grid
 * point (handled by emcSimulation with PMScheme)
 * @param nrPartTypes number of used particle types
 * @param nrContacts number of contacts in the device
 * @param nrCurrSteps number of steps used to track the mean current
 * at each contact (nr. non-transient steps)
 */
template <class T, class DeviceType> class emcSimulationResults {
  static const SizeType Dim = DeviceType::Dimension;

  /// type that holds information of the current for each
  /// particle type (idx 0) at each contact (idx 1)
  typedef std::vector<std::vector<T>> CurrentMeasurement;
  /// type that counts the netto number of particles that
  /// leave / enter the device for each particle type (idx 0)
  /// at each contact (idx 1)
  typedef std::vector<std::vector<int>> NettoParticleCounter;
  typedef emcGrid<T, Dim> GridType;

  const SizeType nrPartTypes, nrContacts, nrCurrSteps;
  const emcSimulationParameter<T, DeviceType> &param;
  GridType currPot, avgPot;
  std::vector<GridType> currConc, avgConc, eField, nrPart;
  SizeType nrCummSteps{0}, nrAvgSteps{0}; /**< counter */
  std::vector<NettoParticleCounter> nettoPart;
  std::vector<CurrentMeasurement> current;
  CurrentMeasurement nettoPartSum;
  std::vector<T> currentFactor;

public:
  emcSimulationResults() = delete;

  emcSimulationResults(const DeviceType &device,
                       const emcSimulationParameter<T, DeviceType> &inParam)
      : nrPartTypes(inParam.getNrParticleTypes()),
        nrContacts(device.getSurface().getNrContacts()),
        nrCurrSteps(inParam.getNrNonTransientSteps()), param(inParam),
        currPot(device.getGridExtent()), avgPot(currPot),
        eField(Dim, GridType(currPot)), avgConc(nrPartTypes, GridType(currPot)),
        currConc(nrPartTypes, GridType(currPot)),
        nrPart(nrPartTypes, GridType(currPot)),
        nettoPart(nrCurrSteps, NettoParticleCounter(
                                   nrPartTypes, std::vector<int>(nrContacts))),
        current(nrCurrSteps,
                CurrentMeasurement(nrPartTypes, std::vector<T>(nrContacts))),
        nettoPartSum(nrPartTypes, std::vector<T>(nrContacts, 0)),
        currentFactor(nrPartTypes) {
    initPotential(device);
    for (const auto &[idxType, partType] : param.particleTypes)
      currentFactor[idxType] =
          param.nrCarriersPerPart * partType->getCharge() / param.stepTime;
  }

  /// @brief updates average characteristics of interest.
  /// Currently this includes updating the average potential and particle
  /// concentration of all included particle types.
  void updateAverageCharacteristics() {
    avgPot += currPot;
    std::transform(avgConc.begin(), avgConc.end(), currConc.begin(),
                   avgConc.begin(),
                   [](GridType &avg, GridType &curr) { return avg + curr; });
    nrAvgSteps++;
  }

  /// @brief updates the current particle conentrations with the currently
  /// stored nr of particles close to each grid point (nrPart)
  /// @param device used device
  void updateCurrentParticleConcentrations(const DeviceType &device) {
    std::array<SizeType, Dim> coord;
    auto extent = device.getGridExtent();
    T factor = 1. / device.getCellVolume();
    for (SizeType idxType = 0; idxType < nrPartTypes; idxType++) {
      coord.fill(0);
      std::transform(nrPart[idxType].begin(), nrPart[idxType].end(),
                     currConc[idxType].begin(), [&](const T &nr) {
                       T conc = device.normalizeDoping(nr) * factor;
                       for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
                         if (coord[idxDim] == 0 ||
                             coord[idxDim] == extent[idxDim] - 1)
                           conc *= 2;
                       }
                       device.advanceCoord(coord);
                       return conc;
                     });
    }
  }

  /// @brief stores the netto number of particles that left / entered
  /// device at this step and updates the average current with the
  /// new netto number of particles. Additionally, it counts the
  /// number of steps that were already recorded.
  /// @param nrRemPart number of removed particles (due to particles that left)
  /// @param nrInjPart number of injected (- deleted) particles (due to
  /// equilibrium conditions at ohmic contacts)
  void updateCurrent(const NettoParticleCounter &nrRemPart,
                     const NettoParticleCounter &nrInjPart) {
    // calculate netto number of particles that left / got injected in this step
    std::transform(nrRemPart.begin(), nrRemPart.end(), nrInjPart.begin(),
                   nettoPart[nrCummSteps].begin(), [](auto a, const auto &b) {
                     std::transform(a.begin(), a.end(), b.begin(), a.begin(),
                                    std::minus<T>());
                     return a;
                   });
    // add netto number of particles to cummulative sum of netto particles
    std::transform(nettoPartSum.begin(), nettoPartSum.end(),
                   nettoPart[nrCummSteps].begin(), nettoPartSum.begin(),
                   [](auto a, const auto &b) {
                     std::transform(a.begin(), a.end(), b.begin(), a.begin(),
                                    std::plus<T>());
                     return a;
                   });
    // calculate averaged current from cummulative sum of netto particles
    for (SizeType idxType = 0; idxType < nrPartTypes; idxType++)
      std::transform(nettoPartSum[idxType].begin(), nettoPartSum[idxType].end(),
                     current[nrCummSteps][idxType].begin(),
                     [this, &idxType](const T &sum) {
                       return sum / (nrCummSteps + 1) * currentFactor[idxType];
                     });
    nrCummSteps++;
  }

  /// @brief writes current results to files, includes current potential, el.
  /// Field and current particle concentration of every particle Type
  void writeCurrentResults(std::string nameSuffix,
                           const DeviceType &device) const {
    if (param.adaptPotentialForWrite)
      writeToFile(currPot, param.namePrefix + "Potential" + nameSuffix,
                  param.adaptPotentialForWrite, device);
    else
      writeToFile(currPot, param.namePrefix + "Potential" + nameSuffix,
                  undoNormalizationPotential, device);

    for (const auto &[idxType, partType] : param.particleTypes)
      writeToFile(currConc[idxType],
                  param.namePrefix + partType->getName() + "Conc" + nameSuffix,
                  undoNormalizationConcentration, device);

    writeToFile(eField[0], param.namePrefix + "EFieldX" + nameSuffix);
    writeToFile(eField[1], param.namePrefix + "EFieldY" + nameSuffix);
    if (Dim == 3)
      writeToFile(eField[2], param.namePrefix + "EFieldZ" + nameSuffix);
  }

  /// @brief write final characteristics of interest of device to
  /// files (includes: avg potential, avg particle concentration
  /// and avg current at each non-transient step)
  /// @param device used device
  void writeFinalResults(const DeviceType &device) const {
    // average the characteristics of interest
    auto calcAvg = [this](const T &val) { return val / nrAvgSteps; };
    GridType finalPot{avgPot.getExtent()};
    std::vector<GridType> finalConc{nrPartTypes, GridType{avgPot.getExtent()}};
    std::transform(avgPot.begin(), avgPot.end(), finalPot.begin(), calcAvg);
    for (SizeType idxType = 0; idxType < nrPartTypes; idxType++)
      std::transform(avgConc[idxType].begin(), avgConc[idxType].end(),
                     finalConc[idxType].begin(), calcAvg);

    // write final characteristics of interest to files
    for (const auto &[idxType, partType] : param.particleTypes) {
      std::string prefix = param.namePrefix + partType->getName();
      writeToFile(nettoPart, current, idxType, param.stepTime,
                  param.transientTime, prefix + "Current");
      writeToFile(finalConc[idxType], prefix + "ConcAvg",
                  undoNormalizationConcentration, device);
    }

    if (param.adaptPotentialForWrite)
      writeToFile(finalPot, param.namePrefix + "PotentialAvg",
                  param.adaptPotentialForWrite, device);
    else
      writeToFile(finalPot, param.namePrefix + "PotentialAvg",
                  undoNormalizationPotential, device);
  }

private:
  /// @brief initializes potential with doping potential in each region which is
  /// calculated with asinh(x) = ln(x + sqrt(x * x + 1))
  void initPotential(const DeviceType &device) {
    typename DeviceType::SizeVec coord;
    for (coord.fill(0); !currPot.isEndCoord(coord); currPot.advanceCoord(coord))
      currPot[coord] =
          std::asinh(0.5 * device.getDopingProfile().getDoping(coord, true));
  }

  template <class, class, class, class, class> friend class emcSimulation;
};

#endif // EMC_SIMULATION_RESULTS_HPP