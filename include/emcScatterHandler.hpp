#ifndef EMC_SCATTER_HANDLER_HPP
#define EMC_SCATTER_HANDLER_HPP

#include <algorithm>
#include <map>
#include <memory>
#include <tuple>

#include <ScatterMechanisms/emcScatterMechanism.hpp>
#include <SurfaceScatterMechanisms/emcSurfaceScatterMechanism.hpp>
#include <emcGrainScatterMechanism.hpp>
#include <emcUtil.hpp>

/*! \brief Class that handles all types of types of scattering, normal
 * scattering mechanisms, grain scattering and surface scattering.
 *
 * Usual Scattering: Class pre-calculates (+ normalizes) scatter
 * tables for each scatter mechanism, it is assumed that the rates of the
 * tables depend on the current valley and region of the particle. That
 * is the reason why scatterTables are calculated individually for each
 * (idxValley, idxRegion)-tuple.
 *
 * Grain Scattering: Class holds the added grain scatter mechanism and
 * has functionality to get the remaining free flight time until next
 * grain scatter event and scatter the particle at the grain.
 *
 * Surface Scattering: Class has the added surface scatter mechanisms
 * for each boundary position. If no surface mechanism is added at
 * the given boundary, specular scattering is performed.
 *
 * @param nrEnergyLevels nr. of energy levels (used in each table)
 * @param maxEnergy maximum energy of particle in the tables [eV]
 * @param tau maximal free flight time per region and valley
 * @param scatterMechanisms contains all added scatterMechanisms
 * @param scatterTables contains scatterTables of each (idxValley,
 * idxRegion)-combination
 * @param idxTableToidxMech relates the index of the table (idx of
 * scatterTables) to the index of the scatterMechanism
 * @param grainScatterMech unique pointer to grain scatter mechanism
 * @param grainTau remaining free flight time until next grain scatter
 * event
 * @param surfaceScatterMechanisms pointer to surface scatter mechanism
 * for each boundary position.
 */

template <class T, class DeviceType, SizeType Dim = DeviceType::Dimension>
class emcScatterHandler {
  typedef emcScatterMechanism<T> ScatterMechanism;
  typedef std::vector<T> ScatterTable;
  typedef std::tuple<SizeType, SizeType> IdxTupleType;

private:
  SizeType nrEnergyLevels;
  T maxEnergy;

  std::map<IdxTupleType, T> tau;
  std::vector<std::unique_ptr<ScatterMechanism>> scatterMechanisms;
  std::map<IdxTupleType, std::vector<ScatterTable>> scatterTables;
  std::map<IdxTupleType, std::vector<SizeType>> idxTableToIdxMech;

  std::unique_ptr<emcGrainScatterMechanism<T>> grainScatterMech = nullptr;
  T grainTau = 1.;

  std::vector<std::unique_ptr<emcSurfaceScatterMechanism<T, DeviceType>>>
      surfaceScatterMechanisms;

  mutable std::uniform_real_distribution<T> dist;
  T dE; // energy step inbetween different levels (in eV)

public:
  emcScatterHandler() : emcScatterHandler(1000, 4.) {}

  emcScatterHandler(SizeType inNrEnergyLevels, T inMaxEnergy)
      : nrEnergyLevels(inNrEnergyLevels), maxEnergy(inMaxEnergy), dist(0., 1.),
        dE(maxEnergy / nrEnergyLevels), surfaceScatterMechanisms(2 * Dim) {}

  /// return tau for specific for (idxValley, idxRegion)-tuple (if no mechanism
  /// is added in that region returns default value)
  T getTau(SizeType idxRegion, SizeType idxValley) const {
    IdxTupleType currIdx = std::make_tuple(idxValley, idxRegion);
    if (tau.find(currIdx) != tau.end())
      return tau.at(currIdx);
    return 2e-15;
  }

  /// return grain tau (if no mechanism is added in that region returns default
  /// value)
  T getGrainTau() const { return grainTau; }

  /// initializes + normalizes scatterTables
  void initScatterTables() {
    for (auto &[idxValleyRegion, tablesOfValleyRegion] : scatterTables) {
      for (SizeType idxTable = 0; idxTable < tablesOfValleyRegion.size();
           idxTable++) {
        SizeType idxMech = idxTableToIdxMech[idxValleyRegion][idxTable];
        auto &currTable = tablesOfValleyRegion[idxTable];
        for (SizeType energyLevel = 0; energyLevel < nrEnergyLevels;
             energyLevel++) {
          currTable.push_back(scatterMechanisms[idxMech]->getScatterRate(
              (energyLevel + 1) * dE, std::get<1>(idxValleyRegion)));
        }
      }
    }

    if (grainScatterMech)
      grainTau = 1. / grainScatterMech->getScatterRate();

    writeTablesToFiles();
    renormalizeTables();
  }

  void setGrainScatterMechanism(
      std::unique_ptr<emcGrainScatterMechanism<T>> &&newMechanism) {
    grainScatterMech.swap(newMechanism);
  }

  template <class DerivedSurfaceScatterMechanism>
  typename std::enable_if<
      std::is_base_of<emcSurfaceScatterMechanism<T, DeviceType>,
                      DerivedSurfaceScatterMechanism>::value>::type
  setSurfaceScatterMechanism(
      std::unique_ptr<DerivedSurfaceScatterMechanism> &&newMechanism,
      emcBoundaryPos boundaryPosition) {
    auto idxSurf = toUnderlying(boundaryPosition);
    newMechanism->setBoundaryPosition(boundaryPosition);
    surfaceScatterMechanisms[idxSurf].reset(newMechanism.release());
  }

  /// add a mechanism to a given valley and all given regions
  template <class DerivedScatterMechanism>
  typename std::enable_if<
      std::is_base_of<ScatterMechanism, DerivedScatterMechanism>::value>::type
  addScatterMechanism(std::unique_ptr<DerivedScatterMechanism> &&mechanism,
                      const std::vector<int> &regions) {
    auto idxValley = mechanism->getIdxValley();
    auto idxMech = scatterMechanisms.size();

    // add scatter mechanism to vector of mechanisms
    scatterMechanisms.push_back(std::move(mechanism));

    // add mapping of scatter table to mechanism
    for (const int idxRegion : regions) {
      IdxTupleType idxValleyRegion = std::make_tuple(idxValley, idxRegion);
      scatterTables[idxValleyRegion].push_back(ScatterTable());
      idxTableToIdxMech[idxValleyRegion].push_back(idxMech);
    }
  }

  /// select scatter mechanism and scatter particle
  void scatterParticle(emcParticle<T> &particle, emcRNG &rng) const {
    IdxTupleType idxValleyRegion =
        std::make_tuple(particle.valley, particle.region);
    auto &tables = scatterTables.at(idxValleyRegion);
    auto &idxMechs = idxTableToIdxMech.at(idxValleyRegion);

    if (!tables.empty()) {
      SizeType energyLevel = getEnergyLevel(particle.energy);
      T rand = dist(rng);
      if (rand > tables.back()[energyLevel])
        return; // self-scattering

      T boundLower{0}, boundUpper;
      for (int idxTable = 0; idxTable < tables.size(); idxTable++) {
        boundUpper = tables[idxTable][energyLevel];
        if (rand >= boundLower && rand < boundUpper) {
          scatterMechanisms[idxMechs[idxTable]]->scatterParticle(particle, rng);
          break;
        }
        boundLower = boundUpper;
      }
    }
  }

  void scatterParticleAtBoundary(emcBoundaryPos boundary,
                                 emcParticle<T> &particle,
                                 const DeviceType &device,
                                 std::array<T, Dim> &pos, emcRNG &rng) const {
    auto idxSurf = toUnderlying(boundary);
    // specular scattering in case no surface scatter mechanism is added
    if (!surfaceScatterMechanisms[idxSurf]) {
      auto maxPos = device.getMaxPos();
      for (SizeType idxDim = 0; idxDim < Dim; idxDim++) {
        if (pos[idxDim] < 0) {
          pos[idxDim] = -pos[idxDim];
          particle.k[idxDim] = -particle.k[idxDim];
        } else if (pos[idxDim] > maxPos[idxDim]) {
          pos[idxDim] = 2 * maxPos[idxDim] - pos[idxDim];
          particle.k[idxDim] = -particle.k[idxDim];
        }
      }
    } else
      surfaceScatterMechanisms[idxSurf]->scatterParticle(particle, pos, rng);
  }

  void scatterParticleAtGrain(emcParticle<T> &particle, emcRNG &rng) const {
    if (grainScatterMech)
      grainScatterMech->scatterParticle(particle, rng);
  }

  void writeTablesToFiles() const {
    std::ofstream os;
    T currEnergy;
    for (auto &[idxValleyRegion, tablesOfValleyRegion] : scatterTables) {
      for (SizeType idxTable = 0; idxTable < tablesOfValleyRegion.size();
           idxTable++) {
        SizeType idxMech = idxTableToIdxMech.at(idxValleyRegion)[idxTable];
        os.open(scatterMechanisms[idxMech]->getName() +
                std::to_string(std::get<1>(idxValleyRegion)) +
                std::to_string(std::get<0>(idxValleyRegion)) +
                "ScatterMechanism.txt");
        currEnergy = dE;
        for (auto &val : tablesOfValleyRegion[idxTable]) {
          os << currEnergy << " " << val << "\n";
          currEnergy += dE;
        }
        os.close();
      }
    }
  }

private:
  SizeType getEnergyLevel(T energy) const {
    SizeType energyLevel = std::floor(energy / dE) - 1;
    if (energyLevel == -1)
      return 0;
    if (energyLevel > nrEnergyLevels - 1)
      return nrEnergyLevels - 1;
    return energyLevel;
  }

  /// helper that normalizes scatterTables and calculates tau for each (valley,
  /// region)-combination
  void renormalizeTables() {
    for (auto &[idxValleyRegion, tables] : scatterTables) {
      SizeType nrOfTables = tables.size();

      if (nrOfTables == 0) // default - free flight time
        tau[idxValleyRegion] = 2e-15;
      else {
        for (SizeType i = 1; i < nrOfTables; i++) {
          // cumm. sum of scatter tables
          std::transform(tables[i].begin(), tables[i].end(),
                         tables[i - 1].begin(), tables[i].begin(),
                         std::plus<T>());
        }

        // find maximal cummulative rate
        T maxCummRate = *std::max_element(tables[nrOfTables - 1].begin(),
                                          tables[nrOfTables - 1].end());

        for (auto &table : tables) {
          std::for_each(table.begin(), table.end(),
                        [&maxCummRate](T &val) { val /= maxCummRate; });
        }

        tau[idxValleyRegion] = 1. / maxCummRate;
      }
    }

    std::cout << "Initialized ScatterHandler ...\n";
    for (const auto &[idxValleyRegion, currTau] : tau) {
      std::cout << "\tidxValley " << std::get<0>(idxValleyRegion);
      std::cout << " idxRegion " << std::get<1>(idxValleyRegion);
      std::cout << ": tau = " << currTau << " s\n";

      // std::cout << "\t  included following mechanisms:\n";
      // for (const auto &idxMech : idxTableToIdxMech[idxValleyRegion]) {
      //   std::cout << "\t\t" << scatterMechanisms[idxMech]->getName() << "\n";
      // }
    }
    if (grainScatterMech)
      std::cout << "\tGrain Scattering tau: " << grainTau << " s\n";
  }
};

#endif