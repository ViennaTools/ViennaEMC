#ifndef BASIC_BULK_PARTICLE_HANDLER_HPP
#define BASIC_BULK_PARTICLE_HANDLER_HPP

#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <ParticleType/emcParticleType.hpp>
#include <emcGrid.hpp>
#include <emcParticleDrift.hpp>
#include <emcParticleInitialization.hpp>
#include <emcUtil.hpp>

/**
 * @brief Handler is used for basic bulk simulations.
 *
 * Doesn't use FMM and assumes periodic boundary conditions
 * at the end of the simulation space. Also the device type in
 * this context is only used to represent the simulation space
 * and the doping within that space.
 *
 * @param device simulated device (in this case only extent of device, doping
 * and spacing are needed)
 * @param idxTypeToPartType map relating the idx of a particleType to the
 * specific particleType
 * @param appliedFieldDir direction of the applied field (normalized)
 * @param appliedField applied electric field
 * @param particles container storing the information of each particle
 * (accessible by index of particle type [idxType] and particle index [idxPart])
 * @param positionsParticles container storing the particle position of each
 * particle (accessible by index of particle type [idxType] and particle index
 * [idxPart])
 * @param rngs vector containing random number generator(s)
 */
template <class T, class DeviceType, SizeType Dim = DeviceType::Dimension>
struct basicBulkParticleHandler {
  typedef emcParticleType<T, DeviceType> ParticleType;
  typedef typename DeviceType::SizeVec SizeVec;
  typedef typename DeviceType::ValueVec ValueVec;
  typedef std::map<SizeType, std::unique_ptr<ParticleType>>
      MapIdxToParticleTypes;

private:
  DeviceType &device; //!< simulated geometry
  MapIdxToParticleTypes
      &idxTypeToPartType;   //!< map coupling idxType and PartType
  ValueVec appliedFieldDir; //!< direction of the el. field (normalized)
  ValueVec appliedField;    //!< strength of electric field

  std::vector<std::vector<emcParticle<T>>> particles;
  std::vector<std::vector<std::array<T, Dim>>> positionsParticles;

  std::uniform_real_distribution<T> distForLog{1e-6, 1.}, dist{0., 1.};
  std::vector<emcRNG> rngs;

public:
  basicBulkParticleHandler() = delete;

  /**
   * @brief Constructs a new object of basicBulkParticleHandler without any
   * applied electric field.
   *
   * @param inDevice simulated device (in this case only extent of device,
   * doping and spacing are required to represent the simulation space).
   * @param inTypes used particle types for the simulation.
   * @param inFieldDirection the direction of the applied field (in case the
   * strength is reset, the field is applied in that direction)
   */
  basicBulkParticleHandler(DeviceType &inDevice, MapIdxToParticleTypes &inTypes,
                           const ValueVec &inFieldDirection)
      : basicBulkParticleHandler(inDevice, inTypes, inFieldDirection, 0) {}

  /**
   * @brief Construct a new object of basicBulkParticleHandler, with the given
   * applied field.
   *
   * @param inDevice simulated device (in this case only extent of device,
   * doping and spacing are required to represent the simulation space).
   * @param inTypes particle types for the simulation.
   * @param inFieldDirection the direction of the applied field
   * @param inFieldStrength field strength
   */
  basicBulkParticleHandler(DeviceType &inDevice, MapIdxToParticleTypes &inTypes,
                           const ValueVec &inFieldDirection, T inFieldStrength)
      : device(inDevice), idxTypeToPartType(inTypes),
        appliedFieldDir(inFieldDirection), particles(inTypes.size()),
        positionsParticles(inTypes.size()) {
    // set electric field
    normalize(appliedFieldDir);
    appliedField = scale(appliedFieldDir, inFieldStrength);

    // initialize scatter tables
    for (const auto &partType : idxTypeToPartType) {
      if (partType.second->isMoved())
        partType.second->initScatterTables();
    }

    // seed random number generator(s)
    rngs.emplace_back(emcRNG(static_cast<long unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count())));
#ifdef _OPENMP
    for (SizeType idxSeed = 1; idxSeed < omp_get_max_threads(); idxSeed++)
      rngs.emplace_back(emcRNG(rngs[0]()));
#endif
  }

  /// @brief function that allows the seeding of the used random number
  /// generator(s).
  /// @param inSeed seed for the rng.
  void setSeed(SizeType inSeed) {
    rngs[0].seed(inSeed);
#ifdef _OPENMP
    for (SizeType idxSeed = 1; idxSeed < omp_get_max_threads(); idxSeed++)
      rngs[idxSeed].seed(rngs[0]());
#endif
  }

  /// @brief resets the strength of the applied field
  /// @param inAppliedFieldStrength new strength of the applied el. field
  void resetAppliedFieldStrength(T inAppliedFieldStrength) {
    appliedField = scale(appliedFieldDir, inAppliedFieldStrength);
  }

  /// @brief generates initial nr of Particles at each coord
  /// NOTE: the number of particles that are used depends on
  /// the used device (that represents the simulation space)
  /// and the used particleType!
  void generateInitialParticles() {
    emcGrid<T, Dim> pot(device.getGridExtent(), 0);
    for (const auto &[idxType, partType] : idxTypeToPartType) {
      SizeVec coord;
      for (coord.fill(0); !device.isEndCoord(coord);
           device.advanceCoord(coord)) {
        auto nrPartToCreate =
            partType->getInitialNrParticles(coord, device, pot);
        while (nrPartToCreate >= 1) {
          addParticle(idxType, coord, rngs[0]);
          nrPartToCreate--;
        }
        if (dist(rngs[0]) < nrPartToCreate)
          addParticle(idxType, coord, rngs[0]);
      }
    }
  }

  SizeType getNrParticles(SizeType idxType) const {
    return positionsParticles[idxType].size();
  }

  ValueVec getAppliedField() const { return appliedField; }

  void printNrParticles() const {
    for (const auto &[idxType, partType] : idxTypeToPartType)
      std::cout << "\t" << positionsParticles[idxType].size() << " "
                << partType->getName() << "\n";
  }

  /**
   * @brief moves particles for a specific amount of time.
   *
   * Motion is based on alternating free-flights and scatter events
   * based on the remaining free-flight time that a particle has left.
   *
   * @param tStep time-step for which the particle moves
   */
  void moveParticles(T tStep) {
    for (const auto &partType : idxTypeToPartType) {
      SizeType idxType = partType.first;
      auto &type = partType.second;
      if (type->isMoved()) {
        auto force = scale(appliedField, type->getCharge());
#pragma omp parallel
        {
          SizeType idxThread = 0;
#ifdef _OPENMP
          idxThread = omp_get_thread_num();
#endif
          auto &currRNG = rngs[idxThread];
#pragma omp for
          for (SizeType idxPart = 0; idxPart < getNrParticles(idxType);
               idxPart++) {
            auto &particle = particles[idxType][idxPart];
            auto &pos = positionsParticles[idxType][idxPart];
            auto valley = type->getValley(particle.valley);
            driftParticle(std::min(particle.tau, tStep), particle, valley, pos,
                          force);
            T tRemaining = tStep - particle.tau;
            while (tRemaining > 0) {
              type->scatterParticle(particle, currRNG);
              T newTau =
                  type->getNewTau(particle.valley, particle.region, currRNG);
              particle.tau += newTau;
              valley = type->getValley(particle.valley);
              driftParticle(std::min(tRemaining, newTau), particle, valley, pos,
                            force);
              tRemaining -= newTau;
            }
            particle.tau -= tStep;

            // handle grain scattering
            particle.grainTau -= tStep;
            if (particle.grainTau <= 0) {
              type->scatterParticleAtGrain(particle, currRNG);
              particle.grainTau = type->getNewGrainTau(currRNG);
            }
          }
        }
      }
    }
  }

  void print(std::string namePrefix, std::string nameSuffix) const {
    for (const auto &[idxType, type] : idxTypeToPartType) {
      std::ofstream os;
      os.open(namePrefix + type->getName() + nameSuffix + ".txt");
      os << device.getMaxPos() << "\n";
      auto nrPart = positionsParticles[idxType].size();
      for (SizeType idxPart = 0; idxPart < nrPart; idxPart++) {
        os << idxPart << " ";
        os << positionsParticles[idxType][idxPart];
        if (type->isMoved()) {
          os << " ";
          os << particles[idxType][idxPart].k << " ";
          os << particles[idxType][idxPart].energy << " ";
          os << particles[idxType][idxPart].subValley << " ";
          os << particles[idxType][idxPart].valley;
        }
        if (idxPart < nrPart - 1)
          os << "\n";
      }
      os.close();
    }
  }

  //! Prints the current drift velocity of each particle.
  void printDriftVelocities(std::ofstream &os) const {
    for (const auto &[idxType, partType] : idxTypeToPartType) {
      if (partType->isMoved()) {
        auto nrPart = positionsParticles[idxType].size();
        for (SizeType idxPart = 0; idxPart < nrPart; idxPart++) {
          auto &part = particles[idxType][idxPart];
          auto valley = partType->getValley(part.valley);
          auto vel = valley->getVelocity(part.k, part.energy, part.subValley);
          auto driftVel = innerProduct(vel, appliedFieldDir);
          os << driftVel;
          if (idxPart < nrPart - 1)
            os << " ";
        }
        os << std::endl;
      }
    }
  }

  //! Prints current velocity of each particle.
  void printVelocities(std::ofstream &os) const {
    for (const auto &[idxType, partType] : idxTypeToPartType) {
      if (partType->isMoved()) {
        auto nrPart = positionsParticles[idxType].size();
        for (SizeType idxPart = 0; idxPart < nrPart; idxPart++) {
          auto &part = particles[idxType][idxPart];
          auto valley = partType->getValley(part.valley);
          auto vel = valley->getVelocity(part.k, part.energy, part.subValley);
          os << vel;
          if (idxPart < nrPart - 1)
            os << " ";
        }
        os << std::endl;
      }
    }
  }

  //! Returns the current valley occupation for each valley of particleType
  //! with given index idxType.
  std::vector<T> getValleyOccupationProbability(SizeType idxType) const {
    auto &partType = idxTypeToPartType[idxType];
    std::vector<T> countPart(partType->getNrValleys(), 0.);
    if (getNrParticles(idxType) != 0) {
      for (auto &part : particles[idxType]) {
        countPart[part.valley]++;
      }
      std::for_each(countPart.begin(), countPart.end(),
                    [&](T &nrPart) { nrPart /= getNrParticles(idxType); });
    }
    return countPart;
  }

  //! Returns the current average energy for each valley of particleType
  //! with given index idxType.
  std::vector<T> getAvgEnergy(SizeType idxType) const {
    auto &partType = idxTypeToPartType[idxType];
    std::vector<T> avgEnergy(partType->getNrValleys(), 0.);
    if (getNrParticles(idxType) != 0) {
      std::vector<SizeType> countPart(partType->getNrValleys(), 0);
      for (auto &part : particles[idxType]) {
        avgEnergy[part.valley] += part.energy;
        countPart[part.valley]++;
      }
      std::transform(countPart.begin(), countPart.end(), avgEnergy.begin(),
                     avgEnergy.begin(), [](SizeType nrPart, T sum) {
                       if (nrPart != 0)
                         return sum / nrPart;
                       else
                         return 0.;
                     });
    }
    return avgEnergy;
  }

  //! Returns the avgDriftVelocity for each for each valley of particleType
  //! with given index idxType.
  std::vector<T> getAvgDriftVelocity(SizeType idxType) const {
    auto &partType = idxTypeToPartType[idxType];
    std::vector<T> avgDriftVel(partType->getNrValleys(), 0.);
    if (getNrParticles(idxType) != 0) {
      std::vector<SizeType> countPart(partType->getNrValleys(), 0.);
      for (const auto &part : particles[idxType]) {
        auto valley = partType->getValley(part.valley);
        auto velocity =
            valley->getVelocity(part.k, part.energy, part.subValley);
        avgDriftVel[part.valley] += innerProduct(velocity, appliedFieldDir);
        countPart[part.valley]++;
      }
      std::transform(countPart.begin(), countPart.end(), avgDriftVel.begin(),
                     avgDriftVel.begin(), [](SizeType nrPart, T sumDrift) {
                       if (nrPart != 0)
                         return sumDrift / nrPart;
                       else
                         return 0.;
                     });
    }
    return avgDriftVel;
  }

  //! Deletes all current particles.
  void deleteParticles() {
    for (const auto &[idxType, _] : idxTypeToPartType) {
      particles[idxType].clear();
      positionsParticles[idxType].clear();
    }
  }

private:
  //! Helper that adds a particle of a specific type near the given coordinate
  void addParticle(SizeType idxType, const SizeVec &coord, emcRNG &rng) {
    positionsParticles[idxType].push_back(initParticlePos(
        coord, device.getGridExtent(), device.getSpacing(), rng));
    auto &partType = idxTypeToPartType.at(idxType);
    if (partType->isMoved()) {
      emcParticle<T> part;
      part = partType->generateInitialParticle(coord, device, rng);
      particles[idxType].push_back(part);
    }
  }

  //! Helper that moves particles and adapts the position if particle
  //! leaves device domain (uses periodic BC)
  void driftParticle(T dt, emcParticle<T> &part,
                     const emcAbstractValley<T> *valley,
                     std::array<T, Dim> &pos, const std::array<T, 3> &force) {
    drift(dt, part, valley, pos, force);
    auto maxPos = device.getMaxPos();
    std::transform(pos.begin(), pos.end(), maxPos.begin(), pos.begin(),
                   [](auto currPos, auto maxPos) {
                     if (currPos < 0)
                       return currPos + maxPos;
                     if (currPos > maxPos)
                       return currPos - maxPos;
                     return currPos;
                   });
  }
};

#endif // BASIC_BULK_PARTICLE_HANDLER_HPP