#ifndef FMM_BULK_PARTICLE_HANDLER_HPP
#define FMM_BULK_PARTICLE_HANDLER_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <unistd.h>

#include <ParticleType/emcParticleType.hpp>
#include <emcDevice.hpp>
#include <emcGrid.hpp>
#include <emcParticleDrift.hpp>
#include <emcParticleInitialization.hpp>
#include <emcUtil.hpp>

// scalFMM-headers
#include <Components/FSimpleLeaf.hpp>
#include <Containers/FOctree.hpp>
#include <Core/FFmmAlgorithmPeriodic.hpp> // use peridodic BC
#include <Kernels/P2P/FP2PParticleContainer.hpp>
#include <Kernels/Rotation/FRotationCell.hpp>

#ifdef USE_CUTOFF_KERNEL
#include <ParticleHandler/FMM/FRotationKernel.hpp>
#else // USE_CUTOFF_KERNEL
#include <Kernels/Rotation/FRotationKernel.hpp>
#endif // USE_CUTOFF_KERNEL

#include <Arranger/FArrangerPeriodic.hpp>
#include <ParticleHandler/FMM/FBasicParticleContainerMover.hpp>

static constexpr unsigned P = 20;

/**
 * @brief Handler is used for bulk simulations with particle-particle
 * interactions.
 *
 * FMM is used to include the real-space particle-particle interactions.
 *
 * @param device simulated device (in this case only extent of device, doping
 * and spacing are needed)
 * @param idxTypeToPartType map relating the idx of a particleType to the
 * specific particleType
 * @param appliedFieldDir direction of the applied field (normalized)
 * @param appliedField applied background Field
 * @param posTree octTree that stores the positions + index + charge of the
 * particles
 * @param kernels class that implements the FMM-functions L2M, M2M, P2P (1 / r)
 * @param algorithm performs FMM algorithm
 * @param arranger class that can rearrange the octree (put the moved particles
 * in the right leaves again)
 */
template <class T, class DeviceType, SizeType Dim = DeviceType::Dimension>
struct FMMBulkParticleHandler {
  typedef emcParticleType<T, DeviceType> ParticleType;
  typedef std::array<SizeType, Dim> SizeVec;
  typedef std::array<T, 3> ValueVec;
  typedef std::map<SizeType, std::unique_ptr<ParticleType>>
      MapIdxToParticleTypes;

  // typedefs of scalFMM library
  typedef FP2PParticleContainer<T, 2> ContainerClass;
  typedef FSimpleLeaf<T, ContainerClass> LeafClass;
  typedef FRotationCell<T, P> CellClass;
  typedef FOctree<T, CellClass, ContainerClass, LeafClass> OctreeClass;
  typedef FRotationKernel<T, CellClass, ContainerClass, P> KernelClass;
  typedef FFmmAlgorithmPeriodic<T, OctreeClass, CellClass, ContainerClass,
                                KernelClass, LeafClass>
      FmmClass;
  typedef FBasicParticleContainerMover<T, OctreeClass, ContainerClass>
      MoverClass;
  typedef FArrangerPeriodic<T, OctreeClass, ContainerClass, MoverClass>
      ArrangerClass;

  DeviceType &device;
  MapIdxToParticleTypes
      &idxTypeToPartType;   //!< map coupling idxType and PartType
  ValueVec appliedFieldDir; //!< direction of the el. field (normalized)
  ValueVec appliedField;    //!< background electric field
  T coulombFactor; // factor to adapt fmm calculations to coulomb results

  // scalFMM classes
  OctreeClass posTree; // tree stroing position and charge of particles
  FmmClass algorithm;  // algorithm for FMM procedure
  std::unique_ptr<KernelClass> kernels; // kernel (1 / r)
  ArrangerClass arranger;               // needed for rearrangement of tree

  // map coupling idxPart and Particle characteristics
  std::map<FSize, emcParticle<T>> idxPartToParticle;
  // map coupling idxPart to idxType
  std::map<FSize, SizeType> idxPartToIdxType;
  std::vector<SizeType> nrParticles; // counts nr of particles of each type
  FSize nextIdx{0};                  // next available index for particle
  emcRNG rng{static_cast<SizeType>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count())};
  std::uniform_real_distribution<T> dist{0., 1.};

public:
  FMMBulkParticleHandler() = delete;

  /**
   * @brief Construct a new object of FMMBulkParticleHandler.
   *
   * NOTE: the Octree requires a cubic simulation space!
   *
   * @param inDevice simulated device (in this case only extent of device,
   * doping and spacing are required to represent the simulation space).
   * @param inTypes particle types for the simulation.
   * @param inFieldDirection the direction of the applied field
   * @param inFieldStrength field strength
   * @param periodicityParam parameter that determines how often simulation box
   * is repeated to represent periodic BC
   */
  FMMBulkParticleHandler(DeviceType &inDevice, MapIdxToParticleTypes &inTypes,
                         const ValueVec &inFieldDirection, T inFieldStrength,
                         SizeType periodicityParam)
      : device(inDevice), appliedFieldDir(inFieldDirection),
        coulombFactor(constants::ke / device.getMaterial().getEpsR()),
        posTree(5, 2, device.getMaxPos()[0],
                {device.getMaxPos()[0] / 2., device.getMaxPos()[0] / 2.,
                 device.getMaxPos()[0] / 2.}),
        algorithm(&posTree, periodicityParam),
        kernels(new KernelClass(algorithm.extendedTreeHeight(),
                                algorithm.extendedBoxWidth(),
                                algorithm.extendedBoxCenter())),
        arranger(&posTree), idxTypeToPartType(inTypes), idxPartToParticle(),
        idxPartToIdxType(), nrParticles(inTypes.size(), 0) {
    algorithm.setKernel(kernels.get());
    // set electric field
    normalize(appliedFieldDir);
    appliedField = scale(appliedFieldDir, inFieldStrength);
    // initialize scatter tables
    for (const auto &partType : idxTypeToPartType) {
      if (partType.second->isMoved())
        partType.second->initScatterTables();
    }
#ifdef USE_CUTOFF_KERNEL
    std::cout << "Using Cutoff-Kernel!\n";
#else  // USE_CUTOFF_KERNEL
    std::cout << "Not Using Cutoff-Kernel!\n";
#endif // USE_CUTOFF_KERNEL
  }

  /// @brief function that allows the seeding of the used random number
  /// generator.
  /// @param inSeed seed for the rng.
  void setSeed(SizeType inSeed) { rng.seed(inSeed); }

  /// @brief generates initial nr of Particles at each coord
  /// NOTE: the number of particles that are used depends on
  /// the used device (that represents the simulation space)
  /// and the used particleType!
  void generateInitialParticles() {
    emcGrid<T, Dim> pot(device.getGridExtent(), 0);
    for (const auto &type : idxTypeToPartType) {
      SizeVec coord;
      for (coord.fill(0); !device.isEndCoord(coord);
           device.advanceCoord(coord)) {
        auto nrPartToCreate =
            type.second->getInitialNrParticles(coord, device, pot);
        while (nrPartToCreate >= 1) {
          addParticle(type.first, coord, rng);
          nrPartToCreate--;
        }
        if (dist(rng) < nrPartToCreate)
          addParticle(type.first, coord, rng);
      }
    }
  }

  SizeType getNrParticles(SizeType idxType) const {
    return nrParticles[idxType];
  }

  ValueVec getAppliedField() const { return appliedField; }

  /// returns the total nr of repeated simulation boxes (used for the calc
  /// of the potential and the force)
  SizeType getNrRepeatedBoxes() const {
    return std::pow(algorithm.theoricalRepetition(), 3);
  }

  /// prints current number of particles for each particleType
  void printNrParticles() const {
    for (const auto &type : idxTypeToPartType)
      std::cout << "\t" << nrParticles[type.first] << " "
                << type.second->getName() << "\n";
  }

  /// @brief resets the forces and potential in the octree
  /// that is used for the calculation of the Coulomb force.
  void resetForcesAndPotential() {
    posTree.forEachCell(resetCell);
    posTree.forEachLeaf(resetLeaf);
  }

  /// @brief calculates the coulomb force at each particle,
  /// uses the fast multipole method.
  void executeFMM() { algorithm.execute(); }

  /// @brief rearrange tree (+ suppress ouput of rearrange function).
  /// redirect stdout buffer to ignore printf in rearrange function
  /// from:
  /// https://stackoverflow.com/questions/13816994/how-to-disable-printf-function
  void rearrangeTree() {
    fpos_t pos;
    fgetpos(stdout, &pos); // save the position in the file stream
    int fd = dup(
        fileno(stdout)); // use the dup() function to create a copy of stdout
    auto tmp = freopen("dummy.txt", "w", stdout); // redirect stdout
    arranger.rearrange();
    fflush(stdout);
    dup2(fd, fileno(stdout)); // restore the stdout
    close(fd);
    clearerr(stdout);
    fsetpos(stdout, &pos); // move to the correct position
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
    ValueVec position, force;
    posTree.forEachLeaf([&](LeafClass *leaf) {
      ContainerClass *tgt = leaf->getTargets();
      for (int idxCont = 0; idxCont < tgt->getNbParticles(); idxCont++) {
        // get particle characteristics
        FSize idxPart = tgt->getPhysicalValues(1, 0)[idxCont];
        SizeType idxType = idxPartToIdxType[idxPart];
        auto partType = idxTypeToPartType[idxType].get();
        if (partType->isMoved()) {
          auto &part = idxPartToParticle[idxPart];
          getPosition(position, tgt, idxCont);
          getForceAtParticle(force, tgt, idxCont, partType->getCharge());
          auto valley = partType->getValley(part.valley);
          // let particle move
          if (part.tau >= tStep) {
            driftParticle(tStep, part, valley, position, force);
          } else {
            driftParticle(part.tau, part, valley, position, force);
            T tRemaining = tStep - part.tau;
            while (tRemaining > 0) {
              partType->scatterParticle(part, rng);
              valley = partType->getValley(part.valley);
              T newTau = partType->getNewTau(part.valley, part.region, rng);
              part.tau += newTau;
              driftParticle(std::min(tRemaining, newTau), part, valley,
                            position, force);
              tRemaining -= newTau;
            }
          }
          part.tau -= tStep;
          // update position
          tgt->getPositions()[0][idxCont] = position[0];
          tgt->getPositions()[1][idxCont] = position[1];
          tgt->getPositions()[2][idxCont] = position[2];

          // handle grain scattering
          part.grainTau -= tStep;
          if (part.grainTau <= 0) {
            partType->scatterParticleAtGrain(part, rng);
            part.grainTau = partType->getNewGrainTau(rng);
          }
        }
      }
    });
  }

  /// creates a file for each particle type and writes the particle
  /// characteristics into that file
  void print(std::string namePrefix, std::string nameSuffix) {
    std::vector<std::ofstream> streams;
    for (const auto &type : idxTypeToPartType) {
      std::ofstream os;
      os.open(namePrefix + type.second->getName() + nameSuffix + ".txt");
      os << device.getMaxPos() << "\n";
      streams.push_back(std::move(os));
    }
    posTree.forEachLeaf([this, &streams](LeafClass *leaf) {
      ContainerClass *src = leaf->getSrc();
      FSize nrPart = src->getNbParticles();
      for (FSize idxCont = 0; idxCont < nrPart; idxCont++) {
        FSize idxPart = src->getPhysicalValues(1, 0)[idxCont];
        SizeType idxType = idxPartToIdxType[idxPart];
        auto &stream = streams[idxType];
        stream << idxPart << " ";
        stream << src->getPositions()[0][idxCont] << " ";
        stream << src->getPositions()[1][idxCont] << " ";
        stream << src->getPositions()[2][idxCont];
        if (idxTypeToPartType.at(idxType)->isMoved()) {
          const auto &currParticle = idxPartToParticle[idxPart];
          stream << " " << currParticle.k << " ";
          stream << currParticle.energy << " ";
          stream << currParticle.subValley << " ";
          stream << currParticle.valley << " ";
          stream << coulombFactor * src->getPotentials()[idxCont] << " ";
          stream << -coulombFactor * src->getForcesX()[idxCont] << " ";
          stream << -coulombFactor * src->getForcesY()[idxCont] << " ";
          stream << -coulombFactor * src->getForcesZ()[idxCont];
        }
        streams[idxType] << "\n";
      }
    });
    for (auto &os : streams)
      os.close();
  }

  //! Returns the current average energy for each valley of particleType
  //! with given index idxType.
  std::vector<T> getAvgEnergy(SizeType idxType) {
    auto &partType = idxTypeToPartType[idxType];
    std::vector<T> avgEnergy(partType->getNrValleys(), 0.);
    if (getNrParticles(idxType) != 0) {
      std::vector<SizeType> countPart(partType->getNrValleys(), 0.);
      posTree.forEachLeaf(
          [this, &avgEnergy, &countPart, &idxType](LeafClass *leaf) {
            ContainerClass *src = leaf->getSrc();
            FSize nrPart = src->getNbParticles();
            for (FSize idxCont = 0; idxCont < nrPart; idxCont++) {
              FSize idxPart = src->getPhysicalValues(1, 0)[idxCont];
              auto currIdxType = idxPartToIdxType[idxPart];
              if (idxType == currIdxType) {
                auto &part = idxPartToParticle[idxPart];
                avgEnergy[part.valley] += part.energy;
                countPart[part.valley]++;
              }
            }
          });
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
  std::vector<T> getAvgDriftVelocity(SizeType idxType) {
    auto &partType = idxTypeToPartType[idxType];
    std::vector<T> avgDriftVel(partType->getNrValleys(), 0.);
    if (getNrParticles(idxType) != 0) {
      std::vector<SizeType> countPart(partType->getNrValleys(), 0.);
      posTree.forEachLeaf([this, &avgDriftVel, &countPart, &idxType,
                           &partType](LeafClass *leaf) {
        ContainerClass *src = leaf->getSrc();
        FSize nrPart = src->getNbParticles();
        for (FSize idxCont = 0; idxCont < nrPart; idxCont++) {
          FSize idxPart = src->getPhysicalValues(1, 0)[idxCont];
          auto currIdxType = idxPartToIdxType[idxPart];
          if (currIdxType == idxType) {
            auto &part = idxPartToParticle[idxPart];
            auto valley = partType->getValley(part.valley);
            auto k = valley->getVelocity(part.k, part.energy, part.subValley);
            avgDriftVel[part.valley] += innerProduct(k, appliedFieldDir);
            countPart[part.valley]++;
          }
        }
      });
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

  //! Returns the current valley occupation for each valley of particleType
  //! with given index idxType.
  std::vector<T> getValleyOccupationProbability(SizeType idxType) {
    auto &partType = idxTypeToPartType[idxType];
    std::vector<T> countPart(partType->getNrValleys(), 0.);
    if (getNrParticles(idxType) != 0) {
      posTree.forEachLeaf([this, &countPart, &idxType](LeafClass *leaf) {
        ContainerClass *src = leaf->getSrc();
        FSize nrPart = src->getNbParticles();
        for (FSize idxCont = 0; idxCont < nrPart; idxCont++) {
          FSize idxPart = src->getPhysicalValues(1, 0)[idxCont];
          auto currIdxType = idxPartToIdxType[idxPart];
          if (idxType == currIdxType) {
            auto &part = idxPartToParticle[idxPart];
            countPart[part.valley]++;
          }
        }
      });
      std::for_each(countPart.begin(), countPart.end(),
                    [&](T &nrPart) { nrPart /= getNrParticles(idxType); });
    }
    return countPart;
  }

private:
  //! Helper that adds a particle of a specific type near the given coordinate
  void addParticle(SizeType idxType, const SizeVec &coord, emcRNG &rng) {
    auto partType = idxTypeToPartType.at(idxType).get();
    auto pos = initParticlePos(coord, device.getGridExtent(),
                               device.getSpacing(), rng);
    posTree.insert({pos[0], pos[1], pos[2]}, partType->getCharge(), nextIdx);
    idxPartToIdxType[nextIdx] = idxType;
    if (partType->isMoved()) {
      emcParticle<T> part;
      part = partType->generateInitialParticle(coord, device, rng);
      idxPartToParticle[nextIdx] = part;
    }
    nrParticles[idxType]++;
    nextIdx++;
  }

  static void resetCell(CellClass *cell) { cell->resetToInitialState(); }

  static void resetLeaf(LeafClass *leaf) {
    ContainerClass *src = leaf->getSrc();
    src->resetForcesAndPotential();
  }

  /// helper to get force of particle in container, force combined
  /// of electric background field and fmm-calculations
  void getForceAtParticle(ValueVec &force, ContainerClass *container,
                          FSize idxCont, T charge) const {
    force[0] = container->getForcesX()[idxCont] * -coulombFactor +
               charge * appliedField[0];
    force[1] = container->getForcesY()[idxCont] * -coulombFactor +
               charge * appliedField[1];
    force[2] = container->getForcesZ()[idxCont] * -coulombFactor +
               charge * appliedField[2];
  }

  //! Helper to get the particle position at idxCont in container
  void getPosition(ValueVec &position, ContainerClass *container,
                   FSize idxCont) const {
    position[0] = container->getPositions()[0][idxCont];
    position[1] = container->getPositions()[1][idxCont];
    position[2] = container->getPositions()[2][idxCont];
  }

  //! Helper that lets particle drift + uses periodic boundary conditions
  //! if particle leaves device
  void driftParticle(T dt, emcParticle<T> &part,
                     const emcAbstractValley<T> *valley, ValueVec &pos,
                     const ValueVec &force) {
    drift(dt, part, valley, pos, force);
    auto deviceMaxPos = device.getMaxPos();
    std::transform(pos.begin(), pos.end(), deviceMaxPos.begin(), pos.begin(),
                   [](auto currPos, auto deviceMaxPos) {
                     if (currPos < 0)
                       return currPos + deviceMaxPos;
                     if (currPos > deviceMaxPos)
                       return currPos - deviceMaxPos;
                     return currPos;
                   });
  }
};

#endif // FMM_BULK_PARTICLE_HANDLER_HPP