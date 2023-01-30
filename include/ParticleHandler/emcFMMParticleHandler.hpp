#ifndef EMC_FMM_PARTICLE_HANDLER_HPP
#define EMC_FMM_PARTICLE_HANDLER_HPP

#include <fstream>

#include <ParticleHandler/emcAbstractParticleHandler.hpp>
#include <emcParticleDrift.hpp>
#include <emcParticleInitialization.hpp>

// scalFMM-headers
#include <Components/FSimpleLeaf.hpp>
#include <Containers/FOctree.hpp>
#include <Containers/FVector.hpp>

// for the openmp headers.
#ifdef _OPENMP
#include "Core/FFmmAlgorithmThread.hpp"
#else
#include <Core/FFmmAlgorithm.hpp>
#endif

#include <Kernels/P2P/FP2PParticleContainer.hpp>

#include <Arranger/FOctreeArranger.hpp>
#include <ParticleHandler/FMM/FBasicParticleContainerMover.hpp>

/// approx. of Rotation-Kernel based on spherical harmonics
/// see: http://users.cecs.anu.edu.au/~Eric.McCreath/students/haigh2011.pdf
#include <Kernels/Rotation/FRotationCell.hpp>
#ifdef USE_CUTOFF_KERNEL
#include <ParticleHandler/FMM/FRotationKernel.hpp>
#else // USE_CUTOFF_KERNEL
#include <Kernels/Rotation/FRotationKernel.hpp>
#endif // USE_CUTOFF_KERNEL

static constexpr unsigned P = 20;

/*! \brief FMM Particle Handler using scalFMM to calculate the particle
 * particle interactions.
 *
 * Particles are stored via their unique index (idxPart). Different
 * characteristics of the particles are stored in different containers
 * (see posTree, idxPartToParticle, idxPartToIdxType) and are
 * accessible via this unique particle index.
 *
 * The movement of the particles is governed by the sum of the given
 * electric field and the potential coming from the particle-particle
 * interactions.
 *
 * Particles that leave / got injected through ohmic contacts are counted
 * with the type NettoParticleCounter and the particle concentration near
 * the ohmic contact is kept constant.
 *
 * @param treeBoxWidth width of the cubic tree extent (= maximal device extent)
 * @param posTree octTree that stores the positions + index + charge of the
 * particles
 * @param kernels class that implements the FMM-functions L2M, M2M, P2P (1 / r)
 * @param algorithm performs FMM algorithm
 * @param arranger class that can rearrange the octree (put the moved particles
 * in the right leaves again)
 * @param idxPartToParticle map that maps the index of a particle (idxPart) to
 * the characteristics of the particle (like waveVec, energy,...)
 * @param idxPartToIdxType map that maps the index of the particle (idxPart) to
 * the index of the particle type (idxType)
 * @param nextIdx index for the next particle (currently always increased when a
 * particle is added)
 */
template <class T, class DeviceType, class PMScheme>
struct emcFMMParticleHandler
    : emcAbstractParticleHandler<T, DeviceType, PMScheme,
                                 DeviceType::Dimension> {
  static_assert(DeviceType::Dimension == 3,
                "FMMParticleContainer can only be used in 3 Dimensions.");

  static const SizeType Dim = DeviceType::Dimension;

  // typedefs of Base-Class
  typedef emcAbstractParticleHandler<T, DeviceType, PMScheme, Dim> Base;
  typedef typename Base::ParticleType ParticleType;
  typedef typename Base::SizeVec SizeVec;
  typedef typename Base::ValueVec ValueVec;
  typedef typename Base::MapIdxToParticleTypes MapIdxToParticleTypes;
  typedef typename Base::NettoParticleCounter NettoParticleCounter;

  // typedefs of scalFMM library
  typedef FP2PParticleContainer<T, 2> ContainerClass;
  typedef FSimpleLeaf<T, ContainerClass> LeafClass;
  typedef FRotationCell<T, P> CellClass;
  typedef FOctree<T, CellClass, ContainerClass, LeafClass> OctreeClass;
  typedef FRotationKernel<T, CellClass, ContainerClass, P> KernelClass;
#ifdef _OPENMP
  typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass,
                              KernelClass, LeafClass>
      FmmClass;
#else
  typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass,
                        LeafClass>
      FmmClass;
#endif
  typedef FBasicParticleContainerMover<T, OctreeClass, ContainerClass>
      MoverClass;
  typedef FOctreeArranger<T, OctreeClass, ContainerClass, MoverClass>
      ArrangerClass;

  ValueVec deviceMaxPos;  //!< maximal position of device
  ValueVec deviceSpacing; //!< spacing of grid of the device
  SizeVec deviceExtent;   //!< extent of the grid of the device
  T treeBoxWidth; //!< width of the cubic tree extent (= maximal device extent)

  OctreeClass
      posTree; //!< octree storing positions + charges + idx of particles
  std::unique_ptr<KernelClass> kernels; //!< class implements P2P, M2M, M2L, ...
  FmmClass algorithm;                   //!< class to perform FMM algorithm
  ArrangerClass arranger;               //!< class to rearrange tree

  /// particle informations which is accessible via idx of the particle
  std::map<FSize, emcParticle<T>> idxPartToParticle; //!< characteristics
  std::map<FSize, SizeType> idxPartToIdxType;        //!< type

  /// helper members
  FSize nextIdx; //!< next available index (always increased)
  std::vector<SizeType> nrParticles; //!< counts nrParticles of all types
  std::vector<emcGrid<T, Dim>> nrPartGrid;
  T coulombFactor; //!< 1 / (4 * pi * eps_r)

public:
  emcFMMParticleHandler() = delete;

  emcFMMParticleHandler(const DeviceType &inDevice, PMScheme &inPMScheme,
                        MapIdxToParticleTypes &inTypes,
                        SizeType inNrCarriersPerParticle, SizeType inSeed)
      : Base(inDevice, inPMScheme, inNrCarriersPerParticle, inTypes, inSeed),
        deviceMaxPos(inDevice.getMaxPos()),
        deviceSpacing(inDevice.getSpacing()),
        deviceExtent(inDevice.getGridExtent()),
        treeBoxWidth(
            *std::max_element(deviceMaxPos.begin(), deviceMaxPos.end())),
        posTree(5, 2, treeBoxWidth,
                {treeBoxWidth / 2., treeBoxWidth / 2., treeBoxWidth / 2.}),
        kernels(new KernelClass(posTree.getHeight(), posTree.getBoxWidth(),
                                posTree.getBoxCenter())),
        algorithm(&posTree, kernels.get()), arranger(&posTree), nextIdx(0),
        nrParticles(inTypes.size(), 0),
        nrPartGrid(inTypes.size(), emcGrid<T, Dim>(deviceExtent, 0)),
        coulombFactor(-constants::ke / inDevice.getMaterial().getEpsR()) {
#ifdef USE_CUTOFF_KERNEL
    std::cout << "Using Cutoff-Kernel!\n";
#else  // USE_CUTOFF_KERNEL
    std::cout << "Not Using Cutoff-Kernel!\n";
#endif // USE_CUTOFF_KERNEL
  }

  bool calcsPartPartInteraction() const { return true; }

  void printNrParticles() const {
    for (const auto &type : Base::idxTypeToPartType)
      std::cout << "\t" << nrParticles[type.first] << " "
                << type.second->getName() << "\n";
  }

  //! \brief Calculates the Force at each Particle and then moves the
  //! particles.
  NettoParticleCounter
  driftScatterParticles(T tStep, std::vector<emcGrid<T, Dim>> &eField) {
    posTree.forEachCell(resetCell);
    posTree.forEachLeaf(resetLeaf);
    rearrangeTree();
    algorithm.execute();
    return moveParticles(tStep, eField);
  }

  //! Counts nr of particles "near" each grid point using given pmScheme
  void assignParticlesToMesh(SizeType idxType,
                             emcGrid<T, Dim> &gridNrParticles) {
    posTree.forEachLeaf([&](LeafClass *leaf) {
      const ContainerClass *src = leaf->getSrc();
      FSize nrPart = src->getNbParticles();
      const auto *pos = src->getPositions();
      const auto *idx = src->getPhysicalValues(1, 0);
      for (FSize idxPart = 0; idxPart < nrPart; idxPart++) {
        if (idxPartToIdxType[idx[idxPart]] == idxType) {
          ValueVec posPart = {pos[0][idxPart], pos[1][idxPart],
                              pos[2][idxPart]};
          Base::pmScheme.assignToMesh(posPart, Base::nrCarriersPerPart,
                                      deviceSpacing, gridNrParticles);
        }
      }
    });
  }

  //! Returns potential at a specific position.
  T getParticlePotential(ValueVec position) {
    T distSum = 0;
    posTree.forEachLeaf([&](LeafClass *leaf) {
      ContainerClass *src = leaf->getSrc();
      const auto *charge = src->getPhysicalValues();
      auto *pos = src->getPositions();
      SizeType nrPart = src->getNbParticles();
      for (FSize idxPart = 0; idxPart < nrPart; idxPart++) {
        ValueVec distVec = {position[0] - pos[0][idxPart],
                            position[1] - pos[1][idxPart],
                            position[2] - pos[2][idxPart]};
        T dist = norm(distVec);
#ifdef USE_CUTOFF_KERNEL
        if (dist < 1e-9)
          dist = 1e-9;
#endif
        distSum += charge[idxPart] / dist;
      }
    });
    return Base::device.normalizeVoltage(-coulombFactor * distSum);
  }

  //! Creates a file for each particleType and writes the particle
  //! characteristics into that file
  void print(std::string namePrefix, std::string nameSuffix) {
    std::vector<std::ofstream> streams;
    for (const auto &type : Base::idxTypeToPartType) {
      std::ofstream os;
      os.open(namePrefix + type.second->getName() + nameSuffix + ".txt");
      os << deviceMaxPos << "\n";
      streams.push_back(std::move(os));
    }
    posTree.forEachLeaf([this, &streams](LeafClass *leaf) {
      ContainerClass *src = leaf->getSrc();
      FSize nrPart = src->getNbParticles();
      const auto *idx = src->getPhysicalValues(1, 0);
      const auto *pos = src->getPositions();
      const T *potential = src->getPotentials();
      for (FSize idxPart = 0; idxPart < nrPart; idxPart++) {
        FSize idxOrigin = idx[idxPart];
        SizeType idxType = idxPartToIdxType[idxOrigin];
        auto &stream = streams[idxType];
        stream << idxOrigin << " ";
        stream << pos[0][idxPart] << " " << pos[1][idxPart] << " "
               << pos[2][idxPart];
        if (Base::idxTypeToPartType.at(idxType)->isMoved()) {
          const auto &currParticle = idxPartToParticle[idxOrigin];
          stream << " " << currParticle.k << " ";
          stream << currParticle.energy << " ";
          stream << currParticle.subValley << " ";
          stream << currParticle.valley << " ";
          stream << coulombFactor * potential[idxPart] << " ";
          stream << coulombFactor * src->getForcesX()[idxPart] << " ";
          stream << coulombFactor * src->getForcesY()[idxPart] << " ";
          stream << coulombFactor * src->getForcesZ()[idxPart];
        }
        streams[idxType] << "\n";
      }
    });
    for (auto &os : streams)
      os.close();
  }

  //! \brief Counts the number of particles near an ohmic contact.
  //! if too many particles are there -> delete excess particles.
  //! if nr particles too low -> inject particles.
  NettoParticleCounter handleOhmicContacts() {
    NettoParticleCounter nettoPart = Base::initNettoParticleCounter();
    auto &surface = Base::device.getSurface();
    SizeVec coord;
    std::for_each(nrPartGrid.begin(), nrPartGrid.end(),
                  [](auto &grid) { grid.fill(0); });
    posTree.forEachLeaf([&](LeafClass *leaf) {
      ContainerClass *src = leaf->getSrc();
      for (int idxCont = 0; idxCont < src->getNbParticles(); idxCont++) {
        FSize idxType = idxPartToIdxType[src->getPhysicalValues(1, 0)[idxCont]];
        if (Base::idxTypeToPartType[idxType]->isInjected()) {
          coord = Base::device.posToCoord(getPosition(src, idxCont));
          if (surface.isOhmicContact(coord)) {
            nrPartGrid[idxType][coord] += Base::nrCarriersPerPart;
            if (nrPartGrid[idxType][coord] > Base::expNrPart[idxType][coord]) {
              removeParticle(src, idxCont);
              nettoPart[idxType][surface.getOhmicContactIdx(coord)]--;
              idxCont--;
            }
          }
        }
      }
    });
    for (auto &type : Base::idxTypeToPartType) {
      if (type.second->isInjected()) {
        auto injPart =
            Base::generateInjectedParticles(type.first, nrPartGrid[type.first]);
        std::transform(nettoPart[type.first].begin(),
                       nettoPart[type.first].end(), injPart.begin(),
                       nettoPart[type.first].begin(), std::plus<T>());
      }
    }
    return nettoPart;
  }

private:
  //! helper that adds a particle of a specific type near the given coordinate
  void addParticle(SizeType idxType, const SizeVec &coord, emcRNG &rng,
                   bool isInitial) {
    const auto &partType = Base::idxTypeToPartType.at(idxType);
    const auto &pos = initParticlePos(coord, deviceExtent, deviceSpacing, rng);
    posTree.insert({pos[0], pos[1], pos[2]}, partType->getCharge(), nextIdx);
    idxPartToIdxType[nextIdx] = idxType;
    if (partType->isMoved()) {
      emcParticle<T> part;
      if (isInitial)
        part = partType->generateInitialParticle(coord, Base::device, rng);
      else
        part = partType->generateInjectedParticle(coord, Base::device, rng);
      idxPartToParticle[nextIdx] = part;
    }
    nrParticles[idxType]++;
    nextIdx++;
  }

  //! helper function that removes particle of type idxType with index idxPart
  void removeParticle(ContainerClass *container, FSize idxCont) {
    FSize idxPart = container->getPhysicalValues(1, 0)[idxCont];
    container->removeParticles(&idxCont, 1);
    FSize idxType = idxPartToIdxType[idxPart];
    if (Base::idxTypeToPartType[idxType]->isMoved())
      idxPartToParticle.erase(idxPart);
    nrParticles[idxType]--;
    idxPartToIdxType.erase(idxPart);
  }

  //! \brief Moves the particles based on the sum of the force steming from the
  //! el. Field and the one from the particle-particle interaction.
  NettoParticleCounter moveParticles(T tStep,
                                     std::vector<emcGrid<T, Dim>> &eField) {
    auto nettoPart = Base::initNettoParticleCounter();
    ValueVec position, frcPoisson, force;
    posTree.forEachLeaf([&](LeafClass *leaf) {
      ContainerClass *tgt = leaf->getTargets();
      for (int idxCont = 0; idxCont < tgt->getNbParticles(); idxCont++) {
        FSize idxPart = tgt->getPhysicalValues(1, 0)[idxCont];
        SizeType idxType = idxPartToIdxType[idxPart];
        auto partType = Base::idxTypeToPartType[idxType].get();
        if (partType->isMoved()) {
          auto &part = idxPartToParticle[idxPart];
          position = getPosition(tgt, idxCont);
          force = Base::pmScheme.interpolateForce(
              eField, position, deviceSpacing,
              tgt->getPhysicalValues()[idxCont]);
          force = getForceAtParticle(tgt, idxCont, force);
          bool removed = false;
          if (part.tau >= tStep) {
            removed = Base::driftParticle(tStep, part, partType, position,
                                          force, this->rngs[0]);
          } else {
            removed = Base::driftParticle(part.tau, part, partType, position,
                                          force, this->rngs[0]);
            T tRemaining = tStep - part.tau;
            while (tRemaining > 0 && !removed) {
              partType->scatterParticle(part, this->rngs[0]);
              T tFreeFlightNew =
                  partType->getNewTau(part.valley, part.region, this->rngs[0]);
              part.tau += tFreeFlightNew;
              // TODO update force
              removed = Base::driftParticle(
                  std::min(tRemaining, tFreeFlightNew), part, partType,
                  position, force, this->rngs[0]);
              tRemaining -= tFreeFlightNew;
            }
          }
          part.tau -= tStep;

          // TODO add grain scattering
          if (removed) {
            removeParticle(tgt, idxCont);
            nettoPart[idxType][Base::device.getSurface().getOhmicContactIdx(
                Base::device.posToCoord(position))]++;
            idxCont--;
          } else {
            tgt->getPositions()[0][idxCont] = position[0];
            tgt->getPositions()[1][idxCont] = position[1];
            tgt->getPositions()[2][idxCont] = position[2];
          }
        }
      }
    });
    return nettoPart;
  }

  static void resetCell(CellClass *cell) { cell->resetToInitialState(); }

  static void resetLeaf(LeafClass *leaf) {
    leaf->getSrc()->resetForcesAndPotential();
  }

  //! helper function that returns the force at the current particle
  ValueVec getForceAtParticle(ContainerClass *container, FSize idxPart,
                              ValueVec &frcPoisson) const {
    ValueVec force;
    force[0] = container->getForcesX()[idxPart] * coulombFactor + frcPoisson[0];
    force[1] = container->getForcesY()[idxPart] * coulombFactor + frcPoisson[1];
    force[2] = container->getForcesZ()[idxPart] * coulombFactor + frcPoisson[2];
    return force;
  }

  //! helper function that returns the position of the current particle
  ValueVec getPosition(ContainerClass *container, FSize idxPart) const {
    ValueVec pos;
    pos[0] = container->getPositions()[0][idxPart];
    pos[1] = container->getPositions()[1][idxPart];
    pos[2] = container->getPositions()[2][idxPart];
    return pos;
  }

  //! rearrange tree (+ suppress output of rearrange function)
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
};

#endif // EMC_FMM_PARTICLE_HANDLER_HPP