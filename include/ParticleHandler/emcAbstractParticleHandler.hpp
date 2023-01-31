#ifndef EMC_ABSTRACT_PARTICLE_HANDLER_HPP
#define EMC_ABSTRACT_PARTICLE_HANDLER_HPP

#include <array>
#include <chrono>
#include <map>
#include <memory>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <ParticleType/emcParticleType.hpp>
#include <emcGrid.hpp>
#include <emcParticle.hpp>
#include <emcUtil.hpp>

/*! \brief Abstract Particle Handler. This class should store and
 * handle all the simulated particles.
 *
 * The derived classes determine how the particles are
 * stored, move and behave at boundaries.
 *
 * @param device device in which the particles are created
 * @param pmScheme Particle - Mesh Scheme
 * @param nrCarriersPerPart nr of carriers that are simulated
 * with one particle (// TODO test with values > 1)
 * @param idxTypeToPartType Map that links the index of a specific
 * particle type to the particle types.
 * @param expNrPart stores the expected nr of particles at ohmic
 * contacts (is needed if particles are injected)
 * @param rngs random number generator(s) (for each thread)
 */
template <class T, class DeviceType, class PMScheme, SizeType Dim>
struct emcAbstractParticleHandler {
  typedef emcParticleType<T, DeviceType> ParticleType;
  typedef typename DeviceType::SizeVec SizeVec;
  typedef typename DeviceType::ValueVec ValueVec;
  //! Type relates each particleType to an unique index (idxType)
  typedef std::map<SizeType, std::unique_ptr<ParticleType>>
      MapIdxToParticleTypes;
  //! type is used to count the number of particles of
  //! each particleType that leave / are injected at each contact
  typedef std::vector<std::vector<int>> NettoParticleCounter;

protected:
  const DeviceType &device;
  PMScheme &pmScheme;
  const SizeType nrCarriersPerPart;
  MapIdxToParticleTypes &idxTypeToPartType;
  std::vector<emcRNG> rngs;

  std::vector<emcGrid<T, Dim>> expNrPart;
  std::uniform_real_distribution<T> dist{0., 1.};

public:
  emcAbstractParticleHandler() = delete;

  /**
   * @brief Construct new emcAbstractParticleHandler object, which assumes
   * that each simulated carrier represents exactly one particle.
   *
   * @param inDevice device which is simulated
   * @param inPMScheme applied particle-mesh scheme
   * @param inIdxToPartTypesMap map which holds all indexed particle types
   * which are used for the simulation
   * @param inSeed seed for the random number generator(s)
   */
  emcAbstractParticleHandler(
      const DeviceType &inDevice, PMScheme &inPMScheme,
      MapIdxToParticleTypes &inIdxToPartTypesMap,
      SizeType inSeed = static_cast<SizeType>(
          std::chrono::high_resolution_clock::now().time_since_epoch().count()))
      : emcAbstractParticleHandler(inDevice, inPMScheme, 1, inIdxToPartTypesMap,
                                   inSeed) {}

  /**
   * @brief Construct new emcAbstractParticleHandler object.
   *
   * @param inDevice device which is simulated
   * @param inPMScheme applied particle-mesh scheme
   * @param inNrCarriersPerPart number of carriers that are represented by
   * one simulated particle (if > 1, superparticle)
   * @param inIdxToPartTypesMap map which holds all indexed particle types
   * which are used for the simulation
   * * @param inSeed seed for the random number generator(s)
   */
  emcAbstractParticleHandler(
      const DeviceType &inDevice, PMScheme &inPMScheme,
      SizeType inNrCarriersPerPart, MapIdxToParticleTypes &inIdxToPartTypesMap,
      SizeType inSeed = static_cast<SizeType>(
          std::chrono::high_resolution_clock::now().time_since_epoch().count()))
      : device(inDevice), pmScheme(inPMScheme),
        nrCarriersPerPart(inNrCarriersPerPart),
        idxTypeToPartType(inIdxToPartTypesMap) {
    calcExpectedNrParticlesAtContact();
    // initialize scatter Handler
    for (const auto &partType : idxTypeToPartType) {
      if (partType.second->isMoved())
        partType.second->initScatterTables();
      // initialize random number generator(s)
      rngs.emplace_back(emcRNG(inSeed));
#ifdef _OPENMP
      for (SizeType idxSeed = 1; idxSeed < omp_get_max_threads(); idxSeed++)
        rngs.emplace_back(emcRNG(rngs[0]()));
#endif
    }
  }

  /*! \brief function that tells if derived class includes the
   * particle - particle interaction in the calculation of
   * the force that is seen by the particles or not.
   *
   * This function is needed because different types if potentials
   * have to be calculated if the particle-particle interaction
   * is included or not (see \file emcSimulation.hpp)
   */
  virtual bool calcsPartPartInteraction() const = 0;

  SizeType getNrParticleTypes() const { return idxTypeToPartType.size(); }

  virtual void printNrParticles() const = 0;

  /*! \brief Generates the intial particles of every particleType.
   *
   * At every grid coordinate the expected number of initial particles
   * are added. If number of expected particles is not integer, a random
   * number is used to determine if an additional particle is added.
   *
   * @param potential equilibrium potential of the grid
   */
  void generateInitialParticles(const emcGrid<T, Dim> &potential) {
    for (const auto &type : idxTypeToPartType) {
      SizeVec coord;
      for (coord.fill(0); !device.isEndCoord(coord);
           device.advanceCoord(coord)) {
        auto nrPartToCreate =
            type.second->getInitialNrParticles(coord, device, potential);
        while (nrPartToCreate >= 1) {
          addParticle(type.first, coord, rngs[0], true);
          nrPartToCreate -= nrCarriersPerPart;
        }
        if (dist(rngs[0]) < nrPartToCreate)
          addParticle(type.first, coord, rngs[0], true);
      }
    }
  }

  /*! \brief Assigns every particle to the given grid.
   *
   * @param idxType index of the required particleType
   * @param gridNrParticles grid that is overwritten, in the end it stores
   * the number of particles that is "near" each grid point
   */
  virtual void assignParticlesToMesh(SizeType idxType,
                                     emcGrid<T, Dim> &gridNrParticles) = 0;

  /*! \brief Moves Particles for the given amount of time using the el. field.
   *
   * @param tStep amount of time for which the particles should move
   * @param eField el. Field calculated on the grid points
   * @return counted number of particles of all particleTpes that leave
   * the device through the different contacts
   */
  virtual NettoParticleCounter
  driftScatterParticles(T tStep, std::vector<emcGrid<T, Dim>> &eField) = 0;

  /*! \brief Handles particles near an ohmic contact.
   *
   * @return counted number of particles of all particleTypes that leave
   * / got injected at different contacts
   */
  virtual NettoParticleCounter handleOhmicContacts() = 0;

  //! prints information on all particles
  virtual void print(std::string namePrefix, std::string nameSuffix) = 0;

  /*! \brief Calculates potential at a specific position.
   *
   * Only interesting in case particle-particle interaction is considered.
   */
  virtual T getParticlePotential(const ValueVec &position) { return 0; }

protected:
  /// \brief Creates one new particle of a specific particleType.
  virtual void addParticle(SizeType idxType, const SizeVec &coord, emcRNG &rng,
                           bool isInitial) = 0;

  /*! \brief Handles the injection of particles at ohmic contacts.
   *
   * Compares the nr of particles in the grid to the nr of particles
   * that are expected near the contact and if the number of actual
   * particles is too low, new particles are injected.
   *
   * @param idxType index of the used particleType
   * @param nrPart grid that stores the nr of particles that
   * are currently present near the contact
   */
  std::vector<int> generateInjectedParticles(SizeType idxType,
                                             emcGrid<T, Dim> &nrPart) {
    SizeVec coord;
    auto &surface = device.getSurface();
    std::vector<int> nrInjPart(surface.getNrContacts(), 0);
    for (coord.fill(0); !nrPart.isEndCoord(coord); nrPart.advanceCoord(coord)) {
      if (surface.isOhmicContact(coord)) {
        T nrDiffEl = expNrPart[idxType][coord] - nrPart[coord];
        while (nrDiffEl > 0) {
          addParticle(idxType, coord, rngs[0], false);
          nrInjPart[surface.getOhmicContactIdx(coord)]++;
          nrDiffEl -= nrCarriersPerPart;
        }
      }
    }
    return nrInjPart;
  }

  /*! \brief Helper function, returns initialized type particle counter.
   *
   * Initializes the NettoParticleCounter with 0, meaning that at the
   * beginning no particle has been injected or deleted at any contact.
   * The netto number of particles that then left / got injected at
   * a contact can then be accessed with the two indices idxType and
   * idxContact.
   */
  NettoParticleCounter initNettoParticleCounter() {
    NettoParticleCounter nettoNrPart(
        idxTypeToPartType.size(),
        std::vector<int>(device.getSurface().getNrContacts(), 0));
    return nettoNrPart;
  }

  /*! \brief Helper function, that moves particle for a specific
   * amount of time and also handles particles that leave through a boundary.
   *
   * @return boolean that tells if particle should be removed or not
   */
  bool driftParticle(T dt, emcParticle<T> &part, ParticleType *particleType,
                     std::array<T, Dim> &pos, const std::array<T, 3> &force,
                     emcRNG &rng) {
    drift(dt, part, particleType->getValley(part.valley), pos, force);
    bool removed =
        handleParticleAtBoundary(part, pos, particleType, device, rng);
    if (!removed) {
      auto coord = device.posToCoord(pos);
      part.region = device.getDopingProfile().getDopingRegionIdx(coord);
    }
    return removed;
  }

  //! \brief Helper function, that returns the current thread index, can
  //! be used for parallelization.
  SizeType getCurrentThreadIdx() {
    SizeType idxThread = 0;
#ifdef _OPENMP
    idxThread = omp_get_thread_num();
#endif
    return idxThread;
  }

private:
  //! \brief Calculates expected number of particles near an ohmic contact.
  void calcExpectedNrParticlesAtContact() {
    emcGrid<T, Dim> tmp(device.getGridExtent(), 0);
    for (const auto &type : idxTypeToPartType) {
      expNrPart.push_back(tmp);
      if (type.second->isInjected()) {
        SizeVec coord;
        for (coord.fill(0); !tmp.isEndCoord(coord); tmp.advanceCoord(coord)) {
          if (device.getSurface().isOhmicContact(coord))
            expNrPart[type.first][coord] =
                type.second->getExpectedNrParticlesAtContact(coord, device);
        }
      }
    }
  }
};

#endif // EMC_ABSTRACT_PARTICLE_HANDLER_HPP