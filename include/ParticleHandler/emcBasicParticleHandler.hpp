#ifndef EMC_BASIC_PARTICLE_HANDLER_HPP
#define EMC_BASIC_PARTICLE_HANDLER_HPP

#ifdef _OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <random>
#include <unistd.h>

#include <ParticleHandler/emcAbstractParticleHandler.hpp>
#include <emcParticleDrift.hpp>
#include <emcParticleInitialization.hpp>

/*! \brief Basic Particle Handler, for case when FMM is not used.
 *
 * Class contains containers for all particles and particle positions
 * of all particle types.
 *
 * @param particles stores all the particles of all particle types;
 * to access a specific particle of a specific type, first use the
 * index of the particle type (idxType) and then the index of the
 * particle (idxPart)
 * @param positionsParticles stores the positions of all particles;
 * again to access a specific particle, first the index of the particle
 * type (idxType) and then the index of the particle (idxPart) has to
 * be used.
 */
template <class T, class DeviceType, class PMScheme>
class emcBasicParticleHandler
    : public emcAbstractParticleHandler<T, DeviceType, PMScheme,
                                        DeviceType::Dimension> {
  static constexpr SizeType Dim = DeviceType::Dimension;

  typedef emcAbstractParticleHandler<T, DeviceType, PMScheme, Dim> Base;
  typedef typename Base::ParticleType ParticleType;
  typedef typename Base::SizeVec SizeVec;
  typedef typename Base::ValueVec ValueVec;
  typedef typename Base::MapIdxToParticleTypes MapIdxToParticleTypes;
  typedef typename Base::NettoParticleCounter NettoParticleCounter;

  std::vector<std::vector<emcParticle<T>>> particles;
  std::vector<std::vector<std::array<T, Dim>>> positionsParticles;

  ValueVec deviceSpacing;
  SizeVec deviceExtent;

public:
  emcBasicParticleHandler() = delete;

  emcBasicParticleHandler(const DeviceType &inDevice, PMScheme &inPMScheme,
                          MapIdxToParticleTypes &inTypes,
                          SizeType inNrCarriersPerParticle, SizeType inSeed)
      : Base(inDevice, inPMScheme, inNrCarriersPerParticle, inTypes, inSeed),
        deviceSpacing(inDevice.getSpacing()),
        deviceExtent(inDevice.getGridExtent()) {
    particles.resize(Base::getNrParticleTypes());
    positionsParticles.resize(Base::getNrParticleTypes());
  }

  bool calcsPartPartInteraction() const { return false; }

  //! Returns nr of particles of a specific particleType.
  SizeType getNrParticles(SizeType idxType) const {
    return positionsParticles[idxType].size();
  }

  void printNrParticles() const {
    for (const auto &type : Base::idxTypeToPartType)
      std::cout << "\t" << positionsParticles[type.first].size() << " "
                << type.second->getName() << "\n";
  }

  //! \brief Moves Particles for the given amount of time based on
  //! the given el. Field.
  NettoParticleCounter
  driftScatterParticles(T tStep, std::vector<emcGrid<T, Dim>> &eField) {
    auto &rngs = Base::rngs;
    auto &surface = Base::device.getSurface();
    auto nettoNrPart = Base::initNettoParticleCounter();

    for (const auto &partType : Base::idxTypeToPartType) {
      SizeType idxType = partType.first;
      auto &type = partType.second;
      if (type->isMoved()) {
        auto charge = type->getCharge();
        std::vector<SizeType> idxToRemove;
#pragma omp parallel
        {
          SizeType idxThread = this->getCurrentThreadIdx();
          auto &currRNG = rngs[idxThread];
#pragma omp for
          for (SizeType idxPart = 0; idxPart < getNrParticles(idxType);
               idxPart++) {
            auto &particle = particles[idxType][idxPart];
            auto &pos = positionsParticles[idxType][idxPart];
            auto force = Base::pmScheme.interpolateForce(eField, pos,
                                                         deviceSpacing, charge);
            bool removed = false;
            if (particle.tau >= tStep) { // drift for whole time-step
              removed = Base::driftParticle(tStep, particle, type.get(), pos,
                                            force, currRNG);
            } else { // drift for remaining time and scatter
              removed = Base::driftParticle(particle.tau, particle, type.get(),
                                            pos, force, currRNG);
              T tRemaining = tStep - particle.tau;
              while (tRemaining > 0 && !removed) {
                type->scatterParticle(particle, currRNG);
                T newTau =
                    type->getNewTau(particle.valley, particle.region, currRNG);
                particle.tau += newTau;
                force = Base::pmScheme.interpolateForce(eField, pos,
                                                        deviceSpacing, charge);
                removed =
                    Base::driftParticle(std::min(tRemaining, newTau), particle,
                                        type.get(), pos, force, currRNG);
                tRemaining -= newTau;
              }
            }
            particle.tau -= tStep;
            if (removed) {
              SizeType idxCont =
                  surface.getOhmicContactIdx(Base::device.posToCoord(pos));
#pragma omp critical
              {
                nettoNrPart[idxType][idxCont]++;
                idxToRemove.push_back(idxPart);
              }
            }
            // TODO handle grain scattering as normal scattering?
            particle.grainTau -= tStep;
            if (particle.grainTau <= 0 && !removed) {
              type->scatterParticleAtGrain(particle, currRNG);
              particle.grainTau = type->getNewGrainTau(currRNG);
            }
          }
        }
        removeParticles(idxType, idxToRemove);
      }
    }
    return nettoNrPart;
  }

  //! Counts nr of particles "near" each grid point using given pmScheme
  void assignParticlesToMesh(SizeType idxType,
                             emcGrid<T, Dim> &gridNrParticles) {
    Base::pmScheme.assignToMesh(positionsParticles[idxType],
                                Base::nrCarriersPerPart, deviceSpacing,
                                gridNrParticles);
  }

  //! \brief Counts the number of particles near an ohmic contact.
  //! if too many particles are there -> delete excess particles.
  //! if nr particles too low -> inject particles.
  NettoParticleCounter handleOhmicContacts() {
    SizeVec coord;
    auto nrRemPart = Base::initNettoParticleCounter();
    auto nrInjPart = Base::initNettoParticleCounter();
    auto &surface = Base::device.getSurface();
    emcGrid<T, Dim> nrPart(deviceExtent, 0);
    for (const auto &[idxType, partType] : Base::idxTypeToPartType) {
      if (partType->isInjected()) {
        nrPart.fill(0);
        for (SizeType idxPart = 0; idxPart < getNrParticles(idxType);
             idxPart++) {
          // count particles at each contact + delete excess particles
          auto &posPart = positionsParticles[idxType][idxPart];
          coord = Base::device.posToCoord(posPart);
          if (surface.isOhmicContact(coord)) {
            if (nrPart[coord] < Base::expNrPart[idxType][coord]) {
              nrPart[coord] += Base::nrCarriersPerPart;
            } else {
              removeParticle(idxType, idxPart);
              idxPart--;
              nrRemPart[idxType][surface.getOhmicContactIdx(coord)]--;
            }
          }
        }
        // generate missing particles
        nrInjPart[idxType] = Base::generateInjectedParticles(idxType, nrPart);

        std::transform(nrRemPart[idxType].begin(), nrRemPart[idxType].end(),
                       nrInjPart[idxType].begin(), nrRemPart[idxType].begin(),
                       std::plus<int>());
      }
    }

    return nrRemPart;
  }

  void print(std::string namePrefix, std::string nameSuffix) {
    for (const auto &type : Base::idxTypeToPartType) {
      std::ofstream os;
      os.open(namePrefix + type.second->getName() + nameSuffix + ".txt");
      os << Base::device.getMaxPos() << "\n";
      auto nrPart = positionsParticles[type.first].size();
      for (SizeType idxPart = 0; idxPart < nrPart; idxPart++) {
        os << idxPart << " ";
        os << positionsParticles[type.first][idxPart];
        if (type.second->isMoved()) {
          os << " ";
          os << particles[type.first][idxPart].k << " ";
          os << particles[type.first][idxPart].energy << " ";
          os << particles[type.first][idxPart].subValley << " ";
          os << particles[type.first][idxPart].valley << " ";
          os << particles[type.first][idxPart].tau;
        }
        if (idxPart < nrPart - 1)
          os << "\n";
      }
      os.close();
    }
  }

private:
  //! helper that adds a particle of a specific type near the given coordinate
  void addParticle(SizeType idxType, const SizeVec &coord, emcRNG &rng,
                   bool isInitial) {
    positionsParticles[idxType].push_back(
        initParticlePos(coord, deviceExtent, deviceSpacing, rng));
    auto &partType = Base::idxTypeToPartType.at(idxType);
    if (partType->isMoved()) {
      emcParticle<T> part;
      if (isInitial)
        part = partType->generateInitialParticle(coord, Base::device, rng);
      else
        part = partType->generateInjectedParticle(coord, Base::device, rng);
      particles[idxType].push_back(part);
    }
  }

  //! helper function that removes single particle of type idxType with index
  //! idxPart
  void removeParticle(SizeType idxType, SizeType idxParticle) {
    positionsParticles[idxType].erase(positionsParticles[idxType].begin() +
                                      idxParticle);
    if (Base::idxTypeToPartType.at(idxType)->isMoved())
      particles[idxType].erase(particles[idxType].begin() + idxParticle);
  }

  //! helper function that removes multiple particles of type idxType with
  //! index idxPart (careful: can't be parallelized!)
  void removeParticles(SizeType idxType, std::vector<SizeType> idxParticles) {
    std::sort(idxParticles.begin(), idxParticles.end());
    bool isMoved = Base::idxTypeToPartType.at(idxType)->isMoved();
    for (SizeType i = 0; i < idxParticles.size(); i++) {
      SizeType currPos = idxParticles[i] - i;
      positionsParticles[idxType].erase(positionsParticles[idxType].begin() +
                                        currPos);
      if (isMoved)
        particles[idxType].erase(particles[idxType].begin() + currPos);
    }
  }
};

#endif // EMC_BASIC_PARTICLE_HANDLER_HPP