#include <PMSchemes/emcNGPScheme.hpp>
#include <ParticleHandler/emcBasicParticleHandler.hpp>
#include <ParticleHandler/emcFMMParticleHandler.hpp>
#include <ParticleType/emcParticleType.hpp>
#include <ValleyTypes/emcParabolicIsotropValley.hpp>
#include <emcDevice.hpp>
#include <emcParticle.hpp>
#include <emcTestAsserts.hpp>

const SizeType Dim = 3;

using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, Dim>;
using PMScheme = emcNGPScheme<NumType, DeviceType>;
using BasicHandler = emcBasicParticleHandler<NumType, DeviceType, PMScheme>;
using FMMHandler = emcFMMParticleHandler<NumType, DeviceType, PMScheme>;

using MapIdxToParticleTypes = FMMHandler::MapIdxToParticleTypes;

/// exemplary moving particle type
template <class T, class DeviceType>
struct emcMovingType : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::SizeVec SizeVec;
  typedef typename DeviceType::ValueVec ValueVec;
  static const SizeType Dim = DeviceType::Dimension;

  std::string getName() const { return "MovingType"; }
  T getMass() const { return constants::me; }
  T getCharge() const { return -constants::q; }
  bool isMoved() const { return true; }
  bool isInjected() const { return false; }

  T getInitialNrParticles(const SizeVec &coord, const DeviceType &device,
                          const emcGrid<T, Dim> &potential) {
    return 1;
  }

  emcParticle<T> generateInitialParticle(const SizeVec &coord,
                                         const DeviceType &device,
                                         emcRNG &rng) {
    emcParticle<T> particle;
    particle.valley = 0;
    particle.subValley = 0;
    particle.k = {0, 0, 0};
    particle.energy = 0;
    particle.tau = 1;
    particle.grainTau = 1;
    return particle;
  }
};

/// exemplary non-moving particle type
template <class T, class DeviceType>
struct emcNonMovingType : public emcParticleType<T, DeviceType> {
  typedef typename DeviceType::SizeVec SizeVec;
  typedef typename DeviceType::ValueVec ValueVec;
  static const SizeType Dim = DeviceType::Dimension;

  std::string getName() const { return "NonMovingType"; }
  T getMass() const { return constants::me; }
  T getCharge() const { return -constants::q; }
  bool isMoved() const { return false; }
  bool isInjected() const { return false; }

  T getInitialNrParticles(const SizeVec &coord, const DeviceType &device,
                          const emcGrid<T, Dim> &potential) {
    return 2;
  }
};

template <class Handler>
void testNrParticlesAtEachCoordinate(Handler &handler, const DeviceType &device,
                                     SizeType idxPartType, NumType expValue) {
  emcGrid<NumType, Dim> nrPart(device.getGridExtent(), 0);
  handler.assignParticlesToMesh(idxPartType, nrPart);
  for (auto &val : nrPart) {
    EMCTEST_ASSERT(val == expValue);
  }
}

using MovingType = emcMovingType<NumType, DeviceType>;
using NonMovingType = emcNonMovingType<NumType, DeviceType>;

int main() {
  // initialize material, device, particle types and handler
  MaterialType testMaterial(11.8, 2329, 1.45e16, 9040, 1.14);
  DeviceType device(testMaterial, {1e-6, 1e-6, 1e-6}, {5e-7, 5e-7, 5e-7});
  device.addConstantDopingRegion({0, 0, 0}, {5e-7, 5e-7, 5e-7}, 1.45e16);

  MapIdxToParticleTypes idxToPartTypes;
  idxToPartTypes[0] = std::make_unique<MovingType>();
  idxToPartTypes[0]->addValley(
      std::make_unique<emcParabolicIsotropValley<NumType>>(1, constants::me,
                                                           1));
  idxToPartTypes[1] = std::make_unique<NonMovingType>();

  PMScheme pmScheme;
  FMMHandler handler(device, pmScheme, idxToPartTypes, 1, 1);
  // BasicHandler handler(device, pmScheme, idxToPartTypes, 1);

  // test particle initialization with NGP scheme -----------------------------
  emcGrid<NumType, Dim> potential(device.getGridExtent(), 0);
  handler.generateInitialParticles(potential);
  testNrParticlesAtEachCoordinate(handler, device, 0, 1);
  testNrParticlesAtEachCoordinate(handler, device, 1, 2);
  // handler.print("testFMMHandler", "Eq");

  // test drift scatter function ----------------------------------------------
  std::vector<emcGrid<NumType, Dim>> field;
  emcGrid<NumType, Dim> exField(device.getGridExtent(), 1e10);
  field.push_back(exField);
  field.push_back(exField);
  field.push_back(exField);

  for (int i = 0; i < 100; i++) {
    handler.driftScatterParticles(1e-15, field);
  }

  // test if non-moving particles moved
  testNrParticlesAtEachCoordinate(handler, device, 1, 2);

  // handler.printNrParticles();
  // handler.print("testFMMHandler", "Final");

  return 0;
}