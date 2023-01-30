#include <ParticleType/emcElectron.hpp>
#include <SurfaceScatterMechanisms/emcConstantSurfaceScatterMechanism.hpp>
#include <SurfaceScatterMechanisms/emcMomentumDependentSurfaceScatterMechanism.hpp>
#include <emcDevice.hpp>
#include <emcParticle.hpp>
#include <emcTestAsserts.hpp>

using NumType = double;
using MaterialType = emcMaterial<NumType>;
using DeviceType = emcDevice<NumType, 3>;
using ConstantSurface = emcConstantSurfaceScatterMechanism<NumType, DeviceType>;
using MovmentDependentSurface =
    emcMomentumDependentSurfaceScatterMechanism<NumType, DeviceType>;

/// small tests of surface scatter mechanisms
int main() {
  emcRNG rng(1);

  // create test - device
  std::array<NumType, 3> maxPos = {1e-6, 1e-6, 1e-6};
  std::array<NumType, 3> gridDelta = {1e-7, 1e-7, 1e-7};
  MaterialType silicon{1, 2329, 1, 9040, 1.15};
  DeviceType device = DeviceType(silicon, maxPos, gridDelta);
  emcElectron<NumType, DeviceType> partType;

  // add surface scatter mechanisms
  partType.setSurfaceScatterMechanism(
      emcBoundaryPos::YMAX, std::make_unique<ConstantSurface>(1, maxPos));
  partType.setSurfaceScatterMechanism(
      emcBoundaryPos::XMIN, std::make_unique<ConstantSurface>(0, maxPos));
  partType.setSurfaceScatterMechanism(
      emcBoundaryPos::XMAX,
      std::make_unique<MovmentDependentSurface>(3e-10, maxPos));

  // create test particle
  emcParticle<NumType> testParticle;
  std::array<NumType, 3> pos, initK, initPos;

  for (SizeType i = 0; i < 1000; i++) {
    // scattering at YMAX -> smooth surface
    pos = {1e-6, 1.1e-6, 1e-6};
    testParticle.k = {1, 1, 1};
    initK = testParticle.k;
    initPos = pos;
    partType.scatterParticleAtBoundary(emcBoundaryPos::YMAX, testParticle,
                                       device, pos, rng);
    EMCTEST_ASSERT(testParticle.k[1] < 0);
    EMCTEST_ASSERT(testParticle.k[0] == initK[0]);
    EMCTEST_ASSERT(testParticle.k[2] == initK[2]);
    EMCTEST_ASSERT(pos[1] == (2 * maxPos[1] - initPos[1]));
    EMCTEST_ASSERT(pos[0] == initPos[0]);
    EMCTEST_ASSERT(pos[2] == initPos[2]);
    EMCTEST_ASSERT_ISCLOSE(norm(testParticle.k), norm(initK),
                           norm(initK) / 1e10);

    // scattering at XMIN -> rough surface
    initPos = {-0.1e-6, 1e-6, 1e-6};
    pos = initPos;
    testParticle.k = {-1, 1, 1};
    initK = testParticle.k;
    partType.scatterParticleAtBoundary(emcBoundaryPos::XMIN, testParticle,
                                       device, pos, rng);
    EMCTEST_ASSERT(testParticle.k[0] > 0);
    EMCTEST_ASSERT(testParticle.k[1] != initK[1]);
    EMCTEST_ASSERT(testParticle.k[2] != initK[2]);
    EMCTEST_ASSERT(pos[0] == -initPos[0]);
    EMCTEST_ASSERT(pos[1] == initPos[1]);
    EMCTEST_ASSERT(pos[2] == initPos[2]);
    EMCTEST_ASSERT_ISCLOSE(norm(testParticle.k), norm(initK),
                           norm(initK) / 1e10);

    // scattering at XMAX -> momentum dependent surface scattering
    initPos = {1.1e-6, 5e-7, 5e-7};
    pos = initPos;
    testParticle.k = {1e10, 1e10, 1e10};
    initK = testParticle.k;
    partType.scatterParticleAtBoundary(emcBoundaryPos::XMAX, testParticle,
                                       device, pos, rng);
    EMCTEST_ASSERT(testParticle.k[0] < 0);
    EMCTEST_ASSERT(testParticle.k[1] != initK[1]);
    EMCTEST_ASSERT(testParticle.k[2] != initK[2]);
    EMCTEST_ASSERT(pos[0] == 2 * maxPos[0] - initPos[0]);
    EMCTEST_ASSERT(pos[1] == initPos[1]);
    EMCTEST_ASSERT(pos[2] == initPos[2]);
    EMCTEST_ASSERT_ISCLOSE(norm(testParticle.k), norm(initK),
                           norm(initK) / 1e10);

    // scatterig at XMIN -> smooth surface
    // (but also adapting position and wave-vector in other dimensions)
    initPos = {-0.1, -0.1, -1};
    pos = initPos;
    initK = {-1, -1, -1};
    testParticle.k = initK;
    partType.scatterParticleAtBoundary(emcBoundaryPos::XMIN, testParticle,
                                       device, pos, rng);
    EMCTEST_ASSERT(testParticle.k[1] > 0);
    EMCTEST_ASSERT(testParticle.k[0] > 0);
    EMCTEST_ASSERT(testParticle.k[2] > 0);
    EMCTEST_ASSERT(pos[1] == -initPos[1]);
    EMCTEST_ASSERT(pos[0] == -initPos[0]);
    EMCTEST_ASSERT(pos[2] == -initPos[2]);
    EMCTEST_ASSERT_ISCLOSE(norm(testParticle.k), norm(initK),
                           norm(initK) / 1e10);

    // scattering at YMIN -> elastic scattering (no added mechanism)
    initPos = {5e-7, -0.1, 5e-7};
    pos = initPos;
    initK = {1, -1, 1};
    testParticle.k = initK;
    partType.scatterParticleAtBoundary(emcBoundaryPos::YMIN, testParticle,
                                       device, pos, rng);
    EMCTEST_ASSERT(testParticle.k[1] > 0);
    EMCTEST_ASSERT(testParticle.k[2] == initK[2]);
    EMCTEST_ASSERT(testParticle.k[0] == initK[0]);
    EMCTEST_ASSERT(pos[1] == -initPos[1]);
    EMCTEST_ASSERT(pos[0] == initPos[0]);
    EMCTEST_ASSERT(pos[2] == initPos[2]);
  }

  return 0;
}