#include <ParticleType/emcElectron.hpp>
#include <ValleyTypes/emcNonParabolicAnistropValley.hpp>
#include <ValleyTypes/emcNonParabolicIsotropValley.hpp>
#include <emcDevice.hpp>
#include <emcParticle.hpp>
#include <emcParticleDrift.hpp>
#include <emcTestAsserts.hpp>
#include <emcUtil.hpp>

#include <algorithm>
#include <math.h>
#include <numeric>

/**
  Tests if the movement of the particle under a given force works as expected.
  \file emcParticleDrift.hpp
*/

static const SizeType Dim = 3;
static const double dT = 1e-12;

using NumType = double;
using MaterialType = emcMaterial<NumType>;
using ValleyType = emcNonParabolicAnisotropValley<NumType>;
using DeviceType = emcDevice<NumType, Dim>;

// initialized particle
void initParticle(emcParticle<NumType> &particle, DeviceType &device,
                  SizeType subValley) {
  particle.valley = 0;
  particle.subValley = subValley;
  particle.energy = 0;
  particle.k = {0, 0, 0};
}

// prints particle (used for debugging)
void printParticle(emcParticle<NumType> &particle,
                   std::array<NumType, Dim> &position) {
  std::cout << "\tk = " << particle.k << "\n";
  std::cout << "\tenergy = " << particle.energy << "\n";
  std::cout << "\tposition = " << position << "\n";
}

// function that calculates expected valley
// expects that force is in ellipse coordinate system
void calcExpectedValues(const std::array<NumType, Dim> &force,
                        const std::array<NumType, Dim> &effMass,
                        std::array<NumType, 3> &expK,
                        std::array<NumType, Dim> &expPos) {
  NumType massCond = std::accumulate(
      effMass.begin(), effMass.end(), 0,
      [](NumType acc, const NumType &val) { return acc + 1. / val; });
  massCond = 3. / massCond;
  std::transform(force.begin(), force.end(), effMass.begin(), expK.begin(),
                 [&](const NumType &f, const NumType &m) {
                   return f * dT * std::sqrt(massCond / m) / constants::hbar;
                 });
  // assumes previous wave-vec is 0
  std::transform(expK.begin(), expK.end(), effMass.begin(), expPos.begin(),
                 [&](const NumType &kNew, const NumType &m) {
                   return kNew / 2. * dT / std::sqrt(m * massCond) *
                          constants::hbar;
                 });
}

// transforms from device coordinate system to given coordinate system
std::array<NumType, Dim> transformToCS(const std::array<NumType, Dim> &dir1,
                                       const std::array<NumType, Dim> &dir2,
                                       const std::array<NumType, Dim> &dir3,
                                       const std::array<NumType, Dim> &vec) {
  std::array<NumType, Dim> result;
  result[0] = dir1[0] * vec[0] + dir1[1] * vec[1] + dir1[2] * vec[2];
  result[0] /= norm(dir1);
  result[1] = dir2[0] * vec[0] + dir2[1] * vec[1] + dir2[2] * vec[2];
  result[1] /= norm(dir2);
  result[2] = dir3[0] * vec[0] + dir3[1] * vec[1] + dir3[2] * vec[2];
  result[2] /= norm(dir3);
  return result;
}

// transforms from given coordinate system back to device coordinate system
std::array<NumType, Dim> transformFromCS(const std::array<NumType, Dim> &dir1,
                                         const std::array<NumType, Dim> &dir2,
                                         const std::array<NumType, Dim> &dir3,
                                         const std::array<NumType, Dim> &vec) {
  std::array<NumType, Dim> result;
  result[0] = dir1[0] * vec[0] / norm(dir1) + dir2[0] * vec[1] / norm(dir2) +
              dir3[0] * vec[2] / norm(dir3);
  result[1] = dir1[1] * vec[0] / norm(dir1) + dir2[1] * vec[1] / norm(dir2) +
              dir3[1] * vec[2] / norm(dir3);
  result[2] = dir1[2] * vec[0] / norm(dir1) + dir2[2] * vec[1] / norm(dir2) +
              dir3[2] * vec[2] / norm(dir3);
  return result;
}

template <class ValleyType> void setValleyCoordSystems(ValleyType &valley) {
  valley.setSubValleyEllipseCoordSystem(0, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
  valley.setSubValleyEllipseCoordSystem(1, {0, 1, 0}, {1, 0, 0}, {0, 0, 1});
  valley.setSubValleyEllipseCoordSystem(2, {0, 0, 1}, {0, 1, 0}, {1, 0, 0});
  valley.setSubValleyEllipseCoordSystem(3, {-1, 1, 1}, {1, 1, 0}, {1, -1, 2});
  valley.setSubValleyEllipseCoordSystem(4, {1, 1, 1}, {-1, 1, 0}, {-1, -1, 2});
  valley.setSubValleyEllipseCoordSystem(5, {-1, -1, 1}, {1, 0, 1}, {-1, 2, 1});
}

int main() {
  emcElectron<NumType, DeviceType> electronType;
  emcRNG rng;

  // test: isotropic material --------------------------------------------
  // regardless of subvalley coordinate system results should be same
  std::array<NumType, 3> effMass = {1, 1, 1};
  emcNonParabolicIsotropValley<NumType> isoValley(1, constants::me, 6, 0.5);

  // create device
  MaterialType isoMaterial{1, 1, 1, 1, 1};
  DeviceType isoDevice{isoMaterial, {1e-6, 1e-6, 1e-6}, {1e-8, 1e-8, 1e-8}};

  // set initial particle characteristics
  emcParticle<NumType> testParticle;
  std::array<NumType, Dim> testPos = {0, 0, 0};
  std::array<NumType, Dim> force = {1e-19, 5e-19, 1e-20};

  // calculates expected values after drift
  std::array<NumType, Dim> expectedK, expectedPos;
  calcExpectedValues(force, effMass, expectedK, expectedPos);

  // compare calculated + expected values
  for (SizeType subValley = 0; subValley < isoValley.getDegeneracyFactor();
       subValley++) {
    std::cout << "Testing subValley " << subValley << " ...\n";
    testPos = {0, 0, 0};
    initParticle(testParticle, isoDevice, subValley);
    drift(dT, testParticle, &isoValley, testPos, force);
    handleParticleAtBoundary(testParticle, testPos, &electronType, isoDevice,
                             rng);

    EMCTEST_ASSERT_ISCLOSE(expectedPos[0], testPos[0], 1e-9);
    EMCTEST_ASSERT_ISCLOSE(expectedPos[1], testPos[1], 1e-9);
    EMCTEST_ASSERT_ISCLOSE(expectedPos[2], testPos[2], 1e-9);
    EMCTEST_ASSERT_ISCLOSE(expectedK[0], testParticle.k[0], 1e-9);
    EMCTEST_ASSERT_ISCLOSE(expectedK[1], testParticle.k[1], 1e-9);
    EMCTEST_ASSERT_ISCLOSE(expectedK[2], testParticle.k[2], 1e-9);
  }

  // test: anisotropic material ------------------------------------------
  // results depend on coordinate system of subvalley

  // create device
  MaterialType anisoMaterial{1, 1, 1, 1, 1};
  effMass = {0.1, 0.5, 1};
  emcNonParabolicAnisotropValley<NumType> anisoValley(effMass, constants::me, 6,
                                                      0.5);
  setValleyCoordSystems(anisoValley);
  DeviceType anisoDevice{anisoMaterial, {1e-6, 1e-6, 1e-6}, {1e-8, 1e-8, 1e-8}};

  force = {1e-17, 1e-17, 1e-17};
  std::array<NumType, 3> dir1, dir2, dir3;
  for (SizeType subValley = 0; subValley < isoValley.getDegeneracyFactor();
       subValley++) {
    std::cout << "Testing subValley " << subValley << " ...\n";
    testPos = {1e-7, 1e-7, 1e-7};
    initParticle(testParticle, anisoDevice, subValley);
    switch (subValley) { // set subValley coordinate system
    case 0:
      dir1 = {1, 0, 0};
      dir2 = {0, 1, 0};
      dir3 = {0, 0, 1};
      break;
    case 1:
      dir1 = {0, 1, 0};
      dir2 = {1, 0, 0};
      dir3 = {0, 0, 1};
      break;
    case 2:
      dir1 = {0, 0, 1};
      dir2 = {0, 1, 0};
      dir3 = {1, 0, 0};
      break;
    case 3:
      dir1 = {-1, 1, 1};
      dir2 = {1, 1, 0};
      dir3 = {1, -1, 2};
      break;
    case 4:
      dir1 = {1, 1, 1};
      dir2 = {-1, 1, 0};
      dir3 = {-1, -1, 2};
      break;
    case 5:
      dir1 = {-1, -1, 1};
      dir2 = {1, 0, 1};
      dir3 = {-1, 2, 1};
      break;
    }

    drift(dT, testParticle, &anisoValley, testPos, force);
    handleParticleAtBoundary(testParticle, testPos, &electronType, anisoDevice,
                             rng);
    // printParticle(testParticle, testPos);

    // calculate expected values
    auto tmpForce = transformToCS(dir1, dir2, dir3, force);
    calcExpectedValues(tmpForce, effMass, expectedK, expectedPos);
    expectedK = transformFromCS(dir1, dir2, dir3, expectedK);
    expectedPos = transformFromCS(dir1, dir2, dir3, expectedPos);
    expectedPos = add(expectedPos, {1e-7, 1e-7, 1e-7});
    NumType expectedEnergy = anisoValley.getEnergy(expectedK);

    // compare expected and calculated values
    EMCTEST_ASSERT_ISCLOSE(expectedK[0], testParticle.k[0], 1e-5);
    EMCTEST_ASSERT_ISCLOSE(expectedK[1], testParticle.k[1], 1e-5);
    EMCTEST_ASSERT_ISCLOSE(expectedK[2], testParticle.k[2], 1e-5);
    EMCTEST_ASSERT_ISCLOSE(expectedPos[0], testPos[0], 1e-10);
    EMCTEST_ASSERT_ISCLOSE(expectedPos[1], testPos[1], 1e-10);
    EMCTEST_ASSERT_ISCLOSE(expectedPos[2], testPos[2], 1e-10);
    EMCTEST_ASSERT_ISCLOSE(expectedEnergy, testParticle.energy, 1e-10);
  }

  return 0;
}