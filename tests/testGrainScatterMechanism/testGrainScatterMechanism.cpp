#include <emcGrainScatterMechanism.hpp>
#include <emcParticle.hpp>
#include <emcTestAsserts.hpp>
#include <emcUtil.hpp>

using NumType = double;

/// tests of grain scatter mechanism
/// \file emcGrainScatterMechanisms.hpp
int main() {
  emcRNG rng(time(NULL));
  emcParticle<NumType> testPart;

  std::vector<std::array<NumType, 3>> initialWaveVectors;
  initialWaveVectors.push_back({1, 0, 0});
  initialWaveVectors.push_back({1, 1, 1});
  initialWaveVectors.push_back({-5, 6, 7});
  initialWaveVectors.push_back({-1, -2, -3});

  for (auto &initK : initialWaveVectors) {
    for (SizeType i = 0; i < 1000; i++) {
      NumType initNorm = norm(initK);
      testPart.k = initK;

      emcGrainScatterMechanism<NumType> grainScatterMech(0.5, 1e-13);
      grainScatterMech.scatterTransmittingParticle(testPart, rng);
      EMCTEST_ASSERT(innerProduct(testPart.k, initK) >= 0);
      EMCTEST_ASSERT_ISCLOSE(norm(testPart.k), initNorm, initNorm / 1e10);

      testPart.k = initK;
      grainScatterMech.scatterReflectingParticle(testPart, rng);
      EMCTEST_ASSERT(innerProduct(testPart.k, initK) <= 0);
      EMCTEST_ASSERT_ISCLOSE(norm(testPart.k), initNorm, initNorm / 1e10);
    }
  }

  return 0;
}