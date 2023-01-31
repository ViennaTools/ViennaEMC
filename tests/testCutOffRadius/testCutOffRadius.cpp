#include <Components/FSimpleLeaf.hpp>
#include <Containers/FOctree.hpp>
#include <Core/FFmmAlgorithm.hpp>
#include <Core/FFmmAlgorithmPeriodic.hpp>
#include <Kernels/P2P/FP2PParticleContainer.hpp>
#include <Kernels/Rotation/FRotationCell.hpp>

#include <ParticleHandler/FMM/FBasicParticleContainerMover.hpp>
// include kernel that uses cutoff approach
#include <ParticleHandler/FMM/FRotationKernel.hpp>

#include <emcConstants.hpp>
#include <emcTestAsserts.hpp>
#include <emcUtil.hpp>

static const FSize P = 20;

using NumType = double;
using ContainerClass = FP2PParticleContainer<NumType, 2>;
using LeafClass = FSimpleLeaf<NumType, ContainerClass>;
using CellClass = FRotationCell<NumType, P>;
using OctreeClass = FOctree<NumType, CellClass, ContainerClass, LeafClass>;
using KernelClass = FRotationKernel<NumType, CellClass, ContainerClass, P>;
using FmmClass = FFmmAlgorithm<OctreeClass, CellClass, ContainerClass,
                               KernelClass, LeafClass>;

// print leaf + particle characteristics
void printLeaf(LeafClass *leaf) {
  auto container = leaf->getSrc();
  FSize nrPart = container->getNbParticles();
  auto pos = container->getPositions();
  NumType *charge = container->getPhysicalValues();
  auto idx = container->getPhysicalValues(0, 1);
  NumType *potential = container->getPotentials();
  NumType *forcesX = container->getForcesX();
  NumType *forcesY = container->getForcesY();
  NumType *forcesZ = container->getForcesZ();
  for (int idxPart = 0; idxPart < nrPart; idxPart++) {
    std::cout << "\t index " << idx[idxPart] << "\n";
    std::cout << "\t\t potential: " << potential[idxPart] << "\n";
    std::cout << "\t\t force: " << forcesX[idxPart] << " ";
    std::cout << forcesY[idxPart] << " ";
    std::cout << forcesZ[idxPart] << "\n";
  }
}

NumType calcExpectedPotential(const std::array<NumType, 3> &pos1,
                              const std::array<NumType, 3> &pos2,
                              NumType physValue = 1.) {
  auto distNorm = norm(subtract(pos1, pos2));
  if (distNorm < 1e-9)
    distNorm = 1e-9;
  return physValue / distNorm;
}

std::array<NumType, 3> calcExpectedForce(const std::array<NumType, 3> &pos1,
                                         const std::array<NumType, 3> &pos2,
                                         NumType physValue = 1.) {
  auto dist = subtract(pos1, pos2);
  auto distNorm = norm(dist);
  if (distNorm < 1e-9) {
    dist = scale(dist, 1e-9 / distNorm);
    distNorm = 1e-9;
  }
  return scale(dist, physValue * physValue / (std::pow(distNorm, 3)));
}

void testParticles(const std::array<NumType, 3> &pos0,
                   const std::array<NumType, 3> &pos1) {
  // create simulation space (octree)
  OctreeClass tree(5, 3, 100e-9, {50e-9, 50e-9, 50e-9});
  KernelClass *kernels = new KernelClass(tree.getHeight(), tree.getBoxWidth(),
                                         tree.getBoxCenter());
  FmmClass algorithm(&tree, kernels);

  // add two particles at given positions
  NumType physValue = 1;
  tree.insert({pos0[0], pos0[1], pos0[2]}, physValue, 0);
  tree.insert({pos1[0], pos1[1], pos1[2]}, physValue, 1);
  algorithm.execute();

  // calculate expected force + potential with cutoff approach
  auto expPot = calcExpectedPotential(pos0, pos1, physValue);
  auto expForce = calcExpectedForce(pos0, pos1, physValue);

  // compare results to expected values
  auto testCalcValues = [&](LeafClass *leaf) {
    auto container = leaf->getSrc();
    FSize nrPart = container->getNbParticles();
    auto idx = container->getPhysicalValues(0, 1);
    NumType *potential = container->getPotentials();
    NumType *forcesX = container->getForcesX();
    NumType *forcesY = container->getForcesY();
    NumType *forcesZ = container->getForcesZ();
    for (int idxPart = 0; idxPart < nrPart; idxPart++) {
      EMCTEST_ASSERT_ISCLOSE(potential[idxPart], expPot,
                             std::fabs(expPot) / 1e10);
      EMCTEST_ASSERT_ISCLOSE(std::fabs(forcesX[idxPart]),
                             std::fabs(expForce[0]),
                             std::fabs(forcesX[idxPart]) / 1e10);
      EMCTEST_ASSERT_ISCLOSE(std::fabs(forcesY[idxPart]),
                             std::fabs(expForce[1]),
                             std::fabs(forcesY[idxPart]) / 1e10);
      EMCTEST_ASSERT_ISCLOSE(std::fabs(forcesZ[idxPart]),
                             std::fabs(expForce[2]),
                             std::fabs(forcesZ[idxPart]) / 1e10);
    }
  };
  tree.forEachLeaf(testCalcValues);
  // tree.forEachLeaf(&printLeaf);
}

/// this file tests the implemented cutoff radius of the Coulomb force at 1nm.
int main() {
  testParticles({50e-9, 50e-9, 50e-9}, {50.1e-9, 50.1e-9, 50.1e-9});
  testParticles({50e-9, 50e-9, 50e-9}, {50.5e-9, 50.1e-9, 50.1e-9});
  testParticles({50e-9, 50e-9, 50e-9}, {52e-9, 53e-9, 54e-9});
  return 0;
}