#include <Components/FSimpleLeaf.hpp>
#include <Containers/FOctree.hpp>
#include <Core/FFmmAlgorithm.hpp>
#include <Kernels/P2P/FP2PParticleContainer.hpp>
#include <Kernels/P2P/FP2PParticleContainerIndexed.hpp>
#include <Kernels/Rotation/FRotationCell.hpp>
#include <Kernels/Rotation/FRotationKernel.hpp>

#include <Arranger/FBasicParticleContainerIndexedMover.hpp>
#include <Arranger/FOctreeArranger.hpp>
#include <ParticleHandler/FMM/FBasicParticleContainerMover.hpp>

#include <emcTestAsserts.hpp>
#include <emcUtil.hpp>

static const FSize P = 20;

// #define INDEXED_MOVER

using NumType = double;

#ifdef INDEXED_MOVER
using ContainerClass = FP2PParticleContainerIndexed<NumType>;
#else
using ContainerClass = FP2PParticleContainer<NumType, 2>;
#endif

using LeafClass = FSimpleLeaf<NumType, ContainerClass>;
using CellClass = FRotationCell<NumType, P>;
using OctreeClass = FOctree<NumType, CellClass, ContainerClass, LeafClass>;
using KernelClass = FRotationKernel<NumType, CellClass, ContainerClass, P>;
using FmmClass = FFmmAlgorithm<OctreeClass, CellClass, ContainerClass,
                               KernelClass, LeafClass>;

#ifdef INDEXED_MOVER
using MoverClass =
    FBasicParticleContainerIndexedMover<NumType, OctreeClass, ContainerClass>;
#else
using MoverClass =
    FBasicParticleContainerMover<NumType, OctreeClass, ContainerClass>;
#endif

using ArrangerClass =
    FOctreeArranger<NumType, OctreeClass, ContainerClass, MoverClass>;

// print leaf + particle characteristics
void printLeaf(LeafClass *leaf) {
  std::cout << "\tNew Leaf:\n";
  auto container = leaf->getSrc();
  FSize nrPart = container->getNbParticles();
  auto pos = container->getPositions();
  NumType *charge = container->getPhysicalValues();
#ifdef INDEXED_MOVER
  auto idx = container->getIndexes();
#else
  auto idx = container->getPhysicalValues(0, 1);
#endif
  NumType *potential = container->getPotentials();
  NumType *forcesX = container->getForcesX();
  NumType *forcesY = container->getForcesY();
  NumType *forcesZ = container->getForcesZ();
  for (int idxPart = 0; idxPart < nrPart; idxPart++) {
    std::cout << "\t\t" << idx[idxPart] << " ";
    std::cout << charge[idxPart] << " ";
    std::cout << pos[0][idxPart] << " ";
    std::cout << pos[1][idxPart] << " ";
    std::cout << pos[2][idxPart] << " ";
    std::cout << potential[idxPart] << "\n";
  }
}

/// this file tests the kernel and the rearrange function of the
/// FMM library (indexed mover + container has problems with
/// updating the idx, can be tested if line 18 is uncommented)
int main() {
  OctreeClass tree(5, 3, 1, {0.5, 0.5, 0.5});
  ArrangerClass arranger(&tree);
  KernelClass *kernels = new KernelClass(tree.getHeight(), tree.getBoxWidth(),
                                         tree.getBoxCenter());
  FmmClass algorithm(&tree, kernels);

  // insert 2 particles in different leafs
  NumType physValue = 1.602e-19;

#ifdef INDEXED_MOVER
  tree.insert({0.9, 0.9, 0.9}, 0, physValue);
  tree.insert({0.5, 0.5, 0.5}, 1, -physValue);
#else
  tree.insert({0.9, 0.9, 0.9}, physValue, 0);
  tree.insert({0.5, 0.5, 0.5}, -physValue, 1);
#endif

  algorithm.execute();

  // test calculation of force + potential
  // test potential calculation
  std::array<NumType, 3> dist = {0.4, 0.4, 0.4};
  auto expPot = physValue / norm(dist);
  auto expForce = physValue * physValue * 0.4 / std::pow(norm(dist), 3);

  std::cout << "Expected Force: " << expForce << "\n";
  std::cout << "Expected Potential: " << expPot << "\n";

  // test calculated values
  auto testCalcValues = [&](LeafClass *leaf) {
    ContainerClass *tgt = leaf->getTargets();
    FSize nrPart = tgt->getNbParticles();
#ifdef INDEXED_MOVER
    auto idx = tgt->getIndexes();
#else
    auto idx = tgt->getPhysicalValues(0, 1);
#endif
    NumType *potential = tgt->getPotentials();
    NumType *forcesX = tgt->getForcesX();
    NumType *forcesY = tgt->getForcesY();
    NumType *forcesZ = tgt->getForcesZ();
    for (int idxPart = 0; idxPart < nrPart; idxPart++) {
      EMCTEST_ASSERT_ISCLOSE(std::abs(potential[idxPart]), expPot, 1e-22);
      EMCTEST_ASSERT_ISCLOSE(std::abs(forcesX[idxPart]), expForce, 1e-40);
      EMCTEST_ASSERT_ISCLOSE(std::abs(forcesY[idxPart]), expForce, 1e-40);
      EMCTEST_ASSERT_ISCLOSE(std::abs(forcesZ[idxPart]), expForce, 1e-40);
    }
  };

  tree.forEachLeaf(testCalcValues);

  std::cout << "Initial Tree: \n";
  tree.forEachLeaf(printLeaf);

  // test moving (see output)
  // move particle with index 0 (to same leave as particle with idx 1)
  tree.forEachLeaf([&](LeafClass *leaf) {
    ContainerClass *src = leaf->getSrc();
    FSize nrPart = src->getNbParticles();
    auto pos = src->getPositions();
#ifdef INDEXED_MOVER
    auto idx = src->getIndexes();
#else
    auto *idx = src->getPhysicalValues(0, 1);
#endif
    for (int idxPart = 0; idxPart < nrPart; idxPart++) {
      if (idx[idxPart] == 0) {
        pos[0][idxPart] = 0.56;
        pos[1][idxPart] = 0.5;
        pos[2][idxPart] = 0.5;
      }
    }
  });

  std::cout << "After Moving Particle 0: (without rearrange) \n";
  tree.forEachLeaf(printLeaf);

  // rearrange tree
  arranger.rearrange();
  std::cout << "After Rearrange: \n";
  tree.forEachLeaf(printLeaf);

  // test if index is updated in the right way
  SizeType counter = 0;
  tree.forEachLeaf([&](LeafClass *leaf) {
    ContainerClass *src = leaf->getSrc();
    FSize nrPart = src->getNbParticles();
    auto pos = src->getPositions();
#ifdef INDEXED_MOVER
    auto idx = src->getIndexes();
#else
    auto *idx = src->getPhysicalValues(0, 1);
#endif
    for (int idxPart = 0; idxPart < nrPart; idxPart++) {
      if (idx[idxPart] == 1)
        counter++;
    }
  });

  EMCTEST_ASSERT(counter == 1);
  return 0;
}