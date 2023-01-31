// adapted on 10.03.2023 by Gollner Laura
// added FReal template parameter for FAbstractMover (line 17)
// and added FParticleType to input parameter (line 37)
// (C) from scalFMM main branch
// include/Arranger/FBasicParticleContainerMover.hpp

#ifndef FBASICPARTICLECONTAINERMOVER_HPP
#define FBASICPARTICLECONTAINERMOVER_HPP

#include <Arranger/FAbstractMover.hpp>

/**
 * This class should be use with the octree arrange to move particles
 * that are stored in a FBasicParticleContainer.
 */
template <class FReal, class OctreeClass, class ContainerClass>
class FBasicParticleContainerMover
    : public FAbstractMover<FReal, OctreeClass, ContainerClass> {
private:
  ContainerClass toStoreRemovedParts;

public:
  FBasicParticleContainerMover() {}

  virtual ~FBasicParticleContainerMover() {}

  /** To get the position of the particle at idx idxPart in leaf lf */
  void getParticlePosition(ContainerClass *lf, const FSize idxPart,
                           FPoint<FReal> *particlePos) {
    (*particlePos) = FPoint<FReal>(lf->getPositions()[0][idxPart],
                                   lf->getPositions()[1][idxPart],
                                   lf->getPositions()[2][idxPart]);
  }

  /** Remove a particle but keep it to reinsert it later*/
  void removeFromLeafAndKeep(ContainerClass *lf,
                             const FPoint<FReal> &particlePos,
                             const FSize idxPart, FParticleType /*type*/) {
    std::array<typename ContainerClass::AttributesClass,
               ContainerClass::NbAttributes>
        particleValues;
    for (int idxAttr = 0; idxAttr < ContainerClass::NbAttributes; ++idxAttr) {
      particleValues[idxAttr] = lf->getAttribute(idxAttr)[idxPart];
    }

    toStoreRemovedParts.push(particlePos, particleValues);

    lf->removeParticles(&idxPart, 1);
  }

  /** Reinsert the previously saved particles */
  void insertAllParticles(OctreeClass *tree) {
    std::array<typename ContainerClass::AttributesClass,
               ContainerClass::NbAttributes>
        particleValues;

    for (FSize idxToInsert = 0;
         idxToInsert < toStoreRemovedParts.getNbParticles(); ++idxToInsert) {
      for (int idxAttr = 0; idxAttr < ContainerClass::NbAttributes; ++idxAttr) {
        particleValues[idxAttr] =
            toStoreRemovedParts.getAttribute(idxAttr)[idxToInsert];
      }
      const FPoint<FReal> particlePos(
          toStoreRemovedParts.getPositions()[0][idxToInsert],
          toStoreRemovedParts.getPositions()[1][idxToInsert],
          toStoreRemovedParts.getPositions()[2][idxToInsert]);

      tree->insert(particlePos, particleValues);
    }

    toStoreRemovedParts.clear();
  }
};

#endif // FBASICPARTICLECONTAINERMOVER_HPP
