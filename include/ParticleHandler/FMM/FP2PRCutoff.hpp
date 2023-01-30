// adapted on 17.03.2023 by Gollner Laura
// adapted to use minimal cut off radius
// (C) from scalFMM main branch include/Kernels/P2P/FP2PR.hpp

#ifndef FP2PR_CUTOFF_HPP
#define FP2PR_CUTOFF_HPP

#include "Utils/FGlobal.hpp"
#include "Utils/FMath.hpp"

/**
 * @brief The FP2PRCutOff namespace.
 * Uses 1 / r - Potential and adapts the potential + force if two particles
 * are too close to each other to the potential / force at a specific cutoff
 * minimal distance.
 */
namespace FP2PRCutoff {
template <class FReal>
inline void MutualParticles(
    const FReal targetX, const FReal targetY, const FReal targetZ,
    const FReal targetPhysicalValue, FReal *targetForceX, FReal *targetForceY,
    FReal *targetForceZ, FReal *targetPotential, const FReal sourceX,
    const FReal sourceY, const FReal sourceZ, const FReal sourcePhysicalValue,
    FReal *sourceForceX, FReal *sourceForceY, FReal *sourceForceZ,
    FReal *sourcePotential) {
  FReal dx = targetX - sourceX;
  FReal dy = targetY - sourceY;
  FReal dz = targetZ - sourceZ;

  FReal inv_square_distance = FReal(1.0) / (dx * dx + dy * dy + dz * dz);
  FReal inv_distance = FMath::Sqrt(inv_square_distance);

  // ################################################################
  // if particles are too close adapt their distance to cutOffDistance
  // of 1 nm (calc force and potential with that distance)
  if (inv_distance > 1e9) {
    auto factor = inv_distance * 1e-9;
    dx *= factor;
    dy *= factor;
    dz *= factor;
    inv_distance = 1e9;
    inv_square_distance = 1e18;
  } // ################################################################

  inv_square_distance *= inv_distance;
  inv_square_distance *= targetPhysicalValue * sourcePhysicalValue;

  dx *= -inv_square_distance;
  dy *= -inv_square_distance;
  dz *= -inv_square_distance;

  *targetForceX += dx;
  *targetForceY += dy;
  *targetForceZ += dz;
  *targetPotential += (inv_distance * sourcePhysicalValue);

  *sourceForceX -= dx;
  *sourceForceY -= dy;
  *sourceForceZ -= dz;
  *sourcePotential += (inv_distance * targetPhysicalValue);
}

template <class FReal>
inline void NonMutualParticles(
    const FReal targetX, const FReal targetY, const FReal targetZ,
    const FReal targetPhysicalValue, FReal *targetForceX, FReal *targetForceY,
    FReal *targetForceZ, FReal *targetPotential, const FReal sourceX,
    const FReal sourceY, const FReal sourceZ, const FReal sourcePhysicalValue) {
  FReal dx = targetX - sourceX;
  FReal dy = targetY - sourceY;
  FReal dz = targetZ - sourceZ;

  FReal inv_square_distance = FReal(1.0) / (dx * dx + dy * dy + dz * dz);
  FReal inv_distance = FMath::Sqrt(inv_square_distance);

  // ################################################################
  // if particles are too close adapt their distance to cutOffDistance
  // of 1 nm (calc force and potential with that distance)
  if (inv_distance > 1e9) {
    auto factor = inv_distance * 1e-9;
    dx *= factor;
    dy *= factor;
    dz *= factor;
    inv_distance = 1e9;
    inv_square_distance = 1e18;
  } // ################################################################

  inv_square_distance *= inv_distance;
  inv_square_distance *= targetPhysicalValue * sourcePhysicalValue;

  // d/dx(1/|x-y|)=-(x-y)/r^3
  dx *= -inv_square_distance;
  dy *= -inv_square_distance;
  dz *= -inv_square_distance;

  *targetForceX += dx;
  *targetForceY += dy;
  *targetForceZ += dz;
  *targetPotential += (inv_distance * sourcePhysicalValue);
}

template <class FReal, class ContainerClass, class ComputeClass,
          int NbFRealInComputeClass>
static void GenericFullMutual(ContainerClass *const FRestrict inTargets,
                              ContainerClass *const inNeighbors[],
                              const int limiteNeighbors) {

  const FSize nbParticlesTargets = inTargets->getNbParticles();
  const FReal *const targetsPhysicalValues = inTargets->getPhysicalValues();
  const FReal *const targetsX = inTargets->getPositions()[0];
  const FReal *const targetsY = inTargets->getPositions()[1];
  const FReal *const targetsZ = inTargets->getPositions()[2];
  FReal *const targetsForcesX = inTargets->getForcesX();
  FReal *const targetsForcesY = inTargets->getForcesY();
  FReal *const targetsForcesZ = inTargets->getForcesZ();
  FReal *const targetsPotentials = inTargets->getPotentials();

  const ComputeClass mOne = ComputeClass(1);

  for (FSize idxNeighbors = 0; idxNeighbors < limiteNeighbors; ++idxNeighbors) {
    if (inNeighbors[idxNeighbors]) {
      const FSize nbParticlesSources =
          inNeighbors[idxNeighbors]->getNbParticles();
      const FReal *const sourcesPhysicalValues =
          inNeighbors[idxNeighbors]->getPhysicalValues();
      const FReal *const sourcesX =
          inNeighbors[idxNeighbors]->getPositions()[0];
      const FReal *const sourcesY =
          inNeighbors[idxNeighbors]->getPositions()[1];
      const FReal *const sourcesZ =
          inNeighbors[idxNeighbors]->getPositions()[2];
      FReal *const sourcesForcesX = inNeighbors[idxNeighbors]->getForcesX();
      FReal *const sourcesForcesY = inNeighbors[idxNeighbors]->getForcesY();
      FReal *const sourcesForcesZ = inNeighbors[idxNeighbors]->getForcesZ();
      FReal *const sourcesPotentials =
          inNeighbors[idxNeighbors]->getPotentials();

      for (FSize idxTarget = 0; idxTarget < nbParticlesTargets; ++idxTarget) {
        FSize idxSource = 0;
        {
          const FSize nbVectorizedInteractions =
              (nbParticlesSources / NbFRealInComputeClass) *
              NbFRealInComputeClass;
          const ComputeClass tx = ComputeClass(targetsX[idxTarget]);
          const ComputeClass ty = ComputeClass(targetsY[idxTarget]);
          const ComputeClass tz = ComputeClass(targetsZ[idxTarget]);
          const ComputeClass tv =
              ComputeClass(targetsPhysicalValues[idxTarget]);
          ComputeClass tfx = ComputeClass::GetZero();
          ComputeClass tfy = ComputeClass::GetZero();
          ComputeClass tfz = ComputeClass::GetZero();
          ComputeClass tpo = ComputeClass::GetZero();

          for (; idxSource < nbVectorizedInteractions;
               idxSource += NbFRealInComputeClass) {
            ComputeClass dx = tx - ComputeClass(&sourcesX[idxSource]);
            ComputeClass dy = ty - ComputeClass(&sourcesY[idxSource]);
            ComputeClass dz = tz - ComputeClass(&sourcesZ[idxSource]);

            ComputeClass inv_square_distance =
                mOne / (dx * dx + dy * dy + dz * dz);

            /*const*/ ComputeClass inv_distance = inv_square_distance.sqrt();

            // #########################################################
            // if particles are too close adapt their distance to cutOffDistance
            // of 1 nm (calc force and potential with that distance)
            auto mask_ab = (inv_distance > ComputeClass(1e9));
            if (mask_ab.isAllTrue()) {
              inv_square_distance = ComputeClass(1e18);
              auto factor = inv_distance * ComputeClass(1e-9);
              dx *= factor;
              dy *= factor;
              dz *= factor;
              inv_distance = ComputeClass(1e9);
            } // #########################################################

            inv_square_distance *= inv_distance;
            inv_square_distance *=
                tv * ComputeClass(&sourcesPhysicalValues[idxSource]);

            dx *= -inv_square_distance;
            dy *= -inv_square_distance;
            dz *= -inv_square_distance;

            tfx += dx;
            tfy += dy;
            tfz += dz;
            tpo +=
                inv_distance * ComputeClass(&sourcesPhysicalValues[idxSource]);

            (ComputeClass(&sourcesForcesX[idxSource]) - dx)
                .storeInArray(&sourcesForcesX[idxSource]);
            (ComputeClass(&sourcesForcesY[idxSource]) - dy)
                .storeInArray(&sourcesForcesY[idxSource]);
            (ComputeClass(&sourcesForcesZ[idxSource]) - dz)
                .storeInArray(&sourcesForcesZ[idxSource]);
            (ComputeClass(&sourcesPotentials[idxSource]) + inv_distance * tv)
                .storeInArray(&sourcesPotentials[idxSource]);
          }

          targetsForcesX[idxTarget] += tfx.horizontalSum();
          targetsForcesY[idxTarget] += tfy.horizontalSum();
          targetsForcesZ[idxTarget] += tfz.horizontalSum();
          targetsPotentials[idxTarget] += tpo.horizontalSum();
        }
        {
          const FReal tx = FReal(targetsX[idxTarget]);
          const FReal ty = FReal(targetsY[idxTarget]);
          const FReal tz = FReal(targetsZ[idxTarget]);
          const FReal tv = FReal(targetsPhysicalValues[idxTarget]);
          FReal tfx = FReal(0.);
          FReal tfy = FReal(0.);
          FReal tfz = FReal(0.);
          FReal tpo = FReal(0.);

          for (; idxSource < nbParticlesSources; idxSource += 1) {
            FReal dx = tx - FReal(sourcesX[idxSource]);
            FReal dy = ty - FReal(sourcesY[idxSource]);
            FReal dz = tz - FReal(sourcesZ[idxSource]);

            FReal inv_square_distance =
                FReal(1) / (dx * dx + dy * dy + dz * dz);

            /*const*/ FReal inv_distance = FMath::Sqrt(inv_square_distance);

            // #########################################################
            // if particles are too close adapt their distance to cutOffDistance
            // of 1 nm (calc force and potential with that distance)
            if (inv_distance > 1e9) {
              inv_square_distance = 1e18;
              auto factor = inv_distance * 1e-9;
              dx *= factor;
              dy *= factor;
              dz *= factor;
              inv_distance = 1e9;
            } // #########################################################

            inv_square_distance *= inv_distance;

            inv_square_distance *= tv * FReal(sourcesPhysicalValues[idxSource]);

            dx *= -inv_square_distance;
            dy *= -inv_square_distance;
            dz *= -inv_square_distance;

            tfx += dx;
            tfy += dy;
            tfz += dz;
            tpo += inv_distance * FReal(sourcesPhysicalValues[idxSource]);

            sourcesForcesX[idxSource] -= dx;
            sourcesForcesY[idxSource] -= dy;
            sourcesForcesZ[idxSource] -= dz;
            sourcesPotentials[idxSource] += inv_distance * tv;
          }

          targetsForcesX[idxTarget] += tfx;
          targetsForcesY[idxTarget] += tfy;
          targetsForcesZ[idxTarget] += tfz;
          targetsPotentials[idxTarget] += tpo;
        }
      }
    }
  }
}

template <class FReal, class ContainerClass, class ComputeClass,
          int NbFRealInComputeClass>
static void GenericInner(ContainerClass *const FRestrict inTargets) {

  const FSize nbParticlesTargets = inTargets->getNbParticles();
  const FReal *const targetsPhysicalValues = inTargets->getPhysicalValues();
  const FReal *const targetsX = inTargets->getPositions()[0];
  const FReal *const targetsY = inTargets->getPositions()[1];
  const FReal *const targetsZ = inTargets->getPositions()[2];
  FReal *const targetsForcesX = inTargets->getForcesX();
  FReal *const targetsForcesY = inTargets->getForcesY();
  FReal *const targetsForcesZ = inTargets->getForcesZ();
  FReal *const targetsPotentials = inTargets->getPotentials();

  const ComputeClass mOne = ComputeClass(1);

  { // In this part, we compute (vectorially) the interaction
    // within the target leaf.

    const FSize nbParticlesSources = nbParticlesTargets;
    const FReal *const sourcesPhysicalValues = targetsPhysicalValues;
    const FReal *const sourcesX = targetsX;
    const FReal *const sourcesY = targetsY;
    const FReal *const sourcesZ = targetsZ;
    FReal *const sourcesForcesX = targetsForcesX;
    FReal *const sourcesForcesY = targetsForcesY;
    FReal *const sourcesForcesZ = targetsForcesZ;
    FReal *const sourcesPotentials = targetsPotentials;

    for (FSize idxTarget = 0; idxTarget < nbParticlesTargets; ++idxTarget) {
      FSize idxSource = idxTarget + 1;
      {
        const FSize nbVectorizedInteractions =
            ((nbParticlesSources - idxSource) / NbFRealInComputeClass) *
                NbFRealInComputeClass +
            idxSource;
        const ComputeClass tx = ComputeClass(targetsX[idxTarget]);
        const ComputeClass ty = ComputeClass(targetsY[idxTarget]);
        const ComputeClass tz = ComputeClass(targetsZ[idxTarget]);
        const ComputeClass tv = ComputeClass(targetsPhysicalValues[idxTarget]);
        ComputeClass tfx = ComputeClass::GetZero();
        ComputeClass tfy = ComputeClass::GetZero();
        ComputeClass tfz = ComputeClass::GetZero();
        ComputeClass tpo = ComputeClass::GetZero();

        for (; idxSource < nbVectorizedInteractions;
             idxSource += NbFRealInComputeClass) {
          ComputeClass dx = tx - ComputeClass(&sourcesX[idxSource]);
          ComputeClass dy = ty - ComputeClass(&sourcesY[idxSource]);
          ComputeClass dz = tz - ComputeClass(&sourcesZ[idxSource]);

          ComputeClass inv_square_distance =
              mOne / (dx * dx + dy * dy + dz * dz);

          /*const*/ ComputeClass inv_distance = inv_square_distance.sqrt();
          // #########################################################
          // if particles are too close adapt their distance to cutOffDistance
          // of 1 nm (calc force and potential with that distance)
          auto mask_ab = (inv_distance > ComputeClass(1e9));
          if (mask_ab.isAllTrue()) {
            inv_square_distance = ComputeClass(1e18);
            auto factor = inv_distance * ComputeClass(1e-9);
            dx *= factor;
            dy *= factor;
            dz *= factor;
            inv_distance = ComputeClass(1e9);
          } // #########################################################

          inv_square_distance *= inv_distance;
          inv_square_distance *=
              tv * ComputeClass(&sourcesPhysicalValues[idxSource]);

          dx *= -inv_square_distance;
          dy *= -inv_square_distance;
          dz *= -inv_square_distance;

          tfx += dx;
          tfy += dy;
          tfz += dz;
          tpo += inv_distance * ComputeClass(&sourcesPhysicalValues[idxSource]);

          (ComputeClass(&sourcesForcesX[idxSource]) - dx)
              .storeInArray(&sourcesForcesX[idxSource]);
          (ComputeClass(&sourcesForcesY[idxSource]) - dy)
              .storeInArray(&sourcesForcesY[idxSource]);
          (ComputeClass(&sourcesForcesZ[idxSource]) - dz)
              .storeInArray(&sourcesForcesZ[idxSource]);
          (ComputeClass(&sourcesPotentials[idxSource]) + inv_distance * tv)
              .storeInArray(&sourcesPotentials[idxSource]);
        }

        targetsForcesX[idxTarget] += tfx.horizontalSum();
        targetsForcesY[idxTarget] += tfy.horizontalSum();
        targetsForcesZ[idxTarget] += tfz.horizontalSum();
        targetsPotentials[idxTarget] += tpo.horizontalSum();
      }
      {
        const FReal tx = FReal(targetsX[idxTarget]);
        const FReal ty = FReal(targetsY[idxTarget]);
        const FReal tz = FReal(targetsZ[idxTarget]);
        const FReal tv = FReal(targetsPhysicalValues[idxTarget]);
        FReal tfx = FReal(0.);
        FReal tfy = FReal(0.);
        FReal tfz = FReal(0.);
        FReal tpo = FReal(0.);

        for (; idxSource < nbParticlesSources; idxSource += 1) {
          FReal dx = tx - FReal(sourcesX[idxSource]);
          FReal dy = ty - FReal(sourcesY[idxSource]);
          FReal dz = tz - FReal(sourcesZ[idxSource]);

          FReal inv_square_distance = FReal(1) / (dx * dx + dy * dy + dz * dz);
          /*const*/ FReal inv_distance = FMath::Sqrt(inv_square_distance);

          // #########################################################
          // if particles are too close adapt their distance to cutOffDistance
          // of 1 nm (calc force and potential with that distance)
          if (inv_distance > 1e9) {
            inv_square_distance = 1e18;
            auto factor = inv_distance * 1e-9;
            dx *= factor;
            dy *= factor;
            dz *= factor;
            inv_distance = 1e9;
          } // #########################################################

          inv_square_distance *= inv_distance;
          inv_square_distance *= tv * FReal(sourcesPhysicalValues[idxSource]);

          dx *= -inv_square_distance;
          dy *= -inv_square_distance;
          dz *= -inv_square_distance;

          tfx += dx;
          tfy += dy;
          tfz += dz;
          tpo += inv_distance * FReal(sourcesPhysicalValues[idxSource]);

          sourcesForcesX[idxSource] -= dx;
          sourcesForcesY[idxSource] -= dy;
          sourcesForcesZ[idxSource] -= dz;
          sourcesPotentials[idxSource] += inv_distance * tv;
        }

        targetsForcesX[idxTarget] += tfx;
        targetsForcesY[idxTarget] += tfy;
        targetsForcesZ[idxTarget] += tfz;
        targetsPotentials[idxTarget] += tpo;
      }
    }
  }
}

template <class FReal, class ContainerClass, class ComputeClass,
          int NbFRealInComputeClass>
static void GenericFullRemote(ContainerClass *const FRestrict inTargets,
                              const ContainerClass *const inNeighbors[],
                              const int limiteNeighbors) {
  const FSize nbParticlesTargets = inTargets->getNbParticles();
  const FReal *const targetsPhysicalValues = inTargets->getPhysicalValues();
  const FReal *const targetsX = inTargets->getPositions()[0];
  const FReal *const targetsY = inTargets->getPositions()[1];
  const FReal *const targetsZ = inTargets->getPositions()[2];
  FReal *const targetsForcesX = inTargets->getForcesX();
  FReal *const targetsForcesY = inTargets->getForcesY();
  FReal *const targetsForcesZ = inTargets->getForcesZ();
  FReal *const targetsPotentials = inTargets->getPotentials();

  const ComputeClass mOne = ComputeClass(1);

  for (FSize idxNeighbors = 0; idxNeighbors < limiteNeighbors; ++idxNeighbors) {
    if (inNeighbors[idxNeighbors]) {
      const FSize nbParticlesSources =
          inNeighbors[idxNeighbors]->getNbParticles();
      const FReal *const sourcesPhysicalValues =
          inNeighbors[idxNeighbors]->getPhysicalValues();
      const FReal *const sourcesX =
          inNeighbors[idxNeighbors]->getPositions()[0];
      const FReal *const sourcesY =
          inNeighbors[idxNeighbors]->getPositions()[1];
      const FReal *const sourcesZ =
          inNeighbors[idxNeighbors]->getPositions()[2];

      for (FSize idxTarget = 0; idxTarget < nbParticlesTargets; ++idxTarget) {
        FSize idxSource = 0;
        {
          const FSize nbVectorizedInteractions =
              (nbParticlesSources / NbFRealInComputeClass) *
              NbFRealInComputeClass;
          const ComputeClass tx = ComputeClass(targetsX[idxTarget]);
          const ComputeClass ty = ComputeClass(targetsY[idxTarget]);
          const ComputeClass tz = ComputeClass(targetsZ[idxTarget]);
          const ComputeClass tv =
              ComputeClass(targetsPhysicalValues[idxTarget]);
          ComputeClass tfx = ComputeClass::GetZero();
          ComputeClass tfy = ComputeClass::GetZero();
          ComputeClass tfz = ComputeClass::GetZero();
          ComputeClass tpo = ComputeClass::GetZero();

          for (; idxSource < nbVectorizedInteractions;
               idxSource += NbFRealInComputeClass) {
            ComputeClass dx = tx - ComputeClass(&sourcesX[idxSource]);
            ComputeClass dy = ty - ComputeClass(&sourcesY[idxSource]);
            ComputeClass dz = tz - ComputeClass(&sourcesZ[idxSource]);

            ComputeClass inv_square_distance =
                mOne / (dx * dx + dy * dy + dz * dz);
            const ComputeClass inv_distance = inv_square_distance.sqrt();

            inv_square_distance *= inv_distance;
            inv_square_distance *=
                tv * ComputeClass(&sourcesPhysicalValues[idxSource]);

            dx *= -inv_square_distance;
            dy *= -inv_square_distance;
            dz *= -inv_square_distance;

            tfx += dx;
            tfy += dy;
            tfz += dz;
            tpo +=
                inv_distance * ComputeClass(&sourcesPhysicalValues[idxSource]);
          }

          targetsForcesX[idxTarget] += tfx.horizontalSum();
          targetsForcesY[idxTarget] += tfy.horizontalSum();
          targetsForcesZ[idxTarget] += tfz.horizontalSum();
          targetsPotentials[idxTarget] += tpo.horizontalSum();
        }
        {
          const FReal tx = FReal(targetsX[idxTarget]);
          const FReal ty = FReal(targetsY[idxTarget]);
          const FReal tz = FReal(targetsZ[idxTarget]);
          const FReal tv = FReal(targetsPhysicalValues[idxTarget]);
          FReal tfx = FReal(0.);
          FReal tfy = FReal(0.);
          FReal tfz = FReal(0.);
          FReal tpo = FReal(0.);

          for (; idxSource < nbParticlesSources; idxSource += 1) {
            FReal dx = tx - FReal(sourcesX[idxSource]);
            FReal dy = ty - FReal(sourcesY[idxSource]);
            FReal dz = tz - FReal(sourcesZ[idxSource]);

            FReal inv_square_distance =
                FReal(1) / (dx * dx + dy * dy + dz * dz);
            const FReal inv_distance = FMath::Sqrt(inv_square_distance);

            inv_square_distance *= inv_distance;
            inv_square_distance *= tv * FReal(sourcesPhysicalValues[idxSource]);

            dx *= -inv_square_distance;
            dy *= -inv_square_distance;
            dz *= -inv_square_distance;

            tfx += dx;
            tfy += dy;
            tfz += dz;
            tpo += inv_distance * FReal(sourcesPhysicalValues[idxSource]);
          }

          targetsForcesX[idxTarget] += tfx;
          targetsForcesY[idxTarget] += tfy;
          targetsForcesZ[idxTarget] += tfz;
          targetsPotentials[idxTarget] += tpo;
        }
      }
    }
  }
}

} // namespace FP2PRCutoff

template <class FReal> struct FP2PRTCutoff {};

#include "InastempCompileConfig.h"

template <> struct FP2PRTCutoff<double> {
  template <class ContainerClass>
  static void FullMutual(ContainerClass *const FRestrict inTargets,
                         ContainerClass *const inNeighbors[],
                         const int limiteNeighbors) {
    FP2PRCutoff::GenericFullMutual<double, ContainerClass, InaVecBestTypeDouble,
                                   InaVecBestTypeDouble::VecLength>(
        inTargets, inNeighbors, limiteNeighbors);
  }

  template <class ContainerClass>
  static void Inner(ContainerClass *const FRestrict inTargets) {
    FP2PRCutoff::GenericInner<double, ContainerClass, InaVecBestTypeDouble,
                              InaVecBestTypeDouble::VecLength>(inTargets);
  }

  template <class ContainerClass>
  static void FullRemote(ContainerClass *const FRestrict inTargets,
                         const ContainerClass *const inNeighbors[],
                         const int limiteNeighbors) {
    FP2PRCutoff::GenericFullRemote<double, ContainerClass, InaVecBestTypeDouble,
                                   InaVecBestTypeDouble::VecLength>(
        inTargets, inNeighbors, limiteNeighbors);
  }
};

template <> struct FP2PRTCutoff<float> {
  template <class ContainerClass>
  static void FullMutual(ContainerClass *const FRestrict inTargets,
                         ContainerClass *const inNeighbors[],
                         const int limiteNeighbors) {
    FP2PRCutoff::GenericFullMutual<float, ContainerClass, InaVecBestTypeFloat,
                                   InaVecBestTypeFloat::VecLength>(
        inTargets, inNeighbors, limiteNeighbors);
  }

  template <class ContainerClass>
  static void Inner(ContainerClass *const FRestrict inTargets) {
    FP2PRCutoff::GenericFullMutual<float, ContainerClass, InaVecBestTypeFloat,
                                   InaVecBestTypeFloat::VecLength>(inTargets);
  }

  template <class ContainerClass>
  static void FullRemote(ContainerClass *const FRestrict inTargets,
                         const ContainerClass *const inNeighbors[],
                         const int limiteNeighbors) {
    FP2PRCutoff::GenericFullRemote<float, ContainerClass, InaVecBestTypeFloat,
                                   InaVecBestTypeFloat::VecLength>(
        inTargets, inNeighbors, limiteNeighbors);
  }
};

#endif // FP2PR_CUTOFF_HPP
