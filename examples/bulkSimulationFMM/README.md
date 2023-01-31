
# Bulk Simulation FMM Example

This example shows how to simulate bulk silicon with an applied background electric field and with using the Fast Multipole Method (FMM) to include the particle-particle interactions.

## General Idea

This simulation performs the same simulation as the [bulkSimulation Example](../bulkSimulation/) (see associated [README](../bulkSimulation/README.md)), the only difference is that here the particle-particle interactions are considered (with the help of the FMM).

With this example it is possible to add multiple moving or non-moving particle types to the simulation and see the effect that they have on each other.
 
### Force Calculation on Particle
Considering the particle-particle interaction in this case just means that the force $F_i$ that is felt by particle $i$ is calculated like that:

![FormulaForceOnParticle](https://latex.codecogs.com/svg.image?\Large&space;F_i&space;=&space;q_i&space;\cdot&space;(E^{b}&space;&plus;&space;E^{pp}_i)&space;)

Here $E^b$ is the constant background electric field, $E^{pp}$ is the electric field that stems from the particle interaction, which is calculated like that for the $i$-th particle:

![FormulaFieldParticleInteraction](https://latex.codecogs.com/svg.image?\Large&space;E^{pp}_i&space;=&space;\sum_{\substack{j=0&space;\\&space;j\neq&space;i}}^n&space;\frac{q_j&space;}{4\pi&space;\epsilon_r&space;\epsilon_0}&space;\frac{\vec{r}_{ij}}{\|\vec{r}_{ij}\|^3})

In this calculation $n$ is the number of particles, $\epsilon_r$ is the dielectric constant of the material, $\epsilon_0$ is the vacuum permittivity and $q_j$ the charge of the j-th particle. Furthermore $\vec{r}_{ij} = \vec{r}_j - \vec{r}_i$ is the distance between the two particle positions.

This sum that is needed for each particle is calculated using the Fast Multipole Method (with the library [scalFMM](https://gitlab.inria.fr/solverstack/ScalFMM)).

### Periodic Boundary Conditions in [scalFMM](https://gitlab.inria.fr/solverstack/ScalFMM)
See Doc Page of scalFMM: [Periodicity in scalFMM](https://solverstack.gitlabpages.inria.fr/ScalFMM/_periodic.html). 

The previously mentioned page states that periodicity in [scalFMM](https://gitlab.inria.fr/solverstack/ScalFMM) is not the usual periodicity, it's just that the simulation box (which has to be cubic and includes all the particles) is repeated for a specific number of times in each direction. The calculated potential and force at each particle in the original box then consists of the interactions with every other particle in the original box and with every particle in the repeated (copied + translated) boxes.

The number of times that  the box is repeated is determined with the parameter *periodicityParam* which can be in (-1, 0, 1, 2, 3,...). In the following table the relation between **periodicityParam** and the number of boxes that are repeated in each direction **nrBoxesDir** and the total number of boxes that are considered in the simulation **nrBoxesTotal**  is shown:

| periodicityParam   | nrBoxesDir|  nrBoxesTotal |
|:----------:|:-------------:|:------:|
| -1 |  3 | 27 |
| 0 |  6 	| 216  |
| 1 | 12 	| 1728  |
| 2 |  24 | 13824 |
| 3 |  48 	| 110592  |

The implementation of this can be seen here: [scalFMM/include/Core/FFMMAlgorithmPeriodic.hpp](https://gitlab.inria.fr/solverstack/ScalFMM/-/blob/master/include/Core/FFmmAlgorithmPeriodic.hpp)

How *periodicityParam* (in this example called *idxLevelAbove*) effects the particle potential can be seen in the following example: [scalFMM/Tests/Kernels/testRotationPeriodicBench.cpp](https://gitlab.inria.fr/solverstack/ScalFMM/-/blob/master/Tests/Kernels/testRotationPeriodicBench.cpp)

### Cut Off Radius Approach

It is also possible to use a cut-off radius for the electric field of the particle-particle interaction, to avoid numerical heating if particles are "too close" to each other. In this case the formula for the particle interaction between two particles $i$ and $j$, whose distance to each other is $\vec{r}_{ij} = \vec{r}_j - \vec{r}_i$ is calculated like before if the particles are further apart than a predefined cutoff radius $r_{cutoff}$ (whose value is set to $1~nm$ in the simulations), if the particles are too closer to each other than this cutoff radius the formula is adapted to:


![FormulaCutoffEField](https://latex.codecogs.com/svg.image?\Large&space;E^{pp-cutoff}(\vec{r}_{ij})&space;=&space;\frac{q_j&space;}{4\pi&space;\epsilon_r&space;\epsilon_0}&space;\frac{\vec{r}_{ij}}{r_{cutoff}^3}&space;\frac{r_{cutoff}}{\|\vec{r}_{ij}\|})

## Customizable Parameter

Most of the customizable simulation parameter are the same as in the [bulkSimulation Example](../bulkSimulation/).

In this simulation one can also add more than one particleType and see the effect of the different particle-types on each other (see emcDonor in example).

Furthermore one can also decide if the cutoff radius approach should be used or not, this can be set by either commenting line 9 (macro *#define USE_CUTOFF_KERNEL*) out or not.

As already mentioned in previous section with *periodicityParam* the handling of the periodic boundary conditions in scalFMM can be altered.

## Simulation Description

The simulation is the same as the one in [bulkSimulation Example](../bulkSimulation/). The only difference is that for the calculation of the particle-particle interaction some additional functions have to be called in each round:

- **executeFMM:** calculates the $E^{pp}$ (the electric field that stems from the particle interactions)
- **resetForcesAndPotential:** resets the calculated electric field
- **rearrangeTree:** the particles are stored in an octree, this function makes sure that each particle is in the right leaf of the octree after the particles were moved.

## Output

Same as [bulkSimulation](../bulkSimulation/) see [README](../bulkSimulation/README.md).

## References

- [scalFMM Library:](https://gitlab.inria.fr/solverstack/ScalFMM) The library that is used for the particle interaction calculations. In the following some of the used classes and some helpful tests that were useful for a better understanding are linked:
	- Kernel:
		- [include/Kernels/Rotation/RotationKernel.hpp:](https://gitlab.inria.fr/solverstack/ScalFMM/-/blob/master/include/Kernels/Rotation/FRotationKernel.hpp) Rotation Kernel of scalFMM using spherical harmonics.
	- Periodic Boundary Conditions: 
		- [Periodicity under scalFMM:](https://solverstack.gitlabpages.inria.fr/ScalFMM/_periodic.html) Doc Page where periodicity in scalFMM is shortly described
		- [Tests/Kernels/testRotationPeriodicBench.cpp:](https://gitlab.inria.fr/solverstack/ScalFMM/-/blob/master/Tests/Kernels/testRotationPeriodicBench.cpp) Test of scalFMM library that shows how to use the rotation kernel with periodic boundary conditions.
		- [Tests/Utils/testFMMAlgorithmPeriodic.cpp:](https://gitlab.inria.fr/solverstack/ScalFMM/-/blob/master/Tests/Utils/testFmmAlgorithmPeriodic.cpp) Test of scalFMM library that gives better insight on how periodicity in scalFMM works.
