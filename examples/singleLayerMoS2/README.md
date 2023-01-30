# Single Layer MoS2 Simulation Example
This example shows how to simulate a free-standing monolayer of molybdenum disulfide (ML-MoS2).

The model of the material is based on different papers that can be found in the following:

- [Phonon-limited mobility in n-type single-layer MoS2 from first principles - Kaasbjerg et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.115317)
- [Intrinsic electrical transport properties of monolayer silicene and MoS2 from first principles - Li et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.115418)
- [Semi-classical transport in MoS2 and MoS2 transistors by a Monte Carlo approach - Pilotto et al.](https://www.sciencedirect.com/science/article/abs/pii/S0038110122000673)

It should be noted that this implementation is still in progress and only the results from the simulations with the parameters from the last listed paper (from Pilotto et al.)
were compared to the results of the paper. 

The description of this example is organized as follows:

- **General Idea of Example**: This section describes the general idea of the example and ideas of what could be tested with that example.
- **Parameter of the Simulation**: This section describes the parameter of the simulation that can be customized.
- **Description of the Simulation Process**: This section describes the simulation process in more detail.
- **Output of the Simulation**: This section describes the output of this simulation and how it could be plotted.
- **References**: This section lists useful links for that example.

## General Idea
The general idea of this example is to simulate the transport characteristics of electrons in ML-MoS2 with various applied electric fields. 

The particle characteristics that are of interest in this kind of simulation are:

- the average drift velocity of the particles (the average velocity in the direction of the applied field) in each valley
- the average energy of the particles in each valley
- the percentage of particles that are present in each valley (valley occupation)

Those average particle characteristics are calculated in each step and stored in seperate files in the end of the simulation.

Different things can be tested when adapting this example, e.g., the effect of different scatter mechanisms on the drift velocity. or the effects of variations in the band structure on the transport characteristics.

The main difference to silicon is that the particle transport in this material happens in two dimensions and that the transport in the third dimension is discarded within the simulation. For this reason a customized particleType *electron2D* is applied. 
 
## Customizable Parameters
The parameters that are used and can be adapted are listed and described in the following.

###  Background Electric Field
The *background electric field* can be set via the parameter *appliedFields* and *appliedFieldDir*. As the name of those parameter states, the first parameter sets the strength (norm) and the second the direction of the electric field. The direction of the electric field doesn't have to be a unit vector, it is normalized later in the code.

### Paper
As mentioned in the introduction, the parameters given in different papers were implemented in this chapter. The selected paper can be implemented with the parameter *selectedPaperForParameter*. The parameter for each given paper are stored in a separate parameter file, information on the assumptions and parameters in the corresponding paper can be found in that file.

### Number of simulated particles
The number of simulated particles can be changed via adapting *maxPos* or *spacing*, as at each grid point a specific number of particles is generated; so if the number of grid points is increased also the number of particles increases.

### Track velocity of particles (for calculation of mobility)
The parameter *sampleVelocityOfParticles* allows the recording of the velocities of all particles which can be used for the calculation of the 
velocity autocorrelation fucntion. If no field is applied this then can be used to calculate the diffusion coefficient / mobility of the material.

##  Simulation Description

Same as described in [Bulk Simulation Example](../bulkSimulation/README.md).

## Output

The results are stored in multiple files:

- **prefix + AvgEnergy (+ parameter).txt**: This file contains at least two columns, the first one represents the time (in seconds) and the following ones represent the mean particle energy (in eV) in each valley.
- **prefix + AvgVelocity (+ parameter).txt**: This file contains at least two columns, the first one represents the time (in seconds) and the following ones represent the mean particle velocity in the direction of the applied electric field (in m/s).
- **prefix + valleyOccupation (+ parameter).txt**: This file contains at least two columns, the first one represents the time (in seconds) and the following ones represent the valley occupation percentage for each valley.
- **prefix + velocityAtDifferentTimes (+ parameter).txt**: This file contains the velocity of all particles at different steps and allows for the calculation of the velocity auto correlation function. It is only created if *sampleVelocityOfParticles* is set to true.

The file [plotSingleLayerMoS2Results.py](plotSingleLayerMoS2Results.py) shows how these resulting files can be plotted. For the usage of this file `emcPlottingFiles` has to be installed, as is described [here](../../README.md). Additionally, the file [calcMobilityFromVACF.py](calcMobilityFromVACF.py) can be used to calculate the mobility from the **prefix + velocityAtDifferentTimes (+ parameter).txt**-file.

## References
The main references are the three listed papers from above, which are repeated here:
- [Phonon-limited mobility in n-type single-layer MoS2 from first principles - Kaasbjerg et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.115317)
- [Intrinsic electrical transport properties of monolayer silicene and MoS2 from first principles - Li et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.115418)
- [Semi-classical transport in MoS2 and MoS2 transistors by a Monte Carlo approach - Pilotto et al.](https://www.sciencedirect.com/science/article/abs/pii/S0038110122000673)


