
# ValleyTypes

Each derived object of the type *ValleyType* should represent one group of valleys of the analytical band structure of a semiconductor. It is assumed that each subvalley contained in that group possesses the same effective masses in its ellipse coordinate system. If only one valley is required the function *getDegeneracyFactor()* should return 1.  Additionally, for anisotropical effective masses the Herring Vogt transformation is used throughout this code.

## Background Information
Valleys stem from the analytical band structure, which uses parabolic and non-parabolic functions to approximate the full band structure of a given material. The idea is that particles close to extrema of the conduction / valence band, behave like free particles with an adapted effective mass *m**.  Each approximated extremum of the bands is then called valley, additionally valleys that possess the same effective masses are usually grouped and also called *valley*.  One example for a group of valleys would be the X-valleys in silicon, which have a degeneracy of 6 (as can be seen [here](https://www.iue.tuwien.ac.at/phd/smirnov/node45.html)).

Further informations on the analytical band structure are given in the following references:

- [Simplified Band Structure Models and Carrier Dynamics](https://nanohub.org/resources/1522/download/simplifiedbandstructurecarrierdynamics_word.pdf): A PDF that describes the ideas behind analytical band structures. 

- [Analytical Bandstructure in Dissertation of Dr. Sergey Smirnov](https://www.iue.tuwien.ac.at/phd/smirnov/node42.html): A description of analytical band structures including the ideas of the Herring-Vogt transformation.

- [Silicon Band Structure Models in Dissertation of Dr. Wilfried Wessner](https://www.iue.tuwien.ac.at/phd/wessner/node31.html) : A short overview on how to approximate the band structure of silicon.

  
## Implementation
*ValleyType* is an abstract class that presents all the functionalities that are needed for the implementation of one specific valley or type of valley. In the following the functions and the ideas behind it will be explained:

### Description of the effective mass:
 The effective mass *m**, which can either be isotropical (direction-independent) or anisotropical is one of the main characteristics of a valley. 

In case the effective mass is isotropical *m** is a scalar value, which is returned by the functions *getEffMassDOS()* and *getEffMassCond()*. Also, the transformation to and out of the ellipse coordinate system is unnecessary as the effective mass in all coordinate systems can be represented by the given scalar value.

In case the effective mass is anisotropical the effective mass is a matrix. This case is a little more complicated. Even though, the mass is a matrix, there still exists one coordinate system, in which this matrix only has non-zero diagonal entries of the form

$$ m^* =  	\begin{pmatrix}
						m_1 & 0 & 0 \\
						0 & m_2 & 0 \\
						0 & 0 & m_3
					\end{pmatrix} \quad .
$$

The coordinate system, which fulfills this condition is then called ellipse coordinate system (ECS) and it has to be noted that this system is not always the same as the device coordinate system (DCS), which is used within the simulation.  The effective mass of each subvalley is then defined uniquely by the determination of the three values $m_1$, $m_2$ and $m_3$ and the specification of the relation between the ECS and the DCS. As already mentioned before, the effective masses of a valley have to be the same for each subvalley. However, the orientation of the ECS of the subvalleys can differ. The conversion of a vector between the ECS and the DCS is handled by the functions *transformToDeviceCoord* and *transformToEllipseCoord*.  Additionally, the anisotropical valleys are handled using the Herring-Vogt Transformation, which is implemented in the code and uses the function *GetVogtTransformationFactor* $v_{HV}$, which is defined as 

$$v_{HV} = \left( \sqrt{\frac{m}{m_1}},  \sqrt{\frac{m}{m_2}},  \sqrt{\frac{m}{m_3}}\right) \quad ,$$ 

where m is either mostly either the DOS effective mass or the conduction effective mass. In case the valley is isotrop $v_{HV}=(1,1,1)$ should be used. Lastly, the conductive effective mass and the density of states (DOS) effective mass, which are returned by *getEffMassDOS()* and *getEffMassCond()* are determined by the geometric and harmonic mean of $m_1$, $m_2$ and $m_3$. 

### Degeneracy of the valley
The degeneracy of the valley, which is retruned by the function *getDegeneracyFactor* is the number of subvalleys that is contained in this group of valleys and should thereby be an integer. As mentioned in the beginning if only one valley should be simulated this fucntion should return 1. 

### Bottom Energy of the valley
The bottom energy of the valley is the energy at extremum of the valley in $eV$.  Currently it is assumed that each subvalley possesses the same bottom energy. This information is needed for specific scatter events between different groups of valleys, if only one group of valley is considered, this value can be set to 0. If multiple group of valleys are considered the energy difference between the bottom energy is the parameter of interest, meaning that a shift in all given energies is not important - as long as it is consistent in all given bottom energies of valleys.

### Non parabolicity factor ($\alpha$)
The non-parabolicity factor (alpha) determines the amount of non-parabolicity of that valley and is given in $eV^{-1}$. This parameter is present in the relation between the energy $E$ and the wave vector $k$ of the given valley

$$\gamma(k) := E (1 + \alpha E) = \frac{\hbar^2 k^2}{ 2 m} \quad ,$$

here an isotropical effective mass is assumed and the relation between those two parameter is usually called dispersion relation. If $\alpha=0$ the approximation of the valley is called parabolic, else it is called non-parabolic. 

### Additional Functionality
The additional functions that are implemented in valleyType, are:

- *getNormWaveVec(energy)*: This function calculates the norm of the wave vector $k$ of a particle based on a given energy using the dispersion relation.
- *getEnergy(waveVector)*: This function calculates the energy of a particle based on its given wavevector. This calculation is again based on the dispresion relation.
- *getGamma(energy)*: This function returns  $\gamma(k) := E (1 + \alpha E)$.
- *getVelocity(waveVector, energy, idxSubValley)*: returns the group velocity $v$ if a particle assigned to a specific valley, which is given by
$$v = \nabla_{k}E(k)/\hbar$$

### Additional References
The implementation of the valleys in this code is based on the description of the valleys in the thesis "Generalized Monte Carlo Tool for Investigating Low-Field and High Field Properties of Materials Using Non-parabolic Band Structure Model" by Raguraj Hathwar, which can be found [here](https://www.researchgate.net/publication/267566065_Generalized_Monte_Carlo_Tool_for_Investigating_Low-Field_and_High_Field_Properties_of_Materials_Using_Non-parabolic_Band_Structure_Model).