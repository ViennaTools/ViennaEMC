import numpy as np
import pandas as pd
from scipy.constants import k, elementary_charge
import matplotlib.pyplot as plt

"""

This file helps with the calculation of the mobility of a material
from the tracked velocities during a simulation (via the velocity
auto correlation function VACF). This calculation is described in
more detail in "The Monte Carlo Method for Semiconductor Device
Simulation" from Jacobini and Lugli (p. 74 ff and p. 143 ff).

To use this file the simulation should be performed with no applied
field and a long simulation time (e.g. 5e-12 s). Additionally, the
parameter nrBetweenOutput should be low (e.g. 100) and the boolean
sampleVelocityOfParticles must be set to true so that the file 
"singleLayerMoS2velocityAtDiffTimes" + parameter + ".txt" is created
and filled.

NOTE: This calculation can take some time. 
"""

### parameter to be adapted #######################################

# path to the file of interest
path = "../../build/"
# name of the file of interest with tracked velocities at different
#steps
fileNames = ["singleLayerMoS2velocityAtDiffTimesE0T300N20808.txt"]

# number of the used steps for the calculation (if set to -1 all 
# possible steps are used)
nrUsedSteps = -1

# number of particles used for the calculation (if set to -1 all
# particles are used)
nrUsedParticles = -1

###################################################################

# integrates the given yValues with the trapezoid rule
def integrateTrapezoid(xMin, xMax, deltaX, yValues):
    xValues = np.arange(xMin, xMax + deltaX / 2, deltaX)
    result = 0
    for i in range(xValues.size - 2):
        result += yValues[i] + yValues[i + 1] 
    return result * deltaX / 2

# returns mobility (in cm^2/Vs) calculated with the Einstein relation
# from the given diffusion coefficient
def getMobility(diffusionCoeff, temperature = 300):
    return diffusionCoeff * elementary_charge / (k * temperature) * 1e4

def readData(file, verbose=True):
    with open(file) as f:     
        dT = float(f.readline().split()[0])       
        data = pd.read_csv(path + "/" + fileName, sep=" ", header=None, skip_blank_lines=True, skiprows=1, engine="python").to_numpy()
    if(verbose):
        print("\t\tnr. used particles in simulation= ", int(float(data.shape[1])/3.))
        print("\t\tnr. used steps recorded velocities in simulation = ", data.shape[0])
    return dT, data

# calculates the mobility from the given file of velocities by calculating via
# the autocorrelation function and the drift coefficient
def calcMobility(path, fileName, nrUsedSteps, nrUsedPart=0, plotMeanVelCorr=True, plotAllVelCorr=False, temperature=300, driftVelocity = 0, degreesOfFreedom = 2):    
    # read file + store data
    print("\tReading Data ....")
    dT, data = readData(path + "/" + fileName)
    print("\t\ttemperature = ", temperature, " K")

    # check input parameter
    print("\tChecking Input Parameter ...")
    if(nrUsedPart <= 0 or nrUsedPart > data.shape[1] / 3.):
        nrUsedPart = int(data.shape[1] / 3)
        print("\t\tWARNING: Nr of used particles for calculation of mobility was invalid. Reset to number of particles used in simulation: ", nrUsedPart)    
    if nrUsedSteps <= 0 or nrUsedSteps > float(data.shape[0]) / 2.:
        nrUsedSteps = int(float(data.shape[0]) / 2.)
        print("\t\tWARNING: Nr of used steps for calculation of mobility was invalid. Reset to half of the number of steps used in simulation: ", nrUsedSteps)
    print("nr of used steps in VACF calculation =" , nrUsedSteps)
    print("nr of used particles in VACF calculation =" , nrUsedPart)

    print("\tCalculating diffusion coefficient ...")
    # calculate velocity auto correlation function (VACF)
    velCorr = np.zeros(nrUsedSteps)
    for idxPart in range(nrUsedPart):
        if(idxPart % (nrUsedPart / 10) == 0):
            print("\t\t Current Particle: ", idxPart, " / ", nrUsedPart)
        # calculate velocity correlation between i * dT and (i-j) * dT
        currVelCorr = np.zeros([nrUsedSteps, nrUsedSteps])
        for i in range(nrUsedSteps):
            for j in range(nrUsedSteps):
                currPart = 3 * idxPart
                currVelCorr[i, j] += data[nrUsedSteps + i , currPart] * data[nrUsedSteps + i - j, currPart]
                currVelCorr[i, j] += data[nrUsedSteps + i , currPart + 1] * data[nrUsedSteps + i - j, currPart + 1]
                # currVelCorr[i, j] += data[nrUsedSteps + i , currPart + 2] * data[nrUsedSteps + i - j, currPart + 2]
        # sum up the velocity correlations with same j
        currVelCorr = currVelCorr.mean(axis=0) - driftVelocity**2
        if plotAllVelCorr:
            plt.plot(np.arange(currVelCorr.size) * dT, currVelCorr, ":")
        velCorr += currVelCorr
    print("\t\t Current Particle: ", nrUsedPart, " / ", nrUsedPart)
    # calculate mean over diffusion coeff + VACF of all particles
    velCorr /= nrUsedPart
    if plotMeanVelCorr:
        plt.plot(np.arange(velCorr.size) * dT, velCorr, color="black", label="mean correlation function")
        plt.ylabel("autocorrelation function [m^2  / s^2]")
        plt.xlabel("time [s]")
        plt.legend()
    # integrate VACF to get diffusion coefficient, return mobility
    diffusionCoeff = integrateTrapezoid(0, nrUsedSteps * dT, dT, velCorr) / degreesOfFreedom
    return getMobility(diffusionCoeff, temperature)



mobility = []
for fileName in fileNames:
    print("Working on ", fileName, "...")
    mobility.append(calcMobility(path, fileName, nrUsedSteps, nrUsedParticles))
    print("\tmobility = ", mobility[-1], " cm^2/Vs")
plt.show()

