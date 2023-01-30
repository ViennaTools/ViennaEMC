from emcPlottingFiles.readResultFile import *
from emcPlottingFiles.plotBulkSimulationResults import *

import matplotlib.pyplot as plt

"""

This file shows the results of the example "bulkSimulationFMM".

Note: To use this file you have to install emcPlottingFiles. For 
instructions on how to do that look at the main README.md.

To use this file, run the example "bulkSimulationFMM", then adapt the
"path" (path to the files).

"""

######## Adaptable Parameters #####################################################

path = "../../build/"           # path to result files (adapt if needed)
prefix = "bulkSimulationFMM"    # file name prefix
suffix = ".txt"                 # file name suffix

# if parameter were used in the filename, set this parameter to true
# see "includeParameterInFileName" in bulkSimulation.cpp file
addedParameterToFileName = False
nrPart = 20808      # nr of particles (adapt if needed)
temp = 300          # temperature [K] (adapt if needed)
appliedField = 1e5  # electric field [V / m] (adapt if needed)

# select the scatter rates that should be plotted (adapt if needed)
plotAcouticRate = True                  # rates from emcAcousticScatterMechanism.hpp
plotZeroOrderIntervalleyRate = True     # rates from emcZeroOrderIntervalleyScatterMechanism.hpp
plotFirstOrderIntervalleyRate = True    # rates from emcFirstOrderIntervalleyScatterMechanism.hpp
plotCoulombRate = False                 # rates from emcCoulombScatterMechanism.hpp

#################################################################################


# create parameter part of filename if needed
param = ""
if addedParameterToFileName:
    param = "E" + str(int(appliedField)) + "T" + str(temp) + "N" + str(nrPart)

# read all result files
resultFiles = []
resultFiles.append(readBulkSimulationAvgFile(path + prefix + "AvgEnergy" + param + ".txt"))
resultFiles.append(readBulkSimulationAvgFile(path + prefix + "AvgDriftVelocity" + param + ".txt"))
resultFiles.append(readBulkSimulationAvgFile(path + prefix + "valleyOccupation" + param + ".txt"))

# plot each result file (time vs. observed parameter)
ylabels = ["avg. energy [eV]", "avg drift velocity [m / s]", "valley occupation [%]"]
titles = ["Avg. Energy Development", "Avg. Drift Velocity Development", "Valley Occupation Development"]
for data, ylabel, title in zip(resultFiles, ylabels, titles):
    plt.figure()
    for valley in range(data.shape[1] - 1):
        plt.plot(data["time"], np.abs(data["valley_" + str(valley)]), label="valley " + str(valley))
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel("time [s]")
    plt.legend()

nrValleys = resultFiles[0].shape[1] - 1
print("Nr. of used valleys = ", nrValleys)

# print valley-wise mean characteristics
print("Avg. characteristics for each valley: ")

print("\t Avg. valley occupation:")
avgOccPerValley = calcValleywiseMean(resultFiles[2])
for idxValley in range(nrValleys):
    print("\t\tvalley ", idxValley, ": ", "{:.2f}".format(avgOccPerValley[idxValley] * 100), " %")

print("\t Avg. DriftVelocity:")
avgDriftVel = calcValleywiseMean(resultFiles[1])
for idxValley in range(nrValleys):
    print("\t\tvalley ", idxValley, ": ", "{:.2f}".format(np.abs(avgDriftVel[idxValley])) , " m / s" )

print("\t Avg. Energy:")
avgEnergy = calcValleywiseMean(resultFiles[0])
for idxValley in range(nrValleys):
    print("\t\tvalley ", idxValley, ": ", "{:.2f}".format(avgEnergy[idxValley] * 1e3), " meV" )

# print mean characteristics (weighted by the valley-occupation)
avgDriftVel = calcValleyOccupationWeightedMean(resultFiles[1], resultFiles[2])
avgEnergy = calcValleyOccupationWeightedMean(resultFiles[0], resultFiles[2])
print("Avg. characteristics weighted by the valley-occupation: ")
print("\tAvg. Drift Velocity = ", "{:.2f}".format(np.abs(avgDriftVel)), " m / s")
print("\tAvg. Energy = ", "{:.2f}".format(avgEnergy * 1e3), " meV")

plt.figure()

# plot scatter mechanisms
if plotAcouticRate:
    acRate = readScatterMechanismFile2(path, "Acoustic", 0, 0)
    plt.semilogy(acRate["energy"], acRate["rate"], label="acoustic", color="red")
if plotZeroOrderIntervalleyRate:
    ivZeroFRate = readScatterMechanismFile2(path, "ZeroInterValleyEmissionF", 0 , 0)
    ivZeroGRate = readScatterMechanismFile2(path, "ZeroInterValleyEmissionG", 0 , 0)
    plt.semilogy(ivZeroFRate["energy"], ivZeroFRate["rate"], label="f-type zero order emission")
    plt.semilogy(ivZeroGRate["energy"], ivZeroGRate["rate"], label="g-type zero order emission")
    
    ivZeroFRate = readScatterMechanismFile2(path, "ZeroInterValleyAbsorptionF", 0 , 0)
    ivZeroGRate = readScatterMechanismFile2(path, "ZeroInterValleyAbsorptionG", 0 , 0)
    plt.semilogy(ivZeroFRate["energy"], ivZeroFRate["rate"], label="f-type zero order absoption")
    plt.semilogy(ivZeroGRate["energy"], ivZeroGRate["rate"], label="g-type zero order absorption")
if plotFirstOrderIntervalleyRate:
    ivFirstFRate = readScatterMechanismFile2(path, "FirstInterValleyEmissionF", 0 , 0)
    ivFirstGRate = readScatterMechanismFile2(path, "FirstInterValleyEmissionG", 0 , 0)
    plt.semilogy(ivFirstFRate["energy"], ivFirstFRate["rate"], label="f-type first order emission")
    plt.semilogy(ivFirstGRate["energy"], ivFirstGRate["rate"], label="g-type first order emission")

    ivFirstFRate = readScatterMechanismFile2(path, "FirstInterValleyAbsorptionF", 0 , 0)
    ivFirstGRate = readScatterMechanismFile2(path, "FirstInterValleyAbsorptionG", 0 , 0)
    plt.semilogy(ivFirstFRate["energy"], ivFirstFRate["rate"], label="f-type first order absorption")
    plt.semilogy(ivFirstGRate["energy"], ivFirstGRate["rate"], label="g-type first order absorption")
if plotCoulombRate:
    coulombRate = readScatterMechanismFile2(path, "Coulomb", 0, 0)
    plt.semilogy(coulombRate["energy"], coulombRate["rate"], label= "coulomb")

plt.title("Scatter Rates")
plt.xlabel("energy [eV]")
plt.ylabel("rate [1 / s]")
plt.legend()
plt.show()