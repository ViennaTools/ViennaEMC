from emcPlottingFiles.readResultFile import *
from emcPlottingFiles.plotBulkSimulationResults import *

from enum import Enum
import matplotlib.pyplot as plt

"""

This file shows the results of the example "singleLayerMoS2".

Note: To use this file you have to install emcPlottingFiles. For 
instructions on how to do that look at the main README.md.

To use this file, run the example "singleLayerMoS2", then adapt the
"path" (path to the files).

"""

# parameter type (like in singleLayerMoS2.cpp)
class PaperType(Enum):
    KAASBJERG = 0
    LI = 1
    PILOTTO = 2

# applied fields type (like in singleLayerMoS2.cpp)
# shows the used applied electric fields for the simulations
# in V / m
class AppliedFieldsType(Enum):
    NO = 0
    LOW = 1
    HIGH = 2
    CUSTOM = 3

nofields = [0]
lowfields = [1e4,  2e4,  3e4,  4e4,  5e4,  6e4,  7e4,  8e4,  9e4,  10e4,
    11e4, 12e4, 13e4, 14e4, 15e4, 16e4, 17e4, 18e4, 19e4, 20e4]
highfields = [1e5,   2e5,   5e5,   10e5,  20e5, 40e5,  60e5,  80e5,  
            100e5, 150e5, 200e5, 250e5, 300e5, 350e5, 400e5]

######## Adaptable Parameters #####################################################
path = "../../build/"       # path to result files (adapt if needed)
prefix = "singleLayerMoS2"  # file name prefix
suffix = ".txt"             # file name suffix

# if parameter were used in the filename, set this parameter to true
# see "includeParameterInFileName" in singleLayerMoS2.cpp file
addedParameterToFileName = True
nrPart = 20808      # nr of particles (adapt if needed)
temp = 300          # temperature [K] (adapt if needed)

# if true results for one single simulation are plotted
plotSingleResult = True
appliedField = 40e5 # electric field [V / m] (adapt if needed)

# if true results of all applied fields of that appliedFieldType
# are plotted
plotResultsOfAllAppliedFields = True
typeFields = AppliedFieldsType.HIGH # used applied field type (adapt if needed)
customfields = [40e5] # custom applied field type (adapt if needed)

plotScatterTables = True
parameterType = PaperType.PILOTTO
#################################################################################


if plotSingleResult:
    # create parameter part of filename if needed
    param = ""
    if addedParameterToFileName:
        param = "E" + str(int(appliedField)) + "T" + str(temp) + "N" + str(nrPart)

    # read all result files
    resultFiles = []
    resultFiles.append(readBulkSimulationAvgFile(path + prefix + "AvgEnergy" + param + ".txt"))
    resultFiles.append(readBulkSimulationAvgFile(path + prefix + "AvgDriftVelocity" + param + ".txt"))
    resultFiles.append(readBulkSimulationAvgFile(path + prefix + "valleyOccupation" + param + ".txt"))

    # plot each result file
    for data, ylabel in zip(resultFiles, ["avg. energy [eV]", "avg drift velocity [m / s]", "valley occupation"]):
        plt.figure()
        for valley in range(data.shape[1] - 1):
            plt.plot(data["time"], np.abs(data["valley_" + str(valley)]), label="valley " + str(valley))
        plt.ylabel(ylabel)
        plt.xlabel("time [s]")
        plt.legend()

    nrValleys = resultFiles[0].shape[1] - 1
    print("Nr. of used Valleys = ", nrValleys)

    # print valley-wise mean characteristics
    print("Avg. characteristics for each valley: ")

    print("\t Avg. valley occupation:")
    avgOccPerValley = calcValleywiseMean(resultFiles[2])
    for idxValley in range(nrValleys):
        print("\t\tvalley ", idxValley, ": ", "{:.2f}".format(avgOccPerValley[idxValley] * 100), "%")

    print("\t Avg. DriftVelocity:")
    avgDriftVel = calcValleywiseMean(resultFiles[1])
    for idxValley in range(nrValleys):
        print("\t\tvalley ", idxValley, ": ", "{:.2f}".format(np.abs(avgDriftVel[idxValley])*1e2) , " cm / s" )

    print("\t Avg. Energy:")
    avgEnergy = calcValleywiseMean(resultFiles[0])
    for idxValley in range(nrValleys):
        print("\t\tvalley ", idxValley, ": ", "{:.2f}".format(avgEnergy[idxValley] * 1e3), " meV" )

    # print mean characteristics (weighted by the valley-occupation)
    avgDriftVel = calcValleyOccupationWeightedMean(resultFiles[1], resultFiles[2])
    avgEnergy = calcValleyOccupationWeightedMean(resultFiles[0], resultFiles[2])
    print("Avg. characteristics weighted by the valley-occupation: ")
    print("\tAvg. Drift Velocity = ", "{:.2f}".format(np.abs(avgDriftVel) * 1e2), " cm / s")
    print("\tAvg. Energy = ", "{:.2f}".format(avgEnergy * 1e3), " meV")

    plt.show()

if plotResultsOfAllAppliedFields:
    avgEnergyValley, avgEnergyMean = [], []
    avgDriftVelValley, avgDriftVelMean = [], []
    avgOccValley = []

    if typeFields == AppliedFieldsType.NO:
        fields = nofields
    elif typeFields == AppliedFieldsType.LOW:
        fields = lowfields
    elif typeFields == AppliedFieldsType.HIGH:
        fields = highfields
    elif typeFields == AppliedFieldsType.CUSTOM:
        fields = customfields
    else:
        print("ERROR: Invalid Fields Type was selected.")
        fields = []

    for field in fields:
        if addedParameterToFileName:
            param = "E" + str(int(field)) + "T" + str(temp) + "N" + str(nrPart)

        dataEnergy = readBulkSimulationAvgFile(path + prefix + "AvgEnergy" + param + ".txt")
        dataDriftVel = readBulkSimulationAvgFile(path + prefix + "AvgDriftVelocity" + param + ".txt")
        dataOccupation = readBulkSimulationAvgFile(path + prefix + "valleyOccupation" + param + ".txt")
        
        avgEnergyMean.append(calcValleyOccupationWeightedMean(dataEnergy,dataOccupation))
        avgEnergyValley.append(calcValleywiseMean(dataEnergy))
        avgDriftVelMean.append(np.abs(calcValleyOccupationWeightedMean(dataDriftVel,dataOccupation) *1e2))
        avgDriftVelValley.append(np.abs(calcValleywiseMean(dataDriftVel) * 1e2))
        avgOccValley.append(calcValleywiseMean(dataOccupation) * 100)

    # convert fields to kV / cm
    fields = [float(i) /  1e5 for i in fields]

    # plot valleyOccupation weighted mean
    fig, ax = plt.subplots(2)
    fig.suptitle('Mean characteristics with different electric fields')
    ax[0].plot(fields , avgDriftVelMean)
    ax[0].set_ylabel("avgDriftVel [cm / s]")
    ax[1].plot(fields, avgEnergyMean)
    ax[1].set_ylabel("avgEnergy [eV]")
    ax[1].set_xlabel("applied el. field [kV / cm]")

    if typeFields == AppliedFieldsType.LOW:
        mobility, err = fitMobility(fields, avgDriftVelMean)
        fittedDrift = [i * mobility for i in fields]
        label = "$\mu$ ~ " + "{:.2f}".format(mobility / 1e3) +  " $cm^2 / Vs$" + ", $\sigma$ ~ " + "{:.2f}".format(err / 1e3)
        ax[0].plot(fields, fittedDrift, ":", label="Fitted "+label)
        ax[0].legend()

    # plot mean per valley
    fig, ax = plt.subplots(3)
    nrValleys = len(avgDriftVelValley[0]) - 1
    fig.suptitle('Mean characteristics of each valley with different electric fields')

    for idxValley in range(nrValleys + 1):
        currAvgEnergy, currAvgDriftVel, currAvgOcc = [], [], []
        for idxField in range(len(fields)):
            currAvgEnergy.append(avgEnergyValley[idxField][idxValley])
            currAvgDriftVel.append(avgDriftVelValley[idxField][idxValley])
            currAvgOcc.append(avgOccValley[idxField][idxValley])
        
        ax[0].plot(fields , currAvgDriftVel, label="Valley " + str(idxValley))
        ax[0].set_ylabel("avg. DriftVel [cm / s]")
        ax[1].plot(fields, currAvgEnergy, label="Valley " + str(idxValley))
        ax[1].set_ylabel("avg. Energy [eV]")
        ax[2].plot(fields, currAvgOcc, label="Valley " + str(idxValley))
        ax[2].set_ylabel("avg. valley Occupation [%]")
        ax[2].set_xlabel("applied el. field [kV / cm]")

    ax[0].legend()
    ax[0].grid()
    ax[1].legend()
    ax[1].grid()
    ax[2].legend()
    ax[2].grid()
plt.show()
    
if plotScatterTables:
    suffix = "ScatterMechanism.txt"

    if parameterType == PaperType.PILOTTO or parameterType == PaperType.LI:
        # K - Valley
        acKValley = readScatterMechanismFile(path + "AcousticSLK00" + suffix)
        acIvToKSlashAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLKAcKToKSlash00" + suffix)
        acIvToQQAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLKAcQToQ00" + suffix)
        acIvToQMAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLKAcMToQ00" + suffix)
        acIvToKSlashEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLKAcKToKSlash00" + suffix)
        acIvToQQEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLKAcQToQ00" + suffix)
        acIvToQMEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLKAcMToQ00" + suffix)

        acK = acKValley["rate"] 
        acK += acIvToKSlashAb["rate"] + acIvToQQAb["rate"] + acIvToQMAb["rate"]
        acK += acIvToKSlashEm["rate"] + acIvToQQEm["rate"] + acIvToQMEm["rate"]

        opIvToKAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLKOpGammaToK00" + suffix)
        opIvToKSlashAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLKOpKToKSlash00" + suffix)
        opIvToQQAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLKOpQToQ00" + suffix)
        opIvToQMAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLKOpMToQ00" + suffix)
        opIvToKEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLKOpGammaToK00" + suffix)
        opIvToKSlashEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLKOpKToKSlash00" + suffix)
        opIvToQQEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLKOpQToQ00" + suffix)
        opIvToQMEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLKOpMToQ00" + suffix)

        opK = opIvToKAb["rate"] + opIvToKSlashAb["rate"] + opIvToQMAb["rate"] + opIvToQQAb["rate"]
        opK += opIvToKEm["rate"] + opIvToKSlashEm["rate"] + opIvToQMEm["rate"] + opIvToQQEm["rate"]

        # Q Valley
        acQValley = readScatterMechanismFile(path + "AcousticSLQ01" + suffix)
        acQIvToQ2Em = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQAcQToQ201" + suffix)
        acQIvToQ3Em = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQAcMToQ301" + suffix)
        acQIvToQ4Em = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQAcKToQ401" + suffix)
        acQIvToKEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQAcQToK01" + suffix)
        acQIvToKSlashEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQAcMToKSlash01" + suffix)
        acQIvToQ2Ab = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQAcQToQ201" + suffix)
        acQIvToQ3Ab = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQAcMToQ301" + suffix)
        acQIvToQ4Ab = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQAcKToQ401" + suffix)
        acQIvToKAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQAcQToK01" + suffix)
        acQIvToKSlashAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQAcMToKSlash01" + suffix)

        acQ = acQValley["rate"]
        acQ += acQIvToQ2Em["rate"] + acQIvToQ3Em["rate"] + acQIvToQ4Em["rate"] + acQIvToKEm["rate"] + acQIvToKSlashEm["rate"]
        acQ += acQIvToQ2Ab["rate"] + acQIvToQ3Ab["rate"] + acQIvToQ4Ab["rate"] + acQIvToKAb["rate"] + acQIvToKSlashAb["rate"]

        opQIvToQEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQOpQToQ101" + suffix)
        opQIvToQ2Em = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQOpQToQ201" + suffix)
        opQIvToQ3Em = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQOpMToQ301" + suffix)
        opQIvToQ4Em = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQOpKToQ401" + suffix)
        opQIvToKEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQOpQToK01" + suffix)
        opQIvToKSlashEm = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLQOpMToKSlash01" + suffix)

        opQIvToQAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQOpQToQ101" + suffix)
        opQIvToQ2Ab = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQOpQToQ201" + suffix)
        opQIvToQ3Ab = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQOpMToQ301" + suffix)
        opQIvToQ4Ab = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQOpKToQ401" + suffix)
        opQIvToKAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQOpQToK01" + suffix)
        opQIvToKSlashAb = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLQOpMToKSlash01" + suffix)

        opQ = opQIvToQ2Em["rate"] + opQIvToQ3Em["rate"] + opQIvToQ4Em["rate"] + opQIvToKEm["rate"] + opQIvToKSlashEm["rate"] + opQIvToQEm["rate"]
        opQ += opQIvToQ2Ab["rate"] + opQIvToQ3Ab["rate"] + opQIvToQ4Ab["rate"] + opQIvToKAb["rate"] + opQIvToKSlashAb["rate"] + opQIvToQAb["rate"] 

        plt.title("Scatter Rates")
        plt.xlabel("energy [eV]")
        plt.ylabel("rate [1/s]")
        plt.semilogy(acKValley["energy"], opK + acK, label="K-Valley")
        plt.semilogy(acQValley["energy"], opQ + acQ, label="Q-Valley")
        plt.legend()
        plt.show()
    
    elif parameterType == PaperType.KAASBJERG:
        acLA = readScatterMechanismFile(path + "AcousticSLLA00" + suffix)
        acTA = readScatterMechanismFile(path + "AcousticSLTA00" + suffix)
        acRate = acLA["rate"] + acTA["rate"]
        plt.semilogy(acLA["energy"], acRate, label="acoustic", color="blue")

        # # ZERO ORDER ODP 
        Iv0AbLO = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLLO00" + suffix)
        Iv0EmLO = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLLO00" + suffix)
        Iv0EmHP = readScatterMechanismFile(path + "ZeroInterValleyEmissionSLHomoPolar00" + suffix)
        Iv0AbHP = readScatterMechanismFile(path + "ZeroInterValleyAbsorptionSLHomoPolar00" + suffix)
        ODP0rate = Iv0AbHP["rate"] + Iv0EmHP["rate"] + Iv0AbLO["rate"] + Iv0EmLO["rate"]
        plt.semilogy(Iv0AbHP["energy"], ODP0rate, label="zero order ODP", color="green")

        # FIRST ORDER ODP
        Iv1AbTA = readScatterMechanismFile(path + "FirstInterValleyAbsorptionSLTA00" + suffix)
        Iv1EmTA = readScatterMechanismFile(path + "FirstInterValleyEmissionSLTA00" + suffix)
        Iv1AbLA = readScatterMechanismFile(path + "FirstInterValleyAbsorptionSLLA00" + suffix)
        Iv1EmLA = readScatterMechanismFile(path + "FirstInterValleyEmissionSLLA00" + suffix)
        Iv1AbTO = readScatterMechanismFile(path + "FirstInterValleyAbsorptionSLTO00" + suffix)
        Iv1EmTO = readScatterMechanismFile(path + "FirstInterValleyEmissionSLTO00" + suffix)
        Iv1AbTOGamma = readScatterMechanismFile(path + "FirstInterValleyAbsorptionSLTOGamma00" + suffix)
        Iv1EmTOGamma = readScatterMechanismFile(path + "FirstInterValleyEmissionSLTOGamma00" + suffix)
        ODP1rate =  Iv1AbTA["rate"] + Iv1AbLA["rate"] + Iv1EmTA["rate"] + Iv1EmLA["rate"]
        ODP1rate += Iv1AbTO["rate"] + Iv1EmTO["rate"] + Iv1AbTOGamma["rate"] + Iv1EmTOGamma["rate"]
        plt.semilogy(Iv1AbTO["energy"], ODP1rate, label="first order ODP", color="orange")

        # FRÖHLICH SCATTERING
        froehlichEm = readScatterMechanismFile(path + "froehlichEmissionSL00" + suffix)
        # plt.semilogy(froehlichEm["energy"], froehlichEm["rate"] , ":", color="grey")
        froehlichAb = readScatterMechanismFile(path + "froehlichAbsorptionSL00" + suffix)
        # plt.semilogy(froehlichAb["energy"], froehlichAb["rate"] , ":", color="grey")

        froehlichRate = froehlichEm["rate"] + froehlichAb["rate"]
        plt.semilogy(froehlichEm["energy"], froehlichRate , label="Fröhlich", color="grey")

        # Total
        totalRate = acRate + ODP0rate + ODP1rate + froehlichRate
        plt.semilogy(acLA["energy"], totalRate, label="Total", color="red")

        plt.title("Scattering rates Kaasbjerg")
        plt.xlim(0, 0.25)
        plt.ylim(1e10, 1e14)
        plt.legend()
        plt.show()
