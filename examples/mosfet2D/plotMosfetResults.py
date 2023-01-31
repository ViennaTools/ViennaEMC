from emcPlottingFiles.plotParticleCharacteristics import *
from emcPlottingFiles.plotCurrentCharacteristics import *
from emcPlottingFiles.plotGrid import *
from emcPlottingFiles.readResultFile import *

"""

This file shows the results of the example "mosfet2D".

Note: To use this file you have to install emcPlottingFiles. For 
instructions on how to do that look at the main README.md.

To use this file, run the example "mosfet2D", then adapt the
"path" (path to the files).

Additionally, some charts only appear ones the previous charts
are closed. 

"""

######## Adaptable Parameters #####################################################
path = "../../build/"       # path to result files (adapt if needed)
prefix = "mosfet"           # file name prefix
suffix = ".txt"             # file name suffix
particleType = "Electrons"  # used particle type name
is3D = False                # determines if 3D simulation is used

includeParameterInFileName = True
Vd = 1 # in V
Vg = 1 # in V

# plots current characteristics of the two pnJunction contacts
plotCurrentCharacteristics = True
labelsOfContacts = ["Substrate", "Source", "Gate", "Drain"] # label of contacts

# plots particle characteristics
plotParticleCharacteristics = True

# plots grid characteristics of the two pnJunction contacts
plotEqGridCharacteristics = True
plotAvgGridCharacteristics = True

# plots scatter tables 
plotScatterTables = True
# select the scatter rates that should be plotted (adapt if needed)
plotAcouticRate = True                  # rates from emcAcousticScatterMechanism.hpp
plotZeroOrderIntervalleyRate = True     # rates from emcZeroOrderIntervalleyScatterMechanism.hpp
plotFirstOrderIntervalleyRate = False   # rates from emcFirstOrderIntervalleyScatterMechanism.hpp
plotCoulombRate = True                  # rates from emcCoulombScatterMechanism.hpp
# select index of the valley and doping region of interest (adapt if needed)
idxValley = 0   # for silicon only one valley is used, possible values: 0
idxRegion = 2   # for mosfet 4 doping regions are added, possible values: 0,1,2,3
#################################################################################

if(includeParameterInFileName):
    prefix += "Vd" + str(int(Vd * 1000)) + "Vg" + str(int(Vg*1000))

if(plotCurrentCharacteristics):
    print("Plotting Current Characteristics ...")
    dataCurrent = readCurrentFile(path + prefix + particleType + "Current" + suffix)
    
    # get nr of contacts (counts all kind of contacts ohmic / gate / schottky)
    nrContacts = getNrContacts(dataCurrent)
    print("\tNr of Contacts: ", nrContacts)

    # print the final current for each contact (here also labels can be used)
    print("\tFinal Mean Current:")
    printFinalCurrent(dataCurrent, contactLabel=labelsOfContacts)

    # plot current and netto nrParticles for each contact
    plotCurrentOverview(dataCurrent, title="Current Overview at each Contact", contactLabel=labelsOfContacts)

    # plot current of single contact (can also be used to compare different results)
    plt.figure()
    plotCurrent(dataCurrent, 1, "Source", "green")
    plotCurrent(dataCurrent, 3, "Drain", "red")
    plt.title("Current of Source and Drain")

    plt.show()


if(plotParticleCharacteristics):
    print("Plotting Particle Characteristics ...")

    # read the final and initial particle data
    partDataEq, maxPos = readParticleFile(path + prefix + particleType + "Eq" + suffix, return_maxPos=True)
    partDataFinal = readParticleFile(path + prefix + particleType + "Final" + suffix)

    # show  initial and final particle positions
    plotParticlePositions(partDataEq, maxPos=maxPos, title="Initial Electron Positions")
    plotParticlePositions(partDataFinal, maxPos=maxPos, title="Final Electron Positions")
    # compare initial and final particle positions
    plotMultipleParticlePositions([partDataEq, partDataFinal], ["Initial", "Final"], ["green", "red"], is3D, maxPos)
    plt.title("Comparison of Initial and Final Electron Positions")
    plt.show()

    # plot histogramm / bar for different particle characteristics to show distribution of them
    plotParticleMeanCharacteristics(partDataEq, title="Initial Particle Characteristic Distributions")
    plotParticleMeanCharacteristics(partDataFinal, title="Final Particle Characteristic Distributions")
    plt.show()

    # plot particle positions with their color showing their energy
    plotParticleCharacteristicAtPositions(partDataFinal, column="energy", maxPos = maxPos, title="Final Spatial Energy distribution")
    # plot particle positions with their color showing their subvalley
    plotParticleCharacteristicAtPositions(partDataFinal, column="subValley", maxPos=maxPos, title="Final Spatial SubValley Distribution")
    plt.show()
    # other characteristics that could be interesting here: "kX", "kY", ...

    # show color of specific directions
    plotDirectionColours()
    # plot initial (random) wave-vector directions
    plotParticleWaveVecDirections2D(partDataEq,maxPos=maxPos, title="Initial Wave-Vector Directions")
    # plot the directions of the final wave-vectors
    plotParticleWaveVecDirections2D(partDataFinal,maxPos=maxPos, title="Final Wave-Vector Directions")
    plt.show()


if plotEqGridCharacteristics:
    print("Plotting Equilibrium (Initial) Grid Characteristics ...")
    # plot average final potential
    dataPotential = readGridFile(path + prefix + "PotentialEq" + suffix)
    plotGrid(dataPotential, "Potential [V]", addColorBar=False)
    plt.tight_layout()

    # plot average electron concentratiom
    dataElConc = readGridFile(path + prefix + particleType + "ConcEq" + suffix)
    plotGrid(dataElConc, r"Conc. [$m^{-3}$]", addColorBar=False)
    plt.tight_layout()

    # plot initial electric field
    dataEFIeldX = readGridFile(path + prefix + "EFieldXEq" + suffix)
    plotGrid(dataEFIeldX, "x-field [V/m]", addColorBar=False)
    plt.tight_layout()

    dataEFIeldX = readGridFile(path + prefix + "EFieldYEq" + suffix)
    plotGrid(dataEFIeldX, "y-field [V/m]", addColorBar=False)
    plt.tight_layout()
    
    plt.show()

if plotAvgGridCharacteristics:
    print("Plotting Average (Final) Grid Characteristics ...")
    # plot average final potential
    dataPotential = readGridFile(path + prefix + "PotentialAvg" + suffix)
    plotGrid(dataPotential, "Potential [V]", addColorBar=False)
    plt.tight_layout()
    # plt.title("Avg. Potential")

    # plot average electron concentratiom
    dataElConc = readGridFile(path + prefix + particleType + "ConcAvg" + suffix)
    plotGrid(dataElConc, r"Conc. [$m^{-3}$]", addColorBar=False)
    plt.tight_layout()    
    plt.show()

if plotScatterTables:
    # plot scatter mechanisms
    if plotAcouticRate:
        acRate = readScatterMechanismFile2(path, "Acoustic", idxValley, idxRegion)
        plt.semilogy(acRate["energy"], acRate["rate"], label="acoustic", color="red")
    if plotZeroOrderIntervalleyRate:
        ivZeroFRate = readScatterMechanismFile2(path, "ZeroInterValleyEmissionF", idxValley, idxRegion)
        ivZeroGRate = readScatterMechanismFile2(path, "ZeroInterValleyEmissionG", idxValley, idxRegion)
        plt.semilogy(ivZeroFRate["energy"], ivZeroFRate["rate"], label="f-type zero order emission")
        plt.semilogy(ivZeroGRate["energy"], ivZeroGRate["rate"], label="g-type zero order emission")
        
        ivZeroFRate = readScatterMechanismFile2(path, "ZeroInterValleyAbsorptionF", idxValley, idxRegion)
        ivZeroGRate = readScatterMechanismFile2(path, "ZeroInterValleyAbsorptionG", idxValley, idxRegion)
        plt.semilogy(ivZeroFRate["energy"], ivZeroFRate["rate"], label="f-type zero order absoption")
        plt.semilogy(ivZeroGRate["energy"], ivZeroGRate["rate"], label="g-type zero order absorption")
    if plotFirstOrderIntervalleyRate:
        ivFirstFRate = readScatterMechanismFile2(path, "FirstInterValleyEmissionF", idxValley, idxRegion)
        ivFirstGRate = readScatterMechanismFile2(path, "FirstInterValleyEmissionG", idxValley, idxRegion)
        plt.semilogy(ivFirstFRate["energy"], ivFirstFRate["rate"], label="f-type first order emission")
        plt.semilogy(ivFirstGRate["energy"], ivFirstGRate["rate"], label="g-type first order emission")

        ivFirstFRate = readScatterMechanismFile2(path, "FirstInterValleyAbsorptionF", idxValley, idxRegion)
        ivFirstGRate = readScatterMechanismFile2(path, "FirstInterValleyAbsorptionG", idxValley, idxRegion)
        plt.semilogy(ivFirstFRate["energy"], ivFirstFRate["rate"], label="f-type first order absorption")
        plt.semilogy(ivFirstGRate["energy"], ivFirstGRate["rate"], label="g-type first order absorption")
    if plotCoulombRate:
        coulombRate = readScatterMechanismFile2(path, "Coulomb", idxValley, idxRegion)
        plt.semilogy(coulombRate["energy"], coulombRate["rate"], label= "coulomb")

    plt.title("Scatter Rates")
    plt.xlabel("energy [eV]")
    plt.ylabel("rate [1 / s]")
    plt.legend()
    plt.show()



    