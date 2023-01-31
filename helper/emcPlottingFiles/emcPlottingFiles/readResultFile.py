import pandas as pd
import numpy as np

"""
This file contains functions that can read in the results-files created by EnsembleMC.
"""


# reads particle file of EnsembleMC, returns dataframe
def readParticleFile(filePath, return_maxPos = False):
    try:
        # read maximal possible Particle possition
        maxPos = pd.read_csv(filePath, sep=" ", header=None, nrows=1).to_numpy()[0]
        if maxPos.shape[0] == 2:
            is3D = False
        else:
            is3D = True

        # read particle data + add header
        data = pd.read_csv(filePath, sep=" ", header=None, skiprows=[0], skip_blank_lines=True)
        if data.shape[1] == 4 or data.shape[1] == 3: # non - moveable particle
            header = ["idx", "posX", "posY"]
        else:
            if data.shape[1] == 9:
                header = ["idx", "posX", "posY", "kX", "kY", "kZ", "energy", "subValley", "valley"]
            else:
                header = ["idx", "posX", "posY", "kX", "kY", "kZ", "energy", "subValley", "valley", "tau"]
        if is3D:
            header.insert(3, "posZ")
        if len(header) < data.shape[1]:
            header.extend(["potential","forceX", "forceY", "forceZ"])
        data.columns = header
    except:
        print("No / Invalid data file was given.")
        exit()
    print("Read particle file from: \t", filePath)
    if not return_maxPos:
        return data
    else:
        return data, maxPos


# reads a file to which a 2D or 3D emcGrid-Object was written
# returns the grid 
def readGridFile(filePath):
    try:
        extent = pd.read_csv(filePath, sep=" ", header=None, nrows=1).to_numpy()
        extent[0] = extent[0][::-1]
        data = pd.read_csv(filePath, sep=" ", header=None, skiprows=[0], skip_blank_lines=True).to_numpy()
        if extent[0].shape[0] == 3:
            data = np.reshape(data, extent[0], 'C')
        return data
    except:
        print("No / Invalid data file was given.")


# reads a file to which a 2D - array of ViennaWD was written to
# returns the grid (in EMC order)
def readGridFileViennaWD(filePath):
    try:
        data = pd.read_csv(filePath, sep=" ", header=None, index_col=3, dtype=np.float64).to_numpy()
        dataShape = [int(np.max(data[:,0])) + 1, int(np.max(data[:,1])) + 1]
        data = data[:,2].reshape(dataShape)
        data = np.transpose(data)
        return data
    except:
        print("No / Invalid data file was given.")

# read ensembleMC current file
def readCurrentFile(filePath):
    try: 
        data = pd.read_csv(filePath, sep=" ", header=None)
        nrContacts = (int)((data.shape[1] - 1) / 2)
        
        #add header
        header = ["time"]
        for idx in range(nrContacts):
            header.append("netParContact" + str(idx))
        for idx in range(nrContacts):
            header.append("currentContact" + str(idx))
        data.columns = header
        print("Read current file from: \t", filePath)
        return data
    except:
        print("No / Invalid data file was given.") 


# read a cur_from_charge.csv file from ViennaWD
def readCurrentFileViennaWD(filePath):
    try: 
        data = pd.read_csv(filePath, sep="\t", header=None)
        header = ["time", "currentContact0", "currentContact1", "netParContact0", "netParContact1",]
        data.columns = header
        data["time"] = data["time"] / 1e12  # change unit to seconds
        data["currentContact1"] *= -1
        print("Read current file from: \t", filePath)
        return data
    except:
        print("No / Invalid data file was given.") 

# read scatterMehanism File (output from scatterHandler)
def readScatterMechanismFile(filePath):
    try: 
        data = pd.read_csv(filePath, sep=" ", header=None, skip_blank_lines=True )
        data.columns = ["energy", "rate"]
        return data
    except:
        print("No / Invalid data file was given.") 

def readScatterMechanismFile2(path, scatterMechName, idxValley, idxRegion):
    filename = scatterMechName + str(idxRegion) + str(idxValley) + "ScatterMechanism.txt"
    return readScatterMechanismFile(path + filename)

# read average file of bulk simulation 
# (usable for valley occupation, avg energy, avg drift velocity)
def readBulkSimulationAvgFile(filePath):
    try: 
        data = pd.read_csv(filePath, sep=" ", header=None, skip_blank_lines=True )
        nrValleys = data.shape[1] - 1
        cols = ["time"]
        for i in range(nrValleys):
            cols.append("valley_" + str(i))
        data.columns = cols
        return data
    except:
        print("No / Invalid data file was given.") 