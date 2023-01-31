import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize

""" 
File contains functionality to plot particle characteristics, read in with readParticleFile().
"""

# plot particle positions of one dataset
def plotParticlePositions(data, maxPos = np.array([]), title=""):        
    if "posZ" in data.columns: # 3d - scatter - plot
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(data["posX"], data["posY"], data["posZ"], alpha=0.5)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        if np.any(maxPos):
            ax.set_xlim(0, maxPos[0])
            ax.set_ylim(0, maxPos[1])
            ax.set_zlim(0, maxPos[2])
    else: # 2d - scatter - plot
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.scatter(data["posX"], data["posY"], s=0.8, alpha=0.5)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        if np.any(maxPos):
            ax.set_xlim(0, maxPos[0])
            ax.set_ylim(0, maxPos[1])
    ax.set_title(title)


# plot particle positions of multiple datasets
# dataList = list of the datasets to plot
# lables = lables of the datasets
# color = colors for the different datasets
# is3D = bool that tells if the data to plot is 3D
def plotMultipleParticlePositions(dataList, labels, color, is3D, maxPos = np.array([])):  
    fig = plt.figure()
    if is3D:
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        if np.any(maxPos):
            ax.set_xlim(0, maxPos[0])
            ax.set_ylim(0, maxPos[1])
            ax.set_zlim(0, maxPos[2])
    else:
        ax = fig.add_subplot(111)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        if np.any(maxPos):
            ax.set_xlim(0, maxPos[0])
            ax.set_ylim(0, maxPos[1])

    for (data, label, c) in zip(dataList, labels, color):      
        if "posZ" in data.columns:
            if not is3D:
                print("Can't plot 2D simulation result in 3D.")
                break
            ax.scatter(data["posX"], data["posY"], data["posZ"], s=0.5, alpha=0.5, label=label)

        else:
            if is3D:
                print("Can't plot 3D simulation result in 2D.")
                break
            ax.scatter(data["posX"], data["posY"], s=0.8, alpha=0.5, label=label, color=c)
        plt.legend()


# plots particle positions and sets the particle color with the value of the given column
def plotParticleCharacteristicAtPositions(data, column="energy", maxPos = np.array([]), title=""):
    if column in data.columns:
        fig = plt.figure()
        if "posZ" in data.columns:
            ax = fig.add_subplot(111, projection='3d')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            if np.any(maxPos):
                ax.set_xlim(0, maxPos[0])
                ax.set_ylim(0, maxPos[1])
                ax.set_zlim(0, maxPos[2])
            scat = ax.scatter(data["posX"], data["posY"], data["posZ"], c=data[column], alpha=0.5)
        else:
            ax = fig.add_subplot(111)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            if np.any(maxPos):
                ax.set_xlim(0, maxPos[0])
                ax.set_ylim(0, maxPos[1])
            scat = ax.scatter(data["posX"], data["posY"], c=data[column], s=0.5, alpha=0.5)
        ax.set_title(title)
        fig.colorbar(scat, shrink=0.5, aspect=5)

    else:
        print("column " + str(column) + "not given in dataset.")


# plots the wave-vector directions at the particle positions (for 2D)
# TODO 3D version of this function
def plotParticleWaveVecDirections2D(data, label = "", maxPos = np.array([]), title="Directions of Wave-Vector"):
    plt.figure()
    # set colormap for specific directions
    ph = np.linspace(0, 2*np.pi, 20)
    u = np.cos(ph)  
    v = np.sin(ph)
    colors = np.arctan2(u, v)
    norm = Normalize()
    norm.autoscale(colors)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.quiver(data["posX"], data["posY"], data["kX"], data["kY"], color = cm.inferno(norm(np.arctan2(data["kX"], data["kY"]))))
    if np.any(maxPos):
        plt.xlim(0, maxPos[0])
        plt.ylim(0, maxPos[1])
    if label:
        plt.legend()
    plt.title(title)


# plots the colours of the specific wave-vector directions (helper for prev. function)
def plotDirectionColours():
    ph = np.linspace(0, 2*np.pi, 20)
    x = np.cos(ph)
    y = np.sin(ph)
    u = np.cos(ph)  
    v = np.sin(ph)
    colors = np.arctan2(u, v)
    norm = Normalize()
    norm.autoscale(colors)
    plt.figure(figsize=(6, 6))
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.quiver(x, y, u, v, color=cm.inferno(norm(colors)), angles='xy',scale_units='xy', scale=1, pivot='mid')


# plots particle mean characteristics
def plotParticleMeanCharacteristics(data, title="Distribution of particle characteristics"):
    fig, axs = plt.subplots(2,3)
    fig.suptitle(title)

    colToPlot = ["kX", "kY", "kZ", "energy", "subValley", "valley"]
    plotType = ["hist", "hist", "hist", "hist", "bar", "bar"]
    idx0, idx1 = 0, 0
    for (col, type) in zip(colToPlot, plotType):
        if type == "hist":
            axs[idx0, idx1].hist(data[col], bins=25)
        else:
            curCol = data[col].to_numpy()
            values, counts = np.unique(curCol, return_counts=True)
            axs[idx0, idx1].set_xticks(values)
            axs[idx0, idx1].bar(values,counts)

        axs[idx0, idx1].set_title(col)
        axs[idx0, idx1].set_ylabel("nr. Particles")
    
        idx1 += 1
        if idx1 == 3:
            idx0 += 1
            idx1 = 0


# TODO plot particle force (if it is given in data)        
# def plotParticleForce(data, maxPos=np.array([]), factor=1e14):
#     if "forceX" in data.columns:
#         #create colormap
#         x = data["posX"].to_numpy()
#         y = data["posY"].to_numpy()
#         z = data["posZ"].to_numpy()
#         u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
#         v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
#         w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
#             np.sin(np.pi * z))

#         c = np.arctan2(v, u)
#         c = (c.ravel() - c.min()) / c.ptp()
#         c = np.concatenate((c, np.repeat(c, 2)))
#         # Colormap
#         c = plt.cm.inferno(c)  

#         fig = plt.figure()
#         ax = fig.add_subplot(111, projection='3d')
#         ax.quiver(data["posX"], data["posY"], data["posZ"], -data["forceX"]*factor, -data["forceY"]*factor, -data["forceZ"]*factor, color=c)
#         ax.set_xlabel('X')
#         ax.set_ylabel('Y')
#         ax.set_zlabel('Z')
#         if np.any(maxPos):
#             ax.set_xlim(0, maxPos[0])
#             ax.set_ylim(0, maxPos[1])
#             ax.set_zlim(0, maxPos[2])
#         ax.set_title("Forces on particles")
#     else:
#         print("Data has no column for force.")
    

    
