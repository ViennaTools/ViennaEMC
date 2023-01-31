import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize

""" 
File contains functionality to plot current characteristics. 
(read current-files in with "readCurrentFile" or "readCurrentFileViennaWD")
"""

# plot current and netto nr Particles of all contacts
# contactLabel = list of labels for each contact
def plotCurrentOverview(data, contactLabel=[], title="Current Overview"):
    nrContacts = getNrContacts(data)
    fig, axs = plt.subplots(2,1)
    fig.suptitle(title)
    if not contactLabel:
        for idxContact in range(nrContacts):
            contactLabel.append("Contact " + str(idxContact))
    if len(contactLabel) != nrContacts:
        print("contactLabel has to have the right length: ", nrContacts)
    else:
        for idxContact in range(nrContacts):
            axs[0].plot(data["time"], data["currentContact" + str(idxContact)] *1000, label = contactLabel[idxContact])
        axs[0].set_ylabel("current [mA]")
        axs[0].set_xlabel("time [s]")
        axs[0].grid()
        axs[0].legend()

        for idxContact in range(nrContacts):
            axs[1].plot(data["time"], data["netParContact" + str(idxContact)], label = contactLabel[idxContact])
        axs[1].set_ylabel("netto nr particles (removed - inserted)")
        axs[1].set_xlabel("time [s]")
        axs[1].legend()

# plot the current of one specific contact
def plotCurrent(data, nrContact, label, color, style = '-'):
    nrContacts = (int)((data.shape[1] - 1) / 2)
    if nrContact < nrContacts:
        plt.plot(data["time"], data["currentContact" + str(nrContact)] *1000, style, color = color, label = label)
    plt.xlabel("time [s]")
    plt.ylabel("current [mA]")
    plt.tight_layout()
    plt.legend()

# helper function, returns nr of contacts
def getNrContacts(data):
    return (int)((data.shape[1] - 1) / 2)

# prints the final mean current of each contact
def printFinalCurrent(data, contactLabel=[]):
    nrContacts = (int)((data.shape[1] - 1) / 2)
    if not contactLabel:
        for idxContact in range(nrContacts):
            contactLabel.append("Contact " + str(idxContact))
    if len(contactLabel) != nrContacts:
        print("contactLabel has to have the right length: ", nrContacts)
    else:
        for idxContact in range(nrContacts):
            print("\t\t" + contactLabel[idxContact] + ":\t\t" + "{:.4f}".format(data["currentContact" + str(idxContact)].iloc[-1] * 1000) + "\tmA")
    

    
