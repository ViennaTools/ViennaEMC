import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

""" 
function calculates valleyOccupation-weighted mean (first calculates
mean of each valley and then sums those mean up, weighted by the valley)

input:
    data              read-in data from outputfile of bulkSimulation
    occupationData    read-in data from valleyOccupation output file
    startAvgTime      values that are recorded after that time are used 
                      for the averaging
"""
def calcValleyOccupationWeightedMean(data, occupationData, startAvgTime=1e-12 ):
    if(data["time"].iat[-1] < startAvgTime):
        print(" PROBLEM: Can't calculate mean: startAvgTime > simulationTime -> Change startAvgTime.")
    data = data[data["time"] > startAvgTime]
    occupationData = occupationData[occupationData["time"] > startAvgTime]

    mean = 0
    for i in range(0, data.shape[1] - 1):
        tmp = data["valley_" + str(i)] * occupationData["valley_" + str(i)]
        mean += tmp.mean()
    return mean

""" 
function calculates mean for each valley

input:
    data              read-in data from outputfile of bulkSimulation
    startAvgTime      values that are recorded after that time are used 
                      for the averaging
"""
# calculates the mean of a characteristic for each valley
def calcValleywiseMean(data, startAvgTime=1e-12 ):
    if(data["time"].iat[-1] < startAvgTime):
        print(" PROBLEM: Can't calculate mean: startAvgTime > simulationTime -> Change startAvgTime.")

    data = data[data["time"] > startAvgTime]
    mean = data.mean(axis=0)
    return mean[1:]


def linearFunc(x, a):
    return a * x

# fits the given fields and dirft velocities to a linear function
# 
# parameter:
#   fields          list of applied electric fields
#   avgDriftVel     list of avgDriftVel for each field
#   maxFieldUsed    maximal electricField that is used for fitting
#   
# return:
#   popt        mobility
#   perr        fitting error of mobility (std-dev)
def fitMobility(fields, avgDriftVel):
    popt, pcov = curve_fit(linearFunc, fields, avgDriftVel)
    perr = np.sqrt(np.diag(pcov)) # std-dev
    return popt[0], perr[0]

