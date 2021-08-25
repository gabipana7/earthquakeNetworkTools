import arviz as az
import matplotlib
import numpy as np
import datetime

from datetime import datetime  
from datetime import timedelta 

import matplotlib.pyplot as plt
import math
import collections

from statsmodels.graphics.tsaplots import plot_acf
from sqlCollectDatabaseWithCubes import sqlCollect
from seismicZones import getCondition
from graphCreation import*


# ---------------------------SETUP AND COLLECTION OF QUAKES------------------------------------#
# Which region do you want to analyze ?
region = input('Input region : Vrancea / Romania / California / Italy / Japan : ')
side = int(input('Insert the cube side lenght 5 / 10 / 20 km : '))

# The side of the cubes that you split the region in 
#side = int(input('Input side of the cube split in km 5 / 10 / 20 : '))

# The condition for the SQL collection
condition,year = getCondition(region)

# Magnitude windows for the condition that collects the database through mySQL
magMin = int(input('Input minimum magnitude: '))
magMax = int(input('Input maximum magnitude: '))
# Add to condition
condition+=f" AND `magnitude`>={magMin} AND `magnitude`<={magMax}"

# Collect the earthquakes
quakes = sqlCollect(condition,side,region,witchCubes=False,energyRelease=True)

# Calculate the total time frame from the first quake to the last quake recorded
timeWindow = max(quakes['date']) - min(quakes['date'])

# Calculate all the intervals between 2 consecutive earthquakes
deltas = [quakes['date'][i+1]-quakes['date'][i] for i in range(len(quakes['date'])-1)]

# Select the biggest time interval
dt = max(deltas)

# Add a day to the biggest time interval
delta = timedelta(days=1) + dt

# Split the total timeframe by this interval in order to get a certain number of equal windows
windows = round(timeWindow/delta)

# Place all the earthquakes in the respective windows together
interval = {}
for i in range(windows):
    interval[i] = [quakes['date'][0] + i*delta,quakes['date'][0] + (i+1)*delta]

# Calculate the energy release of the quakes in each window
quakesList = {}
energyInterval = dict.fromkeys(interval.keys(),0)
for i in range(len(interval)):
    quakesList[i] =list(set(quakes.index[interval[i][0]<= quakes['date']].tolist()).intersection(quakes.index[quakes['date'] < interval[i][1]].tolist()))
    for k in quakesList[i]:
        energyInterval[i] += quakes['energyRelease'][k]

# Append these values 
data = list(energyInterval.values())
values = np.array(data)


# Plot the autocorrelation function for these energies. 
f = plt.figure()
ax = f.add_subplot(111)
plot_acf(values,lags=len(energyInterval)-1,ax=ax)

plt.xlabel(r'k (lags)')
plt.ylabel(r'C(k)')

plt.text(0.7,0.9,f'dt = {dt}\nno. events = {len(quakes["date"])}\n{magMin}<mag<{magMax}',
     transform=ax.transAxes)

plt.title(f'Temporal corelation for {region}',fontsize = 14)
plt.savefig(f'./temporalCorrelation{region}/TemporalCorrelation_{region}_{magMin}<mag<{magMax}.png')
plt.close()


# WITHOUT ZERO to better see plot feature
f = plt.figure()
ax = f.add_subplot(111)
plot_acf(values,lags=len(energyInterval)-1,ax=ax,zero=False)

plt.xlabel(r'k (lags)')
plt.ylabel(r'C(k)')

plt.text(0.7,0.9,f'dt = {dt}\nno. events = {len(quakes["date"])}\n{magMin}<mag<{magMax}',
     transform=ax.transAxes)

plt.title(f'Temporal corelation for {region}',fontsize = 14)
plt.savefig(f'./temporalCorrelation{region}/TemporalCorrelation_{region}_{magMin}<mag<{magMax}_noZero.png')
plt.close()