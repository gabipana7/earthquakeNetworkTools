import arviz as az
import matplotlib
import numpy as np
import datetime

import matplotlib.pyplot as plt
import math
import collections

from sqlCollectDatabaseWithCubes import sqlCollect
from seismicZones import getCondition

# POORLY OPTIMIZED CODE . NEEDS WORK !!!!!!!!!!! 

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


fig = plt.figure(figsize=(8,8))
for T,idx in zip((1,7,14,30),(221,222,223,224)):

    M2={}
    Cm2={}
    distances=[]

    for i in range(len(quakes['energyRelease'])-1):
        try:
            j = quakes.loc[quakes.date >= quakes['date'][i] + datetime.timedelta(days=T)].index[0]
            distance=round(math.sqrt((quakes['xLatitude'][i]-quakes['xLatitude'][j])**2 + 
                                            (quakes['yLongitude'][i]-quakes['yLongitude'][j])**2 +
                                            (quakes['zDepth'][i]-quakes['zDepth'][j])**2))

            if distance in M2.keys():
                M2[distance] += quakes['energyRelease'][i]*quakes['energyRelease'][j]
                Cm2[distance] += 1

            else:
                M2[distance] = quakes['energyRelease'][i]*quakes['energyRelease'][j]
                Cm2[distance] = 1

        except:
            break


    for sub in Cm2.keys():
        M2[sub] = M2[sub] / Cm2[sub]

    ax = fig.add_subplot(idx) 

    ax.scatter( Cm2.keys(), M2.values() )

    #ax.set_xscale('log')
    ax.set_yscale('log') 

    ax.set_xlabel('r')
    ax.set_ylabel('M2')
    ax.set_title(f'{T}days',fontsize = 8)

plt.tight_layout()
plt.suptitle(f'{region} - {side}km - {magMin}<mag<{magMax}',fontsize=10)
plt.savefig(f'M2_{region}_{side}km_{magMin}<mag<{magMax}.png')