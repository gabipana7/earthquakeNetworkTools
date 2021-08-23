# Import libraries
import arviz as az
import pandas as pd
import matplotlib.pyplot as plt
import math
import seaborn as sns
import numpy as np
import matplotlib.image as mpimg

from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cbook import get_sample_data

# Import the collection from database tool
from sqlCollectDatabaseWithCubes import sqlCollect
from seismicZones import getCondition


# Set the style of the plots to a more beautiful format
#az.style.use('arviz-darkgrid')


# ---------------------------SETUP AND COLLECTION OF QUAKES------------------------------------#
# Which region do you want to analyze ?
region = input('Input region : Vrancea / Romania / California / Italy / Japan : ')

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
quakes = sqlCollect(condition,region,side=[],withCubes=False,energyRelease=False)
# ----------------------------------------------------------------------------------------------#


# -------------------------SCALES FOR THE DOTS BASED ON MAG-----------------------------------#
# Make scales for the dots changing with magnitude
# Set scale for the colors of dots
magnitudes=[int(x) for x in quakes['magnitude']]

# Set scale for the sizes of dots
mean = int(quakes['magnitude'].mean())
max1 = int(quakes['magnitude'].max())
min1 = int(quakes['magnitude'].min())
magnitudesScale=[((x-mean)/(max1-min1))*100 for x in quakes['magnitude']]
# ----------------------------------------------------------------------------------------------#


# ----------------------------------------PLOTS------------------------------------------------#
# 3D Scatter plot of earthquakes locations
# Initialize figure and 3D axes
fig = plt.figure(figsize=(10,10))
ax = Axes3D(fig)

# Scatter the earthquakes and use scales to customize the dots
im = ax.scatter(quakes['longitude'], quakes['latitude'],-quakes['depth'],
           s=magnitudesScale ,c=magnitudes, cmap='Greens', marker='o',
           edgecolor='black', linewidth=0.1, alpha=0.7)

# Label the axes
ax.set_xlabel('Longitude', fontsize=16)
ax.set_ylabel('Latitude', fontsize=16)
ax.set_zlabel('Depth', fontsize=16)

# Show the colorbar used on the magnitude scale
cbar = fig.colorbar(im, ax=ax, orientation='horizontal',fraction=0.046, pad=0.04)
# Label the colorbar
cbar.set_label('Magnitude',fontsize=16)


# ----------------------------------2D MAP PROJECTION----------------------------------------------#
# Choose if you wish to include a projection of the region's map as 2D image under the 3D scatterplot
# You require png images with the maps for this to work
withMap = input('Do you wish the scatter to have a projection of the region map ? True / False : ')

if withMap == 'True':
	
	if region == 'Romania':
		img = mpimg.imread('romania.png')
		# Create lists of points ranging from:
		# For x : longitude 19.8 - 30.2
		xx = np.linspace(19.8,30.2,img.shape[0])
		# For y : latitude 43.5941 - 48.4 
		yy = np.linspace(43.5941,48.4,img.shape[1])

		# Create the x component of the grid for the image plot
		x = np.ndarray((img.shape[0],1))
		# Create y component of the grid for the image plot
		y = np.ndarray((1,img.shape[1]))
		# Create the z component (choose at which depth the image is plotted)
		z = -196 * np.ones(x.shape)

	if region == 'Vrancea':
		img = mpimg.imread('vrancea.png')
		# Create lists of points ranging from:
		# For x : longitude 26 - 27
		xx = np.linspace(26,27,img.shape[0])
		# For y : latitude 45 - 46 
		yy = np.linspace(45,46,img.shape[1])

		# Create the x component of the grid for the image plot
		x = np.ndarray((img.shape[0],1))
		# Create y component of the grid for the image plot
		y = np.ndarray((1,img.shape[1]))
		# Create the z component (choose at which depth the image is plotted)
		z = -200 * np.ones(x.shape)

	if region == 'California':
		img = mpimg.imread('california.png')
		# Create lists of points ranging from:
		# For x : longitude min max
		xx = np.linspace(-121.953,-114.026,img.shape[0])
		# For y : latitude min max
		yy = np.linspace(32,37,img.shape[1])

		# Create the x component of the grid for the image plot
		x = np.ndarray((img.shape[0],1))
		# Create y component of the grid for the image plot
		y = np.ndarray((1,img.shape[1]))
		# Create the z component (choose at which depth the image is plotted)
		z = -32 * np.ones(x.shape)

	if region == 'Italy':
		img = mpimg.imread('italy.png')
		# Create lists of points ranging from:
		# For x : longitude min max
		xx = np.linspace(-6.08,36.02,img.shape[0])
		# For y : latitude min max
		yy = np.linspace(30.61,47.998,img.shape[1])

		# Create the x component of the grid for the image plot
		x = np.ndarray((img.shape[0],1))
		# Create y component of the grid for the image plot
		y = np.ndarray((1,img.shape[1]))
		# Create the z component (choose at which depth the image is plotted)
		z = -645  * np.ones(x.shape)

	if region == 'Japan':
		img = mpimg.imread('japan.png')
		# Create lists of points ranging from:
		# For x : longitude min max
		xx = np.linspace(118.90483,156.68133,img.shape[0])
		# For y : latitude min max
		yy = np.linspace(17.40933,50.42683,img.shape[1])

		# Create the x component of the grid for the image plot
		x = np.ndarray((img.shape[0],1))
		# Create y component of the grid for the image plot
		y = np.ndarray((1,img.shape[1]))
		# Create the z component (choose at which depth the image is plotted)
		z = -700  * np.ones(x.shape)


	# Add any other region here by respecting the format:
		# read the image file using mpimg.imread('file.png')
		# Create lists of points in the shape of the image
		# x coordinate
		# y coordinate
		# z coordinate - the depth at which the 2D image will be displayed, preferably under the scatterplot


	# Assign x values (first component) the proper longitude
	for i in range(len(x)):
	    x[i][0] = float(xx[i])

	# Assign y values (second component) the proper latitude
	for j in range(len(y[0])):
	    y[0][j] = float(yy[j])


	# Plot the overlay with the map under the scattered eathquakes 
	ax.plot_surface(x,y,z, facecolors=img, rstride=1, cstride=1)


	# Show and save
	plt.savefig(f'quakes{region}_{magMin}mag{magMax}_from{year}_with2DMap.png')
	plt.show()
# ----------------------------------------------------------------------------------------------#	


else: 
	# Show and save
	plt.savefig(f'quakes{region}_{magMin}mag{magMax}_from{year}.png')
	plt.show()
# ----------------------------------------------------------------------------------------------#	