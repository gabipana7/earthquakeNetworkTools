import mysql.connector
import pandas as pd
import math

# Have the credentials in a separate file in order not to leak them online 
from credentials import host,user,password

def sqlCollect(condition,region,side=[],withCubes=True,energyRelease=False):
	mydb = mysql.connector.connect(
	host=host,
	user=user,
	password=password,
	database="cutremure",
	auth_plugin='mysql_native_password'
	)

	mycursor = mydb.cursor()

	# Call the database with the set conditions
	mycursor.execute(condition)

	# Append the database to the result object array
	result = mycursor.fetchall()

	# Initialize our lists    
	date= []
	latitude = []
	longitude = []
	depth = []
	magnitude = []

	# Append the result of sql into our lists
	for i in range(len(result)):
		date.append(result[i][0])
		latitude.append(result[i][1])
		longitude.append(result[i][2])
		depth.append(result[i][3])
		magnitude.append(result[i][4])
	  
	# Force the results into float type
	latitude = [float(item) for item in latitude]
	longitude = [float(item) for item in longitude]
	depth = [float(item) for item in depth]
	magnitude = [float(item) for item in magnitude]
	
	# If withCubes == TRUE => collect database and create cube indexes
	if withCubes:
		# Create cubes indexes for all the dimensions
		# ------------------------------------------

		# Set min max for the dimensions
		minLat = min(latitude)
		maxLat = max(latitude)
		minLong = min(longitude)
		maxLong = max(longitude)
		minDepth = min(depth)
		maxDepth = max(depth)


		# Find the number of cubes, given the size of the cube in km ( previously inputted )
		# Depending on region, the size of one degree lat/long differs:

		# IN VRANCEA 1 deg Lat = 111km / 1 deg Long = 79km
		if region=='Vrancea' or region=='Romania':
			longitudeInKm = 79
		
		# IN CALI 1 deg Lat = 111km / 1 deg Long = 94km
		if region=='California':
			longitudeInKm = 94

		# IN ITALY 1 deg Lat = 111km / 1 deg Long = 84km
		if region=='Italy':
			longitudeInKm = 84

		# IN JAPAN 1 deg Lat = 111km / 1 deg Long = 91km
		if region=='Japan':
			longitudeInKm = 91

		x = round((maxLat-minLat)*111 / side)
		y = round((maxLong-minLong)*longitudeInKm / side)
		z = round((maxDepth-minDepth) / side)

		# Create cubes indexes for all the dimensions 
		xLatitude = [math.floor((i-minLat)*x/(maxLat-minLat))+1 for i in latitude]
		yLongitude = [math.floor((i-minLong)*y/(maxLong-minLong))+1 for i in longitude]
		zDepth = [math.floor((i-minDepth)*z/(maxDepth-minDepth))+1 for i in depth ]
		cubeIndex=[]
		# and a general cubeindex for graph formation
		for i in range(len(xLatitude)):
			cubeIndex.append(int(str(xLatitude[i])+str(yLongitude[i])+str(zDepth[i])))
			#cubeIndex.append((xLatitude[i]-1)*x+(yLongitude[i]-1)*y+zDepth[i])


		# Input coordinates of every cubeIndex
		cubeLatitude=[round(minLat + (side/111)*(i-1+1/2),4) for i in xLatitude]
		cubeLongitude=[round(minLong + (side/longitudeInKm)*(i-1+1/2),4) for i in yLongitude]
		cubeDepth=[round(minDepth + (side)*(i-1+1/2),4) for i in zDepth]

		# if energyReleas = True, add a column with the energy release converted from magnitude 
		if energyRelease:
			energyRelease=[math.pow(10, 5.24 + 1.44*i) for i in magnitude]

			d = {'date':date,'latitude':latitude,'longitude':longitude,
				'depth':depth,'magnitude':magnitude, 'energyRelease':energyRelease,'x':xLatitude,
				'y':yLongitude,'z':zDepth,'cubeIndex':cubeIndex,
				'cubeLatitude':cubeLatitude,'cubeLongitude':cubeLongitude,
				'cubeDepth':cubeDepth}
			return pd.DataFrame(data=d)


		else:
			# Create DataFrame
			d = {'date':date,'latitude':latitude,'longitude':longitude,
				'depth':depth,'magnitude':magnitude,'x':xLatitude,
				'y':yLongitude,'z':zDepth,'cubeIndex':cubeIndex,
				'cubeLatitude':cubeLatitude,'cubeLongitude':cubeLongitude,
					'cubeDepth':cubeDepth}
			return pd.DataFrame(data=d)
	
	# if withCubes==False , do not create cubes (used for simple scatterplot visualization of quakes)
	else:
	
		if energyRelease:
			energyRelease=[math.pow(10, 5.24 + 1.44*i) for i in magnitude]

			d = {'date':date,'latitude':latitude,'longitude':longitude,'depth':depth,
			     'magnitude':magnitude,'energyRelease':energyRelease}
			return pd.DataFrame(data=d)
		else:
			d = {'date':date,'latitude':latitude,'longitude':longitude,'depth':depth,'magnitude':magnitude}
			return pd.DataFrame(data=d)
