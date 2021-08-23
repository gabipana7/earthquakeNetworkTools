import arviz as az
import matplotlib
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pickle

from graphCreation import*
from sqlCollectDatabaseWithCubes import sqlCollect


# Set the style of the plots to a more beautiful format
az.style.use('arviz-darkgrid')


region = input('Input region : Vrancea / Romania / California / Italy / Japan : ')
edgeWeight = input('Input if connectivity takes into account edge weight : True / False : ')


# VRANCEA
if region=='Vrancea':
	condition=(f"SELECT * FROM romplus WHERE `dateandtime`>='1976-01-01 00:00:00'"
		   f" AND `latitude`>=45.2 AND `latitude`<=46 AND `longitude`>=26" 
		   f" AND `longitude`<=27 AND `depth`>=50 AND `depth`<=200")
	
# ROMANIA
if region=='Romania':
	condition=f"SELECT * FROM romplus WHERE `dateandtime`>='1976-01-01 00:00:00'"

# CALIFORNIA
if region=='California':
	condition=(f"SELECT * FROM california WHERE `dateandtime`>='1984-01-01 00:00:00'"
				f" AND (`magtype` LIKE 'l' OR `magtype` LIKE 'w')")

# ITALY
if region=='Italy':
	condition=f"SELECT * FROM italy WHERE `depth`>=0"

# JAPAN
if region=='Japan':
	condition=(f"SELECT * FROM japan WHERE `dateandtime`>='1992-01-01 00:00:00'"
				f" AND `latitude`>0 AND `longitude`>0 AND `depth`>0"
				f" AND `magtype` LIKE 'V'")



# Power law used to fit the connectivity distribution
def power_law(x, a, b):
	return a*np.power(x, -b)

connectivityData={}
connectivityData['5km']={}
connectivityData['10km']={}

for (magMin,magMax) in ((1,10),(2,10),(1,3),(3,7),(2,4),(4,7)):
	for side in (5,10):
		
		# Magnitude windows for the condition that collects the database through mySQL
		condition+=f" AND `magnitude`>={magMin} AND `magnitude`<={magMax}"
		
		# Collect the database and create the graph
		quakes = sqlCollect(condition,side,region)
		
		# Process graph and its connectivity with edge weight
		if edgeWeight=='True':
			G = graphCreation(quakes,withEdgeWeight=True)
		
			# Degree
			degree_dict = dict(G.degree(G.nodes(),weight='weight'))
			connectivity=[]
			for item in degree_dict.values():
				connectivity.append(item)

		else:
			G = graphCreation(quakes)
			
			# Degree
			degree_dict = dict(G.degree(G.nodes()))
			connectivity=[]
			for item in degree_dict.values():
				connectivity.append(item)
		

		# Histogram of connectivity
		hist, bins = np.histogram(connectivity,bins=round(math.sqrt(len(connectivity))));
		
		# Create the data k,y from the hist and bins
		# k = bins centers . First create empty array of the length of hist
		k = np.zeros_like(hist)
		# Append to it the centers of the bins
		for i in range(1,len(bins)):
			k[i-1]=((bins[i]+bins[i-1])/2)
		
		# Check for zeros in the hist list and cut both k and y where the first zero occurs
		for i in range(len(hist)):
			if hist[i]==0:
				y=np.array(hist[:i])
				k=np.array(k[:i])
				break
			# If there is no zeros, make y=hist ( full data is taken into acoount)
			else:
				y=hist
		

		# Use try because you may not have enough data for a fit
		try:
			# Renormalize the ydata to 1 => P_k
			P_k = np.array([float(i)/sum(y) for i in y])
			# Compute the mean degree of G : k_mean
			k_mean = np.dot(k,P_k)
			# Compute the remaining degree distribution q_k
			q = np.array([(k[i+1]*P_k[i+1])/k_mean for i in range(len(k)-1)])
			# Compute the entropy
			H_q = -np.dot(q,np.log(q))


			# Compute the power_law fit to our data 
			pars, cov = curve_fit(f=power_law,xdata=k,ydata=P_k,maxfev=5000)
			# Compute the chi_squared goodness of fit = sum( ( observed - expected )^2 / expected )
			chi_squared = np.sum((P_k-power_law(k,*pars))**2/power_law(k,*pars))


			# APPEND RESULTS TO DATA
			connectivityData[f'{side}km'][f'{magMin}<magnitude<{magMax}']=[np.round(pars[1],4),np.round(chi_squared,4)]


			if side == 5:
				# Plot the results
				fig, ax = plt.subplots(1,2,figsize=(15,5))
				
				# The data, scattered
				ax[0].scatter(k,P_k)
				ax[0].set_xscale('log')
				ax[0].set_yscale('log')
				# The fit
				ax[0].plot(k,power_law(k,*pars),
						   label=f'$\gamma$ = {np.round(pars[1],4)}\n$\chi^2$ = {np.round(chi_squared,4)}',
						   color='red')

				# Legend : gamma coefficient of fit and chi_squared goodness of fit
				ax[0].legend(loc='upper right',fontsize=16,frameon=True)

				# Title of connectivity distribution ( data + fit )
				ax[0].set_title('cube size = 5 km ')
				ax[0].set_xlabel('connectivity k')
				ax[0].set_ylabel('P(k)')

			if side == 10 :
				# Connectivity distribution ( data + fit)
				# The data, scattered
				ax[1].scatter(k,P_k)
				ax[1].set_xscale('log')
				ax[1].set_yscale('log')
				# The fit
				ax[1].plot(k,power_law(k,*pars),
						   label=f'$\gamma$ = {np.round(pars[1],4)}\n$\chi^2$ = {np.round(chi_squared,4)}',
						   color='red')


				# Legend : gamma coefficient of fit and chi_squared goodness of fit
				ax[1].legend(loc='upper right',fontsize=16,frameon=True)

				# Title of connectivity distribution ( data + fit )
				ax[1].set_title('cube size = 10 km')
				ax[1].set_xlabel('connectivity k')
				ax[1].set_ylabel('P(k)')


				# Put a title and save the figure, depending on edgeWeight and the values of the different variables
				if edgeWeight=="True":
					plt.savefig(f'./connectivity{region}/connectivity{region}Weighted_{magMin}magnitude{magMax}.png');
					plt.close()

				else:
					plt.savefig(f'./connectivity{region}/connectivity{region}_{magMin}magnitude{magMax}.png');
					plt.close()

		# If there is not enough data to visualize, say so
		except:
			if side == 5:
				# Plot the results
				fig, ax = plt.subplots(1,2,figsize=(15,5))

				ax[0].text(0.5, 0.5, 'Not enough data for a regression', fontsize=16, 
					   ha='center', va='center', transform=ax[0].transAxes)

				ax[0].set_title('cube size = 5 km ')
				ax[0].set_xlabel('connectivity k')
				ax[0].set_ylabel('P(k)')
				
				
			if side == 10 :
				ax[1].text(0.5, 0.5, 'Not enough data for a regression', fontsize=16, 
					   ha='center', va='center', transform=ax[1].transAxes)

				ax[1].set_title('cube size = 10 km')
				ax[1].set_xlabel('connectivity k')
				ax[1].set_ylabel('P(k)')

			
				# Put a title and save the figure, depending on edgeWeight and the values of the different variables
				if edgeWeight=="True":
					plt.savefig(f'./connectivity{region}/connectivity{region}Weighted_{magMin}magnitude{magMax}.png');
					plt.close()

				else:
					plt.savefig(f'./connectivity{region}/connectivity{region}_{magMin}magnitude{magMax}.png');
					plt.close()


		# Extract the magnitude restrictions from the condition 
		condition = condition.replace(f" AND `magnitude`>={magMin} AND `magnitude`<={magMax}", '')


# Export the dictionary containing the results : gamma, chi squared and shannon entropy
if edgeWeight=="True":
	pickle.dump(connectivityData, open(f"./connectivity{region}/connectivity{region}DataWeighted.p", "wb"))
else:
	pickle.dump(connectivityData, open(f"./connectivity{region}/connectivity{region}Data.p", "wb"))
