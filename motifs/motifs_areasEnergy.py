import arviz as az
import matplotlib
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from motifs_areasVolumesEnergy import*
from graphCreation import*
from sqlCollectDatabaseWithCubes import sqlCollect

# Set the style of the plots to a more beautiful format

region = input('Input region : Vrancea / Romania / California / Italy / Japan : ')


# VRANCEA
if region=='Vrancea':
	condition=(f"SELECT * FROM romplus WHERE `dateandtime`>='1976-01-01 00:00:00'"
		   f" AND `latitude`>=44.9 AND `latitude`<=46.2 AND `longitude`>=25.5"
		   f" AND `longitude`<=27.4 AND `depth`>=50 AND `depth`<=200")

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



# SQUARES VOLUMES
motif = 'triangles'


for mag in (1,2,3):
    
    # Magnitude windows for the condition that collects the database through mySQL
    condition+=f" AND `magnitude`>={mag}"
    

    for side in (5,10):
        # Collect the database and create the graph
        quakes = sqlCollect(condition,side,region,energyRelease=True)
        G = graphCreation(quakes)
        
        motifNodes,energyMotif = totalMeanEnergyMotif(region,side,mag,motif,G,quakes)
        
        areasWeightTotalMag,areasWeightMeanMag = areasInTriangles(motifNodes,energyMotif,G,quakes)
        
        # TOTAL ENERGY PLOTS
        hist, bins = np.histogram(areasWeightTotalMag,bins=round(math.sqrt(len(areasWeightTotalMag))))
        
        # Create the x as hist with zeros, force into floats ! 
        x = np.zeros_like(hist.astype(float))
        for i in range(1,len(bins)):
            x[i-1]=((bins[i]+bins[i-1])/2)

        for i in range(len(hist)):
            if hist[i]==0:
                y=np.array(hist[:i])
                x=np.array(x[:i])
                break
            else:
                y=hist

        y_norm = [float(i)/sum(y) for i in y]

        pars, cov = curve_fit(f=power_law,xdata=x,ydata=y_norm,maxfev=5000)

        # Compute the chi_squared goodness of fit = sum( ( observed - expected )^2 / expected )
        chi_squared = np.sum((y_norm-power_law(x,*pars))**2/power_law(x,*pars))
        
        if side == 5:
            # Plot the results
            fig, ax = plt.subplots(1,2,figsize=(15,5))

            # The data, scattered
            ax[0].scatter(x,y_norm)
            ax[0].set_xscale('log')
            ax[0].set_yscale('log')
            # The fit
            ax[0].plot(x,power_law(x,*pars),
                       label=f'$\gamma$ = {np.round(pars[1],4)}\n$\chi^2$ = {np.round(chi_squared,4)}',
                       color='red')


            # Legend : gamma coefficient of fit and chi_squared goodness of fit
            ax[0].legend(loc='upper right',fontsize=16,frameon=True)

            # Title of connectivity distribution ( data + fit )
            ax[0].set_title('cube size = 5 km ')
            ax[0].set_xlabel(r'$S_{TE}$')
            ax[0].set_ylabel(r'$P(S_{TE})$')

        if side == 10 :
            # Connectivity distribution ( data + fit)
            # The data, scattered
            ax[1].scatter(x,y_norm)
            ax[1].set_xscale('log')
            ax[1].set_yscale('log')
            # The fit
            ax[1].plot(x,power_law(x,*pars),
                       label=f'$\gamma$ = {np.round(pars[1],4)}\n$\chi^2$ = {np.round(chi_squared,4)}',
                       color='red')


            # Legend : gamma coefficient of fit and chi_squared goodness of fit
            ax[1].legend(loc='upper right',fontsize=16,frameon=True)

            # Title of connectivity distribution ( data + fit )
            ax[1].set_title('cube size = 10 km')
            ax[1].set_xlabel(r'$S_{TE}$')
            ax[1].set_ylabel(r'$P(S_{TE})$')
            
            #plt.suptitle(f'Surface distributions in {region} - magnitude>{mag}',fontsize=18)
            fig.savefig(f'./motifs{region}Statistics/quakes{region}_totalEnergy_{mag}mag_trianglesAreas.png')
           
        
        # MEAN ENERGY PLOTS 
        #fig.clear()
        hist, bins = np.histogram(areasWeightMeanMag,bins=round(math.sqrt(len(areasWeightMeanMag))))

        # Create the x as hist with zeros, force into floats ! 
        x = np.zeros_like(hist.astype(float))
        for i in range(1,len(bins)):
            x[i-1]=((bins[i]+bins[i-1])/2)

        for i in range(len(hist)):
            if hist[i]==0:
                y=np.array(hist[:i])
                x=np.array(x[:i])
                break
            else:
                y=hist

        y_norm = [float(i)/sum(y) for i in y]

        pars, cov = curve_fit(f=power_law,xdata=x,ydata=y_norm,maxfev=5000)

        # Compute the chi_squared goodness of fit = sum( ( observed - expected )^2 / expected )
        chi_squared = np.sum((y_norm-power_law(x,*pars))**2/power_law(x,*pars))

        if side == 5:
            # Plot the results
            fig2, ax2 = plt.subplots(1,2,figsize=(15,5))

            # The data, scattered
            ax2[0].scatter(x,y_norm)
            ax2[0].set_xscale('log')
            ax2[0].set_yscale('log')
            # The fit
            ax2[0].plot(x,power_law(x,*pars),
                       label=f'$\gamma$ = {np.round(pars[1],4)}\n$\chi^2$ = {np.round(chi_squared,4)}',
                       color='red')


            # Legend : gamma coefficient of fit and chi_squared goodness of fit
            ax2[0].legend(loc='upper right',fontsize=16,frameon=True)

            # Title of connectivity distribution ( data + fit )
            ax2[0].set_title('cube size = 5 km ')
            ax2[0].set_xlabel(r'$S_{ME}$')
            ax2[0].set_ylabel(r'$P(S_{ME})$')

        if side == 10 :
            # Connectivity distribution ( data + fit)
            # The data, scattered
            ax2[1].scatter(x,y_norm)
            ax2[1].set_xscale('log')
            ax2[1].set_yscale('log')
            # The fit
            ax2[1].plot(x,power_law(x,*pars),
                       label=f'$\gamma$ = {np.round(pars[1],4)}\n$\chi^2$ = {np.round(chi_squared,4)}',
                       color='red')


            # Legend : gamma coefficient of fit and chi_squared goodness of fit
            ax2[1].legend(loc='upper right',fontsize=16,frameon=True)

            # Title of connectivity distribution ( data + fit )
            ax2[1].set_title('cube size = 10 km')
            ax2[1].set_xlabel(r'$S_{ME}$')
            ax2[1].set_ylabel(r'$P(S_{ME})$')

            #plt.suptitle(f'Surface distributions in {region} - magnitude>{mag}',fontsize=18)
            fig2.savefig(f'./motifs{region}Statistics/quakes{region}_meanEnergy_{mag}mag_trianglesAreas.png')

    # Extract the magnitude restrictions from the condition 
    condition = condition.replace(f" AND `magnitude`>={mag}", '')