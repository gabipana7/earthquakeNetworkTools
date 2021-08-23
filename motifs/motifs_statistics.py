import ast
import numpy as np
import datetime
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
from scipy.optimize import curve_fit
import math

# Extract motifs from the motif files depending on region, cube side, magnitude filter and the specific motif
def motifOpen(region='Vrancea', side=5, mag=0, motif='triangles'):
	fileMotif = open(f'quakes{region}_{side}km_{mag}mag_{motif}.txt', 'r')
	linesMotif = fileMotif.readlines()

	# Properly evaluate the Lines to get the Lists
	motifNodes=[]
	for item in linesMotif:
		motifNodes.append(ast.literal_eval(item))
		 
	if motif == 'triangles':

		# Create the Graph containing all the triangles
		H = nx.Graph()
		for item in motifNodes:
			H.add_edge(int(item[0]),int(item[1]))
			H.add_edge(int(item[1]),int(item[2]))
			H.add_edge(int(item[0]),int(item[2]))

	
	if motif == 'squares':
		# Create graph containing all the squares   
		H = nx.Graph()
		for item in motifNodes:
			H.add_edge(int(item[0]),int(item[1]))
			H.add_edge(int(item[1]),int(item[2]))
			H.add_edge(int(item[2]),int(item[3]))
			H.add_edge(int(item[0]),int(item[3]))

	
	if motif == 'pentagons':
		H = nx.Graph()
		for item in motifNodes:
			H.add_edge(int(item[0]),int(item[1]))
			H.add_edge(int(item[1]),int(item[2]))
			H.add_edge(int(item[2]),int(item[3]))
			H.add_edge(int(item[3]),int(item[4]))
			H.add_edge(int(item[0]),int(item[4]))

	return(H)
	

# Mean magnitudes in motifs
def meanMagnitude(graphMotifs,originalG,quakesDataFrame):
	# Unpack the attribute quake_index from original graph to triangles graph 
	quakesInMotifs=[]
	for node in graphMotifs.nodes():
		if node in originalG.nodes():
			quakesInMotifs.append(originalG.nodes[node]['quake_index'])

			
	# unpack the list of lists obtained => A list of all the quakes present in the pentagons
	quakesInMotifs = [item for sublist in quakesInMotifs for item in sublist]


	# Calculate mean magnitude in Pentagons
	magnitudes=0
	for item in quakesInMotifs:
		magnitudes += quakesDataFrame['magnitude'][int(item)]

	meanMagnitude = magnitudes/(len(quakesInMotifs))
	return(quakesInMotifs,meanMagnitude)


# Extract the top X quakes and see in how many motifs of a certain type they exist in 
def topXwhere(quakesDataFrame, quakesInMotifs,X=50,):
	# Extract the top X quakes by magnitude ( high to low )
	qV2 = quakesDataFrame.sort_values('magnitude', ascending=False)[:X]

	isInMotif=0
	for index, row in qV2.iterrows():
		if str(index) in quakesInMotifs:
			isInMotif +=1

	return(isInMotif)


# Mean time between events in cells in different motifs
def meanTime(quakesDataFrame,originalG,graphMotifs):
	# Need to unpack the quakes in each node 
	allQuakes = nx.get_node_attributes(originalG,'quake_index')
	timeDiff = {}

	# Use this to calculate 
	for node in allQuakes:
	    if len(allQuakes[node])==1:
	        if allQuakes[node][0]=='0':
	            continue
	        else:
	            timeDiff[node] = quakesDataFrame['date'][int(allQuakes[node][0])] - quakesDataFrame['date'][int(allQuakes[node][0])-1]
	    if len(allQuakes[node])==2:
	        timeDiff[node] = quakesDataFrame['date'][int(allQuakes[node][1])] - quakesDataFrame['date'][int(allQuakes[node][0])]
	    else:
	        diff_list=[]
	        for x,y in zip(allQuakes[node][0::], allQuakes[node][1::]):
	            #print(x,y)
	            diff_list.append(quakesDataFrame['date'][int(y)]-quakesDataFrame['date'][int(x)])
	            timeDiff[node]= sum(diff_list,datetime.timedelta()) / len(diff_list)

	# Set the time difference attribute for each node
	nx.set_node_attributes(originalG, timeDiff, 'timeDiff')

	# Time difference in TRIANGLES:
	timeDiffMotifs=[]
	for node in graphMotifs.nodes():
	    if node in originalG.nodes():
	        timeDiffMotifs.append(originalG.nodes[node]['timeDiff'])

	# Need to convert in seconds, use dataframe because the timedifference can  be too big to work in nanoseconds
	#do this with dataframe
	df1 = pd.DataFrame(timeDiffMotifs)
	# convert in seconds
	df1[0] = df1[0].dt.total_seconds()
	seconds1 = df1[0].sum() / len(df1[0])

	return(str(datetime.timedelta(seconds=seconds1)))



# DEGREE MEASURES FOR MOTIFS

def power_law(x, a, b):
    return a*np.power(x, b)

def connectivity(graphMotifs):

	# Motifs connectivity
	degree_dict_motifs = dict(graphMotifs.degree(graphMotifs.nodes()))
	nx.set_node_attributes(graphMotifs, degree_dict_motifs, 'degree')

	sorted_degree_motifs = sorted(degree_dict_motifs.items(), 
	                       key=itemgetter(1), reverse=True)

	motifList=[]
	for item in sorted_degree_motifs:
	    motifList.append(item[1])
	    
	hist, bins = np.histogram(motifList,bins=round(math.sqrt(len(motifList))));

	exceptionControl = True
	try:
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

		y_norm = (y - np.min(y))/np.ptp(y)

		pars, cov = curve_fit(f=power_law, xdata=x, ydata=y_norm, maxfev=5000)
		return(motifList,x,y_norm,pars,exceptionControl)
		
		
	except:
		#print('Not enough data to process a power law regression for motifs in:\n')
		exceptionControl = False
		motifList = False
		x = False
		y_norm = False
		pars = False
		return(motifList,x,y_norm,pars,exceptionControl)
		
		
