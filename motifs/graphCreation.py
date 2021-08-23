import datetime
import pandas as pd
import networkx as nx


# Takes the quakes DataFrame and returns the quakes network in networkx graph (G) with attributes:
#(cube dimensions and index of quakes with same cubeIndex)
def graphCreation(quakes,withEdgeWeight=False):

	# Create graph -  QUAKES ONE BY ONE
	G = nx.Graph()
	# Iterate through the earthquakes checking for magnitude condition
	for i in range(len(quakes['cubeIndex'])):
		# Create a node for this earthquake
		G.add_node(quakes['cubeIndex'][i])
		# check the following earthquakes chronologically
		for j in range(i+1, len(quakes['cubeIndex'])):
			# I want only the first occurence (so break at the end)
			# create a node for the target quake index
			G.add_node(quakes['cubeIndex'][j])
			
			# Check if you need the edge weights
			if withEdgeWeight==True:
				# WITH EDGE WEIGHT
				# check to see if there is an edge between them
				if G.has_edge(quakes['cubeIndex'][i], quakes['cubeIndex'][j]):
					# we added this one before, just increase the weight by one
					G[quakes['cubeIndex'][i]][quakes['cubeIndex'][j]]['weight'] += 1
					break

				else:
					# new edge. add with weight=1
					G.add_edge(quakes['cubeIndex'][i], quakes['cubeIndex'][j], weight=1)
					break
			else:
				# WITHOUT EDGE WEIGHT
				G.add_edge(quakes['cubeIndex'][i], quakes['cubeIndex'][j])
				break

	# Setup network attributes as dictionaries
	# Only dictionaries are supported by networkx
	quake_index = {}
	quake_zDepth = {}
	quake_yLongitude = {}
	quake_xLatitude = {}

	# Turn the dataframe into a dictionary, the indexing quantity
	#is the index of the quake ( the row value )
	for index, row in quakes.iterrows():
		
		# Create a dictionary with cubeIndex : [ list of quake indexes ] (the quakes in that certain cube)
		# if the cube has already been indexed
		if row['cubeIndex'] in quake_index.keys():
			# just add that event in the dictionary
			quake_index[row['cubeIndex']].append(str(index))
		# if not yet indexed
		else:
			# create that dictionary entry
			quake_index[row['cubeIndex']]=[]
			# and add that event in the dictionary
			quake_index[row['cubeIndex']].append(str(index))
			
		# Create dictionaries for each dimension ( cubeIndex : dimension)
		quake_zDepth[row['cubeIndex']] = str(row['zDepth'])
		quake_yLongitude[row['cubeIndex']] = str(row['yLongitude'])
		quake_xLatitude[row['cubeIndex']] = str(row['xLatitude'])


	nx.set_node_attributes(G, quake_zDepth, 'quake_zDepth')
	nx.set_node_attributes(G, quake_yLongitude, 'quake_yLongitude')
	nx.set_node_attributes(G, quake_xLatitude, 'quake_xLatitude')
	nx.set_node_attributes(G, quake_index, 'quake_index')


	# Relabel the nodes
	nodeList=[]
	for n in G.nodes():
		nodeList.append(n)

	# Create the mapping : dict with {G node value} : 
	#new value ( from 1 to n = len of G nodes )    
	mapping = {}
	for i in range(len(nodeList)):
		# i+1 to create from 1 to n ( not from 0 ) 
		mapping[nodeList[i]] = i+1

	# Create new graph with relabeled nodes
	G = nx.relabel_nodes(G, mapping)

	return(G)
