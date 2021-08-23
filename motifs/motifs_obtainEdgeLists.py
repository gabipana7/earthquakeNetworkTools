import pandas as pd
import networkx as nx


from sqlCollectDatabaseWithCubes import sqlCollect


# THIS CODE OBTAINS THE EDGELISTS AS TXT USED FOR MOTIF DETECTION IN THE WEB APP NEMOMAP
# https://bioresearch.css.uwb.edu/biores/nemo/
# AND THE GEXF FOR MAGNITUDES > 4 FOR VISUALIZATION IN PARAVIEW

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


# Select the cube side length in km: ( 5 / 10 km)
for side in (5,10):
	# Select desired magnitude filter 
	for mag in (1,2,3,4):
		# Magnitude windows for the condition that collects the database through mySQL
		condition+=f" AND `magnitude`>={mag}"
		# Collect the database and create the graph
		quakes = sqlCollect(condition,side,region)

		# Create graph
		G = nx.Graph()
		for i in range(len(quakes['cubeIndex'])):
			# Create a node for this earthquake
			G.add_node(quakes['cubeIndex'][i])
			# check the following earthquakes chronologically
			for j in range(i+1, len(quakes['cubeIndex'])):
				# I want only the first occurence ( so use break )
				# create a node for the target quake index
				G.add_node(quakes['cubeIndex'][j])
				# Add edges
				# ADD EDGE WEIGHT ONLY FOR MAG == 4 ( visualization purposes )
				if mag == 4:
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
				

		# If magnitude filter = 4 , export gexf for visualization in paraview ( after converting with vtk )
		if mag == 4:
			# Setup network attributes as dictionaries ( dimensions for visualization purposes )
			# Only dictionaries are supported by networkx
			quake_zDepth = {}
			quake_yLongitude = {}
			quake_xLatitude = {}

			# Iterate through the rows of our database
			for index, row in quakes.iterrows():
				# Create dictionaries for each dimension ( cubeIndex : dimension) , as cubeIndex was used to create G
				quake_zDepth[row['cubeIndex']] = row['zDepth']
				quake_yLongitude[row['cubeIndex']] = row['yLongitude']
				quake_xLatitude[row['cubeIndex']] = row['xLatitude']

			# Append attributes to graph
			nx.set_node_attributes(G, quake_zDepth, 'quake_zDepth')
			nx.set_node_attributes(G, quake_yLongitude, 'quake_yLongitude')
			nx.set_node_attributes(G, quake_xLatitude, 'quake_xLatitude')

			# Write gexf 
			# Relabel the nodes
			nodeList=[]
			for n in G.nodes():
				nodeList.append(n)
			# Create the mapping : dict with {G node value} : 
			#new value ( from 1 to n = len of G nodes )    
			mapping = {}
			for i in range(len(nodeList)):
				# i+1 to create from 1 to n ( not from 0 ) idk why it does not want from 0 
				mapping[nodeList[i]] = i+1
			# Create new graph with relabeled nodes
			G = nx.relabel_nodes(G, mapping)
			
			nx.write_gexf(G,f'./motifs{region}Visualization/quakes{region}_{side}km_{mag}mag.gexf')

		if mag != 4:
			# Relabel the nodes
			nodeList=[]
			for n in G.nodes():
				nodeList.append(n)
			# Create the mapping : dict with {G node value} : 
			#new value ( from 1 to n = len of G nodes )    
			mapping = {}
			for i in range(len(nodeList)):
				# i+1 to create from 1 to n ( not from 0 ) idk why it does not want from 0 
				mapping[nodeList[i]] = i+1
			# Create new graph with relabeled nodes
			G = nx.relabel_nodes(G, mapping)

		nx.write_edgelist(G, f'quakes{region}_{side}km_{mag}mag.txt', data=False)


		# Extract the magnitude restrictions from the condition
		condition = condition.replace(f" AND `magnitude`>={mag}", '')

