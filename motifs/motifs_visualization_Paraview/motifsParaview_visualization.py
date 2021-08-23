import vtk
import networkx as nx

import ast
import numpy as np

from writeNodesEdges import writeObjects
from writeMotifs import writeObjectsMotifs


# ---------------------------SETUP AND COLLECTION OF QUAKES------------------------------------#
# Which region do you want to analyze ?
region = input('Input region : Vrancea / Romania / California / Italy / Japan : ')

# The side of the cubes that you split the region in 
side = int(input('Input side of the cube split in km 5 / 10 / 20 - recommended = 5 : '))

# Select desired magnitude threshold
mag = int(input('Select desired magnitude threshold - recommended = 4 : '))

G = nx.read_gexf(f'quakes{region}_{side}km_{mag}mag.gexf')



#---------------------------------------TRIANGLES-----------------------------------------#

# Use readlines() to open the 
fileTriangle = open(f'quakes{region}_{side}km_{mag}mag_triangles.txt', 'r')
linesTriangle = fileTriangle.readlines()

# Properly evaluate the Lines to get the Lists
triangleNodes=[]
for item in linesTriangle:
    triangleNodes.append(ast.literal_eval(item))
    
# Graph containing triangles only
H = nx.Graph()
for item in triangleNodes:
    H.add_edge(int(item[0]),int(item[1]))
    H.add_edge(int(item[1]),int(item[2]))
    H.add_edge(int(item[0]),int(item[2]))
    
# Set the triangle attribute = 0 to each edge
nx.set_edge_attributes(G, 0, name='triangle')

# Iterate through our triangle only network edges
for (u,v) in H.edges():
    # Assign to original network edges the attribute triangle = 1
    G[u][v]['triangle'] = 1


# Get attributes that go into VTK function

# X Y Z coords
lat =[]
long = []
depth =[]
for n in G.nodes():
    lat.append(int(G.nodes[n]['quake_xLatitude']))
    long.append(int(G.nodes[n]['quake_yLongitude']))
    depth.append(int(G.nodes[n]['quake_zDepth']))

minLat = min(lat)
maxLat = max(lat)
minLong = min(long)
maxLong = max(long)
minDepth = min(depth)
maxDepth = max(depth)
maxDimension = max(maxLat,maxLong,maxDepth)


coords=[]
for n in G.nodes():
    coords.append([ np.float32(round((int(G.nodes[n]['quake_xLatitude']) - minLat)*(maxLat/maxDimension)/(maxLat-minLat),3)),
                    np.float32(round((int(G.nodes[n]['quake_yLongitude']) - minLong)*(maxLong/maxDimension)/(maxLong-minLong),3)),
                    np.float32(round((int(G.nodes[n]['quake_zDepth']) - minDepth)*(maxDepth/maxDimension)/(maxDepth-minDepth),3))])
  
# Degree of nodes edges    
degree = [d for n, d in G.degree()]

# Weight of edges
weights = []
for (i,j) in G.edges():
    weights.append(G[i][j]['weight'])
    
# Triangle quality of edges
triangles = []
for (i,j) in G.edges():
    triangles.append(G[i][j]['triangle'])


# Write the VTK file that goes in Paraview
# The network
writeObjects(nodeCoords=coords,
             edges=G.edges(),
             scalar=degree, name='degree',
             scalar2=weights, name2='weight',
             escalar2=triangles, ename2='triangle',
#             nodeLabel=nodeLabel,
             fileout=f'network{region}Motifs_{side}km_{mag}mag_triangles')

# Only the motifs
writeObjectsMotifs(nodeCoords=coords,
             motifCoords=triangleNodes,
             fileout=f'network{region}Motifs_{side}km_{mag}mag_trianglesOnly')




#---------------------------------------SQUARES-----------------------------------------#

# Use readlines() to open the motifs txt document
fileSquare = open(f'quakes{region}_{side}km_{mag}mag_squares.txt', 'r')
linesSquare = fileSquare.readlines()

# Properly evaluate the Lines to get the Lists
squareNodes=[]
for item in linesSquare:
    squareNodes.append(ast.literal_eval(item))

# Graph containing squares only
J = nx.Graph()
for item in squareNodes:
    J.add_edge(int(item[0]),int(item[1]))
    J.add_edge(int(item[1]),int(item[2]))
    J.add_edge(int(item[2]),int(item[3]))
    J.add_edge(int(item[0]),int(item[3]))

# Set the square attribute = 0 to each edge
nx.set_edge_attributes(G, 0, name='square')


# Iterate through our squares only network edges
for (u,v) in J.edges():
    # Assign to original network edges the attribute square = 1
    G[u][v]['square'] = 1


# Get attributes that go into VTK function

# X Y Z coords
lat =[]
long = []
depth =[]
for n in G.nodes():
    lat.append(int(G.nodes[n]['quake_xLatitude']))
    long.append(int(G.nodes[n]['quake_yLongitude']))
    depth.append(int(G.nodes[n]['quake_zDepth']))

minLat = min(lat)
maxLat = max(lat)
minLong = min(long)
maxLong = max(long)
minDepth = min(depth)
maxDepth = max(depth)
maxDimension = max(maxLat,maxLong,maxDepth)

coords=[]
for n in G.nodes():
    coords.append([ np.float32(round((int(G.nodes[n]['quake_xLatitude']) - minLat)*(maxLat/maxDimension)/(maxLat-minLat),3)),
                    np.float32(round((int(G.nodes[n]['quake_yLongitude']) - minLong)*(maxLong/maxDimension)/(maxLong-minLong),3)),
                    np.float32(round((int(G.nodes[n]['quake_zDepth']) - minDepth)*(maxDepth/maxDimension)/(maxDepth-minDepth),3))])

# Degree of nodes edges    
degree = [d for n, d in G.degree()]

# Weight of edges
weights = []
for (i,j) in G.edges():
    weights.append(G[i][j]['weight'])
    
# Square quality of edges
squares = []
for (i,j) in G.edges():
    squares.append(G[i][j]['square'])


# Write the VTK file that goes in Paraview
# The network
writeObjects(coords, edges=G.edges(),
             scalar=degree, name='degree',
             scalar2=weights, name2='weight',
             escalar2=squares, ename2='square',
             fileout=f'network{region}Motifs_{side}km_{mag}mag_squares')

# Only the motifs
writeObjectsMotifs(nodeCoords=coords,
             motifCoords=squareNodes,
             fileout=f'network{region}Motifs_{side}km_{mag}mag_squaresOnly')