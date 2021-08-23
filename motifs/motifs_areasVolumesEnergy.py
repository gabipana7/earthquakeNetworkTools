# MOTIFS SURFACING
import ast
from geopy.distance import distance
from geopy.point import Point
import math
import numpy as np
import matplotlib


#-------------------------MOTIFS FILE PROCESSING--------------------------------#

#-------------------------TOTAL / MEAN  ENERGY--------------------------------#


# Calculates the total and mean energy per motif
def totalMeanEnergyMotif(region,side,mag,motif,originalG,quakesDataFrame):
    # Open the motifs file
    fileMotif = open(f'quakes{region}_{side}km_{mag}mag_{motif}.txt', 'r')
    linesMotif = fileMotif.readlines()

    # Properly evaluate the Lines to get the Lists
    motifNodes=[]
    for item in linesMotif:
        motifNodes.append(ast.literal_eval(item))


    # Dictionary for setting motif : 
    # {Motif explicitly by nodes : [total energy in motif, mean energy in Motif]}
    energyMotif = {}

    # Iterate through the Motifs in our list that contains the Motifs as a list of it's nodes
    for motifs in motifNodes:
        
        # Initialize the lists ( 3 elements ) containing the total and mean energy in Motif ( per nodes )
        energyInMotif = []
        meanEnergyInMotif = []
        
        # Iterate through the nodes of the motif
        for i in range(len(motifs)):
            
            # List the quakes in the node
            quakesInNode = originalG.nodes[int(motifs[i])]['quake_index']
            # Set the energy of that node to 0
            energyInNode = 0
            
            # Iterathe through the quakes of that node
            for quake in quakesInNode:
                # Raise the total energy accumulated
                energyInNode += quakesDataFrame['energyRelease'][int(quake)]
            # Append the total and mean energy of that node to the lists for our Motif
            energyInMotif.append(energyInNode)
            meanEnergyInMotif.append(energyInNode/len(quakesInNode))
        
        # After computing for all the nodes of a Motif append to the dictionary the values of total and mean energy
        if motif=='triangles':
            energyMotif[str(motifs)]=[sum(energyInMotif),sum(meanEnergyInMotif)/3]
        if motif=='squares':
            energyMotif[str(motifs)]=[sum(energyInMotif),sum(meanEnergyInMotif)/4]

    return(motifNodes,energyMotif)



#-------------------------TRIANGLES--------------------------------#

# Calculate the areas per motif
# Uses the magnitudesMotif to calculate area/totalMag and area/meanMag per motif
def areasInTriangles(motifNodes,energyMotif,originalG,quakesDataFrame):
    
    # Initialize lists for areas
    areas=[]
    areasWeightTotalMag=[]
    areasWeightMeanMag=[]

    # Initialize the dimensions for the points of the trianlge : x=lat,y=long,z=depth 
    x=[0,0,0]
    y=[0,0,0]
    z=[0,0,0]

    # Iterate through the motifs 
    for triangle in motifNodes:
        for i in range(3):
            # For each point of the triangle you have 3 dimensions - lat, long, depth
            # Select from the quakesDataFrame the dimension column based on the index of any quake in our motif node
            #(all the quakes in that motif will have the same cube, i.e. the same dimensions, so here [0] means first quake)
            x[i] = quakesDataFrame['cubeLatitude'][int(originalG.nodes[int(triangle[i])]['quake_index'][0])]
            y[i] = quakesDataFrame['cubeLongitude'][int(originalG.nodes[int(triangle[i])]['quake_index'][0])]
            z[i] = quakesDataFrame['cubeDepth'][int(originalG.nodes[int(triangle[i])]['quake_index'][0])]
        
        # The 3 points of the triangle in GEOPY Point
        X=Point(x[0],y[0],z[0])
        Y=Point(x[1],y[1],z[1])
        Z=Point(x[2],y[2],z[2])
        
        # Calculate the sides of the trianlge using GEOPY distance in km
        # First calculate the flat distance using GEOPY ( IT DOES NOT SUPPORT ALTITUDE)
        flat_distance_a = distance(X[:2], Y[:2]).km 
        flat_distance_b = distance(Y[:2], Z[:2]).km
        flat_distance_c = distance(X[:2], Z[:2]).km
        # Introduce the altitude by doing an euclidian distance
        a = math.sqrt(flat_distance_a**2 + (X[2] - Y[2])**2)
        b = math.sqrt(flat_distance_b**2 + (Y[2] - Z[2])**2)
        c = math.sqrt(flat_distance_c**2 + (X[2] - Z[2])**2)

        # calculate semiperimeter
        sp = (a+b+c)/2
        
        # Use Heron's formula Area = sqrt(semiperimeter(sp-a)(sp-b)(sp-c))
        # We have the areas of the triangle, append them to lists ! 

        A = math.sqrt(abs(sp*round(sp-a,4)*round(sp-b,4)*round(sp-c,4)))
        if A < 5:
            continue
        else:

            AWeightTotalMag = A / energyMotif[str(triangle)][0]
            areasWeightTotalMag.append(AWeightTotalMag)
            
            AWeightMeanMag = A / energyMotif[str(triangle)][1]
            areasWeightMeanMag.append(AWeightMeanMag)
            
            
            

    return(areasWeightTotalMag,areasWeightMeanMag)



#-------------------------SQUARES--------------------------------#

# Calculate the volumes per motif
# Uses the magnitudesMotif to calculate volume/totalMag and volume/meanMag per motif
def volumesInSquares(motifNodes,energyMotif,originalG,quakesDataFrame):
    # Initialize lists for volumes
    volumes=[]
    volumesWeightTotalMag=[]
    volumesWeightMeanMag=[]

    # Initialize the dimensions for the points of the trianlge : x=lat,y=long,z=depth 
    x=[0,0,0,0]
    y=[0,0,0,0]
    z=[0,0,0,0]

    # Iterate through the motifs 
    for square in motifNodes:
        for i in range(4):
            # For each point of the square you have 3 dimensions - lat, long, depth
            x[i] = quakesDataFrame['cubeLatitude'][int(originalG.nodes[int(square[i])]['quake_index'][0])]
            y[i] = quakesDataFrame['cubeLongitude'][int(originalG.nodes[int(square[i])]['quake_index'][0])]
            z[i] = quakesDataFrame['cubeDepth'][int(originalG.nodes[int(square[i])]['quake_index'][0])]
        
        # The 4 points of the square in GEOPY Point
        X=Point(x[0],y[0],z[0])
        Y=Point(x[1],y[1],z[1])
        Z=Point(x[2],y[2],z[2])
        K=Point(x[3],y[3],z[3])
        
        # Calculate the sides of the tetrahedron using GEOPY distance in km
        # First calculate the flat distance using GEOPY ( IT DOES NOT SUPPORT ALTITUDE)
        flat_distance_W = distance(X[:2], Y[:2]).km 
        flat_distance_V = distance(Y[:2], Z[:2]).km
        flat_distance_U = distance(Z[:2], X[:2]).km
        flat_distance_u = distance(Y[:2], K[:2]).km
        flat_distance_v = distance(X[:2], K[:2]).km
        flat_distance_w = distance(Z[:2], K[:2]).km
        # Introduce the altitude by doing an euclidian distance
        W = math.sqrt(flat_distance_W**2 + (X[2] - Y[2])**2)
        V = math.sqrt(flat_distance_V**2 + (Y[2] - Z[2])**2)
        U = math.sqrt(flat_distance_U**2 + (Z[2] - X[2])**2)
        u = math.sqrt(flat_distance_u**2 + (Y[2] - K[2])**2)
        v = math.sqrt(flat_distance_v**2 + (X[2] - K[2])**2)
        w = math.sqrt(flat_distance_w**2 + (Z[2] - K[2])**2)

        # calculate elements that go into elements that go into Heron formula
        A = (w-U+v)*(U+v+w)
        B = (u-V+w)*(V+w+u)
        C = (v-W+u)*(W+u+v)
        
        a = (U-v+w)*(v-w+U)
        b = (V-w+u)*(w-u+V)
        c = (W-u+v)*(u-v+W)
        
        # elements that go into Heron formula
        p = math.sqrt(abs(a*B*C))
        q = math.sqrt(abs(b*C*A))
        r = math.sqrt(abs(c*A*B))
        s = math.sqrt(abs(a*b*c))
        

        # Use Heron's formula for tetrahedrons
        # We have the volume of tetrahedrons, append it to lists !

        Volume = math.sqrt(abs((-p+q+r+s)*(p-q+r+s)*(p+q-r+s)*(p+q+r-s)))/(192*u*v*w)
        if Volume < 5:
            continue
        else:

            VolumeWeightTotalMag = Volume / energyMotif[str(square)][0]
            volumesWeightTotalMag.append(VolumeWeightTotalMag)
            
            VolumeWeightMeanMag = Volume / energyMotif[str(square)][1]
            volumesWeightMeanMag.append(VolumeWeightMeanMag)

         
    return(volumesWeightTotalMag,volumesWeightMeanMag)
        

#-------------------------POWER LAW FOR REGRESSION--------------------------------#

# Define the function to do the regression with 
def power_law(x, a, b):
    return a*np.power(x, -b)
