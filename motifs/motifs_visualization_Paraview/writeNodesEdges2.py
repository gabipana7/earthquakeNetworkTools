import vtk

def writeObjects(nodeCoords,
                 motifCoords=[],
                 edges = [],
                 scalar = [], name = '', power = 1,
                 scalar2 = [], name2 = '', power2 = 1,
                 escalar = [], ename = '', epower = 1,
                 escalar2 = [], ename2 = '', epower2 =1,
                 nodeLabel = [],
                 method = 'vtkPolyData',
                 fileout = 'test'):
    """
    Store points and/or graphs as vtkPolyData or vtkUnstructuredGrid.
    Required argument:
    - nodeCoords is a list of node coordinates in the format [x,y,z]
    Optional arguments:
    - edges is a list of edges in the format [nodeID1,nodeID2]
    - scalar/scalar2 is the list of scalars for each node
    - name/name2 is the scalar's name
    - power/power2 = 1 for r~scalars, 0.333 for V~scalars
    - nodeLabel is a list of node labels
    - method = 'vtkPolyData' or 'vtkUnstructuredGrid'
    - fileout is the output file name (will be given .vtp or .vtu extension)
    """

    # INITIALIZE THE NODES 
    points = vtk.vtkPoints()
#   vertices = vtk.vtkCellArray()
    for node in nodeCoords:
#        point_id = points.InsertNextPoint(node)
        points.InsertNextPoint(node)

#        vertices.InsertNextCell(1)
#        vertices.InsertCellPoint(point_id)

    if motifCoords:
        # INITIALIZE TRIANGLES
        if len(motifCoords[0]) == 3:
            triangles = vtk.vtkCellArray()
    #        triangles.Allocate(len(motifCoords))

            # Initialize the points for the triangle

            for motif in motifCoords:
                triangle = vtk.vtkTriangle()
                
                # Initialize the keys
                
                triangle.GetPointIds().SetId(0,int(motif[0]))
                triangle.GetPointIds().SetId(1,int(motif[1]))
                triangle.GetPointIds().SetId(2,int(motif[2]))
            
                # Append the newly created triangle to the triangles
                triangles.InsertNextCell(triangle)  


        # INITIALIZE TETRAHEDRONS
        if len(motifCoords[0]) == 5:
            tetrahedrons = vtk.vtkCellArray()

            # Initialize the poins of tetrahedron

            for motif in motifCoords:
                tetra = vtk.vtkTetra()
                    
                tetra.GetPointIds().SetId(0,int(motif[0]))  
                tetra.GetPointIds().SetId(1,int(motif[1])) 
                tetra.GetPointIds().SetId(2,int(motif[2])) 
                tetra.GetPointIds().SetId(3,int(motif[3]))

                # Append the newly created tetra to the tetrahedrons
                tetrahedrons.InsertNextCell(tetra)


        # INITIALIZE TETRAHEDRONS BY HAND
        if len(motifCoords[0]) == 4:
            tetrahedrons = vtk.vtkCellArray()

            # Initialize the poins of tetrahedron

            for motif in motifCoords:
                triangle = vtk.vtkTriangle()
                    
                triangle.GetPointIds().SetId(0,int(motif[0]))  
                triangle.GetPointIds().SetId(1,int(motif[1])) 
                triangle.GetPointIds().SetId(2,int(motif[2])) 
                tetrahedrons.InsertNextCell(triangle)

                triangle.GetPointIds().SetId(0,int(motif[3]))  
                triangle.GetPointIds().SetId(1,int(motif[0])) 
                triangle.GetPointIds().SetId(2,int(motif[1])) 
                tetrahedrons.InsertNextCell(triangle)

                triangle.GetPointIds().SetId(0,int(motif[3]))  
                triangle.GetPointIds().SetId(1,int(motif[1])) 
                triangle.GetPointIds().SetId(2,int(motif[2])) 
                tetrahedrons.InsertNextCell(triangle)

                triangle.GetPointIds().SetId(0,int(motif[3]))  
                triangle.GetPointIds().SetId(1,int(motif[0])) 
                triangle.GetPointIds().SetId(2,int(motif[2])) 
                tetrahedrons.InsertNextCell(triangle)



    # INITIALIZE THE EDGES
    if edges:
        lines = vtk.vtkCellArray()
        lines.Allocate(len(edges))
        for edge in edges:
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0,edge[0])  
            line.GetPointIds().SetId(1,edge[1])  # line from point edge[0] to point edge[1]
            lines.InsertNextCell(line)
  
    if scalar:
        attribute = vtk.vtkFloatArray()
        attribute.SetNumberOfComponents(1)
        attribute.SetName(name)
        attribute.SetNumberOfTuples(len(scalar))
        for i, j in enumerate(scalar):   # i becomes 0,1,2,..., and j runs through scalars
            attribute.SetValue(i,j**power)

    if scalar2:
        attribute2 = vtk.vtkFloatArray()
        attribute2.SetNumberOfComponents(1)
        attribute2.SetName(name2)
        attribute2.SetNumberOfTuples(len(scalar2))
        for i, j in enumerate(scalar2):   # i becomes 0,1,2,..., and j runs through scalar2
            attribute2.SetValue(i,j**power2)
            
    if escalar:
        eattribute = vtk.vtkFloatArray()
        eattribute.SetNumberOfComponents(1)
        eattribute.SetName(ename)
        eattribute.SetNumberOfTuples(len(edges))
        for i, j in enumerate(escalar):   # i becomes 0,1,2,..., and j runs through scalar2
            eattribute.SetValue(i,j**epower)

    if escalar2:
        eattribute2 = vtk.vtkFloatArray()
        eattribute2.SetNumberOfComponents(1)
        eattribute2.SetName(ename2)
        eattribute2.SetNumberOfTuples(len(edges))
        for i, j in enumerate(escalar2):  # i becomes 0,1,2,..., and j runs through escalars
            eattribute2.SetValue(i, j**epower2)


    if nodeLabel:
        label = vtk.vtkStringArray()
        label.SetName('tag')
        label.SetNumberOfValues(len(nodeLabel))
        for i, j in enumerate(nodeLabel):   # i becomes 0,1,2,..., and j runs through scalar
            label.SetValue(i,j)

    if method == 'vtkPolyData':
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
#        polydata.SetVerts(vertices)

        if motifCoords:
            if len(motifCoords[0]) == 3:
                polydata.SetPolys(triangles)
            if len(motifCoords[0]) == 4:
                polydata.SetPolys(tetrahedrons)
            
        if edges:
            polydata.SetLines(lines)
        if scalar:
            polydata.GetPointData().AddArray(attribute)
        if scalar2:
            polydata.GetPointData().AddArray(attribute2)
        if escalar:
            polydata.GetCellData().AddArray(eattribute)
        if escalar2:
            polydata.GetCellData().AddArray(eattribute2)
        if nodeLabel:
            polydata.GetPointData().AddArray(label)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fileout+'.vtp')
        writer.SetInputData(polydata)
        writer.Write()


    elif method == 'vtkUnstructuredGrid':
        # caution: ParaView's Tube filter does not work on vtkUnstructuredGrid
        grid = vtk.vtkUnstructuredGrid()
        grid.SetPoints(points)

        if motifCoords:
#            if len(motifCoords[0]) == 3:
#                grid.SetCells(triangles)
            if len(motifCoords[0]) == 4:
                grid.SetCells(vtk.VTK_TETRA,tetrahedrons)


        if edges:
            grid.SetCells(vtk.VTK_LINE, line)
        if scalar:
            grid.GetPointData().AddArray(attribute)
        if scalar2:
            grid.GetPointData().AddArray(attribute2)
        if escalar:
            grid.GetCellData().AddArray(eattribute)
        if escalar2:
            grid.GetCellData().AddArray(eattribute2)
        if nodeLabel:
            grid.GetPointData().AddArray(label)
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(fileout+'.vtu')
        writer.SetInputData(grid)
        writer.Write()