'''
Created on Jan 19, 2019

@author: pierre
'''
import os
import pandas as pd
import vtk

# from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy
# https://www.programcreek.com/python/example/108192/vtk.util.numpy_support.vtk_to_numpy


def concatenatePolyData(polyDataList, removeDupliatePoints):
    print("concatenate the inputs...")
    # Append the meshes
    appendFilter = vtk.vtkAppendPolyData()
    for polyData in polyDataList:
        appendFilter.AddInputData(polyData)
    appendFilter.Update()
    # Remove any duplicate points.
    if removeDupliatePoints:
        print("removing the dupliated points...")
        cleanFilter = vtk.vtkCleanPolyData()
        cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
        cleanFilter.Update()
        return cleanFilter.GetOutput()
    else:
        return appendFilter.GetOutput()

def appendNormal(polyData, verbose, autoOrient):
    if verbose:
        print("create normals ...")
    # create normals
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(polyData)
    normals.ComputeCellNormalsOn()
    normals.SetSplitting(0)
    if not autoOrient:
        normals.AutoOrientNormalsOn()
        normals.SetAutoOrientNormals(True)
        print("AutoOrientNormalsOn and SetAutoOrientNormals(True)")
    normals.ConsistencyOn()
    normals.Update()
    polyData = normals.GetOutput()
    return polyData

def getDataFrameAeroData(polyData, config):
    Narray = polyData.GetCellData().GetNumberOfArrays()
    available_vars = []
    for i in range(0, Narray):
        nameVar = polyData.GetCellData().GetArrayName(i)
        available_vars.append(nameVar)
    checkOK = True
    for key in config.varDict.keys():
        if not config.varDict[key] in available_vars:
            print("could not find:  %s" % config.varDict[key])
            checkOK = False
    if not checkOK:
        print("Could not find at least on of the following:")
        for key in config.varDict.keys():
            print("%30s : %-30s" % (key, config.varDict[key]))
        print("Here is what is available:")
        for i in range(0, Narray):
            nameVar = polyData.GetCellData().GetArrayName(i)
            print("%3d : %s" % (i, polyData.GetCellData().GetArrayName(i)))
        exit(1)
    else:
        print("Use the following data:")
        for key in config.varDict.keys():
            print("%30s : %-30s" % (key, config.varDict[key]))
    dict_new = {}
    if config.IsPressureIntegrationPossible():
        cps = polyData.GetCellData().GetArray(config.varDict["Cp"])
        cps = vtk_to_numpy(cps)
        dict_new['Cp'] = cps
    if config.IsViscousIntegrationPossible():
        cfxs = polyData.GetCellData().GetArray(config.varDict["Cfx"])
        cfys = polyData.GetCellData().GetArray(config.varDict["Cfy"])
        cfzs = polyData.GetCellData().GetArray(config.varDict["Cfz"])
        cfxs = vtk_to_numpy(cfxs)
        cfys = vtk_to_numpy(cfys)
        cfzs = vtk_to_numpy(cfzs)
        dict_new['Cfx'] = cfxs
        dict_new['Cfy'] = cfys
        dict_new['Cfz'] = cfzs
    if config.IsHeatFluxIntegrationPossible():
        Qs = polyData.GetCellData().GetArray(config.varDict["q"])
        Qs = vtk_to_numpy(Qs)
        dict_new['q'] = Qs
    if config.IsChemHeatFluxIntegrationPossible():
        Qchems = polyData.GetCellData().GetArray(config.varDict["qd"])
        Qchems = vtk_to_numpy(Qchems)
        dict_new['qd'] = Qchems
    df = pd.DataFrame(dict_new)
    return df


'''
recover  Normals and center from cells.
'''


def getDataFrameGeo(polyData, verbose, reverse_normal):
    #
    if verbose:
        print("recover triangle centers...")
    normals = polyData.GetCellData().GetArray("Normals")
    normals = vtk_to_numpy(normals)
    if not reverse_normal:
        normals *= -1.
    #
    cells = polyData.GetPolys()
    nbOfCells = cells.GetNumberOfCells()
    centers_x = []
    centers_y = []
    centers_z = []
    areas = []
    for i in range(0, nbOfCells):  # TODO this is slow. Remove that loop and
        tria = polyData.GetCell(i)
        area = tria.ComputeArea()
        center = [0., 0., 0.]
        # TODO : test if tria is a triangle. If not, handle other cases
        tria.TriangleCenter(tria.GetPoints().GetPoint(
            0), tria.GetPoints().GetPoint(1), tria.GetPoints().GetPoint(2), center)
        areas.append(area)
        centers_x.append(center[0])
        centers_y.append(center[1])
        centers_z.append(center[2])
    dict_new = {
        'area': areas,
        'nx': normals[:, 0],
        'ny': normals[:, 1],
        'nz': normals[:, 2],
        'cx': centers_x,
        'cy': centers_y,
        'cz': centers_z,
    }
    df = pd.DataFrame(dict_new)
    return df

def getPolyDataCroppedFromData(polyData, list_variableToKeepForWriting):
    Narray = polyData.GetCellData().GetNumberOfArrays()
    nameVars = []
    for i in range(0, Narray):
        nameVar = polyData.GetCellData().GetArrayName(i)
        nameVars.append(nameVar)
    nameVars.sort()
    # polyData.GetCellData().Dele
    for nameVar in nameVars:
        if nameVar in list_variableToKeepForWriting:
            # keep it
            print("keep : %s" % nameVar)
            pass
        else:
            polyData.GetCellData().RemoveArray(nameVar)
            # remove it
            pass
    return polyData

def appendFrictionVector(polyData, config):
    cfxs = polyData.GetCellData().GetArray(config.varDict["Cfx"])
    cfys = polyData.GetCellData().GetArray(config.varDict["Cfy"])
    cfzs = polyData.GetCellData().GetArray(config.varDict["Cfz"])
    cfxs = vtk_to_numpy(cfxs)
    cfys = vtk_to_numpy(cfys)
    cfzs = vtk_to_numpy(cfzs)
    cells = polyData.GetPolys()
    nbOfCells = cells.GetNumberOfCells( ) 
    pcoords = vtk.vtkDoubleArray()
    pcoords.SetNumberOfComponents(3)
    pcoords.SetName("vectCf")
    for i in range (0 , nbOfCells) : # so slow !!! , numpize this ! 
        pcoords.InsertNextTuple3( cfxs[i] , cfys[i] , cfzs[i]  ) ; 
    polyData.GetCellData().AddArray(pcoords) 
    return polyData
    
def ExtractDataFromTheClosestCellCenter( polyData ,  cvsFilePath) :
    # load df 
    if not os.path.isfile(cvsFilePath):
        raise Exception('File {:s} does not exist'.format(cvsFilePath))
    df = pd.read_csv(cvsFilePath )
    # prepare data to be used from polydata
    # var 
    Narray = polyData.GetCellData().GetNumberOfArrays()
    available_vars = []
    for i in range(0, Narray ):
        nameVar = polyData.GetCellData().GetArrayName(i)
        available_vars.append(nameVar)
    # CoG 
    cells = polyData.GetPolys()
    nbOfCells = cells.GetNumberOfCells()
    centers_x = []
    centers_y = []
    centers_z = []
    for i in range ( 0, nbOfCells ) :  # TODO this is slow. Remove that loop and 
        tria = polyData.GetCell(i)
        center = [0.,0.,0.]
        # TODO : test if tria is a triangle. If not, handle other cases 
        tria.TriangleCenter(tria.GetPoints().GetPoint(0), tria.GetPoints().GetPoint(1), tria.GetPoints().GetPoint(2), center)
        centers_x.append(center[0])
        centers_y.append(center[1])
        centers_z.append(center[2])
    # look up for closest point based on cell CoG
    kDTree = vtk.vtkKdTree()
    points = vtk.vtkPoints()
    for i in range ( 0, nbOfCells )  : # TODO : this is slow, use numpy array instead 
        points.InsertNextPoint(centers_x[i], centers_y[i] , centers_z[i])
    kDTree.BuildLocatorFromPoints(points)
    # TODO check header variables : existence of x y z 
    ids = []
    points = []
    for index, row in df.iterrows():
        x = row['x [m]']
        y = row['y [m]']
        z = row['z [m]']
        point = [x , y , z ]
        points.append(point)
        closestPointDist =  vtk.reference(0.0) 
        testPoint = [ x ,y ,z ] 
        id0 = kDTree.FindClosestPoint(testPoint, closestPointDist)
        ids.append(id0)
        print("closest point for point : (%8.8f , %8.8f, %8.8f) is point id : %d   Closest distance is : %s "%(x,y,z, id,  closestPointDist ) )
    # recover data based on the cell id.
    datas = [] 
    dict_new = {}
    for var in available_vars : 
        vtkdata = polyData.GetCellData().GetArray(var)
        npdata = vtk_to_numpy(vtkdata)
        croppedData = npdata[ids] 
        # print(croppedData.shape)
        # get ride of vect 
        if len(croppedData.shape ) == 1  :
            print( "append variable : %s"% var )
            datas.append(croppedData)
            dict_new[var ] = croppedData
    # write data 
    df_to_append = pd.DataFrame(dict_new)
    df_appened = pd.concat([df, df_to_append ] , axis = 1 ) 
    filename = "interpolation.csv"
    df_appened.to_csv(filename ) 
    print("wrote interpolated data info : %s"%filename)
        
    vtkPoints = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()
    for i , point in enumerate( points ) :
        pid = [0]
        pid[0] = vtkPoints.InsertNextPoint(  point)
        vertices.InsertNextCell( 1, pid)
    polyDataProbeLocation = vtk.vtkPolyData()
    polyDataProbeLocation.SetPoints(vtkPoints)
    polyDataProbeLocation.SetVerts(vertices)
    return polyDataProbeLocation
        
    
