'''
Created on Jan 19, 2019

@author: pierre
'''
import os
import pandas as pd
import vtk
import vtkConvert 

# from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy
# https://www.programcreek.com/python/example/108192/vtk.util.numpy_support.vtk_to_numpy

def errorNonManifoldSurface():
    raise ValueError("Surface is not manifold")

def IsManifold(polyData ):
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.SetInputData(polyData)
    featureEdges.BoundaryEdgesOn() # <--- 
    featureEdges.FeatureEdgesOff()
    featureEdges.ManifoldEdgesOff()
    featureEdges.NonManifoldEdgesOn() # <--- 
    featureEdges.Update()
    poly = featureEdges.GetOutput()
    N = poly.GetNumberOfCells()
    if N == 0 : 
        return True
    else :
        print(poly)
    return False 

def writePolyData(polyData, relativePath ):
    writer = vtk.vtkPolyDataWriter()
    #writer.SetFileTypeToBinary()
    writer.SetInputData(polyData)
    writer.SetFileName(relativePath)
    writer.Write()
    print('wrote vtk data file : %s'%relativePath)


def concatenatePolyData(polyDataList, removeDupliatePoints):
    print("concatenate the inputs...")
    # Append the meshes
    #appendFilter = vtk.vtkAppendPolyData()
    appendFilter = vtk.vtkAppendFilter()
    for polyData in polyDataList:
        appendFilter.AddInputData( polyData ) 
    #appendFilter.MergePointsOn()
    appendFilter.Update()
    ug = appendFilter.GetOutput()
    polyMerge = vtkConvert.uGrid2PpolyData(ug)
    # Remove any duplicate points.
    if removeDupliatePoints:
        print("removing the dupliated points... (read vtkCleanPolyData for details...)")
        cleanFilter = vtk.vtkCleanPolyData()
        cleanFilter.SetInputData(polyMerge)
        cleanFilter.Update()
        return cleanFilter.GetOutput()
    return polyMerge

def appendNormal(polyData, verbose, autoOrient):
    if autoOrient:
        cleanFilter = vtk.vtkCleanPolyData()
        cleanFilter.SetInputData(polyData) 
        #cleanFilter.SetTolerance(1e-7)
        #cleanFilter.PointMergingOn()
        cleanFilter.Update()
        poly2 = cleanFilter.GetOutput()
        if not IsManifold(poly2) :
            raise RuntimeError('Error in concatenatePolyData') from errorNonManifoldSurface(); 
        polyData = poly2 
    if verbose:
        print("create normals ...")
    # create normals
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(polyData)
    normals.ComputeCellNormalsOn()
    normals.ComputePointNormalsOn()
    normals.SetSplitting(180)
    if autoOrient:
        normals.AutoOrientNormalsOn()
        normals.SetAutoOrientNormals(True)
        normals.ConsistencyOn()
        print("AutoOrientNormalsOn and SetAutoOrientNormals(True)")
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
    normalExts = polyData.GetCellData().GetArray("Normals")
    normalExts = vtk_to_numpy(normalExts)
    if not reverse_normal:
        normalExts *= -1.
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
        'nx': normalExts[:, 0],
        'ny': normalExts[:, 1],
        'nz': normalExts[:, 2],
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
    vtkStringArray = vtk.vtkStringArray()
    vtkStringArray.SetNumberOfValues(df.shape[0])
    vtkStringArray.SetName("labels")
    for index, row in df.iterrows():
        x = row['x [m]']
        y = row['y [m]']
        z = row['z [m]']
        tag = row['label']
        vtkStringArray.SetValue(index, tag)
        point = [x , y , z ]
        points.append(point)
        closestPointDist =  vtk.reference(0.0) 
        testPoint = [ x ,y ,z ] 
        id0 = kDTree.FindClosestPoint(testPoint, closestPointDist)
        ids.append(id0)
        print("%5s closest point for point : (%8.8f , %8.8f, %8.8f) is point id : %d Closest distance is : %s "%(tag, x,y,z, id0,  closestPointDist ) )
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
    filename2 = "interpolation-transposed.csv"
    df_appened.T.to_csv(filename2 ) 
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
    polyDataProbeLocation.GetPointData().AddArray(vtkStringArray)
    return polyDataProbeLocation
        
    
