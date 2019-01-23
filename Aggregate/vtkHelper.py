'''
Created on Jan 19, 2019

@author: pierre
'''
import os 
import pandas as pd  
import vtk 
#from Aggregate.geometry.Point import Point

from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
#https://www.programcreek.com/python/example/108192/vtk.util.numpy_support.vtk_to_numpy

def concatenatePolyData( polyDataList, removeDupliatePoints):
    print("concatenate the inputs...")
    #Append the meshes
    appendFilter = vtk.vtkAppendPolyData()
    for polyData in polyDataList : 
        appendFilter.AddInputData(polyData)
    appendFilter.Update()
    # Remove any duplicate points.
    if removeDupliatePoints :
        print("removing the dupliated points...")
        cleanFilter = vtk.vtkCleanPolyData() 
        cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
        cleanFilter.Update()
        return cleanFilter.GetOutput()
    else :
        return appendFilter.GetOutput()

def appendNormal(polyData, verbose, autoOrient ):
    if verbose :
        print("create normals ...")
    # create normals 
    normals =  vtk.vtkPolyDataNormals()
    normals.SetInputData(polyData)
    normals.ComputeCellNormalsOn()
    normals.SetSplitting(0)
    if  not autoOrient :
        normals.AutoOrientNormalsOn()
    normals.ConsistencyOn()
    normals.Update()
    polyData = normals.GetOutput()
    return polyData

def getDataFrameAeroData (polyData, config ):
    Narray = polyData.GetCellData().GetNumberOfArrays()
    available_vars = []
    for i in range(0, Narray ):
        nameVar = polyData.GetCellData().GetArrayName(i)
        available_vars.append(nameVar)
    checkOK = True 
    for key in config.varDict.keys():
        if not config.varDict[key] in available_vars :
            print("could not find:  %s"%config.varDict[key])
            checkOK = False 
    if not checkOK : 
        print("Could not find at least on of the following:" )
        for key in config.varDict.keys():
            print("%30s : %-30s"%(key,config.varDict[key] ))
        print("Here is what is available:" )
        for i in range(0, Narray ):
            nameVar = polyData.GetCellData().GetArrayName(i)
            print("%3d : %s"%(i, polyData.GetCellData().GetArrayName(i)))
        exit(1)
    else : 
        print("Use the following data:" )
        for key in config.varDict.keys():
            print("%30s : %-30s"%(key,config.varDict[key] ))
    dict_new = {}
    if config.IsPressureIntegrationPossible() :
        cps  = polyData.GetCellData().GetArray(config.varDict["Cp"])
        cps = vtk_to_numpy(cps)
        dict_new['Cp'] =  cps 
    if config.IsViscousIntegrationPossible() :
        cfxs = polyData.GetCellData().GetArray(config.varDict["Cfx"])
        cfys = polyData.GetCellData().GetArray(config.varDict["Cfy"])
        cfzs = polyData.GetCellData().GetArray(config.varDict["Cfz"])
        cfxs = vtk_to_numpy(cfxs)
        cfys = vtk_to_numpy(cfys)
        cfzs = vtk_to_numpy(cfzs)
        dict_new['Cfx'] = cfxs 
        dict_new['Cfy'] = cfys 
        dict_new['Cfz'] = cfzs
    if config.IsHeatFluxIntegrationPossible() : 
        Qs = polyData.GetCellData().GetArray(config.varDict["q"])
        Qs = vtk_to_numpy(Qs)
        dict_new['q'] = Qs
    if config.IsChemHeatFluxIntegrationPossible() : 
        Qchems = polyData.GetCellData().GetArray(config.varDict["qd"])
        Qchems = vtk_to_numpy(Qchems)
        dict_new['qd'] = Qchems
    df = pd.DataFrame(dict_new) 
    return df 

'''
recover  Normals and center from cells.
'''
def getDataFrameGeo (polyData, verbose, reverse_normal):
    # 
    if verbose :
        print("recover triangle centers...")
    normals = polyData.GetCellData().GetArray("Normals")
    normals = vtk_to_numpy(normals)
    if  not reverse_normal :
        normals *= -1. 
    # 
    cells = polyData.GetPolys()
    nbOfCells = cells.GetNumberOfCells()
    centers_x = []
    centers_y = []
    centers_z = []
    areas = []
    for i in range ( 0, nbOfCells ) :  # TODO this is slow. Remove that loop and 
        tria = polyData.GetCell(i)
        area = tria.ComputeArea()
        center = [0.,0.,0.]
        # TODO : test if tria is a triangle. If not, handle other cases 
        tria.TriangleCenter(tria.GetPoints().GetPoint(0), tria.GetPoints().GetPoint(1), tria.GetPoints().GetPoint(2), center)
        areas.append(area)
        centers_x.append(center[0])
        centers_y.append(center[1])
        centers_z.append(center[2])
    dict_new = {
        'area': areas , 
        'nx': normals[:,0] , 
        'ny': normals[:,1] , 
        'nz': normals[:,2] , 
        'cx': centers_x , 
        'cy': centers_y , 
        'cz': centers_z , 
    }
    df = pd.DataFrame(dict_new) 
    return df 

def getPolyDataByLoadingFile( absolutePathName, fileExtention  ):
        # test if file exists
        if not os.path.isfile(absolutePathName):
            raise Exception('File {:s} does not exist'.format(absolutePathName))
        # Get extension
        basename = os.path.basename(absolutePathName)
        if fileExtention == "": 
            splitname = basename.split(".")
            if len(splitname) < 2:
                print("Could not detect the file extension of {:s}".format(absolutePathName) ) 
                exit()
            fileExtention = splitname[-1]
        else :
            print("FORCED using file extension: {:s}".format(absolutePathName) ) 
        # Select reader
        if fileExtention == 'ply':
            reader = vtk.vtkPLYReader()
        elif fileExtention ==  'stl':
            reader = vtk.vtkSTLReader()
        elif fileExtention == 'vtk':
            reader = vtk.vtkPolyDataReader()
        elif fileExtention == 'vtu':
            #reader = vtk.vtkUnstructuredGridReader()
            reader = vtk.vtkXMLUnstructuredGridReader()
            #reader = vtk.vtkXMLReader()
        elif fileExtention == 'pvd' or fileExtention == 'pvtu':
            reader = vtk.vtkXMLPUnstructuredGridReader()
        elif fileExtention == 'vtp':
            reader = vtk.vtkPolyDataMapper()
        else:
            raise Exception('Filetype ({:s}) must be either "ply", "stl", "vtk", "vtu", "vtp" '.format(fileExtention))
        print("Reading: {:s}".format(absolutePathName) )
        # Load file
        reader.SetFileName(absolutePathName)
        reader.Update()
        output = reader.GetOutput()
        # uGrid -> polydata
        dsSurfaceFilt = vtk.vtkDataSetSurfaceFilter() 
        dsSurfaceFilt.SetInputData(output)
        dsSurfaceFilt.Update()
        vtkPolyData = dsSurfaceFilt.GetOutput()
        return vtkPolyData
        

        