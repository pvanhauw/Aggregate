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
        print("AutoOrientNormalsOn")
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


def getPolyDataByLoadingFile(absolutePathName, fileExtention=""):
        # test if file exists
        if not os.path.isfile(absolutePathName):
            raise Exception(
                'File {:s} does not exist'.format(absolutePathName))
        # Get extension
        basename = os.path.basename(absolutePathName)
        if fileExtention == "":
            splitname = basename.split(".")
            if len(splitname) < 2:
                print("Could not detect the file extension of {:s}".format(
                    absolutePathName))
                exit()
            fileExtention = splitname[-1]
        else:
            print("FORCED using file extension: {:s}".format(absolutePathName))
        # Select reader
        if fileExtention == 'ply':
            reader = vtk.vtkPLYReader()
        elif fileExtention == 'stl':
            reader = vtk.vtkSTLReader()
        elif fileExtention == 'vtk':
            reader = vtk.vtkPolyDataReader()
        elif fileExtention == 'vtu':
            # reader = vtk.vtkUnstructuredGridReader()
            reader = vtk.vtkXMLUnstructuredGridReader()
            # reader = vtk.vtkXMLReader()
        elif fileExtention == 'pvd' or fileExtention == 'pvtu':
            reader = vtk.vtkXMLPUnstructuredGridReader()
        elif fileExtention == 'vtp':
            reader = vtk.vtkPolyDataMapper()
        else:
            raise Exception(
                'Filetype ({:s}) must be either "ply", "stl", "vtk", "vtu", "vtp" '.format(fileExtention))
        print("Reading: {:s}".format(absolutePathName))
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


def writePolyData(polyData, relativePath, outputFormat):
    relativePathBase = relativePath.split(".")
    relativePath = "%s.%s" % (relativePathBase[:-1].Join(), outputFormat)
    if outputFormat == "vtp":
        writer = vtk.vtkXMLUnstructuredGridWriter()
    else:
        raise Exception(
            'Output type ({:s}) unknowed. Type must be either "ply", "stl", "vtk", "vtu", "vtp", "tria" '.format(outputFormat))

    writer.SetFileName(relativePath)
    writer.SetDataModeToAscii()
    writer.Update()


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

# Apply a vtkTransform on a vtkPolyData and return a vtkPolyData
def ApplyTransformationOnPolyData(polyData , transform): 
    transformFilter = vtk.vtkTransformPolyDataFilter() 
    transformFilter.SetTransform(transform  )
    transformFilter.SetInputData(0, polyData )  
    transformFilter.Update() 
    return transformFilter.GetOutput()

# Apply a vtkTransform on a vtkPolyData and return a vtkPolyData
def Transform(polyData , args ): 
    if not args.translate == [0., 0., 0.] :
        transform = vtk.vtkTransform()
        tx = args.translate[0]
        ty = args.translate[1]
        tz = args.translate[2]
        transform.Translate( tx,  ty , tz  ) 
        print("translating of : [%s, %s %s] ..."%(tx,ty,tz))
        polyData = ApplyTransformationOnPolyData(polyData, transform)
    if not args.rotate == [1., 0., 0., 0., 1., 0., 0., 0., 1.]  :
        transform = vtk.vtkTransform()
        matrix4x4 = vtk.vtkMatrix4x4() 
        # http://www.glprogramming.com/red/appendixf.html#name1
        # row, col, val
        matrix4x4.SetElement(0,0, args.rotate[0] ) 
        matrix4x4.SetElement(0,1, args.rotate[1] ) 
        matrix4x4.SetElement(0,2, args.rotate[2] ) 
        matrix4x4.SetElement(1,0, args.rotate[3] )  
        matrix4x4.SetElement(1,1, args.rotate[4] ) 
        matrix4x4.SetElement(1,2, args.rotate[5] ) 
        matrix4x4.SetElement(2,0, args.rotate[6] )  
        matrix4x4.SetElement(2,1, args.rotate[7] ) 
        matrix4x4.SetElement(2,2, args.rotate[8] ) 
        # T = (0,0,0) : no translation 
        matrix4x4.SetElement(0,3, 0. )  
        matrix4x4.SetElement(1,3, 0. )  
        matrix4x4.SetElement(2,3, 0. )  
        matrix4x4.SetElement(3,0, 0. )  
        matrix4x4.SetElement(3,1, 0. )  
        matrix4x4.SetElement(3,2, 0. )  
        matrix4x4.SetElement(3,3, 1 )  
        transform.SetMatrix(matrix4x4) 
        print("Apply : X = R * X with \n    [%s, %s %s]\nR = [%s, %s %s]\n    [%s, %s %s]"%(
            args.rotate[0],  args.rotate[1],  args.rotate[2],  
            args.rotate[3],  args.rotate[4],  args.rotate[5],  
            args.rotate[6],  args.rotate[7],  args.rotate[8])) 
        polyData = ApplyTransformationOnPolyData(polyData, transform)

    if not args.rx == 0.  : 
        transform = vtk.vtkTransform()
        transform.RotateX(args.rx) 
        print("rotation around ox of %s [deg] ..."%(args.rx))
        polyData = ApplyTransformationOnPolyData(polyData, transform)
    if not args.ry == 0.  : 
        transform = vtk.vtkTransform()
        transform.RotateY(args.ry) 
        print("rotation around oy of %s [deg] ..."%(args.ry))
        polyData = ApplyTransformationOnPolyData(polyData, transform)
    if not args.rz == 0.  : 
        transform = vtk.vtkTransform()
        transform.RotateX(args.rz) 
        print("rotation around oz of %s [deg] ..."%(args.rz))
        polyData = ApplyTransformationOnPolyData(polyData, transform)


        
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
    for index, row in df.iterrows():
        x = row['x [m]']
        y = row['y [m]']
        z = row['z [m]']
        points = [x , y , z ]
        closestPointDist =  vtk.reference(0.0) 
        testPoint = [ x ,y ,z ] 
        id = kDTree.FindClosestPoint(testPoint, closestPointDist)
        ids.append(id)
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
            print( "append varaible : %s"% var )
            datas.append(croppedData)
            dict_new[var ] = croppedData
    # write data 
    df_to_append = pd.DataFrame(dict_new)
    df_appened = pd.concat([df, df_to_append ] , axis = 1 ) 
    filename = "interpolation.csv"
    df_appened.to_csv(filename ) 
    print("wrote interpolated data info : %s"%filename)
    


    
    
    
        
