'''
Created on Jan 19, 2019

@author: pierre
'''
import os
import vtk
import vtkConvert

# from vtk.util.numpy_support import numpy_to_vtk
from vtk.util.numpy_support import vtk_to_numpy
# https://www.programcreek.com/python/example/108192/vtk.util.numpy_support.vtk_to_numpy



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
    vtkPolyData =  vtkConvert.uGrid2PpolyData(output)
    return vtkPolyData


def writePolyData(polyData, relativePath, outputFormat = None):
    defaultExtension = "vtk"
    relativePathBase = relativePath.split(".")
    if len( relativePath ) > 1 :
        extension = relativePathBase[:-1]
    else :
        # no extension found
        extension = defaultExtension
    # override extension 
    if not outputFormat == None :
        extension = defaultExtension
    # rename the file with the correct extension if possible
    if len( relativePath ) > 1 :
        relativePathUsed = "%s.%s" % ("".join(relativePathBase[:-1]), outputFormat)
    else :
        # if not append one 
        relativePathUsed = "%s.%s" % (relativePath, outputFormat)
    
    print('writting: %s'%relativePathUsed)
    if outputFormat == "vtk":
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileTypeToBinary()
        '''
    elif outputFormat == "vtu":
        writer = vtk.vtkXMLUnstructuredGridWriter()
        uGrid = vtkConvert.polyData2Ugrid(polyData)
        '''
    else:
        raise Exception(
            #'Output type ({:s}) unknowed. Type must be either "ply", "stl", "vtk", "vtu", "vtp", "tria" '.format(outputFormat))
            'Output type ({:s}) unknowed. Type must be either "stl", "vtk", "vtu", "vtp", "tria" '.format(outputFormat))
    writer.SetInputData(polyData)
    writer.SetFileName(relativePathUsed)
    #writer.SetDataModeToAscii()
    writer.Write()


        
