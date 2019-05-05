'''
Created on Jan 19, 2019

@author: pierre
'''
import os
import vtk


def uGrid2PpolyData(uGrid) :
    dsSurfaceFilt = vtk.vtkDataSetSurfaceFilter()
    dsSurfaceFilt.SetInputData(uGrid)
    dsSurfaceFilt.Update()
    vtkPolyData = dsSurfaceFilt.GetOutput()
    return vtkPolyData

def triangulate(polyData)  :
    triangulateFilter = vtk.vtkTriangleFilter()
    triangulateFilter.SetInputData(polyData)
    triangulateFilter.Update()
    polyData = triangulateFilter.GetOutput()
    return polyData

def getUgridFromPolyData(polyData):
    appendFilter = vtk.vtkAppendFilter() 
    appendFilter.AddInputData(polyData)
    appendFilter.Update()
    ug = vtk.vtkUnstructuredGrid()
    ug.ShallowCopy(appendFilter.GetOutput())
    return ug
