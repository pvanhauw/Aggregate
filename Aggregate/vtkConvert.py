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


        
