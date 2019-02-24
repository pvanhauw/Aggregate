'''
Created on Jan 19, 2019

@author: pierre
'''
import vtk 

    
def PrintDataArrays(polyData  ):
    Narray = polyData.GetCellData().GetNumberOfArrays()
    nameVars = [] 
    print("Available data")
    for i in range(0, Narray ):
        nameVar = polyData.GetCellData().GetArrayName(i)
        nameVars.append(nameVar)
    nameVars.sort()
    for nameVar in nameVars : 
        print("%15s %s"%(" ", nameVar))

def setVariableToDisplay(config , polyData , variableToDisplay ):
    display_variable = False 
    Narray = polyData.GetCellData().GetNumberOfArrays()
    nameVars = [] 
    VARIABLE = None
    print("Available data")
    for i in range(0, Narray ):
        nameVar = polyData.GetCellData().GetArrayName(i)
        nameVars.append(nameVar)
    nameVars.sort()
    for nameVar in nameVars : 
        if variableToDisplay == nameVar :
            VARIABLE = variableToDisplay
            display_variable = True
            print("%15s %s"%("SELECTED -> ", nameVar))
        else : 
            print("%15s %s"%(" ", nameVar))
    if not display_variable :
        if config.DoIntegrate():
            VARIABLE = "direction [-]"
            display_variable = True 
            print("%15s %s"%("SELECTED -> ", nameVar))
    return display_variable , VARIABLE 
    
def RenderAndInteracte(input , config , variableToDisplay ) : 
    display_variable , VARIABLE = setVariableToDisplay(config , input , variableToDisplay )
    input.GetCellData().SetActiveScalars(VARIABLE)
    mapMesh = vtk.vtkDataSetMapper()
    #mapMesh = vtk.vtkPolyDataMapper()
    mapMesh.SetInputData(input)
    if display_variable :
        mapMesh.SetScalarRange(input.GetCellData().GetArray(VARIABLE).GetRange())
    mapMesh.SetScalarModeToUseCellData()

    meshActor = vtk.vtkActor()
    meshActor.SetMapper(mapMesh)
    meshActor.GetProperty().SetRepresentationToWireframe()

    # Create the rendering window, renderer, and interactive renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
   
    # Add the actors to the renderer, set the background and size
    ren.AddActor(meshActor)
    ren.SetBackground(0, 0, 0)
    renWin.SetSize(1650, 1050)
 
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1)

    iren.Initialize()
    renWin.Render()
    iren.Start()
        
        
        