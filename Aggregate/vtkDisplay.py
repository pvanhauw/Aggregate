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

def lookUpTable( lutNum , lookupStyle , reverse = False ):
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(lutNum) 
    ctf = vtk.vtkColorTransferFunction() 
    # https://github.com/Kitware/ParaView/blob/ece600f179965e6c6b93de1c2e3d69f3fa2adbef/Wrapping/Python/paraview/_colorMaps.py
    cl = [] 
    if lookupStyle == "Rainbow Desaturated" :
        ctf.SetColorSpaceToHSV() 
        cl.append([0.278431372549 ,0.278431372549 ,0.858823529412 ])
        cl.append([0.0 ,0.0 ,0.360784313725 ])
        cl.append([0.0 ,1.0 ,1.0 ])
        cl.append([0.0 ,0.501960784314 ,0.0 ])
        cl.append([1.0 ,1.0 ,0.0 ])
        cl.append([1.0 ,0.380392156863 ,0.0 ])
        cl.append([0.419607843137 ,0.0 ,0.0 ])
        cl.append([0.878431372549 ,0.301960784314 ,0.301960784314 ])
    elif lookupStyle == "Rainbow Blended White": 
        ctf.SetColorSpaceToRGB () 
        cl.append([1.0 ,1.0 ,1.0 ])
        cl.append([0.0 ,0.0 ,1.0 ])
        cl.append([0.0 ,1.0 ,1.0 ])
        cl.append([0.0 ,1.0 ,0.0 ])
        cl.append([1.0 ,1.0 ,0.0 ])
        cl.append([1.0 ,0.0 ,0.0 ])
        cl.append([0.878431372549 ,0.0 ,1.0 ])
    elif lookupStyle == "erdc_iceFire_H":
        ctf.SetColorSpaceToLab()
        cl.append([4.05432e-07 ,0 ,5.90122e-06])
        cl.append([0 ,0.120401 ,0.302675])
        cl.append([0 ,0.216583 ,0.524574])
        cl.append([0.0552475 ,0.345025 ,0.6595])
        cl.append([0.128047 ,0.492588 ,0.720288])
        cl.append([0.188955 ,0.641309 ,0.792092])
        cl.append([0.327673 ,0.784935 ,0.873434])
        cl.append([0.60824 ,0.892164 ,0.935547])
        cl.append([0.881371 ,0.912178 ,0.818099])
        cl.append([0.951407 ,0.835621 ,0.449279])
        cl.append([0.904481 ,0.690489 ,0])
        cl.append([0.85407 ,0.510864 ,0])
        cl.append([0.777093 ,0.33018 ,0.00088199])
        cl.append([0.672862 ,0.139087 ,0.00269398])
        cl.append([0.508815 ,0 ,0])
        cl.append([0.299417 ,0.000366289 ,0.000547829])
        cl.append([0.0157519 ,0.00332021 ,4.55569e-08])
    else :
        print("lookupStyle undefined")
        exit(1)
    vv = [float(xx)/float(len(cl)) for xx in range(len(cl))] 
    if reverse :
        vv.reverse() 
    for pt, color in zip(vv, cl): 
        ctf.AddRGBPoint(pt, color[0], color[1], color[2]) 
    for ii, ss in enumerate([float(xx)/float(lutNum) for xx in range(lutNum)]): 
        cc = ctf.GetColor(ss) 
        lut.SetTableValue(ii, cc[0], cc[1], cc[2], 1.0) 
    return lut 
    
def RenderAndInteracte(input , config , variableToDisplay ) : 
    display_variable , VARIABLE = setVariableToDisplay(config , input , variableToDisplay )
    input.GetCellData().SetActiveScalars(VARIABLE)
    mapMesh = vtk.vtkDataSetMapper()
    #mapMesh = vtk.vtkPolyDataMapper()
    mapMesh.SetInputData(input)
    if display_variable :
        mapMesh.SetScalarRange(input.GetCellData().GetArray(VARIABLE).GetRange())
    mapMesh.SetScalarModeToUseCellData()
    
    lut = lookUpTable( lookupStyle = "Rainbow Blended White" ,lutNum = 32) 
    minvalue = input.GetCellData().GetArray(VARIABLE).GetRange()[0]
    maxvalue = input.GetCellData().GetArray(VARIABLE).GetRange()[1]
    lut.SetRange(minvalue, maxvalue) 
    lut.Build()
    mapMesh.SetLookupTable(lut)

    meshActor = vtk.vtkActor()
    meshActor.SetMapper(mapMesh)
    #meshActor.GetProperty().SetRepresentationToWireframe()

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
 
    # create the scalar_bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetOrientationToHorizontal()
    scalar_bar.SetLookupTable(lut)
    
    # create the scalar_bar_widget
    scalar_bar_widget = vtk.vtkScalarBarWidget()
    scalar_bar_widget.SetInteractor(iren)
    scalar_bar_widget.SetScalarBarActor(scalar_bar)
    scalar_bar_widget.On()
    scalar_bar.SetTitle(VARIABLE)

    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1)

    iren.Initialize()
    renWin.Render()
    iren.Start()
        
        
        