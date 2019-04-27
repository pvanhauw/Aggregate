'''
Created on Jan 19, 2019

@author: pierre
'''
import os 
import pandas as pd  
import vtk 
#from Aggregate.geometry.Point import Point
import numpy as np
import vtkTransform

from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
#https://www.programcreek.com/python/example/108192/vtk.util.numpy_support.vtk_to_numpy

import vtkHelper

# best practice based on enum
# https://stackoverflow.com/questions/702834/whats-the-common-practice-for-enums-in-python
class CellOrientationXYZ:
    X_MINUS, X_PLUS, Y_MINUS, Y_PLUS, Z_MINUS, Z_PLUS = list( range (1 , 1 + 6)) 
    labels = {}
    labels[X_MINUS] = "minusX"
    labels[X_PLUS]  = "plusX"
    labels[Y_MINUS] = "minusY"
    labels[Y_PLUS]  = "plusY"
    labels[Z_MINUS] = "minusZ"
    labels[Z_PLUS]  = "plusZ"

def concatenatePolyData( polyDataList, removeDupliatePoints):
    print("concatenate the inputs...")
    #Append the meshes
    appendFilter = vtk.vtkAppendFilter()
    appendFilter.MergePointsOn()
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

def computeCoG(df):
    w =  df["area"].sum()
    x = df["cx"] * df["area"] / w 
    y = df["cy"] * df["area"] / w 
    z = df["cz"] * df["area"] / w 
    print("Cog coordinate (area weighted): %f %f %f"%(x.sum(),y.sum(),z.sum()))
    x = df["cx"].mean()
    y = df["cy"].mean()
    z = df["cz"].mean()
    print("Cog coordinate                : %f %f %f"%(x,y,z))

def ComputeIntegrationQOIsAllCells(df, config):
    # recover matrix for aeroframe transformation 
    R = vtkTransform.GetRotationMatrixFromObjetFrameIntoAeroFrame(config.alpha_deg, config.beta_deg, config.aeroframeAlphaBetaConvention)
    xE1, yE1, zE1, xE2, yE2, zE2, xE3, yE3, zE3 = R[0,0], R[1,0], R[2,0], R[0,1], R[1,1], R[2,1], R[0,2], R[1,2], R[2,2]
    # https://realpython.com/fast-flexible-pandas/ 
    # PRESSURE 
    if config.IsPressureIntegrationPossible() :
        if config.verbose :
            print("pressure ...") 
        df["CxP"] = df["Cp"] * df["nx"] * df["area"] / config.Sref
        df["CyP"] = df["Cp"] * df["ny"] * df["area"] / config.Sref
        df["CzP"] = df["Cp"] * df["nz"] * df["area"] / config.Sref
        df["ClP"] = (df["cy"] * df["CzP"] - df["cz"] * df["CyP"] ) / config.Sref / config.Lref
        df["CmP"] = (df["cz"] * df["CxP"] - df["cx"] * df["CzP"] ) / config.Sref / config.Lref
        df["CnP"] = (df["cx"] * df["CyP"] - df["cy"] * df["CxP"] ) / config.Sref / config.Lref
        df["CDP"] = (df["CxP"] * xE1 + df["CyP"] * yE1 + df["CzP"]* zE1 ) / config.Sref
        df["CQP"] = (df["CxP"] * xE2 + df["CyP"] * yE2 + df["CzP"]* zE2 ) / config.Sref
        df["CLP"] = (df["CxP"] * xE3 + df["CyP"] * yE3 + df["CzP"]* zE3 ) / config.Sref
    # VISCOUS 
    if config.IsViscousIntegrationPossible() :
        if config.verbose :
            print("Viscous ...") 
        df["CxV"] = df["Cfx"] * df["area"] / config.Sref
        df["CyV"] = df["Cfy"] * df["area"] / config.Sref
        df["CzV"] = df["Cfz"] * df["area"] / config.Sref
        df["ClV"] = (df["cy"] * df["Cfz"] - df["cz"] * df["Cfy"] ) * df["area"] / config.Sref / config.Lref
        df["CmV"] = (df["cz"] * df["Cfx"] - df["cx"] * df["Cfz"] ) * df["area"] / config.Sref / config.Lref
        df["CnV"] = (df["cx"] * df["Cfy"] - df["cy"] * df["Cfx"] ) * df["area"] / config.Sref / config.Lref
        df["CDV"] = df["CxV"] * xE1 + df["CyV"] * yE1 + df["CzV"]* zE1 / config.Sref
        df["CQV"] = df["CxV"] * xE2 + df["CyV"] * yE2 + df["CzV"]* zE2 / config.Sref
        df["CLV"] = df["CxV"] * xE3 + df["CyV"] * yE3 + df["CzV"]* zE3 / config.Sref
    # TOTAL
    if config.IsPressureIntegrationPossible() or config.IsViscousIntegrationPossible() :
        for j in config.infixes : 
            varTot= "%s%s%s"%("C",j,"")
            varP  = "%s%s%s"%("C",j,"P")
            varV  = "%s%s%s"%("C",j,"V")
            df[varTot] = 0 
            if config.IsPressureIntegrationPossible() :
                df[varTot] += df[varP] 
            if config.IsViscousIntegrationPossible() :
                df[varTot] += df[varV] 
    # ENERGY 
    if config.IsHeatFluxIntegrationPossible() : 
        if config.verbose :
            print("Flux ...") 
        df["Q"] = df["q"] * df["area"]
    if config.IsChemHeatFluxIntegrationPossible() : 
        if config.verbose :
            print("Flux chem ...") 
        df["Qd"] = df["qd"] * df["area"]
    return df 
            
def getOrientation(df):
    dotX = df["nx"] 
    dotY = df["ny"] 
    dotZ = df["nz"] 
    getOrientationXYZ_vectorized = np.vectorize(getOrientationXYZ)
    df_orientation= getOrientationXYZ_vectorized(dotX, dotY, dotZ)
    return df_orientation

def getOrientationXYZ(dotX, dotY, dotZ):
    if (np.abs(dotX) >= np.abs(dotY)) and (np.abs(dotX) >= np.abs(dotZ)) :
        if dotX > 0 :
            ret = CellOrientationXYZ.X_MINUS 
        else :
            ret = CellOrientationXYZ.X_PLUS 
    elif (np.abs(dotY) >= np.abs(dotX)) and (np.abs(dotY) >= np.abs(dotZ)) :
        if dotY > 0 :
            ret = CellOrientationXYZ.Y_MINUS 
        else :
            ret = CellOrientationXYZ.Y_PLUS 
    else :
        if dotZ > 0 :
            ret = CellOrientationXYZ.Z_MINUS 
        else :
            ret = CellOrientationXYZ.Z_PLUS 
    return ret 

def returnSnappedDirectionSerie2(df , var  ):
    if df["direction"] == df["select_direction"] :
        return df[var]
    else :
        return 0.

def returnSnappedDirectionSerie(df , var ,  direction ):
    df_copy = df[[var, "direction"]]
    snapDirection = df_copy["direction"].isin([direction])
    df_copy.loc[~snapDirection, var] =  0. 
    return df_copy

def IntegrateOverCells(df, qois): 

    # loop over all cell for all QOIs, for all direction
    df_tmps = []
    for direction in CellOrientationXYZ.labels.keys() :
        direction_lab = CellOrientationXYZ.labels[direction]
        data = []
        for qoi in qois  :
            df2  = returnSnappedDirectionSerie(df , qoi ,  direction ) 
            valsum = df2[qoi].sum()
            #print( df[qoi].sum())
            data.append(valsum)
        df_tmp = pd.DataFrame( [data], columns = qois , index = [direction_lab]  )
        df_tmps.append(df_tmp)
    # total : 
    data = []
    for qoi in qois  :
        valsum = df[qoi].sum()
        data.append(valsum)
    df_tmp = pd.DataFrame( [data], columns = qois , index = ["total"]  )
    df_tmps.append(df_tmp)
    df_mat = pd.concat(df_tmps, axis = 0 ) 
    return df_mat 

def ConcateRow(df_integration ):
    header = []
    data = []
    for index in df_integration.index : 
        for col in list(df_integration) :
            newCol =  "%s%s"%(index, col)
            val = df_integration.loc[index,col]
            header.append(newCol)
            data.append(val)
    return pd.DataFrame( [data], columns = header, index = [0]  )

def translateCenters(df, xcog, ycog, zcog):
    computeCoG(df)
    # account for CoG 
    df["cx"] = df["cx"].apply(lambda x : x + xcog)
    df["cy"] = df["cy"].apply(lambda x : x + ycog)
    df["cz"] = df["cz"].apply(lambda x : x + zcog)
    return df 

def plotIntegration(df, config):
    variable_to_plots = []
    if config.IsPressureIntegrationPossible() or config.IsViscousIntegrationPossible() :
        for j in config.infixes : 
            varTot= "%s%s%s"%("C",j,"")
            varP  = "%s%s%s"%("C",j,"P")
            varV  = "%s%s%s"%("C",j,"V")
            variable_to_plots.append(varTot)
            if config.IsPressureIntegrationPossible() :
                variable_to_plots.append(varP)
            if config.IsViscousIntegrationPossible() :
                variable_to_plots.append(varV)
    if len(variable_to_plots) == 0 :
        return 0 
    # A/R 
    height= 1050
    width= 1680 
    scale_ratio = 100
    import matplotlib.pyplot as plt
    fig, axarr = plt.subplots(2,  figsize=(width / scale_ratio ,height/ scale_ratio),  gridspec_kw = {'height_ratios':[4, 6 ] } )
    # 
    YS = [variable_to_plots ] 
    vals = list(CellOrientationXYZ.labels.values())
    i = 0
    for j, Y in  enumerate(YS) : 
        Y = YS[j] 
        ax = axarr[i]
        df_display = df[Y].T["total"]
        ax = df_display.plot.bar(stacked=True, ax = ax , grid = True  )
        i += 1 
        ax = axarr[i]
        df_display = df[Y].T[vals]
        ax = df_display.plot.bar(stacked=True, ax = ax , grid = True  )
    fig.tight_layout()
    fig = ax.get_figure() 
    fig.suptitle( config.outPutName) 
    filename = "%s.png"%(config.outPutName) 
    fig.savefig(filename ) 
    print("wrote plot in : %s"% filename)  
    plt.close('all') 

def integrate(polyData, config ):        
    # recover data for integration
    df_aero = getDataFrameAeroData (polyData, config )
    df_geo = vtkHelper.getDataFrameGeo (polyData, config.verbose, config.reverse_normal ) 
    df_geo = translateCenters(df_geo, - config.xcog, - config.ycog, - config.zcog)
    df_geo['direction'] = getOrientation(df_geo)
    df = pd.concat([df_aero, df_geo], axis = 1 )
    # append to polydata
    direction = numpy_to_vtk(df["direction"].values)
    direction.SetName("direction [-]")
    polyData.GetCellData().AddArray(direction)
    # integration
    df = ComputeIntegrationQOIsAllCells(df, config)
    df_integration = IntegrateOverCells(df, config.getQois())
    # print 
    print("Compute the following quantities: ")
    for x in list(df_integration) :
        print(x)
    print(df_integration)
    # save in csv line style 
    file_name = os.path.join( os.getcwd() , config.outPutName+".csv" )
    print("wrote integration data in : %s"% file_name)
    df_integration_lin = ConcateRow(df_integration   )
    df_integration_lin.to_csv(file_name, sep=';', encoding='utf-8')
    # save in xls column style
    file_name = os.path.join( os.getcwd() , config.outPutName+".xls" )
    print("wrote integration data in : %s"% file_name)
    df_integration.to_excel(file_name,  encoding='utf-8', sheet_name = "integration")
    
    # plot 
    plotIntegration(df_integration, config)
    
    if not config.no_write_vtk :
        vtk_file_name = os.path.join( os.getcwd() , config.outPutName+".vtk" )
        writer = vtk.vtkPolyDataWriter()
        writer.SetInputData(polyData)
        writer.SetFileName(vtk_file_name )
        writer.Write()
        del writer 

        
        