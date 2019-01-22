'''
Created on Jan 19, 2019

@author: pierre
'''
import os 
import pandas as pd  
import argparse
from argparse import RawTextHelpFormatter
import vtk 
#from geometry.Point import Point
import numpy as np

from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from numpy import poly
#https://www.programcreek.com/python/example/108192/vtk.util.numpy_support.vtk_to_numpy

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
            
class Config(object):
    # necessary varaible to carried out different types of integration 
    pressureVars    = ["Cp"]
    ViscousVars     = ["Cfx", "Cfy","Cfz"]
    HeatVars        = ["q"]
    HeatChemVars    = ["qd"]
    
    # Cx, Cy, Cz coefficient 
    coord_BON = ["x","y","z"] 
    # CD (drag) , CQ (sideforce), CL (lift)
    coord_AERO = ["D","Q","L"]
    # Cl, Cm, Cn momentum 
    coord_Moments = ["l","m","n"] 
    infixes = coord_BON + coord_AERO + coord_Moments 
    
    def __init__(self, args, varDict) : 
        self.xcog = args.cog[0]
        self.ycog = args.cog[1] 
        self.zcog = args.cog[2]
        self.alpha_deg = args.alpha_deg
        self.beta_deg = args.beta_deg
        self.autoOrient = args.autoOrient
        self.reverse_normal = args.reverse_normal
        self.reverse_aero_convention = args.reverse_aero_convention
        self.varDict = varDict
        self.verbose = args.verbose
        self.outPutName = args.outPutName
        self.Sref = args.Sref 
        self.Lref = args.Lref 
    
    def IsPressureIntegrationPossible(self):
        possible = True 
        for var in self.pressureVars :
            if not var in  self.varDict.keys():
                possible = False 
                break 
        return possible
    
    def IsViscousIntegrationPossible(self):
        possible = True 
        for var in self.ViscousVars :
            if not var in  self.varDict.keys():
                possible = False 
                break 
        return possible
    
    def IsHeatFluxIntegrationPossible(self):
        possible = True 
        for var in self.HeatVars :
            if not var in  self.varDict.keys():
                possible = False 
                break 
        return possible

    def IsChemHeatFluxIntegrationPossible(self):
        possible = True 
        for var in self.HeatChemVars :
            if not var in  self.varDict.keys():
                possible = False 
                break 
        return possible

    def DoIntegrate(self):
        if self.IsChemHeatFluxIntegrationPossible() or self.IsHeatFluxIntegrationPossible() or self.IsViscousIntegrationPossible() or self.IsPressureIntegrationPossible() :
            return True 
        else :
            return False 
    '''
    IJK 
        I = "C" always
        J = x,y,z,D,S,L,m,n,l
        K = P,V,'' ('' -> all) 
    '''
    def getQois(self) :
        Ks = []
        QOIs = [] 
        if self.IsPressureIntegrationPossible() :
            Ks.append("P")
        if self.IsViscousIntegrationPossible() :
            Ks.append("V")
        if self.IsPressureIntegrationPossible() or self.IsViscousIntegrationPossible() :
            Ks.append("")
        for j in self.infixes : 
            for k in Ks :
                Cs = "C%s%s"%(j,k)
                QOIs.append(Cs)
        if self.IsHeatFluxIntegrationPossible() : 
            QOIs.append("Q")
        if self.IsChemHeatFluxIntegrationPossible() : 
            QOIs.append("Qd")
        QOIs.append("area")
        return QOIs 

    def getE123xyz(self ) : 
        import math 
        alpha = math.radians(self.alpha_deg)
        beta = math.radians(self.beta_deg)
        if self.reverse_aero_convention:
            xE1 =  math.cos( beta) * math.cos(alpha)
            yE1 = -math.sin( beta )
            zE1 =  math.cos(beta) * math.sin(alpha) 
            xE2 = -math.cos( alpha) * math.sin(beta)
            yE2 = +math.cos( beta)  
            zE2 =  math.sin(alpha) * math.sin(beta) 
            xE3 = -math.sin( alpha)
            yE3 = 0 
            zE3 = +math.cos( alpha)
        else :
            if (beta == 0 ) :
                xE1 =  math.cos( beta) * math.cos(alpha)
                yE1 =  -math.sin( beta) * math.cos(alpha)
                zE1 =  math.sin(alpha)           
                xE2 = +math.sin( beta)
                yE2 = +math.cos( beta)
                zE2 =  0 
                xE3 = -math.cos( beta) * math.sin(alpha)
                yE3 = +math.sin( beta) * math.sin(alpha)
                zE3 =  math.cos(alpha) 
            else  :
                print ("convention std with beta non null not handled. Exit") 
                exit(32) 
        return xE1, yE1, zE1, xE2, yE2, zE2, xE3, yE3, zE3 
    
def main(): 
    parser = argparse.ArgumentParser(description='Read surface file and integrate pressure, friction and other (eg flux)\nElement can be concatenated together provided they have the similar data', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--list_input", help="list of files" ,  nargs='+', type = str , default = [] ,   required = True  )  
    parser.add_argument("-f", "--forceFormat" , help="force the reader" , type = str , default = "") 
    parser.add_argument("-o", "--outPutName" , help="default name for output" , type = str , default = "output") 
    parser.add_argument("-Cp", "--Cp"     , help="integration variable value for pressure coefficient" , type = str , default = "") 
    parser.add_argument("-Cfx",  help="integration variable for Cfx" , type = str , default = "") 
    parser.add_argument("-Cfy",  help="integration variable for Cfy" , type = str , default = "") 
    parser.add_argument("-Cfz",  help="integration variable for Cfz" , type = str , default = "") 
    parser.add_argument("-q",  help="integration variable for per-surface-quantity (eg flux)" , type = str , default = "") 
    parser.add_argument("-qd", help="integration variable2 for per-surface quantity (eg flux chem)" , type = str , default = "") 
    parser.add_argument("-alpha", "--alpha_deg"   , help="angle of attack" , type = float , default = 0 , required = False  ) 
    parser.add_argument("-beta", "--beta_deg"   , help="angle of side slip" , type = float , default = 0 , required = False  ) 
    parser.add_argument("-cog", nargs=3, metavar=('xcog', 'ycog', 'zcog'), help="coordinnates of the center of gravity for moment computation", type=float, default=[0., 0., 0.] , required = False )
    parser.add_argument("-Sref", help="reference surface" , type = float , default = 1. ) 
    parser.add_argument("-Lref", help="reference length" , type = float , default = 1. ) 
    parser.add_argument("-reverse_normal",  help="reverse normal (multiply all by -1)",action="store_true") 
    parser.add_argument("-reverse_aero_convention",  help="TODO",action="store_true") 
    parser.add_argument("-autoOrient", "--autoOrient", help="used autoorient feature from VTK ",action="store_false") 
    parser.add_argument("-c", "--concatenate", help="concatenate all the files, either to integrate or for using the openGL_GUI",action="store_true") 
    parser.add_argument("-t", "--translate" , nargs=3, metavar=('tx', 'ty', 'tz'), help="translate by (tx, ty, tz)", type=float, default=[0., 0., 0.] , required = False )
    parser.add_argument("-r", "--rotate", nargs=9, metavar=('r11', 'r12', 'r13','r21', 'r22', 'r23','r31', 'r32', 'r33',), 
                        help="Apply R * X, with [r11 r12 r13, r21 r22 r23, r31 r32 r33]", type=float, default=[0., 0., 0., 0., 0., 0., 0., 0., 0.] , required = False )
    parser.add_argument("-v", "--verbose" , help="extra output" ,  action="store_true") 
    parser.add_argument("-gl", "--openGL_GUI", help="launch openGL window to vizualize your data",action="store_true") 
    args = parser.parse_args()
    list_input = args.list_input
    cwd= os.getcwd()
    varDict = {}
    if args.Cp != "" :
        varDict["Cp"] = args.Cp
    if args.Cfx != "" :
        varDict["Cfx"] = args.Cfx
    if args.Cfy != "" :
        varDict["Cfy"] = args.Cfy
    if args.Cfz != "" :
        varDict["Cfz"] = args.Cfz
    if args.q != "" :
        varDict["q"] = args.q
    if args.qd != "" :
        varDict["qd"] = args.qd
    print("used autoOrient : %s"%args.autoOrient )
    print("used reverse_normal : %s"%args.reverse_normal )
    config = Config( args,  varDict  )
    
    if config.DoIntegrate():
        VARIABLE = "direction [-]"
        display_variable = True 
    else : 
       #VARIABLE = "Heat Flux: Net [J/s]"
       #VARIABLE = "Temperature translation [K]"
       #VARIABLE = "Cp [-]"
       VARIABLE = ""
       display_variable = False
 
    polyDataList = []
    for relativePath in list_input :
        absolutePath = os.path.join(cwd, relativePath)
        polyData = getPolyDataByLoadingFile(absolutePath , args.forceFormat , config.verbose )
        # create and recover normals 
        polyData = appendNormal(polyData, config )
        if args.concatenate :
            polyDataList.append(polyData)
        else :
            if config.DoIntegrate():
                integrate(polyData, config)
            if args.openGL_GUI :
                Launch(polyData , VARIABLE, display_variable )
    if args.concatenate :
        removeDupliatePoints = False
        polyDataConcatenated = concatenatePolyData( polyDataList, removeDupliatePoints)
        if config.DoIntegrate():
            integrate(polyDataConcatenated, config)
        if args.openGL_GUI :
            Launch(polyDataConcatenated , VARIABLE, display_variable )
    
    if not args.translate == [0., 0., 0.] :
        pass # TODO
    if not args.translate == [0., 0., 0., 0., 0., 0., 0., 0., 0.]  :
        pass # TODO

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


def appendNormal(polyData, config ):
    if config.verbose :
        print("create normals ...")
    # create normals 
    normals =  vtk.vtkPolyDataNormals()
    normals.SetInputData(polyData)
    normals.ComputeCellNormalsOn()
    normals.SetSplitting(0)
    if  not config.autoOrient :
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
            print("counod not find:  %s"%config.varDict[key])
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

'''
recover  Normals and center from cells.
'''
def getDataFrameGeo (polyData, config ):
    # 
    if config.verbose :
        print("recover triangle centers...")
    normals = polyData.GetCellData().GetArray("Normals")
    normals = vtk_to_numpy(normals)
    if  not config.reverse_normal :
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
    computeCoG(df)
    # account for CoG 
    df["cx"] = df["cx"].apply(lambda x : x - config.xcog)
    df["cy"] = df["cy"].apply(lambda x : x - config.ycog)
    df["cz"] = df["cz"].apply(lambda x : x - config.zcog)
    df['direction'] = 0.0
    df['direction'] = getOrientation(df, config)
    return df 

def ComputeIntegrationQOIsAllCells(df, config):
    # recover matrix for aeroframe transformation 
    xE1, yE1, zE1, xE2, yE2, zE2, xE3, yE3, zE3 = config.getE123xyz()
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
            
def getOrientation(df, config):
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

def IntegrateOverCells(df, config): 
    # QOis
    qois = config.getQois()
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

def ConcateRow(df_integration , config):
    header = []
    data = []
    for index in df_integration.index : 
        for col in list(df_integration) :
            newCol =  "%s%s"%(index, col)
            val = df_integration.loc[index,col]
            header.append(newCol)
            data.append(val)
    return pd.DataFrame( [data], columns = header, index = [0]  )

def integrate(polyData, config ):        
    # recover data for integration
    df_aero = getDataFrameAeroData (polyData, config )
    df_geo = getDataFrameGeo (polyData, config )
    df = pd.concat([df_aero, df_geo], axis = 1 )
    # append to polydata
    direction = numpy_to_vtk(df["direction"].values)
    direction.SetName("direction [-]")
    polyData.GetCellData().AddArray(direction);
    # integration
    df = ComputeIntegrationQOIsAllCells(df, config)
    df_integration = IntegrateOverCells(df, config)
    # print 
    print("Compute the following quantities: ")
    for x in list(df_integration) :
        print(x)
    print(df_integration)
    # save 
    file_name = os.path.join( os.getcwd() , config.outPutName+".csv" )
    # csv 
    #df_integration.to_csv(file_name, sep=';', encoding='utf-8')
    # csv lin style 
    df_integration_lin = ConcateRow(df_integration , config )
    df_integration_lin.to_csv(file_name, sep=';', encoding='utf-8')
    
    vtk_file_name = os.path.join( os.getcwd() , config.outPutName+".vtk" )
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polyData)
    writer.SetFileName(vtk_file_name )
    writer.Write()
    del writer 

def getPolyDataByLoadingFile( absolutePathName, fileExtention,  verbose = False ):
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
        
def Launch(input , VARIABLE, display_variable) : 
    input.GetCellData().SetActiveScalars(VARIABLE)
    mapMesh = vtk.vtkDataSetMapper()
    #mapMesh = vtk.vtkPolyDataMapper()
    mapMesh.SetInputData(input)
    #mapMesh.SetInputConnection(input)
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
        
main()
        
        