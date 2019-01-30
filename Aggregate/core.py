'''
Created on Jan 19, 2019

@author: pierre
'''
import os 
import argparse
from argparse import RawTextHelpFormatter
import vtk 
#from Aggregate.geometry.Point import Point

import vtkDisplay as vtkDisplay
import vtkHelper
import integrate
#https://www.programcreek.com/python/example/108192/vtk.util.numpy_support.vtk_to_numpy
            
class Config(object):
    # necessary variables to carried out different types of integration.
    pressureVars    = ["Cp"]
    ViscousVars     = ["Cfx", "Cfy", "Cfz"]
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
        self.no_write_vtk = args.no_write_vtk
    
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
    parser.add_argument("-k", "--list_variableToKeepForWriting", help="list of variable to write in the output vtk file. Everything is written if nothing is specified" ,  nargs='+', type = str , default = [] ,   required = False  )  
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
                        help="Apply R * X, with [r11 r12 r13, r21 r22 r23, r31 r32 r33]", type=float, default=[1., 0., 0., 0., 1., 0., 0., 0., 1.] , required = False )
    parser.add_argument("-v", "--verbose" , help="extra output" ,  action="store_true") 
    parser.add_argument("-nw", "--no_write_vtk" , help="do not write vtk output when concatening" ,  action="store_true") 
    parser.add_argument("-glvar", "--variableToDisplay" , help="variable that will be displayed when using the rendering  window" ,   type = str , default = "") 
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
    
    polyDataList = []
    for relativePath in list_input :
        absolutePath = os.path.join(cwd, relativePath)
        polyData = vtkHelper.getPolyDataByLoadingFile(absolutePath , args.forceFormat )
        # create and recover normals 
        polyData = vtkHelper.appendNormal(polyData, config.verbose, config.autoOrient )
        if args.concatenate :
            polyDataList.append(polyData)
        else :
            if config.DoIntegrate():
                integrate.integrate(polyData, config)
            if args.openGL_GUI :
                vtkDisplay.RenderAndInteracte(polyData ,config , args.variableToDisplay )
    if args.concatenate :
        removeDupliatePoints = False
        polyDataConcatenated = vtkHelper.concatenatePolyData( polyDataList, removeDupliatePoints)
        if config.DoIntegrate():
            integrate.integrate(polyDataConcatenated, config)
    
        if not args.translate == [0., 0., 0.] :
            pass # TODO
        if not args.translate == [1., 0., 0., 0., 1., 0., 0., 0., 1.]  :
            pass # TODO
        if not config.no_write_vtk:
            vtk_file_name = os.path.join( os.getcwd() , "%s-%s.vtk"%(config.outPutName, "concat") )
            print('writting: %s'%vtk_file_name)
            if len(args.list_variableToKeepForWriting ):
                polyDataConcatenated = vtkHelper.getPolyDataCroppedFromData(polyDataConcatenated, args.list_variableToKeepForWriting )
            writer = vtk.vtkPolyDataWriter()
            writer.SetInputData(polyDataConcatenated)
            writer.SetFileName(vtk_file_name )
            writer.SetFileTypeToBinary()
            writer.Write()
        if args.openGL_GUI :
            vtkDisplay.RenderAndInteracte(polyDataConcatenated , config , args.variableToDisplay )
  
main()
        
        