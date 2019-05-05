'''
Created on Jan 19, 2019

@author: pierre
'''
import os 
import argparse
import vtkDisplay
import vtkTransform
import vtkConvert
import vtkIO 
import vtkHelper
import integrate
import cProfile
import pstats
import io
import AeroFrameAlphaBetaConvention
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
        self.autoOrient = not args.noAutoOrient
        self.reverse_normal = args.reverse_normal
        self.varDict = varDict
        self.verbose = args.verbose
        self.outPutName = args.outPutName
        self.Sref = args.Sref 
        self.Lref = args.Lref 
        self.no_write_vtk = args.no_write_vtk
        # defaut choice for SRF -> ARF (Solid Reference Frame -> Aero Reference Frame)
        self.aeroframeAlphaBetaConvention = AeroFrameAlphaBetaConvention.AeroFrameAlphaBetaConvention.SO3_minusAlpha_minusBeta
    
    def IsPressureIntegrationPossible(self):
        possible = True 
        for var in self.pressureVars :
            if var not in  self.varDict.keys():
                possible = False 
                break 
        return possible
    
    def IsViscousIntegrationPossible(self):
        possible = True 
        for var in self.ViscousVars :
            if  var not in  self.varDict.keys():
                possible = False 
                break 
        return possible
    
    def IsHeatFluxIntegrationPossible(self):
        possible = True 
        for var in self.HeatVars :
            if var not in  self.varDict.keys():
                possible = False 
                break 
        return possible

    def IsChemHeatFluxIntegrationPossible(self):
        possible = True 
        for var in self.HeatChemVars :
            if var not in  self.varDict.keys():
                possible = False 
                break 
        return possible

    def DoIntegrate(self):
        return bool(self.IsChemHeatFluxIntegrationPossible() or self.IsHeatFluxIntegrationPossible() or self.IsViscousIntegrationPossible() or self.IsPressureIntegrationPossible()) 
        
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

def profile(fnc):
    
    """A decorator that uses cProfile to profile a function"""
    def inner(*args, **kwargs):
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return inner

#@profile
def main(): 
    parser = argparse.ArgumentParser(description='Read surface file and integrate pressure, friction and other (eg flux)\nElement can be concatenated together provided they have the similar data', formatter_class=argparse.RawTextHelpFormatter)
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
    parser.add_argument("-noAutoOrient", "--noAutoOrient", help="disable used autoorient feature from VTK ",action="store_true") 
    parser.add_argument("-t", "--translate" , nargs=3, metavar=('tx', 'ty', 'tz'), help="translate by (tx, ty, tz)", type=float, default=[0., 0., 0.] , required = False )
    parser.add_argument("-r", "--rotate", nargs=9, metavar=('r11', 'r12', 'r13','r21', 'r22', 'r23','r31', 'r32', 'r33',), 
                        help="Apply R * X, with [r11 r12 r13, r21 r22 r23, r31 r32 r33]", type=float, default=[1., 0., 0., 0., 1., 0., 0., 0., 1.] , required = False )
    parser.add_argument("-rx", help="rotation around ox of angle arg" , type = float , default = 0. ) 
    parser.add_argument("-ry", help="rotation around oy of angle arg" , type = float , default = 0. ) 
    parser.add_argument("-rz", help="rotation around oz of angle arg" , type = float , default = 0. ) 
    parser.add_argument("-nw", "--no_write_vtk" , help="do not write vtk output when concatening" ,  action="store_true") 
    parser.add_argument("-glvar", "--variableToDisplay" , help="variable that will be displayed when using the rendering  window" ,   type = str , default = "") 
    parser.add_argument("-gl", "--openGL_GUI", help="launch openGL window to vizualize your data",action="store_true") 
    parser.add_argument("-e", "--extractFromCVS", help="csv file containing (x [m],y [m],z [m]) data that are used as coordinnates to extract data from the input.\nThe data that is the closest is used", default = "" ,required = False ) 
    parser.add_argument("-v", "--verbose" , help="extra output" ,  action="store_true") 
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
    config = Config( args,  varDict  )
    print("used autoOrient : %s"%config.autoOrient )
    print("used reverse_normal : %s"%config.reverse_normal )
    
    # concatenate data while computing normals individually for each files 
    polyDataList = []
    for relativePath in list_input :
        absolutePath = os.path.join(cwd, relativePath)
        polyData = vtkIO.getPolyDataByLoadingFile(absolutePath , args.forceFormat )
        # triangulate 
        polyData = vtkConvert.triangulate(polyData)
        vtkDisplay.PrintDataArrays(polyData)
        # create and recover normals 
        #if config.DoIntegrate():
        polyData = vtkHelper.appendNormal(polyData, config.verbose, config.autoOrient )
        polyDataList.append(polyData)
    # 
    removeDupliatePoints = False 
    polyDataConcatenated = vtkHelper.concatenatePolyData( polyDataList, removeDupliatePoints)
    #
    polyDataConcatenated = vtkTransform.Transform(polyDataConcatenated , args)
    #
    if config.DoIntegrate():
        integrate.integrate(polyDataConcatenated, config)
    # 
    if not config.no_write_vtk:
        vtk_file_name = os.path.join( os.getcwd() , "%s.vtp"%(config.outPutName) )
        if len(args.list_variableToKeepForWriting ):
            polyDataConcatenated = vtkHelper.getPolyDataCroppedFromData(polyDataConcatenated, args.list_variableToKeepForWriting )
        if config.IsViscousIntegrationPossible() :
            polyDataConcatenated = vtkHelper.appendFrictionVector(polyDataConcatenated, config  )
        vtkIO.writePolyData(polyDataConcatenated, vtk_file_name, "vtk")
    #        
    if not args.extractFromCVS == "" :
        polyDataProbeLocation = vtkHelper.ExtractDataFromTheClosestCellCenter( polyDataConcatenated ,  args.extractFromCVS)
    else : 
        polyDataProbeLocation = None
    # 
    if args.openGL_GUI :
        vtkDisplay.RenderAndInteracte(polyDataConcatenated , config , args.variableToDisplay , polyDataProbeLocation)
  

main()
        
        
