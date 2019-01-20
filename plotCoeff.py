import argparse
import os , sys 
from argparse import RawTextHelpFormatter

import pandas as pd 
# alternatively, add these to your python path
executable = os.environ["ATDB_HOME"]
sys.path.insert(0, executable+'/src/Python')
sys.path.insert(0, executable+'/src/Scripts-matlibplot')

'''
Plot data 
'''
def main(): 
    parser = argparse.ArgumentParser(description='Read csv file. Plot data', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--inputFiles", nargs='+', help="tecplot file intput",  required = True ) 
    parser.add_argument("-t", "--tag", help="tag prefix for output" ,  type = str , default = "") 
    parser.add_argument("-v", "--verbose", help="display bebug",action="store_true") 
    #
    args = parser.parse_args()
    fileRelativePaths   = args.inputFiles 
    tag = args.tag
    verbose = args.verbose
    #
    cwd = os.getcwd() 
    for fileRelativePath in fileRelativePaths : 
        fileAbsoluePath = os.path.join(cwd, fileRelativePath) 
        if os.path.exists(fileAbsoluePath) :
            plot(fileAbsoluePath,  verbose, tag ) 
        else :
            print("file: %s does not exist"%fileAbsoluePath) 

def plot(fileAbsoluePath , verbose,  tag, ) :
    baseName = os.path.basename(fileAbsoluePath)
    #df = pd.read_csv(fileAbsoluePath, skiprows=[0], header=0, sep = ';')
    df = pd.read_csv(fileAbsoluePath,  sep = ';')
    header = list(df) 
    if verbose :
        print (header )
    if len(header) < 2 : 
        print("not enough varaible in %s "%(baseName))
        for i , var in enumerate(header) : 
            print("%3d var:    %s"%(i, var))
        exit(2) 
    if verbose : 
        print("variables are :")
        for i , var in enumerate(header) : 
            print("%3d var:    %s"%(i, var))
    #
    df_tmps = []
    newHeader = getCoeffs () 
    for suff in getSuffixes ()   : # Xm , Ym, ... AL
        data = []
        for coeff in getCoeffs ()   : #  Cx, drag, ... Qd, Q , area 
            #df_tmp = dfT.apply( lambda x : x[0] if suff in  x.name else 0.  )  
            data.append(df.loc[0,"%s%s"%(suff, coeff)]) 
        df_tmp = pd.DataFrame( [data], columns = newHeader , index = [suff]  )
        df_tmps.append(df_tmp)
    df_mat = pd.concat(df_tmps, axis = 0 )
    if verbose :
        print(df_mat)
    YVARS = getCoeffs()
    #YVARS_C = [x for x in YVARS if x not in ['Q', 'Qd' , 'area',"Mcl","Mcm","Mcn" ] ] 
    YVARS_CXYZ  = ["Cx","CxP", "CxV"] 
    YVARS_CXYZ += ["Cy","CyP", "CyV"] 
    YVARS_CXYZ += ["Cz","CzP", "CzV"] 
    YVARS_CDSL =  ["CD","CDP", "CDV"] 
    YVARS_CDSL += ["CQ","CQP", "CQV"] 
    YVARS_CDSL += ["CL","CLP", "CLV"] 
    YVARS_Q = ['Q'] # , 'Qd'] 
    YVARS_A = ['area'] 
    YVARS_M = ["Cl","ClP","ClV" ]
    YVARS_M += ["Cm","CmP","CmV" ]
    YVARS_M += ["Cn","CnP","CnV" ]
    
    XVARS = getSuffixes ()
    XVARS_XYZ = [x for x in XVARS if x not in ['int', 'ext', 'line', 'total' ] ] 
    XVARS_IO  = ['int', 'ext', 'line' ]
    XVARS_ALL = ['total' ] 
    
    # res
    height= 1050
    width= 1680 
    doPlotly = False 
    if doPlotly :
        pass 
    else :
        import matplotlib.pyplot as plt
        fig, axarr = plt.subplots(2, 5, figsize=(width/ 100. ,height/ 100.),  gridspec_kw = {'height_ratios':[4, 6 ] , 'width_ratios':[9,9,3,3, 9] } )
        #fig, ax = plt.subplots( figsize=() )
        # 
        YS = [YVARS_CXYZ ,  YVARS_CDSL  , YVARS_Q, YVARS_A , YVARS_M  ] 
        titles = ["Cx/y/z", "Cdrag/side/lift", "HeatFlux [W]", "area" , "Cl/m/n"]
        for j in range(len(YS)) : 
            i = 0 
            Y = YS[j] 
            ax = axarr[i][j]
            df_display = df_mat[Y].T[XVARS_ALL]
            ax = df_display.plot.bar(stacked=True, ax = ax , grid = True   );
            #
            #i += 1 
            #ax = axarr[i][j]
            #df_display = df_mat[Y].T[XVARS_IO]
            #ax = df_display.plot.bar(stacked=True, ax = ax , grid = True , title = title   );
            #
            i += 1 
            ax = axarr[i][j]
            df_display = df_mat[Y].T[XVARS_XYZ]
            ax = df_display.plot.bar(stacked=True, ax = ax , grid = True , title = titles [j]   );
        fig.tight_layout()
        fig = ax.get_figure() 
        fig.suptitle( "Sref = 1 [m2]"  )
        filename = "%s%s.png"%(tag, baseName.replace(".csv","")) 
        fig.savefig(filename ) 
        print("wrote in : %s"% filename)  
        plt.close('all') 

def getSuffixes ():
    suffixes = []
    suffixes.append("plusX")
    suffixes.append("minusX")
    suffixes.append("plusY")
    suffixes.append("minusY")
    suffixes.append("plusZ")
    suffixes.append("minusZ")
    #suffixes.append("int")
    #suffixes.append("ext")
    #suffixes.append("line")
    suffixes.append("total")
    return suffixes

'''
def getSuffixes2 ():
    suffixes = []
    suffixes.append("Xp")
    suffixes.append("Xm")
    suffixes.append("YP")
    suffixes.append("Ym")
    suffixes.append("Zp")
    suffixes.append("Zm")
    #suffixes.append("inside")
    #suffixes.append("outside")
    #suffixes.append("other")
    #suffixes.append("All")
    return suffixes
'''

def getCoeffs(parseQdiff = False):
    coeffs = []
    coeffs.append("CxP")
    coeffs.append("CyP")
    coeffs.append("CzP")
    coeffs.append("CxV")
    coeffs.append("CyV")
    coeffs.append("CzV")
    coeffs.append("Cx")
    coeffs.append("Cy")
    coeffs.append("Cz")
    coeffs.append("CDP")
    coeffs.append("CQP")
    coeffs.append("CLP")
    coeffs.append("CDV")
    coeffs.append("CQV")
    coeffs.append("CLV")
    coeffs.append("CD")
    coeffs.append("CQ")
    coeffs.append("CL")
    coeffs.append("ClP")
    coeffs.append("CmP")
    coeffs.append("CnP")
    coeffs.append("ClV")
    coeffs.append("CmV")
    coeffs.append("CnV")
    coeffs.append("Cl")
    coeffs.append("Cm")
    coeffs.append("Cn")
    coeffs.append("Q")
    if  parseQdiff  :
        coeffs.append("Qd")
    coeffs.append("area")
    return coeffs

main() 


