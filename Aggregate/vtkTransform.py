'''
Created on Jan 19, 2019

@author: pierre
'''
import math 
import numpy as np 
import vtk
import AeroFrameAlphaBetaConvention

def concatenatePolyData(polyDataList, removeDupliatePoints):
    print("concatenate the inputs...")
    # Append the meshes
    appendFilter = vtk.vtkAppendPolyData()
    for polyData in polyDataList:
        appendFilter.AddInputData(polyData)
    appendFilter.Update()
    # Remove any duplicate points.
    if removeDupliatePoints:
        print("removing the dupliated points...")
        cleanFilter = vtk.vtkCleanPolyData()
        cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
        cleanFilter.Update()
        return cleanFilter.GetOutput()
    else:
        return appendFilter.GetOutput()

# Apply a vtkTransform on a vtkPolyData and return a vtkPolyData
def ApplyTransformationOnPolyData(polyData , transform): 
    transformFilter = vtk.vtkTransformPolyDataFilter() 
    transformFilter.SetTransform(transform  )
    transformFilter.SetInputData(0, polyData )  
    transformFilter.Update() 
    return transformFilter.GetOutput()

# Apply a vtkTransform on a vtkPolyData and return a vtkPolyData
def Transform(polyData , args ): 
    if not args.translate == [0., 0., 0.] :
        transform = vtk.vtkTransform()
        tx = args.translate[0]
        ty = args.translate[1]
        tz = args.translate[2]
        transform.Translate( tx,  ty , tz  ) 
        print("translating of : [%s, %s %s] ..."%(tx,ty,tz))
        polyData = ApplyTransformationOnPolyData(polyData, transform)
    if not args.rotate == [1., 0., 0., 0., 1., 0., 0., 0., 1.]  :
        transform = vtk.vtkTransform()
        matrix4x4 = vtk.vtkMatrix4x4() 
        # http://www.glprogramming.com/red/appendixf.html#name1
        # row, col, val
        matrix4x4.SetElement(0,0, args.rotate[0] ) 
        matrix4x4.SetElement(0,1, args.rotate[1] ) 
        matrix4x4.SetElement(0,2, args.rotate[2] ) 
        matrix4x4.SetElement(1,0, args.rotate[3] )  
        matrix4x4.SetElement(1,1, args.rotate[4] ) 
        matrix4x4.SetElement(1,2, args.rotate[5] ) 
        matrix4x4.SetElement(2,0, args.rotate[6] )  
        matrix4x4.SetElement(2,1, args.rotate[7] ) 
        matrix4x4.SetElement(2,2, args.rotate[8] ) 
        # T = (0,0,0) : no translation 
        matrix4x4.SetElement(0,3, 0. )  
        matrix4x4.SetElement(1,3, 0. )  
        matrix4x4.SetElement(2,3, 0. )  
        matrix4x4.SetElement(3,0, 0. )  
        matrix4x4.SetElement(3,1, 0. )  
        matrix4x4.SetElement(3,2, 0. )  
        matrix4x4.SetElement(3,3, 1. )  
        transform.SetMatrix(matrix4x4) 
        print("Apply : X = R * X with \n    [%s, %s %s]\nR = [%s, %s %s]\n    [%s, %s %s]"%(
            args.rotate[0],  args.rotate[1],  args.rotate[2],  
            args.rotate[3],  args.rotate[4],  args.rotate[5],  
            args.rotate[6],  args.rotate[7],  args.rotate[8])) 
        polyData = ApplyTransformationOnPolyData(polyData, transform)

    if not args.rx == 0.  : 
        transform = vtk.vtkTransform()
        transform.RotateX(args.rx) 
        print("rotation around ox of %s [deg] ..."%(args.rx))
        polyData = ApplyTransformationOnPolyData(polyData, transform)
    if not args.ry == 0.  : 
        transform = vtk.vtkTransform()
        transform.RotateY(args.ry) 
        print("rotation around oy of %s [deg] ..."%(args.ry))
        polyData = ApplyTransformationOnPolyData(polyData, transform)
    if not args.rz == 0.  : 
        transform = vtk.vtkTransform()
        transform.RotateX(args.rz) 
        print("rotation around oz of %s [deg] ..."%(args.rz))
        polyData = ApplyTransformationOnPolyData(polyData, transform)

    return polyData
        

# https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
# https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula 
def rotation_matrix(axis, theta_rad):
    """ 
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta_rad radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta_rad / 2.0)
    b, c, d = -axis * math.sin(theta_rad / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d 
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d 
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])



'''
alpha and beta define the rotation to go from the Solide reference frame to the aero refence frame 
through the parameter alphaFirst 
1) first  rotation is Ry(-alpha)
2) second rotation is applied after the first one: Rz(-beta)  
'''
def GetForceVectorInAeroFrame( forceVector, alpha_deg, beta_deg, aeroFrameAlphaBetaConvention ) :
    R = GetRotationMatrixFromObjetFrameIntoAeroFrame(alpha_deg, beta_deg, aeroFrameAlphaBetaConvention)
    return np.dot(R, forceVector) 

'''
alpha and beta define angle applied for 2 rotations, one being carried out after the other, the order can be specied 
through the parameter alphaFirst 
1) first  rotation is Ry(-alpha)
2) second rotation is applied after the first one: Rz(-beta)  
'''
def GetRotationMatrixFromObjetFrameIntoAeroFrame(alpha_deg, beta_deg, aeroFrameAlphaBetaConvention):
    # rotation of -alpha around oy 
    Ry = rotation_matrix( [0, 1, 0], math.radians(-alpha_deg ))
    # rotation of -beta around oz
    Rzm = rotation_matrix( [0, 0, 1], math.radians(-beta_deg ))
    # rotation of -beta around oz
    Rzp = rotation_matrix( [0, 0, 1], math.radians(beta_deg ))
    if aeroFrameAlphaBetaConvention == AeroFrameAlphaBetaConvention.AeroFrameAlphaBetaConvention.SO3_minusBeta_minusAlpha :
        R = np.matmul(Rzm, Ry) 
    elif aeroFrameAlphaBetaConvention == AeroFrameAlphaBetaConvention.AeroFrameAlphaBetaConvention.SO3_plusBeta_minusAlpha :
        R = np.matmul(Rzp, Ry) 
    else : 
        print("Undefined ot not implemented aeroFrameAlphaBetaConvention: %s"%aeroFrameAlphaBetaConvention)
        exit(2) 
    return R 






