import unittest 
import vtkTransform
import AeroFrameAlphaBetaConvention

class Test_AeroFrame(unittest.TestCase):
    def test_plusBetaMinusALphaTransform(self):
        forceVector = [1.,0, 0]
        alpha_deg = 90.
        beta_deg = 90.
        aeroFrameAlphaBetaConvension = AeroFrameAlphaBetaConvention.AeroFrameAlphaBetaConvention.SO3_plusBeta_minusAlpha
        res = vtkTransform.GetForceVectorInAeroFrame( forceVector, alpha_deg, beta_deg, aeroFrameAlphaBetaConvension)
        expectedRes =  [ 0, 0, 1 ]
        #print(expectedRes)
        #print(res)
        for i in [0, 1, 2] :
            self.assertAlmostEqual(res[i], expectedRes [i]  )

    def test_minusBetaMinusAlphaTransform(self):
        forceVector = [1.,0, 0]
        alpha_deg = 90.
        beta_deg = 90.
        aeroFrameAlphaBetaConvension = AeroFrameAlphaBetaConvention.AeroFrameAlphaBetaConvention.SO3_minusBeta_minusAlpha
        res = vtkTransform.GetForceVectorInAeroFrame( forceVector, alpha_deg, beta_deg, aeroFrameAlphaBetaConvension)
        expectedRes =  [ 0, 0 , 1 ]
        #print(expectedRes)
        #print(res)
        for i in [0, 1, 2] :
            self.assertAlmostEqual(res[i], expectedRes [i]  )

if __name__ == '__main__':
    unittest.main()
