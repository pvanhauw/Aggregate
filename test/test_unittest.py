import unittest 
import vtkTransform

class Test_AeroFrame(unittest.TestCase):
    def test_AlphaBetaTransform(self):
        forceVector = [1.,0, 0]
        alpha_deg = 90.
        beta_deg = 90.
        alpha_first = True 
        res = vtkTransform.GetForceVectorInAeroFrame( forceVector, alpha_deg, beta_deg, alpha_first)
        expectedRes =  [ 0, 0, 1 ]
        #print(expectedRes)
        #print(res)
        for i in [0, 1, 2] :
            self.assertAlmostEqual(res[i], expectedRes [i]  )

    def test_BetaAlphaTransform(self):
        forceVector = [1.,0, 0]
        alpha_deg = 90.
        beta_deg = 90.
        alpha_first = False
        res = vtkTransform.GetForceVectorInAeroFrame( forceVector, alpha_deg, beta_deg, alpha_first)
        expectedRes =  [ 0, 1, 0 ]
        #print(expectedRes)
        #print(res)
        for i in [0, 1, 2] :
            self.assertAlmostEqual(res[i], expectedRes [i]  )



if __name__ == '__main__':
    unittest.main()
