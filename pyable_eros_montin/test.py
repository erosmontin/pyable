import unittest
import SimpleITK as sitk
import numpy as np
import imaginable as ima
class TestImaginable(unittest.TestCase):
    def setUp(self):
        # Create a test image
        self.one3D_size = (11, 11, 9)
        self.one3D = np.ones(self.one3D_size)
        self.labelmappable = ima.LabelMapable()
        self.labelmappable.setImageFromNumpy(self.one3D)
        self.imaginable = ima.Imaginable()
        self.imaginable.setImageFromNumpy(self.one3D)
        self.roiable = ima.Roiable()
        self.roiable.setImageFromNumpy(self.one3D)
        
    def test_Labelmappable_getCenterOfGravityIndex(self):
        center_gravity = self.labelmappable.getCenterOfGravityIndex()
        self.assertTrue(np.array_equal(center_gravity, np.array(self.one3D_size) // 2))
        

    def test_Roiable_getCenterOfGravityIndex(self):
        center_gravity = self.roiable.getCenterOfGravityIndex()
        self.assertTrue(np.array_equal(center_gravity, np.array(self.one3D_size) // 2))

    
    def test_Labelmappable_getCenteroid(self):
        center= self.labelmappable.getCenteroidIndex()
        self.assertTrue(np.array_equal(center, np.array(self.one3D_size) // 2))

    def test_Roiable_getCenteroid(self):
        center= self.roiable.getCenteroidIndex()
        self.assertTrue(np.array_equal(center, np.array(self.one3D_size) // 2))
                
if __name__ == '__main__':
    print("\n\n\n")
    print("***************************DA TEST************************************")
    print("\n\n\n")
    unittest.main(verbosity=2)