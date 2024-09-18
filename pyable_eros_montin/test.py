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

    
    def test_Labelmappable_getCentroid(self):
        center= self.labelmappable.getCentroidIndex()
        self.assertTrue(np.array_equal(center, np.array(self.one3D_size) // 2))

    def test_Roiable_getCentroid(self):
        center= self.roiable.getCentroidIndex()
        self.assertTrue(np.array_equal(center, np.array(self.one3D_size) // 2))
    
    def test_getBoundingBox(self):
        "Test for getBoundingBox method"
        bounding_box = self.imaginable.getBoundingBox()
        print(bounding_box )
        self.assertTrue(np.array_equal(bounding_box, np.array(((0, 0, 0), [f-1 for f in self.one3D_size]))))
if __name__ == '__main__':
    print("\n\n\n")
    print("***************************DA TEST************************************")
    print("\n\n\n")
    unittest.main(verbosity=2)