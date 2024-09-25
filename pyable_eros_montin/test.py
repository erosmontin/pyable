import unittest
import SimpleITK as sitk
import numpy as np
import imaginable as ima
class TestImaginable(unittest.TestCase):
    def setUp(self):
        # Create a test image
        self.one3D_size = (11, 11, 9)
        self.one3D = np.ones(self.one3D_size)
        self.one3D_origin = (-20,-10,-5)
        self.labelmappable = ima.LabelMapable()
        self.labelmappable.setImageFromNumpy(self.one3D)
        self.labelmappable.setImageOrigin(self.one3D_origin)
        self.imaginable = ima.Imaginable()
        self.imaginable.setImageFromNumpy(self.one3D)
        self.imaginable.setImageOrigin(self.one3D_origin)
        self.roiable = ima.Roiable()
        self.roiable.setImageFromNumpy(self.one3D)
        self.roiable.setImageOrigin(self.one3D_origin)
        
    def test_Labelmappable_getCenterOfGravityIndex(self):
        center_gravity = self.labelmappable.getCenterOfGravityIndex()
        self.assertTrue(np.array_equal(center_gravity, [s//2 for s in self.one3D.shape]))
        

    def test_Roiable_getCenterOfGravityIndex(self):
        center_gravity = self.roiable.getCenterOfGravityIndex()
        self.assertTrue(np.array_equal(center_gravity, [s//2 for s in self.one3D.shape]))

    
    def test_Labelmappable_getCentroid(self):
        center= self.labelmappable.getCentroidIndex()
        self.assertTrue(np.array_equal(center, [s//2 for s in self.one3D.shape]))

    def test_Roiable_getCentroid(self):
        center= self.roiable.getCentroidIndex()
        self.assertTrue(np.array_equal(center, [s//2 for s in self.one3D.shape]))
    
    def test_getBoundingBox(self):
        "Test for getBoundingBox method"
        bounding_box = self.imaginable.getBoundingBox()
        print(bounding_box )
        self.assertTrue(np.array_equal(bounding_box, np.array(((0, 0, 0), [f-1 for f in self.one3D_size]))))
    
    def test_padImage(self):
        "Test for padImage method"
        lower_padding = (1, 1, 1)
        upper_padding = (1, 1, 1)
        padding_value = 0
        self.imaginable.padImage(lower_padding, upper_padding, padding_value)
        padded_image = self.imaginable.getImage()
        # Assuming the original image is of size (10, 10, 10), after padding, it should be (12, 12, 12)
        OUTPUT_SIZE = self.one3D_size + np.array(lower_padding) + np.array(upper_padding)
        self.assertEqual(padded_image.GetSize(), tuple(OUTPUT_SIZE))
if __name__ == '__main__':
    print("\n\n\n")
    print("***************************DA TEST************************************")
    print("\n\n\n")
    unittest.main(verbosity=2)