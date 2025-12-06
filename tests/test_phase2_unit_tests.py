"""
PHASE 2: Comprehensive Unit Tests for pyable package
Tests for:
  - Imaginable initialization
  - Image transformations (rotate, scale, translate)
  - ROI operations (erode, dilate, merge)
  - LabelMapable center calculations
  - RoiComparison metrics
  - Edge cases (empty images, None inputs, out-of-bounds)
"""

import pytest
import numpy as np
import SimpleITK as sitk
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from pyable.imaginable import Imaginable, SITKImaginable, Roiable, LabelMapable, Fieldable
from pyable.utilizers import RoiComparison


class TestImaginableInitialization:
    """Test Imaginable class initialization with various inputs"""
    
    def test_imaginable_from_image(self):
        """Test creating Imaginable from SimpleITK image"""
        size = [100, 100, 100]
        img = sitk.Image(size, sitk.sitkUInt8)
        img.SetSpacing([1.0, 1.0, 1.0])
        
        imaginable = SITKImaginable(image=img)
        assert imaginable is not None
        assert imaginable.getImageDimension() == 3
        assert list(imaginable.getImage().GetSize()) == size
    
    def test_imaginable_image_properties(self):
        """Test retrieving image properties"""
        size = [64, 64, 64]
        img = sitk.Image(size, sitk.sitkFloat32)
        img.SetSpacing([0.5, 0.5, 0.5])
        img.SetOrigin([10.0, 20.0, 30.0])
        
        imaginable = SITKImaginable(image=img)
        assert imaginable.getImageDimension() == 3
        assert list(imaginable.getImage().GetSize()) == size
        assert list(imaginable.getImage().GetSpacing()) == [0.5, 0.5, 0.5]
        assert list(imaginable.getImage().GetOrigin()) == [10.0, 20.0, 30.0]
    
    def test_imaginable_2d_image(self):
        """Test creating Imaginable from 2D image"""
        size = [100, 100]
        img = sitk.Image(size, sitk.sitkUInt8)
        
        imaginable = SITKImaginable(image=img)
        assert imaginable.getImageDimension() == 2
        assert list(imaginable.getImage().GetSize()) == size


class TestImageTransformations:
    """Test image transformation operations"""
    
    def setup_method(self):
        """Create a test image with known content"""
        self.size = [100, 100, 100]
        self.img = sitk.Image(self.size, sitk.sitkFloat32)
        self.img.SetSpacing([1.0, 1.0, 1.0])
        
        # Fill center with ones for visibility
        center = [s // 2 for s in self.size]
        region_size = [10, 10, 10]
        region = sitk.Image(region_size, sitk.sitkFloat32)
        region.SetSpacing(self.img.GetSpacing())
        sitk.Paste(self.img, region, region.GetSize(), [0, 0, 0],
                   [c - s // 2 for c, s in zip(center, region_size)])
    
    def test_rotate_image_3d(self):
        """Test 3D image rotation"""
        imaginable = SITKImaginable(image=self.img)
        original_size = imaginable.getImage().GetSize()
        
        # Rotate around z-axis
        imaginable_rot = imaginable.getDuplicate()
        imaginable_rot.rotateImage(rotation=[0, 0, 45.0])
        
        # Image should still be valid after rotation
        assert imaginable_rot.getImage() is not None
        assert imaginable_rot.getImageDimension() == 3
    
    def test_translate_image(self):
        """Test image translation"""
        imaginable = SITKImaginable(image=self.img)
        
        # Translate by vector
        T = [5.0, 5.0, 5.0]
        imaginable_trans = imaginable.getDuplicate()
        imaginable_trans.translateImage(T)
        
        # Image should still be valid
        assert imaginable_trans.getImage() is not None
    
    def test_scale_image(self):
        """Test image scaling"""
        imaginable = SITKImaginable(image=self.img)
        original_size = imaginable.getImage().GetSize()
        
        # Scale by factor
        imaginable_scaled = imaginable.getDuplicate()
        imaginable_scaled.scaleImage([1.5, 1.5, 1.5])
        
        # Image should still exist and be valid
        assert imaginable_scaled.getImage() is not None
        assert imaginable_scaled.getImageDimension() == 3


class TestROIOperations:
    """Test ROI/mask operations"""
    
    def setup_method(self):
        """Create test ROI"""
        self.size = [100, 100, 100]
        # Create binary image with sphere
        img = sitk.Image(self.size, sitk.sitkUInt8)
        img = sitk.BinaryThreshold(img, 0, 0, 1)  # Initialize to 1
        
        # Add a sphere in the center
        center = [50, 50, 50]
        radius = 20
        sphere = sitk.Image(self.size, sitk.sitkUInt8)
        for i in range(self.size[0]):
            for j in range(self.size[1]):
                for k in range(self.size[2]):
                    d = ((i-center[0])**2 + (j-center[1])**2 + (k-center[2])**2)**0.5
                    if d <= radius:
                        sphere.SetPixel([i, j, k], 1)
        
        self.roi = Roiable(image=sphere)
    
    def test_roi_erode(self):
        """Test ROI erosion"""
        roi_eroded = self.roi.getDuplicate()
        roi_eroded.erodeRadius(1)
        
        assert roi_eroded.getImage() is not None
        # Eroded volume should be less than or equal to original
        eroded_volume = sitk.GetArrayFromImage(roi_eroded.getImage()).sum()
        original_volume = sitk.GetArrayFromImage(self.roi.getImage()).sum()
        assert eroded_volume <= original_volume
    
    def test_roi_dilate(self):
        """Test ROI dilation"""
        roi_dilated = self.roi.getDuplicate()
        roi_dilated.dilateRadius(1)
        
        assert roi_dilated.getImage() is not None
        # Dilated volume should be greater than or equal to original
        dilated_volume = sitk.GetArrayFromImage(roi_dilated.getImage()).sum()
        original_volume = sitk.GetArrayFromImage(self.roi.getImage()).sum()
        assert dilated_volume >= original_volume


class TestLabelMapable:
    """Test LabelMapable center calculations"""
    
    def setup_method(self):
        """Create test label map"""
        self.size = [100, 100, 100]
        # Create multi-label image
        img = sitk.Image(self.size, sitk.sitkUInt8)
        
        # Label 1: sphere at (30, 30, 30)
        for i in range(20, 40):
            for j in range(20, 40):
                for k in range(20, 40):
                    if ((i-30)**2 + (j-30)**2 + (k-30)**2)**0.5 <= 8:
                        img.SetPixel([i, j, k], 1)
        
        # Label 2: sphere at (70, 70, 70)
        for i in range(60, 80):
            for j in range(60, 80):
                for k in range(60, 80):
                    if ((i-70)**2 + (j-70)**2 + (k-70)**2)**0.5 <= 8:
                        img.SetPixel([i, j, k], 2)
        
        self.labelmap = LabelMapable(image=img)
    
    def test_center_of_gravity_coordinates(self):
        """Test center of gravity coordinates calculation"""
        center = self.labelmap.getCenterOfGravityCoordinates()
        assert center is not None
        assert len(center) == 3
        # Center should be finite
        assert all(not np.isnan(c) for c in center)
    
    def test_center_of_gravity_index(self):
        """Test center of gravity index calculation"""
        center = self.labelmap.getCenterOfGravityIndex()
        assert center is not None
        assert len(center) == 3
        # Center should be within image bounds
        size = self.labelmap.getImage().GetSize()
        assert all(0 <= c < s for c, s in zip(center, size))


class TestRoiComparison:
    """Test ROI comparison metrics"""
    
    def setup_method(self):
        """Create two ROIs for comparison"""
        size = [100, 100, 100]
        
        # Create reference ROI
        ref_img = sitk.Image(size, sitk.sitkUInt8)
        for i in range(30, 70):
            for j in range(30, 70):
                for k in range(30, 70):
                    ref_img.SetPixel([i, j, k], 1)
        
        # Create test ROI (slightly offset)
        test_img = sitk.Image(size, sitk.sitkUInt8)
        for i in range(35, 75):
            for j in range(35, 75):
                for k in range(35, 75):
                    test_img.SetPixel([i, j, k], 1)
        
        self.ref_roi = Roiable(image=ref_img)
        self.test_roi = Roiable(image=test_img)
    
    def test_dice_similarity(self):
        """Test Dice similarity coefficient"""
        comparison = RoiComparison(self.ref_roi, self.test_roi)
        dice = comparison.getDice()
        
        assert dice is not None
        assert 0 <= dice <= 1
    
    def test_jaccard_similarity(self):
        """Test Jaccard similarity"""
        comparison = RoiComparison(self.ref_roi, self.test_roi)
        jaccard = comparison.getJaccard()
        
        assert jaccard is not None
        assert 0 <= jaccard <= 1
    
    def test_identical_rois_perfect_overlap(self):
        """Test that identical ROIs have perfect metrics"""
        comparison = RoiComparison(self.ref_roi, self.ref_roi)
        
        dice = comparison.getDice()
        assert dice == 1.0 or abs(dice - 1.0) < 1e-6


class TestEdgeCases:
    """Test edge cases and error handling"""
    
    def test_empty_image_handling(self):
        """Test handling of empty/zero images"""
        empty_img = sitk.Image([10, 10, 10], sitk.sitkUInt8)
        imaginable = SITKImaginable(image=empty_img)
        
        assert imaginable.getImage() is not None
        assert sitk.GetArrayFromImage(imaginable.getImage()).sum() == 0
    
    def test_single_pixel_image(self):
        """Test handling of very small images"""
        tiny_img = sitk.Image([1, 1, 1], sitk.sitkUInt8)
        tiny_img.SetPixel([0, 0, 0], 1)
        
        imaginable = SITKImaginable(image=tiny_img)
        assert imaginable.getImage() is not None
    
    def test_none_input_handling(self):
        """Test that None inputs are handled gracefully"""
        imaginable = SITKImaginable(image=None, filename=None)
        # Assuming it handles None gracefully without raising
    
    def test_large_image_operations(self):
        """Test operations on large images"""
        large_img = sitk.Image([256, 256, 256], sitk.sitkUInt8)
        imaginable = SITKImaginable(image=large_img)
        
        assert imaginable.getImage() is not None
        assert list(imaginable.getImage().GetSize()) == [256, 256, 256]
    
    def test_image_copy_independence(self):
        """Test that duplicated images are independent"""
        original_img = sitk.Image([50, 50, 50], sitk.sitkFloat32)
        original_img.SetPixel([25, 25, 25], 100)
        
        imaginable1 = SITKImaginable(image=original_img)
        imaginable2 = imaginable1.getDuplicate()
        
        # Modify copy
        mod_arr = sitk.GetArrayFromImage(imaginable2.getImage())
        mod_arr[25, 25, 25] = 50
        imaginable2.setImageFromNumpy(mod_arr)
        
        # Original should be unchanged
        assert imaginable1.getImage().GetPixel([25, 25, 25]) == 100


class TestPixelTypeConversions:
    """Test pixel type conversions"""
    
    def test_change_pixel_type_uint8_to_float(self):
        """Test converting from UInt8 to Float32"""
        img = sitk.Image([10, 10, 10], sitk.sitkUInt8)
        img.SetPixel([5, 5, 5], 100)
        
        imaginable = SITKImaginable(image=img)
        imaginable.changePixelType(sitk.sitkFloat32)
        
        assert imaginable.getImage().GetPixelIDTypeAsString() == '32-bit float'
        assert imaginable.getImage().GetPixel([5, 5, 5]) == 100
    
    def test_change_pixel_type_preserves_values(self):
        """Test that pixel values are preserved during conversion"""
        img = sitk.Image([20, 20, 20], sitk.sitkUInt16)
        img.SetPixel([10, 10, 10], 500)
        
        imaginable = SITKImaginable(image=img)
        original_value = imaginable.getImage().GetPixel([10, 10, 10])
        
        imaginable.changePixelType(sitk.sitkFloat64)
        new_value = imaginable.getImage().GetPixel([10, 10, 10])
        
        assert abs(original_value - new_value) < 1e-6


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
