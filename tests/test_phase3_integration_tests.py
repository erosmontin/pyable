"""
PHASE 3: Integration Tests for pyable package
Tests for:
  - Full pipeline: load → transform → save
  - VTK conversions
  - Image overlays
"""

import pytest
import numpy as np
import SimpleITK as sitk
import sys
import os
import tempfile

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from pyable.imaginable import Imaginable, SITKImaginable, ROIable, LabelMapable
from pyable.meshable import vtk2sitk, sitk2vtk

try:
    import vtk
    VTK_AVAILABLE = True
except ImportError:
    VTK_AVAILABLE = False


class TestFullPipeline:
    """Test complete image processing pipeline"""
    
    def setup_method(self):
        """Create test data"""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test image
        self.size = [100, 100, 100]
        self.img = sitk.Image(self.size, sitk.sitkUInt8)
        
        # Add data to image
        for i in range(30, 70):
            for j in range(30, 70):
                for k in range(30, 70):
                    self.img.SetPixel([i, j, k], 100)
    
    def teardown_method(self):
        """Clean up temp files"""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_load_create_transform_save(self):
        """Test pipeline: create → transform → save"""
        # Step 1: Create
        imaginable = SITKImaginable(image=self.img)
        assert imaginable.getImage() is not None
        
        # Step 2: Transform (scale)
        imaginable_scaled = imaginable.getDuplicate()
        imaginable_scaled.scaleImage([0.5, 0.5, 0.5])
        
        # Step 3: Save
        output_path = os.path.join(self.temp_dir, "test_output.nii.gz")
        sitk.WriteImage(imaginable_scaled.getImage(), output_path)
        
        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 0
        
        # Step 4: Reload and verify
        loaded_img = sitk.ReadImage(output_path)
        assert loaded_img is not None
    
    def test_multi_step_transformation(self):
        """Test multiple transformations in sequence"""
        imaginable = SITKImaginable(image=self.img)
        
        # Chain transformations
        imaginable.rotateImage(angle=15.0)
        original_volume = sitk.GetArrayFromImage(imaginable.getImage()).sum()
        
        imaginable_scaled = imaginable.getDuplicate()
        imaginable_scaled.scaleImage([1.1, 1.1, 1.1])
        
        imaginable_translated = imaginable_scaled.getDuplicate()
        imaginable_translated.translateImage([5.0, 5.0, 5.0])
        
        # All should be valid
        assert imaginable.getImage() is not None
        assert imaginable_scaled.getImage() is not None
        assert imaginable_translated.getImage() is not None
    
    def test_image_arithmetic_operations(self):
        """Test arithmetic operations on images"""
        imaginable1 = SITKImaginable(image=self.img)
        imaginable2 = SITKImaginable(image=self.img)
        
        # Addition
        img_sum = imaginable1.addImage(imaginable2)
        assert img_sum is not None
        
        # Subtraction
        img_diff = imaginable1.getDuplicate()
        img_diff.subtractImage(imaginable2)
        assert img_diff is not None
        
        # Multiplication
        img_mult = imaginable1.getDuplicate()
        img_mult.multiplyImage(2.0)
        assert img_mult is not None


@pytest.mark.skipif(not VTK_AVAILABLE, reason="VTK not available")
class TestVTKConversions:
    """Test VTK to SimpleITK conversions"""
    
    def setup_method(self):
        """Create test VTK image"""
        if not VTK_AVAILABLE:
            pytest.skip("VTK not available")
        
        # Create VTK image
        self.vtk_image = vtk.vtkImageData()
        self.vtk_image.SetDimensions(11, 11, 11)
        self.vtk_image.SetSpacing(1.0, 1.0, 1.0)
        self.vtk_image.SetOrigin(0.0, 0.0, 0.0)
        
        # Fill with test data
        scalars = vtk.vtkImageData()
        scalars.SetDimensions(11, 11, 11)
    
    def test_vtk_to_sitk_conversion(self):
        """Test converting VTK image to SimpleITK"""
        if not VTK_AVAILABLE:
            pytest.skip("VTK not available")
        
        # Create VTK image
        vtk_img = vtk.vtkImageData()
        vtk_img.SetDimensions(11, 11, 11)
        vtk_img.SetSpacing(1.0, 1.0, 1.0)
        
        # Add scalars
        scalars = vtk.vtkUnsignedCharArray()
        scalars.SetNumberOfComponents(1)
        scalars.SetNumberOfTuples(11 * 11 * 11)
        for i in range(11 * 11 * 11):
            scalars.SetValue(i, 100)
        
        vtk_img.GetPointData().SetScalars(scalars)
        
        # Convert
        sitk_img = vtk2sitk(vtk_img)
        
        assert sitk_img is not None
        assert sitk_img.GetDimension() == 3
    
    def test_sitk_to_vtk_conversion(self):
        """Test converting SimpleITK image to VTK"""
        if not VTK_AVAILABLE:
            pytest.skip("VTK not available")
        
        # Create SimpleITK image
        sitk_img = sitk.Image([10, 10, 10], sitk.sitkUInt8)
        sitk_img.SetSpacing([1.0, 1.0, 1.0])
        sitk_img.SetOrigin([0.0, 0.0, 0.0])
        
        # Convert
        vtk_img = sitk2vtk(sitk_img)
        
        assert vtk_img is not None
        assert vtk_img.GetNumberOfPoints() > 0
    
    @pytest.mark.skipif(not VTK_AVAILABLE, reason="VTK not available")
    def test_roundtrip_conversion(self):
        """Test SITK → VTK → SITK roundtrip"""
        if not VTK_AVAILABLE:
            pytest.skip("VTK not available")
        
        # Create original
        original = sitk.Image([20, 20, 20], sitk.sitkUInt8)
        original.SetSpacing([1.0, 1.0, 1.0])
        
        # Fill with test pattern
        arr = sitk.GetArrayFromImage(original)
        arr[5:15, 5:15, 5:15] = 100
        original = sitk.GetImageFromArray(arr)
        original.SetSpacing([1.0, 1.0, 1.0])
        
        # SITK → VTK
        vtk_img = sitk2vtk(original)
        
        # VTK → SITK
        recovered = vtk2sitk(vtk_img)
        
        assert recovered is not None
        assert recovered.GetDimension() == original.GetDimension()


class TestImageOverlays:
    """Test image overlay operations"""
    
    def setup_method(self):
        """Create test images for overlay"""
        size = [100, 100, 100]
        
        # Create base image
        self.base_img = sitk.Image(size, sitk.sitkUInt8)
        for i in range(20, 80):
            for j in range(20, 80):
                for k in range(20, 80):
                    self.base_img.SetPixel([i, j, k], 100)
        
        # Create overlay image (mask)
        self.overlay_img = sitk.Image(size, sitk.sitkUInt8)
        for i in range(30, 70):
            for j in range(30, 70):
                for k in range(30, 70):
                    self.overlay_img.SetPixel([i, j, k], 1)
    
    def test_mask_application(self):
        """Test applying mask to image"""
        base = SITKImaginable(image=self.base_img)
        mask = ROIable(image=self.overlay_img)
        
        # Apply mask
        masked = base.getDuplicate()
        masked_arr = sitk.GetArrayFromImage(masked.getImage())
        mask_arr = sitk.GetArrayFromImage(mask.getImage())
        
        masked_arr[mask_arr == 0] = 0
        masked.setImageFromNumpy(masked_arr)
        
        assert masked.getImage() is not None
    
    def test_multi_label_overlay(self):
        """Test overlay with multiple labels"""
        size = [100, 100, 100]
        
        # Create multi-label image
        multi_label = sitk.Image(size, sitk.sitkUInt8)
        
        # Label 1
        for i in range(20, 50):
            for j in range(20, 50):
                for k in range(20, 50):
                    multi_label.SetPixel([i, j, k], 1)
        
        # Label 2
        for i in range(50, 80):
            for j in range(50, 80):
                for k in range(50, 80):
                    multi_label.SetPixel([i, j, k], 2)
        
        labelmap = LabelMapable(image=multi_label)
        img = labelmap.getImage()
        
        assert img is not None
        # Verify multiple labels exist
        arr = sitk.GetArrayFromImage(img)
        assert 1 in arr
        assert 2 in arr


class TestImageResamplingAndAlignment:
    """Test image resampling and alignment"""
    
    def test_resample_on_target_image(self):
        """Test resampling one image to match another's grid"""
        # Create reference image
        ref = sitk.Image([50, 50, 50], sitk.sitkUInt8)
        ref.SetSpacing([2.0, 2.0, 2.0])
        ref.SetOrigin([10.0, 10.0, 10.0])
        
        # Create moving image
        mov = sitk.Image([100, 100, 100], sitk.sitkUInt8)
        mov.SetSpacing([1.0, 1.0, 1.0])
        mov.SetOrigin([0.0, 0.0, 0.0])
        
        # Resample
        moving_img = SITKImaginable(image=mov)
        moving_img.resampleOnTargetImage(ref)
        
        resampled = moving_img.getImage()
        assert resampled is not None
        assert list(resampled.GetSpacing()) == list(ref.GetSpacing())
    
    def test_crop_image(self):
        """Test image cropping"""
        img = sitk.Image([100, 100, 100], sitk.sitkUInt8)
        
        # Fill entire image
        arr = sitk.GetArrayFromImage(img)
        arr[:] = 100
        img = sitk.GetImageFromArray(arr)
        
        imaginable = SITKImaginable(image=img)
        
        # Crop to center
        lower = [25, 25, 25]
        upper = [75, 75, 75]
        cropped = imaginable.cropImage(lowerB=lower, upperB=upper)
        
        assert cropped is not None
        cropped_size = list(cropped.getImage().GetSize())
        expected_size = [u - l for l, u in zip(lower, upper)]
        assert cropped_size == expected_size
    
    def test_image_pad(self):
        """Test image padding"""
        img = sitk.Image([50, 50, 50], sitk.sitkUInt8)
        
        imaginable = SITKImaginable(image=img)
        padded = imaginable.getPaddedImage([10, 10, 10])
        
        assert padded is not None
        padded_size = list(padded.getImage().GetSize())
        expected_size = [s + 2*p for s, p in zip([50, 50, 50], [10, 10, 10])]
        assert padded_size == expected_size


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
