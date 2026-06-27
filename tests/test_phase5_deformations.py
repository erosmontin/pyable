"""
Unit tests for deformation and registration functionality in pyable.

Tests cover:
- Basic deformation operations
- Multiple transform types
- Geometry management
- ROI-aware warping
- Label preservation
"""

import unittest
import tempfile
import numpy as np
from pathlib import Path
import sys
import SimpleITK as sitk

# Ensure the workspace copy of pyable is imported before any installed package.
sys.path.insert(0, str(Path(__file__).parent.parent.resolve()))

from pyable import SITKImaginable, Roiable, LabelMapable
from pyable import deformations


class TestDeformationFieldInitialization(unittest.TestCase):
    """Test displacement field creation and initialization."""

    def setUp(self):
        """Create test image."""
        self.size = (64, 64, 64)
        self.origin = (0, 0, 0)
        self.spacing = (1, 1, 1)
        self.direction = (1, 0, 0, 0, 1, 0, 0, 0, 1)
        
        # Create test image
        self.test_image = sitk.Image(self.size, sitk.sitkFloat32)
        self.test_image.SetOrigin(self.origin)
        self.test_image.SetSpacing(self.spacing)
        self.test_image.SetDirection(self.direction)

    def test_initialize_deformation_field_3d(self):
        """Test creating 3D displacement field."""
        df, size, origin, spacing, direction = deformations.initialize_deformation_field(
            self.test_image, vector_components=3
        )
        
        self.assertEqual(df.GetSize(), self.size)
        self.assertEqual(df.GetOrigin(), self.origin)
        self.assertEqual(df.GetSpacing(), self.spacing)
        # Verify it's a vector image by checking it has components
        self.assertTrue(df.GetNumberOfComponentsPerPixel() > 1)

    def test_initialize_deformation_field_2d(self):
        """Test creating 2D displacement field."""
        size_2d = (64, 64)
        test_image_2d = sitk.Image(size_2d, sitk.sitkFloat32)
        
        df, size, origin, spacing, direction = deformations.initialize_deformation_field(
            test_image_2d, vector_components=2
        )
        
        self.assertEqual(len(df.GetSize()), 2)

    def test_initialize_from_file(self):
        """Test displacement field initialization from image file."""

        with tempfile.TemporaryDirectory() as temporary_directory:
            image_path = (
                Path(temporary_directory)
                / "test_displacement_field.nii.gz"
            )

            sitk.WriteImage(
                self.test_image,
                str(image_path),
            )

            df, size, origin, spacing, direction = (
                deformations.initialize_deformation_field(
                    str(image_path)
                )
            )

            self.assertEqual(df.GetSize(), self.size)
            self.assertEqual(size, self.size)
            self.assertEqual(origin, self.origin)
            self.assertEqual(spacing, self.spacing)
            self.assertEqual(direction, self.direction)


class TestTransformToDisplacementField(unittest.TestCase):
    """Test transform to displacement field conversion."""

    def setUp(self):
        """Create test transforms."""
        self.size = (64, 64, 64)
        self.origin = (0, 0, 0)
        self.spacing = (1, 1, 1)
        self.direction = (1, 0, 0, 0, 1, 0, 0, 0, 1)
        
        # Create identity transform
        self.identity_transform = sitk.Transform(3, sitk.sitkIdentity)

    def test_convert_identity_to_displacement_field(self):
        """Test converting identity transform to displacement field."""
        df = deformations.transform_to_displacement_field(
            self.identity_transform,
            self.size,
            self.origin,
            self.spacing,
            self.direction
        )
        
        self.assertEqual(df.GetSize(), self.size)
        
        # Check that all displacements are zero (identity)
        array = sitk.GetArrayFromImage(df)
        # All displacement vectors should be near zero
        self.assertTrue(np.allclose(array, 0, atol=1e-6))

    def test_convert_translation_to_displacement_field(self):
        """Test converting translation transform."""
        # Create translation of (5, 5, 5) mm
        tfm = sitk.TranslationTransform(3, [5, 5, 5])
        
        df = deformations.transform_to_displacement_field(
            tfm,
            self.size,
            self.origin,
            self.spacing,
            self.direction
        )
        
        self.assertEqual(df.GetSize(), self.size)
        
        # Displacement should be constant (5, 5, 5)
        array = sitk.GetArrayFromImage(df)
        # Check that each voxel has approximately (5, 5, 5) displacement
        # array shape is (z, y, x, components) so take mean across spatial dims
        if len(array.shape) == 4:  # 3D with components
            mean_displacement = np.mean(array[:, :, :, :], axis=(0, 1, 2))
        else:
            mean_displacement = np.mean(array, axis=(0, 1, 2) if len(array.shape) > 3 else (0, 1))
        
        np.testing.assert_array_almost_equal(mean_displacement, [5, 5, 5], decimal=0)


class TestDeformationApplication(unittest.TestCase):
    """Test applying deformation to images."""

    def setUp(self):
        """Create test data."""
        self.size = (64, 64, 64)
        
        # Create test image with gradient
        self.test_image = sitk.Image(self.size, sitk.sitkFloat32)
        array = np.arange(64**3, dtype=np.float32).reshape(self.size)
        self.test_image = sitk.GetImageFromArray(array)
        
        # Create identity displacement field
        self.df, _, _, _, _ = deformations.initialize_deformation_field(self.test_image)

    def test_apply_identity_displacement_field(self):
        """Test applying identity displacement field."""
        warped = deformations.apply_deformation_field(
            self.test_image,
            self.df,
            interpolator='linear'
        )
        
        # Should be very close to original
        self.assertEqual(warped.GetSize(), self.test_image.GetSize())

    def test_apply_displacement_field_with_target(self):
        """Test applying displacement field with target geometry."""
        target_image = sitk.Image((32, 32, 32), sitk.sitkFloat32)
        
        warped = deformations.apply_deformation_field(
            self.test_image,
            self.df,
            target_image=target_image,
            interpolator='linear'
        )
        
        self.assertEqual(warped.GetSize(), target_image.GetSize())


class TestImageableDeformationMethods(unittest.TestCase):
    """Test deformation methods on Imaginable objects."""

    def setUp(self):
        """Create test data."""
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # Create test image
        size = (64, 64, 64)
        array = np.random.rand(size[0], size[1], size[2]).astype(np.float32)
        image = sitk.GetImageFromArray(array)
        
        self.image_path = Path(self.temp_dir.name) / 'test_image.nii.gz'
        sitk.WriteImage(image, str(self.image_path))
        
        # Create identity transform
        self.transform = sitk.Transform(3, sitk.sitkIdentity)
        self.transform_path = Path(self.temp_dir.name) / 'identity.tfm'
        sitk.WriteTransform(self.transform, str(self.transform_path))

    def tearDown(self):
        """Clean up."""
        self.temp_dir.cleanup()

    def test_imaginable_apply_transform(self):
        """Test applyTransform method on Imaginable."""
        img = SITKImaginable(str(self.image_path))
        result = img.applyTransform(str(self.transform_path), interpolator='linear')
        
        # Should return self for chaining
        self.assertIsInstance(result, SITKImaginable)
        
        # Image should still have same size
        self.assertEqual(list(img.getImageSize()), [64, 64, 64])

    def test_imaginable_method_chaining(self):
        """Test method chaining with deformation methods."""
        img = SITKImaginable(str(self.image_path))
        
        # Should support method chaining
        result = img.applyTransform(str(self.transform_path)).cast('uint8')
        
        self.assertIsInstance(result, SITKImaginable)

    def test_imaginable_align_geometry(self):
        """Test alignGeometry method."""
        img1 = SITKImaginable(str(self.image_path))
        img2 = SITKImaginable(str(self.image_path))
        
        # Modify geometry of img2
        modified = img2.getImage()
        modified.SetOrigin((10, 10, 10))
        img2.setImage(modified)
        
        # Align to img1
        img2.alignGeometry(img1.getImage())
        
        # Should now have img1's geometry
        self.assertEqual(img2.getImage().GetOrigin(), img1.getImage().GetOrigin())


class TestROIDeformation(unittest.TestCase):
    """Test deformation of ROI/segmentation masks."""

    def setUp(self):
        """Create test ROI."""
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # Create test ROI (binary mask)
        size = (64, 64, 64)
        center = tuple(s // 2 for s in size)
        radius = 10
        
        array = np.zeros(size, dtype=np.uint8)
        # Create sphere
        x, y, z = np.ogrid[0:size[0], 0:size[1], 0:size[2]]
        mask = (x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2 <= radius**2
        array[mask] = 1
        
        roi_image = sitk.GetImageFromArray(array)
        
        self.roi_path = Path(self.temp_dir.name) / 'roi.nii.gz'
        sitk.WriteImage(roi_image, str(self.roi_path))
        
        # Create identity transform
        self.transform = sitk.Transform(3, sitk.sitkIdentity)
        self.transform_path = Path(self.temp_dir.name) / 'identity.tfm'
        sitk.WriteTransform(self.transform, str(self.transform_path))

    def tearDown(self):
        """Clean up."""
        self.temp_dir.cleanup()

    def test_roiable_warp_roi(self):
        """Test warpROI method preserves labels."""
        roi = Roiable(str(self.roi_path))
        
        # Create displacement field
        df, _, _, _, _ = deformations.initialize_deformation_field(roi.getImage())
        df_path = Path(self.temp_dir.name) / 'deform.mha'
        sitk.WriteImage(df, str(df_path))
        
        # Warp ROI
        roi.warpROI(str(df_path))
        
        # Check labels are preserved (should only have 0 and 1)
        array = sitk.GetArrayFromImage(roi.getImage())
        unique_values = np.unique(array)
        self.assertTrue(np.all(np.isin(unique_values, [0, 1])))

    def test_roiable_apply_transform_to_roi(self):
        """Test applyTransformToROI method."""
        roi = Roiable(str(self.roi_path))
        
        # Apply transform
        roi.applyTransformToROI(str(self.transform_path))
        
        # Check labels are preserved
        array = sitk.GetArrayFromImage(roi.getImage())
        unique_values = np.unique(array)
        self.assertTrue(np.all(np.isin(unique_values, [0, 1])))


class TestLabelMapDeformation(unittest.TestCase):
    """Test deformation of multi-label maps."""

    def setUp(self):
        """Create test label map."""
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # Create test label map (3 labels)
        size = (64, 64, 64)
        array = np.zeros(size, dtype=np.uint8)
        
        # Label 1: left hemisphere
        array[:size[0]//2, :, :] = 1
        # Label 2: right hemisphere
        array[size[0]//2:, :, :] = 2
        # Label 3: center sphere
        center = tuple(s // 2 for s in size)
        x, y, z = np.ogrid[0:size[0], 0:size[1], 0:size[2]]
        mask = (x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2 <= 5**2
        array[mask] = 3
        
        labelmap = sitk.GetImageFromArray(array)
        
        self.labelmap_path = Path(self.temp_dir.name) / 'labels.nii.gz'
        sitk.WriteImage(labelmap, str(self.labelmap_path))
        
        # Create identity transform
        self.transform = sitk.Transform(3, sitk.sitkIdentity)
        self.transform_path = Path(self.temp_dir.name) / 'identity.tfm'
        sitk.WriteTransform(self.transform, str(self.transform_path))

    def tearDown(self):
        """Clean up."""
        self.temp_dir.cleanup()

    def test_labelmapable_warp_labelmap(self):
        """Test warpLabelMap preserves all label values."""
        labels = LabelMapable(str(self.labelmap_path))
        
        # Create displacement field
        df, _, _, _, _ = deformations.initialize_deformation_field(labels.getImage())
        df_path = Path(self.temp_dir.name) / 'deform.mha'
        sitk.WriteImage(df, str(df_path))
        
        # Get original labels
        original_labels = set(np.unique(sitk.GetArrayFromImage(labels.getImage())))
        
        # Warp label map
        labels.warpLabelMap(str(df_path))
        
        # Check all original labels are preserved
        warped_labels = set(np.unique(sitk.GetArrayFromImage(labels.getImage())))
        self.assertEqual(original_labels, warped_labels)

    def test_labelmapable_apply_transform(self):
        """Test applyTransformToLabelMap preserves labels."""
        labels = LabelMapable(str(self.labelmap_path))
        
        # Get original labels
        original_labels = set(np.unique(sitk.GetArrayFromImage(labels.getImage())))
        
        # Apply transform
        labels.applyTransformToLabelMap(str(self.transform_path))
        
        # Check all labels preserved
        warped_labels = set(np.unique(sitk.GetArrayFromImage(labels.getImage())))
        self.assertEqual(original_labels, warped_labels)


class TestGeometryAlignment(unittest.TestCase):
    """Test geometry management and alignment."""

    def setUp(self):
        """Create test images with different geometries."""
        self.temp_dir = tempfile.TemporaryDirectory()
        
        size = (64, 64, 64)
        
        # Create reference image
        ref = sitk.Image(size, sitk.sitkFloat32)
        ref.SetOrigin((0, 0, 0))
        ref.SetSpacing((1, 1, 1))
        self.ref_path = Path(self.temp_dir.name) / 'reference.nii.gz'
        sitk.WriteImage(ref, str(self.ref_path))
        
        # Create moving image with different geometry
        mov = sitk.Image(size, sitk.sitkFloat32)
        mov.SetOrigin((10, 10, 10))
        mov.SetSpacing((2, 2, 2))
        self.mov_path = Path(self.temp_dir.name) / 'moving.nii.gz'
        sitk.WriteImage(mov, str(self.mov_path))

    def tearDown(self):
        """Clean up."""
        self.temp_dir.cleanup()

    def test_align_geometry(self):
        """Test align_geometry function."""
        ref = sitk.ReadImage(str(self.ref_path))
        mov = sitk.ReadImage(str(self.mov_path))
        
        # Align
        aligned = deformations.align_geometry(mov, ref)
        
        # Check geometry matches
        self.assertEqual(aligned.GetOrigin(), ref.GetOrigin())
        self.assertEqual(aligned.GetSpacing(), ref.GetSpacing())
        self.assertEqual(aligned.GetDirection(), ref.GetDirection())

    def test_imaginable_align_geometry(self):
        """Test alignGeometry method on Imaginable."""
        ref = SITKImaginable(str(self.ref_path))
        mov = SITKImaginable(str(self.mov_path))
        
        # Align
        mov.alignGeometry(ref.getImage())
        
        # Check geometry matches
        self.assertEqual(
            mov.getImage().GetOrigin(),
            ref.getImage().GetOrigin()
        )


class TestInterpolationMethods(unittest.TestCase):
    """Test different interpolation methods."""

    def setUp(self):
        """Create test data."""
        self.size = (64, 64, 64)
        self.test_image = sitk.Image(self.size, sitk.sitkFloat32)
        
        # Create displacement field
        self.df, _, _, _, _ = deformations.initialize_deformation_field(self.test_image)

    def test_linear_interpolation(self):
        """Test linear interpolation."""
        warped = deformations.apply_deformation_field(
            self.test_image,
            self.df,
            interpolator='linear'
        )
        self.assertIsNotNone(warped)

    def test_nearest_interpolation(self):
        """Test nearest-neighbor interpolation."""
        warped = deformations.apply_deformation_field(
            self.test_image,
            self.df,
            interpolator='nearest'
        )
        self.assertIsNotNone(warped)

    def test_bspline_interpolation(self):
        """Test B-spline interpolation."""
        warped = deformations.apply_deformation_field(
            self.test_image,
            self.df,
            interpolator='bspline'
        )
        self.assertIsNotNone(warped)


if __name__ == '__main__':
    unittest.main()
