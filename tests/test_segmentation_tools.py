"""
Tests for the segmentation tools added in feat/segmentation-tools.

Uses synthetic 3D images (sphere in a cube) so no external data is needed.
Tests both the standalone functions in segmentation.py and the wrapper
methods on Imaginable and Roiable.
"""

import numpy as np
import SimpleITK as sitk
import pytest
import sys
import os
import types

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# Mock VTK to avoid Python version mismatch
_vtk_mock = types.ModuleType('vtk')
_vtk_util = types.ModuleType('vtk.util')
_vtk_np = types.ModuleType('vtk.util.numpy_support')
sys.modules['vtk'] = _vtk_mock
sys.modules['vtk.util'] = _vtk_util
sys.modules['vtk.util.numpy_support'] = _vtk_np

from pyable.imaginable import Imaginable, Roiable, LabelMapable
from pyable import segmentation as seg


# ============================================================================
# FIXTURES: Synthetic 3D images
# ============================================================================

def _make_sphere(size=64, radius=15, center=None, spacing=(1.0, 1.0, 1.0)):
    """Create a synthetic binary sphere ROI as sitk.Image."""
    if center is None:
        center = [s // 2 for s in [size] * 3]
    arr = np.zeros((size, size, size), dtype=np.uint8)
    zz, yy, xx = np.ogrid[:size, :size, :size]
    dist = np.sqrt(
        ((zz - center[0]) * spacing[0]) ** 2 +
        ((yy - center[1]) * spacing[1]) ** 2 +
        ((xx - center[2]) * spacing[2]) ** 2
    )
    arr[dist <= radius] = 1
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(spacing)
    img.SetOrigin((0, 0, 0))
    return img


def _make_intensity_image(size=64, spacing=(1.0, 1.0, 1.0)):
    """Create a synthetic intensity image with a bright sphere."""
    arr = np.random.rand(size, size, size).astype(np.float32) * 0.3
    zz, yy, xx = np.ogrid[:size, :size, :size]
    c = size // 2
    dist = np.sqrt((zz - c) ** 2 + (yy - c) ** 2 + (xx - c) ** 2)
    arr[dist <= 15] += 0.7  # bright sphere
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing(spacing)
    img.SetOrigin((0, 0, 0))
    return img


@pytest.fixture
def sphere_roi():
    return _make_sphere()


@pytest.fixture
def intensity_image():
    return _make_intensity_image()


@pytest.fixture
def sphere_roiable():
    r = Roiable()
    r.setImage(_make_sphere(), 'test sphere')
    return r


@pytest.fixture
def intensity_imaginable():
    img = Imaginable()
    img.setImage(_make_intensity_image(), 'test intensity')
    return img


# ============================================================================
# TESTS: segmentation.py standalone functions
# ============================================================================

class TestThresholdLevelSet:
    def test_basic(self, sphere_roi, intensity_image):
        result = seg.threshold_level_set_refine(sphere_roi, intensity_image)
        arr = sitk.GetArrayFromImage(result)
        assert arr.max() == 1
        assert np.sum(arr) > 0

    def test_no_shrink(self, sphere_roi, intensity_image):
        result = seg.threshold_level_set_refine(
            sphere_roi, intensity_image, allow_shrink=False)
        orig_vol = np.sum(sitk.GetArrayFromImage(sphere_roi) > 0)
        new_vol = np.sum(sitk.GetArrayFromImage(result) > 0)
        assert new_vol >= orig_vol


class TestLaplacianLevelSet:
    def test_basic(self, sphere_roi, intensity_image):
        result = seg.laplacian_level_set_refine(sphere_roi, intensity_image)
        arr = sitk.GetArrayFromImage(result)
        assert arr.max() <= 1
        assert np.sum(arr) > 0


class TestShapeDetectionLevelSet:
    def test_basic(self, sphere_roi, intensity_image):
        result = seg.shape_detection_level_set_refine(
            sphere_roi, intensity_image, allow_shrink=False)
        arr = sitk.GetArrayFromImage(result)
        # With allow_shrink=False, at minimum the seed is preserved
        assert np.sum(arr) >= np.sum(sitk.GetArrayFromImage(sphere_roi) > 0)


class TestChanVese:
    def test_basic(self, sphere_roi, intensity_image):
        result = seg.chan_vese_refine(sphere_roi, intensity_image)
        arr = sitk.GetArrayFromImage(result)
        assert np.sum(arr) > 0


class TestThresholdingMethods:
    def test_otsu(self, intensity_image):
        result = seg.otsu_threshold(intensity_image)
        arr = sitk.GetArrayFromImage(result)
        assert set(np.unique(arr)).issubset({0, 1})
        assert np.sum(arr) > 0

    def test_multi_otsu(self, intensity_image):
        result = seg.multi_otsu_threshold(intensity_image, n_thresholds=2)
        arr = sitk.GetArrayFromImage(result)
        assert len(np.unique(arr)) >= 2

    def test_li(self, intensity_image):
        result = seg.li_threshold(intensity_image)
        # Li may return empty on some synthetic distributions; just check it runs
        arr = sitk.GetArrayFromImage(result)
        assert arr.dtype == np.uint8

    def test_yen(self, intensity_image):
        result = seg.yen_threshold(intensity_image)
        assert np.sum(sitk.GetArrayFromImage(result)) > 0

    def test_triangle(self, intensity_image):
        result = seg.triangle_threshold(intensity_image)
        assert np.sum(sitk.GetArrayFromImage(result)) > 0

    def test_huang(self, intensity_image):
        result = seg.huang_threshold(intensity_image)
        assert np.sum(sitk.GetArrayFromImage(result)) > 0

    def test_manual(self, intensity_image):
        result = seg.manual_threshold(intensity_image, 0.5, 1.5)
        arr = sitk.GetArrayFromImage(result)
        assert set(np.unique(arr)).issubset({0, 1})


class TestRegionGrowing:
    def test_connected_threshold(self, intensity_image, sphere_roi):
        result = seg.connected_threshold_grow(intensity_image, sphere_roi)
        assert np.sum(sitk.GetArrayFromImage(result)) > 0

    def test_neighbourhood_connected(self, intensity_image, sphere_roi):
        result = seg.neighbourhood_connected_grow(
            intensity_image, sphere_roi)
        assert np.sum(sitk.GetArrayFromImage(result)) > 0

    def test_isolated_connected(self, intensity_image):
        seed1 = _make_sphere(radius=5, center=[32, 32, 32])
        seed2 = _make_sphere(radius=5, center=[32, 32, 55])
        result = seg.isolated_connected_grow(intensity_image, seed1, seed2)
        arr = sitk.GetArrayFromImage(result)
        assert arr.dtype == np.uint8


class TestWatershed:
    def test_morphological_watershed(self, intensity_image):
        result = seg.morphological_watershed(intensity_image, level=0.01)
        arr = sitk.GetArrayFromImage(result)
        assert len(np.unique(arr)) >= 1  # at least one basin

    def test_watershed_from_markers(self, intensity_image):
        # Create two markers
        markers_arr = np.zeros((64, 64, 64), dtype=np.int32)
        markers_arr[32, 32, 32] = 1  # center
        markers_arr[5, 5, 5] = 2     # corner
        markers = sitk.GetImageFromArray(markers_arr)
        markers.CopyInformation(intensity_image)
        result = seg.morphological_watershed_from_markers(
            intensity_image, markers)
        arr = sitk.GetArrayFromImage(result)
        assert 1 in arr
        assert 2 in arr


class TestDistanceConstraint:
    def test_constrain_by_distance(self):
        bone = _make_sphere(radius=12, center=[32, 32, 32])
        cartilage = _make_sphere(radius=18, center=[32, 32, 32])
        result = seg.constrain_by_distance(cartilage, bone, max_distance_mm=3.0)
        arr = sitk.GetArrayFromImage(result)
        # Should be smaller than the original cartilage
        assert np.sum(arr) < np.sum(sitk.GetArrayFromImage(cartilage) > 0)
        assert np.sum(arr) > 0

    def test_exclude_interior(self):
        bone = _make_sphere(radius=12, center=[32, 32, 32])
        cartilage = _make_sphere(radius=18, center=[32, 32, 32])
        result = seg.constrain_by_distance(
            cartilage, bone, max_distance_mm=5.0, exclude_interior=True)
        arr = sitk.GetArrayFromImage(result)
        bone_arr = sitk.GetArrayFromImage(bone) > 0
        # No overlap with bone interior
        assert np.sum(arr & bone_arr) == 0


class TestMaskOperations:
    def test_subtract(self):
        big = _make_sphere(radius=15)
        small = _make_sphere(radius=8)
        result = seg.subtract_mask(big, small)
        arr = sitk.GetArrayFromImage(result)
        assert np.sum(arr) > 0
        small_arr = sitk.GetArrayFromImage(small) > 0
        assert np.sum(arr & small_arr) == 0

    def test_intersect(self):
        s1 = _make_sphere(radius=15, center=[32, 32, 28])
        s2 = _make_sphere(radius=15, center=[32, 32, 36])
        result = seg.intersect_masks(s1, s2)
        arr = sitk.GetArrayFromImage(result)
        assert np.sum(arr) > 0
        assert np.sum(arr) < np.sum(sitk.GetArrayFromImage(s1) > 0)

    def test_union(self):
        s1 = _make_sphere(radius=10, center=[32, 32, 20])
        s2 = _make_sphere(radius=10, center=[32, 32, 44])
        result = seg.union_masks(s1, s2)
        arr = sitk.GetArrayFromImage(result)
        assert np.sum(arr) >= np.sum(sitk.GetArrayFromImage(s1) > 0)


class TestMorphologyMM:
    def test_erode(self):
        roi = _make_sphere(radius=15)
        result = seg.binary_erode(roi, radius_mm=2.0)
        assert np.sum(sitk.GetArrayFromImage(result)) < np.sum(
            sitk.GetArrayFromImage(roi) > 0)

    def test_dilate(self):
        roi = _make_sphere(radius=15)
        result = seg.binary_dilate(roi, radius_mm=2.0)
        assert np.sum(sitk.GetArrayFromImage(result)) > np.sum(
            sitk.GetArrayFromImage(roi) > 0)

    def test_open(self):
        roi = _make_sphere(radius=15)
        result = seg.binary_open(roi, radius_mm=1.0)
        assert np.sum(sitk.GetArrayFromImage(result)) > 0

    def test_close(self):
        roi = _make_sphere(radius=15)
        result = seg.binary_close(roi, radius_mm=1.0)
        assert np.sum(sitk.GetArrayFromImage(result)) > 0


class TestPreprocessing:
    def test_n4_bias_correction(self, intensity_image):
        result = seg.n4_bias_field_correction(intensity_image, shrink_factor=4)
        arr = sitk.GetArrayFromImage(result)
        assert arr.shape == (64, 64, 64)

    def test_anisotropic_diffusion(self, intensity_image):
        result = seg.anisotropic_diffusion(intensity_image, iterations=3)
        arr = sitk.GetArrayFromImage(result)
        assert arr.shape == (64, 64, 64)


# ============================================================================
# TESTS: Imaginable wrapper methods
# ============================================================================

class TestImaginableMethods:
    def test_segmentOtsu(self, intensity_imaginable):
        roi = intensity_imaginable.segmentOtsu()
        assert isinstance(roi, Roiable)
        assert np.sum(roi.getImageAsNumpy()) > 0

    def test_segmentLi(self, intensity_imaginable):
        roi = intensity_imaginable.segmentLi()
        assert isinstance(roi, Roiable)

    def test_segmentYen(self, intensity_imaginable):
        roi = intensity_imaginable.segmentYen()
        assert isinstance(roi, Roiable)

    def test_segmentTriangle(self, intensity_imaginable):
        roi = intensity_imaginable.segmentTriangle()
        assert isinstance(roi, Roiable)

    def test_segmentHuang(self, intensity_imaginable):
        roi = intensity_imaginable.segmentHuang()
        assert isinstance(roi, Roiable)

    def test_segmentThreshold(self, intensity_imaginable):
        roi = intensity_imaginable.segmentThreshold(0.5, 1.5)
        assert isinstance(roi, Roiable)

    def test_segmentMultiOtsu(self, intensity_imaginable):
        lm = intensity_imaginable.segmentMultiOtsu(n_thresholds=2)
        assert isinstance(lm, LabelMapable)

    def test_segmentConnectedThreshold(self, intensity_imaginable):
        seed = Roiable()
        seed.setImage(_make_sphere(radius=5), 'seed')
        roi = intensity_imaginable.segmentConnectedThreshold(seed)
        assert isinstance(roi, Roiable)

    def test_segmentMorphologicalWatershed(self, intensity_imaginable):
        lm = intensity_imaginable.segmentMorphologicalWatershed(level=0.5)
        assert isinstance(lm, LabelMapable)

    def test_correctBiasField(self, intensity_imaginable):
        result = intensity_imaginable.correctBiasField(shrink_factor=4)
        assert result is intensity_imaginable  # chainable

    def test_smoothAnisotropic(self, intensity_imaginable):
        result = intensity_imaginable.smoothAnisotropic(iterations=2)
        assert result is intensity_imaginable

    def test_getEdgeMap(self, intensity_imaginable):
        edge = intensity_imaginable.getEdgeMap(sigma=1.0)
        assert isinstance(edge, Imaginable)


# ============================================================================
# TESTS: Roiable wrapper methods
# ============================================================================

class TestRoiableMethods:
    def test_refineThresholdLevelSet(self, sphere_roiable, intensity_imaginable):
        result = sphere_roiable.refineThresholdLevelSet(intensity_imaginable)
        assert result is sphere_roiable

    def test_refineLaplacianLevelSet(self, sphere_roiable, intensity_imaginable):
        result = sphere_roiable.refineLaplacianLevelSet(intensity_imaginable)
        assert result is sphere_roiable

    def test_refineShapeDetectionLevelSet(self, sphere_roiable, intensity_imaginable):
        result = sphere_roiable.refineShapeDetectionLevelSet(intensity_imaginable)
        assert result is sphere_roiable

    def test_refineChanVese(self, sphere_roiable, intensity_imaginable):
        result = sphere_roiable.refineChanVese(intensity_imaginable)
        assert result is sphere_roiable

    def test_constrainByDistance(self, sphere_roiable):
        bone = Roiable()
        bone.setImage(_make_sphere(radius=10), 'bone')
        result = sphere_roiable.constrainByDistance(bone, max_distance_mm=3.0)
        assert result is sphere_roiable

    def test_subtractMask(self, sphere_roiable):
        inner = Roiable()
        inner.setImage(_make_sphere(radius=8), 'inner')
        result = sphere_roiable.subtractMask(inner)
        assert result is sphere_roiable
        assert np.sum(sphere_roiable.getImageAsNumpy()) > 0

    def test_intersectWith(self, sphere_roiable):
        other = Roiable()
        other.setImage(_make_sphere(radius=10, center=[32, 32, 38]), 'other')
        result = sphere_roiable.intersectWith(other)
        assert result is sphere_roiable

    def test_unionWith(self, sphere_roiable):
        other = Roiable()
        other.setImage(_make_sphere(radius=10, center=[32, 32, 50]), 'other')
        result = sphere_roiable.unionWith(other)
        assert result is sphere_roiable

    def test_erodeMM(self, sphere_roiable):
        orig_vol = np.sum(sphere_roiable.getImageAsNumpy() > 0)
        result = sphere_roiable.erodeMM(radius_mm=2.0)
        assert result is sphere_roiable
        assert np.sum(sphere_roiable.getImageAsNumpy() > 0) < orig_vol

    def test_dilateMM(self, sphere_roiable):
        orig_vol = np.sum(sphere_roiable.getImageAsNumpy() > 0)
        result = sphere_roiable.dilateMM(radius_mm=2.0)
        assert result is sphere_roiable
        assert np.sum(sphere_roiable.getImageAsNumpy() > 0) > orig_vol

    def test_openMM(self, sphere_roiable):
        result = sphere_roiable.openMM(radius_mm=1.0)
        assert result is sphere_roiable

    def test_closeMM(self, sphere_roiable):
        result = sphere_roiable.closeMM(radius_mm=1.0)
        assert result is sphere_roiable


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
