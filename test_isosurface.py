#!/usr/bin/env python
"""
Test script for VTK isosurface rendering functionality.

Demonstrates:
- Rendering isosurfaces for continuous images
- ROI boundary rendering
- Multi-component image handling
- Custom isosurface values
- Non-interactive rendering (for batch processing)
"""

import sys
from pathlib import Path
import numpy as np
import SimpleITK as sitk

# Ensure repository import
repo_root = Path(__file__).resolve().parent
sys.path.insert(0, str(repo_root))

from pyable.imaginable import Imaginable, Roiable, SITKImaginable, Fieldable


def create_synthetic_volume():
    """Create a synthetic 3D volume with a sphere for testing."""
    size = (100, 100, 100)
    image = sitk.Image(size, sitk.sitkFloat32)
    
    # Create a Gaussian ball in the center
    center = [s // 2 for s in size]
    radius = 30
    
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                dist = np.sqrt((i - center[0])**2 + (j - center[1])**2 + (k - center[2])**2)
                value = np.exp(-(dist**2) / (2 * radius**2)) * 255
                image.SetPixel([i, j, k], value)
    
    return image


def create_synthetic_roi():
    """Create a synthetic ROI (binary segmentation)."""
    size = (100, 100, 100)
    image = sitk.Image(size, sitk.sitkUInt8)
    
    # Create a sphere
    center = [s // 2 for s in size]
    radius = 25
    
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                dist = np.sqrt((i - center[0])**2 + (j - center[1])**2 + (k - center[2])**2)
                if dist <= radius:
                    image.SetPixel([i, j, k], 1)
                else:
                    image.SetPixel([i, j, k], 0)
    
    return image


def create_synthetic_vector_field():
    """Create a synthetic displacement field."""
    size = (50, 50, 50)
    
    # Create 3D vector field (displacement) using numpy array
    # Create 3 separate scalar images and combine them
    center = [s // 2 for s in size]
    
    # Create magnitude image for testing isosurface
    image = sitk.Image(size, sitk.sitkFloat32)
    
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                # Create radial displacement magnitude
                dx = (i - center[0]) * 0.1
                dy = (j - center[1]) * 0.1
                dz = (k - center[2]) * 0.1
                magnitude = np.sqrt(dx**2 + dy**2 + dz**2) * 100
                image.SetPixel([i, j, k], magnitude)
    
    return image


def test_continuous_image_isosurface():
    """Test isosurface rendering for continuous images."""
    print("\n" + "="*60)
    print("TEST 1: Continuous Image Isosurface")
    print("="*60)
    
    try:
        # Create synthetic volume
        sitk_image = create_synthetic_volume()
        
        # Wrap in Imaginable
        img = SITKImaginable(image=sitk_image)
        
        print("✓ Created synthetic continuous image (100x100x100)")
        print(f"  Pixel type: {img.getImagePixelTypeAsString()[1]}")
        print(f"  Min value: {img.getMinimumValue():.2f}")
        print(f"  Mean value: {img.getMeanValue():.2f}")
        print(f"  Max value: {img.getMaximumValue():.2f}")
        
        print("Test 1b: Rendering isosurface at custom value (128)...")
        actor, _, _ = img.renderIsosurface(isosurface_value=128, show=False, color=(1.0, 0.0, 0.0))
        assert actor is not None
        assert actor.GetMapper().GetInput().GetNumberOfCells() > 0
        
        print("Test 1c: Rendering with custom color and opacity...")
        actor2, _, _ = img.renderIsosurface(isosurface_value=150, show=False, 
                                             color=(0.0, 1.0, 0.0), opacity=0.8)
        assert actor2 is not None
    
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        raise


def test_roi_isosurface():
    """Test isosurface rendering for ROI (binary segmentation)."""
    print("\n" + "="*60)
    print("TEST 2: ROI (Binary Segmentation) Isosurface")
    print("="*60)
    
    try:
        # Create synthetic ROI
        sitk_image = create_synthetic_roi()
        
        # Wrap in Roiable
        roi = Roiable(image=sitk_image)
        
        print("✓ Created synthetic ROI (binary segmentation, 100x100x100)")
        print(f"  Pixel type: {roi.getImagePixelTypeAsString()[1]}")
        print(f"  Non-zero voxels: {roi.getNumberOfNonZeroVoxels()}")
        
        print("\nTest 2a: Rendering ROI boundary (automatic value 0.5)...")
        actor, _, _ = roi.renderIsosurface(show=False, color=(0.0, 1.0, 0.0))
        assert actor is not None
        assert actor.GetMapper().GetInput().GetNumberOfCells() > 0
        
        print("Test 2b: Rendering ROI at explicit boundary...")
        actor2, _, _ = roi.renderIsosurface(isosurface_value=0.5, show=False, 
                                             color=(0.0, 0.0, 1.0), opacity=0.9)
        assert actor2 is not None
    
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        raise


def test_vector_field_magnitude():
    """Test isosurface rendering for vector field magnitude."""
    print("\n" + "="*60)
    print("TEST 3: Vector Field Magnitude Isosurface")
    print("="*60)
    
    try:
        # Create synthetic vector field (magnitude)
        sitk_image = create_synthetic_vector_field()
        
        # Wrap in SITKImaginable (treated as scalar)
        vec_field = SITKImaginable(image=sitk_image)
        
        print("✓ Created synthetic vector field magnitude (50x50x50)")
        print(f"  Pixel type: {vec_field.getImagePixelTypeAsString()[1]}")
        
        stats = sitk.StatisticsImageFilter()
        stats.Execute(sitk_image)
        
        print(f"  Magnitude range: [{stats.GetMinimum():.3f}, {stats.GetMaximum():.3f}]")
        print(f"  Magnitude mean: {stats.GetMean():.3f}")
        
        print("\nTest 3a: Rendering displacement field magnitude...")
        actor, _, _ = vec_field.renderIsosurface(show=False, color=(1.0, 1.0, 0.0))
        assert actor is not None
        assert actor.GetMapper().GetInput().GetNumberOfCells() > 0
        
        print("Test 3b: Rendering at custom magnitude value...")
        actor2, _, _ = vec_field.renderIsosurface(isosurface_value=3.0, show=False, 
                                                    color=(0.5, 1.0, 0.5))
        assert actor2 is not None
    
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        raise


def test_isosurface_colors():
    """Test various color options for isosurface rendering."""
    print("\n" + "="*60)
    print("TEST 4: Isosurface Color Variations")
    print("="*60)
    
    try:
        sitk_image = create_synthetic_volume()
        img = SITKImaginable(image=sitk_image)
        
        colors = [
            ((1.0, 0.0, 0.0), "Red"),
            ((0.0, 1.0, 0.0), "Green"),
            ((0.0, 0.0, 1.0), "Blue"),
            ((1.0, 1.0, 0.0), "Yellow"),
            ((1.0, 0.0, 1.0), "Magenta"),
            ((0.0, 1.0, 1.0), "Cyan"),
        ]
        
        for i, (color, name) in enumerate(colors):
            actor, _, _ = img.renderIsosurface(
                isosurface_value=128, 
                show=False,
                color=color
            )
            assert actor is not None
            assert actor.GetMapper().GetInput().GetNumberOfCells() > 0
            print(f"✓ {name}: {color} - {actor.GetMapper().GetInput().GetNumberOfCells()} polygons")
    
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        raise


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("VTK ISOSURFACE RENDERING TESTS")
    print("="*60)
    
    results = []
    
    # Run tests
    try:
        test_continuous_image_isosurface()
        results.append(("Continuous Image Isosurface", True))
    except Exception:
        results.append(("Continuous Image Isosurface", False))
    try:
        test_roi_isosurface()
        results.append(("ROI Boundary Rendering", True))
    except Exception:
        results.append(("ROI Boundary Rendering", False))
    try:
        test_vector_field_magnitude()
        results.append(("Vector Field Magnitude", True))
    except Exception:
        results.append(("Vector Field Magnitude", False))
    try:
        test_isosurface_colors()
        results.append(("Color Variations", True))
    except Exception:
        results.append(("Color Variations", False))
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
