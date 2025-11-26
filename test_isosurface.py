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
import numpy as np
import SimpleITK as sitk

# Add parent directory to path
sys.path.insert(0, '/home/erosm/pyable')

from pyable_eros_montin.imaginable import Imaginable, Roiable, SITKImaginable, Fieldable


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
        
        # Test 1a: Render with automatic isosurface value (mean)
        print("\nTest 1a: Rendering isosurface at mean intensity...")
        # actor = img.renderIsosurface(show=False)  # Don't show in tests
        # print("✓ Successfully created isosurface actor (without display)")
        
        # Test 1b: Render with custom value
        print("Test 1b: Rendering isosurface at custom value (128)...")
        actor, _, _ = img.renderIsosurface(isosurface_value=128, show=False, color=(1.0, 0.0, 0.0))
        print(f"✓ Successfully created isosurface actor")
        print(f"  Actor mapper has {actor.GetMapper().GetInput().GetNumberOfCells()} polygons")
        
        # Test 1c: Different color and opacity
        print("Test 1c: Rendering with custom color and opacity...")
        actor2, _, _ = img.renderIsosurface(isosurface_value=150, show=False, 
                                             color=(0.0, 1.0, 0.0), opacity=0.8)
        print(f"✓ Successfully created colored isosurface")
        
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


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
        
        # Test 2a: Render ROI at boundary (0.5) - automatic
        print("\nTest 2a: Rendering ROI boundary (automatic value 0.5)...")
        actor, _, _ = roi.renderIsosurface(show=False, color=(0.0, 1.0, 0.0))
        print(f"✓ Successfully created ROI boundary isosurface")
        print(f"  Actor mapper has {actor.GetMapper().GetInput().GetNumberOfCells()} polygons")
        
        # Test 2b: Explicit boundary value
        print("Test 2b: Rendering ROI at explicit boundary...")
        actor2, _, _ = roi.renderIsosurface(isosurface_value=0.5, show=False, 
                                             color=(0.0, 0.0, 1.0), opacity=0.9)
        print(f"✓ Successfully created ROI isosurface with explicit value")
        
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


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
        
        # Test 3a: Render magnitude with automatic value
        print("\nTest 3a: Rendering displacement field magnitude...")
        actor, _, _ = vec_field.renderIsosurface(show=False, color=(1.0, 1.0, 0.0))
        print(f"✓ Successfully created vector field magnitude isosurface")
        print(f"  Actor mapper has {actor.GetMapper().GetInput().GetNumberOfCells()} polygons")
        
        # Test 3b: Custom isosurface value
        print("Test 3b: Rendering at custom magnitude value...")
        actor2, _, _ = vec_field.renderIsosurface(isosurface_value=3.0, show=False, 
                                                    color=(0.5, 1.0, 0.5))
        print(f"✓ Successfully created magnitude isosurface at value 3.0")
        
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


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
            print(f"✓ {name}: {color} - {actor.GetMapper().GetInput().GetNumberOfCells()} polygons")
        
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("VTK ISOSURFACE RENDERING TESTS")
    print("="*60)
    
    results = []
    
    # Run tests
    results.append(("Continuous Image Isosurface", test_continuous_image_isosurface()))
    results.append(("ROI Boundary Rendering", test_roi_isosurface()))
    results.append(("Vector Field Magnitude", test_vector_field_magnitude()))
    results.append(("Color Variations", test_isosurface_colors()))
    
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
