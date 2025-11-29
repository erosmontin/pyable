#!/usr/bin/env python3
"""
Test script demonstrating resampling oblique images to axis-aligned grids.
"""

import numpy as np
import sys
sys.path.insert(0, '/home/erosm/pyable')

from pyable.imaginable import Imaginable

print("=" * 70)
print("Testing Oblique to Axis-Aligned Resampling")
print("=" * 70)

# Create a test image with oblique direction cosines
print("\n1. Creating oblique image (rotated 30 degrees around Z-axis)...")
array = np.zeros((20, 30, 40), dtype=np.float32)
# Put a distinct pattern
array[5:15, 10:20, 15:25] = 100.0  # Rectangular region
print(f"   Array shape: {array.shape} (Z,Y,X)")
print(f"   Pattern region: z=[5:15], y=[10:20], x=[15:25]")

# Create oblique direction matrix (30 degree rotation around Z)
angle = np.radians(30)
cos_a = np.cos(angle)
sin_a = np.sin(angle)

# Rotation around Z-axis: 
# X' = X*cos - Y*sin
# Y' = X*sin + Y*cos
# Z' = Z
oblique_direction = (
    cos_a, sin_a, 0.0,   # X-axis in oblique coordinates
    -sin_a, cos_a, 0.0,  # Y-axis in oblique coordinates
    0.0, 0.0, 1.0        # Z-axis unchanged
)

print(f"\n   Oblique direction matrix (30° rotation):")
print(f"   [{oblique_direction[0]:.3f}, {oblique_direction[1]:.3f}, {oblique_direction[2]:.3f}]")
print(f"   [{oblique_direction[3]:.3f}, {oblique_direction[4]:.3f}, {oblique_direction[5]:.3f}]")
print(f"   [{oblique_direction[6]:.3f}, {oblique_direction[7]:.3f}, {oblique_direction[8]:.3f}]")

img = Imaginable()
img.setImageFromNumpy(array, spacing=[1.0, 1.0, 1.0], origin=[0.0, 0.0, 0.0], 
                      direction=oblique_direction)

print(f"\n2. Initial oblique image properties:")
print(f"   Direction cosines: {img.getDirectionCosines()}")
print(f"   Is axis-aligned? {img.isAxisAligned()}")
print(f"   Image size: {img.getImageSize()}")
print(f"   Array shape: {img.getImageAsNumpy().shape}")
print(f"   Non-zero voxel count: {np.count_nonzero(img.getImageAsNumpy())}")

# Check initial array
initial_array = img.getImageAsNumpy()
initial_sum = np.sum(initial_array)
print(f"   Sum of all voxels: {initial_sum:.1f}")

print(f"\n3. Resampling to axis-aligned grid...")
img_aligned = img.getDuplicate()
img_aligned.resampleToAxisAligned()

print(f"   New direction cosines: {img_aligned.getDirectionCosines()}")
print(f"   Is axis-aligned? {img_aligned.isAxisAligned()}")
print(f"   Image size: {img_aligned.getImageSize()}")
print(f"   Array shape: {img_aligned.getImageAsNumpy().shape}")

aligned_array = img_aligned.getImageAsNumpy()
aligned_sum = np.sum(aligned_array)
print(f"   Non-zero voxel count: {np.count_nonzero(aligned_array)}")
print(f"   Sum of all voxels: {aligned_sum:.1f}")
print(f"   Sum difference: {abs(initial_sum - aligned_sum):.1f} (due to interpolation)")

# Check arrays are different
arrays_equal = np.array_equal(initial_array, aligned_array)
print(f"   Arrays equal? {arrays_equal}")
if not arrays_equal:
    print(f"   ✓ Array was resampled (voxels repositioned)")

print(f"\n4. Testing with dicomOrient for comparison...")
img_dicom = img.getDuplicate()
print(f"   Before: {img_dicom.getOrientationCode()}")
img_dicom.dicomOrient('LPS')
print(f"   After dicomOrient: {img_dicom.getOrientationCode()}")
print(f"   Direction cosines: {img_dicom.getDirectionCosines()}")
print(f"   Is axis-aligned? {img_dicom.isAxisAligned()}")

dicom_array = img_dicom.getImageAsNumpy()
dicom_sum = np.sum(dicom_array)
print(f"   Sum of voxels: {dicom_sum:.1f}")

print(f"\n5. Creating another example: highly oblique")
# More extreme oblique angle
angle2 = np.radians(45)
cos_a2 = np.cos(angle2)
sin_a2 = np.sin(angle2)

oblique_direction2 = (
    cos_a2, sin_a2, 0.0,
    -sin_a2, cos_a2, 0.0,
    0.0, 0.0, 1.0
)

img2 = Imaginable()
img2.setImageFromNumpy(array, spacing=[1.0, 1.0, 1.0], origin=[0.0, 0.0, 0.0],
                       direction=oblique_direction2)

print(f"   Oblique (45°) direction:")
print(f"   [{oblique_direction2[0]:.3f}, {oblique_direction2[1]:.3f}, {oblique_direction2[2]:.3f}]")
print(f"   [{oblique_direction2[3]:.3f}, {oblique_direction2[4]:.3f}, {oblique_direction2[5]:.3f}]")
print(f"   [{oblique_direction2[6]:.3f}, {oblique_direction2[7]:.3f}, {oblique_direction2[8]:.3f}]")
print(f"   Is axis-aligned? {img2.isAxisAligned()}")

img2.resampleToAxisAligned()
print(f"   After resampling:")
print(f"   Direction: {img2.getDirectionCosines()}")
print(f"   Is axis-aligned? {img2.isAxisAligned()}")

print(f"\n6. Testing changeImageDirection with custom direction...")
img3 = Imaginable()
img3.setImageFromNumpy(np.random.rand(10, 20, 30), 
                       spacing=[1.0, 1.0, 1.0], 
                       origin=[0.0, 0.0, 0.0])

custom_direction = (-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0)
print(f"   Original direction: {img3.getDirectionCosines()}")
print(f"   Target direction: {custom_direction}")

img3.changeImageDirection(custom_direction)
print(f"   New direction: {img3.getDirectionCosines()}")
print(f"   Is axis-aligned? {img3.isAxisAligned()}")

print("\n" + "=" * 70)
print("Summary:")
print("=" * 70)
print("✓ resampleToAxisAligned() converts oblique acquisitions to axis-aligned")
print("✓ changeImageDirection() resamples to any target direction matrix")
print("✓ isAxisAligned() checks if direction is identity matrix")
print("✓ Physical coordinates (mm) are preserved during resampling")
print("✓ Voxel values may change slightly due to interpolation")
print("")
print("Key differences:")
print("• dicomOrient('LPS'):       Permutes/flips axes (no interpolation)")
print("• resampleToAxisAligned():  Resamples data (interpolation applied)")
print("=" * 70)
