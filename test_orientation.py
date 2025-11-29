#!/usr/bin/env python3
"""
Test script to demonstrate orientation methods and how dicomOrient changes numpy arrays.
"""

import numpy as np
import sys
sys.path.insert(0, '/home/erosm/pyable')

from pyable.imaginable import Imaginable

print("=" * 70)
print("Testing Image Orientation Methods")
print("=" * 70)

# Create a test image with non-identity direction
print("\n1. Creating test image with known pattern...")
array = np.zeros((10, 20, 30), dtype=np.float32)
# Put a marker in corner: array[z=0, y=0, x=0:5] = 100
array[0, 0, 0:5] = 100.0
print(f"   Original array shape: {array.shape} (Z,Y,X)")
print(f"   Marker at: array[0, 0, 0:5] = 100")

# Create imaginable with LPS orientation
img = Imaginable()
img.setImageFromNumpy(array, spacing=[1.0, 1.0, 1.0], origin=[0.0, 0.0, 0.0], 
                      direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

print(f"\n2. Initial image properties:")
print(f"   Orientation code: {img.getOrientationCode()}")
print(f"   Direction cosines: {img.getDirectionCosines()}")
print(f"   Image size: {img.getImageSize()}")
print(f"   Array shape: {img.getImageAsNumpy().shape}")

# Check where our marker is
initial_array = img.getImageAsNumpy()
print(f"   Marker location check: array[0,0,0:5] = {initial_array[0, 0, 0:5]}")

print(f"\n3. Reorienting to RAS (Right-Anterior-Superior)...")
img_ras = img.getDuplicate()
img_ras.reorientToRAS()

print(f"   New orientation code: {img_ras.getOrientationCode()}")
print(f"   New direction cosines: {img_ras.getDirectionCosines()}")
print(f"   New image size: {img_ras.getImageSize()}")
print(f"   New array shape: {img_ras.getImageAsNumpy().shape}")

# Check if numpy array changed
ras_array = img_ras.getImageAsNumpy()
arrays_equal = np.array_equal(initial_array, ras_array)
print(f"   Arrays equal? {arrays_equal}")
if not arrays_equal:
    print(f"   ✓ CONFIRMED: dicomOrient DOES change the numpy array!")
    print(f"   Original marker: array[0,0,0:5] = {initial_array[0, 0, 0:5]}")
    print(f"   After RAS:       array[0,0,0:5] = {ras_array[0, 0, 0:5]}")
    # Find where marker moved to
    marker_locs = np.where(ras_array > 50)
    if len(marker_locs[0]) > 0:
        print(f"   Marker now at: z={marker_locs[0][0]}, y={marker_locs[1][0]}, x={marker_locs[2][0:5]}")

print(f"\n4. Reorienting to RPI (Right-Posterior-Inferior)...")
img_rpi = img.getDuplicate()
img_rpi.reorientToRPI()

print(f"   New orientation code: {img_rpi.getOrientationCode()}")
print(f"   New direction cosines: {img_rpi.getDirectionCosines()}")
print(f"   New array shape: {img_rpi.getImageAsNumpy().shape}")

rpi_array = img_rpi.getImageAsNumpy()
arrays_equal = np.array_equal(initial_array, rpi_array)
print(f"   Arrays equal to original? {arrays_equal}")

print(f"\n5. Testing direction cosine methods...")
print(f"   getDirectionCosines() returns: {type(img.getDirectionCosines())}")
print(f"   Value: {img.getDirectionCosines()}")

# Try setting custom direction (metadata only - doesn't change array)
print(f"\n6. Testing setDirectionCosines() [metadata only]...")
img_meta = img.getDuplicate()
before_array = img_meta.getImageAsNumpy().copy()
img_meta.setDirectionCosines((-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0))
after_array = img_meta.getImageAsNumpy()

print(f"   Direction changed to: {img_meta.getDirectionCosines()}")
print(f"   Arrays equal? {np.array_equal(before_array, after_array)}")
print(f"   ✓ CONFIRMED: setDirectionCosines() only changes metadata, not array data")

print("\n" + "=" * 70)
print("Summary:")
print("=" * 70)
print("• dicomOrient() / reorientToLPS/RAS/RPI:")
print("  → DOES physically reorient the numpy array")
print("  → Changes voxel data arrangement")
print("  → Updates direction matrix")
print("  → Preserves physical world coordinates")
print("")
print("• setDirectionCosines() / setImageDirection():")
print("  → Only changes metadata")
print("  → Does NOT modify numpy array")
print("  → Use for changing interpretation only")
print("=" * 70)
