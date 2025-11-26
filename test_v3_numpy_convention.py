"""
Test script for PyAble v3 - verifies the new numpy (Z,Y,X) convention
"""
import numpy as np
import SimpleITK as sitk
import sys
sys.path.insert(0, '/home/erosm/pyable')

from pyable_eros_montin import SITKImaginable

print("=" * 70)
print("PyAble v3 - Testing New (Z,Y,X) Convention")
print("=" * 70 + "\n")

# Test 1: getImageAsNumpy now returns (Z,Y,X)
print("Test 1: getImageAsNumpy() returns (Z,Y,X) ordering...")
arr_zyx = np.random.rand(10, 20, 30).astype(np.float32)
img_sitk = sitk.GetImageFromArray(arr_zyx)
img = SITKImaginable(image=img_sitk)

arr = img.getImageAsNumpy()
assert arr.shape == (10, 20, 30), f"Expected (10,20,30), got {arr.shape}"
assert np.allclose(arr, arr_zyx), "Content mismatch"
print(f"  ✓ getImageAsNumpy() shape: {arr.shape} (Z,Y,X)")

# Test 2: All ZYX methods are consistent
arr_zyx2 = img.getImageAsNumpyZYX()
arr_pytorch = img.getImageAsNumpyForPyTorch()
assert np.array_equal(arr, arr_zyx2), "getImageAsNumpyZYX() mismatch"
assert np.array_equal(arr, arr_pytorch), "getImageAsNumpyForPyTorch() mismatch"
print("  ✓ All ZYX methods return same array")

# Test 3: XYZ deprecated method
arr_xyz = img.getImageAsNumpyXYZ()
assert arr_xyz.shape == (30, 20, 10), f"XYZ shape should be (30,20,10), got {arr_xyz.shape}"
print(f"  ✓ getImageAsNumpyXYZ() shape: {arr_xyz.shape} (X,Y,Z) - deprecated")

# Test 4: Round-trip setImageFromNumpy/getImageAsNumpy
print("\nTest 2: Round-trip numpy conversion...")
arr_test = np.random.rand(8, 12, 16).astype(np.float32)
img2 = SITKImaginable()
img2.setImageFromNumpy(arr_test)
arr_back = img2.getImageAsNumpy()
assert arr_back.shape == (8, 12, 16), f"Shape mismatch: {arr_back.shape}"
assert np.allclose(arr_test, arr_back), "Content mismatch"
print(f"  ✓ setImageFromNumpy/getImageAsNumpy round-trip works")
print(f"  ✓ Shape preserved: {arr_test.shape} -> {arr_back.shape}")

# Test 5: getBoundingBox returns (k,j,i)
print("\nTest 3: getBoundingBox() returns (k,j,i) indices...")
arr_bbox = np.zeros((20, 30, 40), dtype=np.float32)
arr_bbox[5:15, 10:20, 15:35] = 1.0
img3 = SITKImaginable()
img_sitk3 = sitk.GetImageFromArray(arr_bbox)
img3.setImage(img_sitk3)

bbox = img3.getBoundingBox(exclude=[0])
if bbox is not None:
    (k_min, j_min, i_min), (k_max, j_max, i_max) = bbox
    print(f"  ✓ Bounding box: k=[{k_min},{k_max}], j=[{j_min},{j_max}], i=[{i_min},{i_max}]")
    assert k_min == 5 and k_max == 14, f"Z bounds wrong"
    assert j_min == 10 and j_max == 19, f"Y bounds wrong"  
    assert i_min == 15 and i_max == 34, f"X bounds wrong"
    print("  ✓ Indices match expected (k,j,i) = (z,y,x) ordering")

# Test 6: Coordinate conversions
print("\nTest 4: Coordinate conversion methods...")
arr4 = np.ones((10, 20, 30), dtype=np.float32)
img_sitk4 = sitk.GetImageFromArray(arr4)
img_sitk4.SetSpacing([2.0, 1.5, 1.0])  # x, y, z spacing
img_sitk4.SetOrigin([10.0, 20.0, 30.0])  # x, y, z origin
img4 = SITKImaginable(image=img_sitk4)

kji = (5, 10, 15)  # array index (z, y, x)
xyz = img4.getPhysicalPointFromArrayIndex(kji)
expected_xyz = (40.0, 35.0, 35.0)  # 10+15*2, 20+10*1.5, 30+5*1
assert np.allclose(xyz, expected_xyz, atol=0.01), f"Expected {expected_xyz}, got {xyz}"
print(f"  ✓ Array index {kji} → physical {xyz} mm")

kji_back = img4.getArrayIndexFromPhysicalPoint(expected_xyz)
assert kji_back == (5, 10, 15), f"Expected (5,10,15), got {kji_back}"
print(f"  ✓ Physical {expected_xyz} mm → array index {kji_back}")

print("\n" + "=" * 70)
print("✓ ALL TESTS PASSED! PyAble v3 is using correct (Z,Y,X) convention")
print("=" * 70)
