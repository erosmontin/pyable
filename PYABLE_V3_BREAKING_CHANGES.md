# PyAble v3 - Breaking Changes Summary

## Overview

PyAble v3 introduces **breaking changes** to make the library consistent with standard numpy and PyTorch conventions. The main change is that **all numpy arrays now use (Z, Y, X) ordering** by default.

## What Changed

### 1. `getImageAsNumpy()` - BREAKING CHANGE ⚠️

**v2 (old):**
```python
arr = img.getImageAsNumpy()  # Returned (X, Y, Z) - NON-STANDARD
```

**v3 (new):**
```python
arr = img.getImageAsNumpy()  # Returns (Z, Y, X) - STANDARD numpy/PyTorch format
```

- Now returns arrays in standard (Z, Y, X) ordering
- For 3D: (Z, Y, X) = (depth/slices, height/rows, width/cols)
- For 2D: (Y, X) = (height/rows, width/cols)
- **Migration:** If you need old behavior, use `getImageAsNumpyXYZ()` (deprecated)

### 2. `setImageFromNumpy()` - BREAKING CHANGE ⚠️

**v2 (old):**
```python
arr = np.random.rand(100, 200, 300)  # (X, Y, Z)
img.setImageFromNumpy(arr)
```

**v3 (new):**
```python
arr = np.random.rand(100, 200, 300)  # (Z, Y, X)
img.setImageFromNumpy(arr)
```

- Now expects arrays in standard (Z, Y, X) ordering
- **Migration:** If you have (X,Y,Z) arrays, use `setImageFromNumpyXYZ()` (deprecated)

### 3. `getBoundingBox()` - BREAKING CHANGE ⚠️

**v2 (old):**
```python
bbox = img.getBoundingBox()  # Returned ((x_min, y_min, z_min), (x_max, y_max, z_max))
```

**v3 (new):**
```python
bbox = img.getBoundingBox()  # Returns ((k_min, j_min, i_min), (k_max, j_max, i_max))
                              # where (k, j, i) = (z, y, x)
```

- Now returns indices in numpy array ordering: (k, j, i) = (z, y, x)
- Directly usable with numpy arrays: `arr[k_min:k_max+1, j_min:j_max+1, i_min:i_max+1]`

### 4. New Methods Added

#### Array Access (Non-Breaking Additions)
- `getImageAsNumpyZYX()` - Alias for `getImageAsNumpy()` with explicit name
- `getImageAsNumpyForPyTorch()` - Explicit alias for PyTorch users
- `getImageAsNumpyXYZ()` - **DEPRECATED** Returns old (X,Y,Z) format

#### Image Object Access
- `getITKImage()` - Explicit alias for `getImage()`, returns SimpleITK.Image
- `getVTKImage()` - Convert to VTK format for 3D rendering

#### Coordinate Conversions
- `getPhysicalPointFromArrayIndex(kji)` - Array index (z,y,x) → physical (x,y,z) mm
- `getArrayIndexFromPhysicalPoint(xyz)` - Physical (x,y,z) mm → array index (z,y,x)
- `getPhysicalPointFromITKIndex(ijk)` - ITK index (i,j,k) → physical (x,y,z) mm
- `getITKIndexFromPhysicalPoint(xyz)` - Physical (x,y,z) mm → ITK index (i,j,k)

#### Transformation Helpers
- `createTranslationMM(tx_mm, ty_mm, tz_mm)` - Build translation vectors
- `createRotationDegrees(rx_deg, ry_deg, rz_deg)` - Build rotation vectors
- `createAffineMatrix2D(angle_deg, scale_x, scale_y, shear)` - Build 2D affine matrices
- `createAffineMatrix3D(rotation, scale)` - Build 3D affine matrices

#### Set from Numpy (Non-Breaking Additions)
- `setImageFromNumpyZYX()` - Alias for `setImageFromNumpy()`
- `setImageFromNumpyXYZ()` - **DEPRECATED** Accepts old (X,Y,Z) format

### 5. Enhanced Documentation

All transformation methods now have comprehensive docstrings clarifying:
- Coordinate systems (physical space in mm vs array indices)
- Axis definitions (X=left-right, Y=anterior-posterior, Z=superior-inferior)
- Examples with explicit axis names

Methods updated:
- `translateImage()` - Documented that T is [tx_mm, ty_mm, tz_mm] in physical space
- `rotateImage()` - Documented rotation axes and degrees
- `transformImageAffine()` - Documented that matrix operates in physical space
- `transformFromRegistration()` - Clarified registration transform usage

## Migration Guide

### Quick Fix for Existing Code

If your code breaks, the quickest fix is to use the deprecated methods:

**Replace:**
```python
arr = img.getImageAsNumpy()  # Was (X,Y,Z) in v2
```

**With:**
```python
arr = img.getImageAsNumpyXYZ()  # Still (X,Y,Z) but deprecated
```

**Replace:**
```python
img.setImageFromNumpy(arr_xyz)  # Was expecting (X,Y,Z) in v2
```

**With:**
```python
img.setImageFromNumpyXYZ(arr_xyz)  # Still accepts (X,Y,Z) but deprecated
```

### Proper Migration (Recommended)

**1. Update array creation/usage:**

```python
# OLD v2:
arr = img.getImageAsNumpy()  # (X, Y, Z)
value = arr[x, y, z]  # Confusing!

# NEW v3:
arr = img.getImageAsNumpy()  # (Z, Y, X)  
k, j, i = z, y, x  # Clear mapping
value = arr[k, j, i]  # Intuitive!
```

**2. Update array setting:**

```python
# OLD v2:
arr_xyz = np.random.rand(100, 200, 300)  # X, Y, Z
img.setImageFromNumpy(arr_xyz)

# NEW v3:
arr_zyx = np.random.rand(300, 200, 100)  # Z, Y, X
img.setImageFromNumpy(arr_zyx)

# Or transpose if you have XYZ data:
arr_xyz = np.random.rand(100, 200, 300)
arr_zyx = np.transpose(arr_xyz, (2, 1, 0))
img.setImageFromNumpy(arr_zyx)
```

**3. Update bounding box usage:**

```python
# OLD v2:
(x_min, y_min, z_min), (x_max, y_max, z_max) = img.getBoundingBox()

# NEW v3:
(k_min, j_min, i_min), (k_max, j_max, i_max) = img.getBoundingBox()
# Where k=z, j=y, i=x
arr = img.getImageAsNumpy()
cropped = arr[k_min:k_max+1, j_min:j_max+1, i_min:i_max+1]
```

**4. Use coordinate conversion methods:**

```python
# NEW v3 - explicit coordinate conversions
arr = img.getImageAsNumpy()  # (Z, Y, X)
k, j, i = 50, 60, 70  # Array indices

# Convert to physical space
xyz_mm = img.getPhysicalPointFromArrayIndex((k, j, i))
print(f"Voxel ({k},{j},{i}) is at {xyz_mm} mm")

# Convert from physical space
xyz_mm = (100.0, 150.0, 80.0)
k, j, i = img.getArrayIndexFromPhysicalPoint(xyz_mm)
value = arr[k, j, i]
```

## PyTorch Integration

**v3 makes PyTorch integration seamless:**

```python
import torch
from pyable_eros_montin import SITKImaginable

# Load image
img = SITKImaginable('scan.nii.gz')

# Get as numpy - now in correct format!
arr = img.getImageAsNumpy()  # (Z, Y, X) - ready for PyTorch!
# OR explicitly:
arr = img.getImageAsNumpyForPyTorch()

# Convert to PyTorch tensor
tensor = torch.from_numpy(arr).float()  # (D, H, W)
tensor = tensor.unsqueeze(0).unsqueeze(0)  # (1, 1, D, H, W)

# Process with model
output = model(tensor)

# Convert back
output_np = output.squeeze().cpu().numpy()  # (Z, Y, X)

# Create result image
result = SITKImaginable()
result.setImageFromNumpy(output_np)  # Expects (Z,Y,X) - perfect!
result_itk = result.getITKImage()
result_itk.CopyInformation(img.getITKImage())
result.write('output.nii.gz')
```

## Coordinate Systems Reference

### Three Coordinate Systems

**1. Array/Numpy Space**
- Ordering: (k, j, i) = (z, y, x)
- Units: Integer indices starting from 0
- Used by: numpy arrays, PyTorch tensors
- Methods: `getImageAsNumpy()`, `getBoundingBox()`

**2. Physical/World Space**
- Ordering: (x, y, z)
- Units: Millimeters (float)
- Used by: All transformations, physical measurements
- Methods: `translateImage()`, `rotateImage()`, `getPhysicalPointFromArrayIndex()`

**3. ITK Index Space**
- Ordering: (i, j, k) = (x, y, z)
- Units: Integer indices starting from 0
- Used by: SimpleITK internal representation
- Methods: `getCoordinatesFromIndex()` (deprecated), `getITKIndexFromPhysicalPoint()`

### Conversion Table

| From → To | Method |
|-----------|--------|
| Array (k,j,i) → Physical (x,y,z) | `getPhysicalPointFromArrayIndex()` |
| Physical (x,y,z) → Array (k,j,i) | `getArrayIndexFromPhysicalPoint()` |
| ITK (i,j,k) → Physical (x,y,z) | `getPhysicalPointFromITKIndex()` |
| Physical (x,y,z) → ITK (i,j,k) | `getITKIndexFromPhysicalPoint()` |

## Testing Your Migration

Run this test to verify your code works with v3:

```python
import numpy as np
from pyable_eros_montin import SITKImaginable

# Create test image
arr_zyx = np.random.rand(10, 20, 30).astype(np.float32)
img = SITKImaginable()
img.setImageFromNumpy(arr_zyx)

# Verify round-trip
arr_back = img.getImageAsNumpy()
assert arr_back.shape == (10, 20, 30), "Shape should be (Z,Y,X)"
assert np.allclose(arr_zyx, arr_back), "Content should match"

print("✓ Your code is compatible with PyAble v3!")
```

## Benefits of v3

1. **Standard Conventions**: Follows numpy/PyTorch standards - no more confusion
2. **Clearer Code**: Array indexing is intuitive: `arr[z, y, x]`
3. **PyTorch Ready**: Direct compatibility with medical imaging PyTorch models
4. **Better Documentation**: All methods clearly document coordinate systems
5. **Explicit Conversions**: New methods make coordinate conversions obvious
6. **Future Proof**: Sets foundation for better deep learning integration

## Orientation Methods - Important Clarification

### Understanding Image Orientation in PyAble v3

PyAble v3 adds clear methods for handling anatomical orientations (LPS, RAS, RPI, etc.). **Important distinction:**

#### Methods that CHANGE the numpy array:

```python
# These physically reorient the image data
img.dicomOrient('LPS')      # Reorient to LPS (DICOM standard)
img.reorientToLPS()         # Convenience: same as dicomOrient('LPS')
img.reorientToRAS()         # Reorient to RAS (NIfTI/neuroimaging)
img.reorientToRPI()         # Reorient to RPI
```

**What happens:**
- ✓ Numpy array IS physically rearranged (voxels move)
- ✓ Direction cosines are updated (usually to near-identity)
- ✓ Physical world coordinates preserved (same anatomy, same locations)
- ✓ Image size may change if axes are permuted

#### Methods that ONLY change metadata:

```python
# These only change interpretation, NOT the data
img.setDirectionCosines((-1, 0, 0, 0, 1, 0, 0, 0, 1))
img.setImageDirection((-1, 0, 0, 0, 1, 0, 0, 0, 1))
```

**What happens:**
- ✓ Direction matrix updated
- ✗ Numpy array stays IDENTICAL (no data movement)
- ✗ Only changes how voxel indices map to world coordinates

#### Query current orientation:

```python
# Get current anatomical orientation
orientation = img.getOrientationCode()      # Returns 'LPS', 'RAS', etc.
direction = img.getDirectionCosines()       # Returns direction matrix tuple
```

#### Example:

```python
# Create test image
img = Imaginable()
img.setImageFromNumpy(np.random.rand(10, 20, 30))

# Check current orientation
print(img.getOrientationCode())  # 'LPS' (default)

# Physically reorient to RAS - THIS CHANGES THE ARRAY
img.reorientToRAS()
arr_ras = img.getImageAsNumpy()  # Different voxel arrangement!

# vs. just changing metadata - ARRAY STAYS SAME
img2 = Imaginable()
img2.setImageFromNumpy(np.random.rand(10, 20, 30))
arr_before = img2.getImageAsNumpy().copy()
img2.setDirectionCosines((-1, 0, 0, 0, -1, 0, 0, 0, 1))
arr_after = img2.getImageAsNumpy()
# arr_before == arr_after  (identical!)
```

**When to use which:**
- Use `dicomOrient()`/`reorientToLPS()` when you need voxels physically arranged in LPS/RAS order
- Use `setDirectionCosines()` when you only need to change metadata/interpretation
- For deep learning: usually want `reorientToLPS()` or `reorientToRAS()` for consistent array layouts

## Backward Compatibility

All v2 behavior is available through deprecated methods:
- `getImageAsNumpyXYZ()` - Returns (X,Y,Z) like old `getImageAsNumpy()`
- `setImageFromNumpyXYZ()` - Accepts (X,Y,Z) like old `setImageFromNumpy()`
- Old coordinate methods still work but are marked deprecated

These will be removed in a future version, so migrate when possible.

## Questions?

- Transformations still in physical space (mm)? **Yes!** Only numpy arrays changed.
- Do I need to change my registration code? **No!** Transformations unchanged.
- Will my old scripts break? **Possibly**, but deprecated methods provide quick fix.
- Is this worth updating? **Yes!** v3 is much clearer and PyTorch-ready.

