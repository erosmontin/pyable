# PyAble v3 - Image Orientation Guide

## Understanding Medical Image Orientations

Medical images have **anatomical orientations** that define how voxel data relates to the patient's body. PyAble v3 provides clear methods to handle these orientations correctly.

## Coordinate Systems Overview

### 1. Array/Numpy Indices (k, j, i)
- Integer indices starting from 0
- In PyAble v3: **(Z, Y, X)** ordering
- Example: `array[k, j, i]` where k=slice, j=row, i=column

### 2. Physical/World Coordinates (x, y, z)
- Float coordinates in millimeters
- Defined by origin, spacing, and direction matrix
- Maps to anatomical directions (Left/Right, Anterior/Posterior, Superior/Inferior)

### 3. Anatomical Orientation Codes
Three-letter codes describing which anatomical direction each axis increases toward:
- **L/R**: Left / Right
- **P/A**: Posterior / Anterior
- **S/I**: Superior / Inferior

Common codes:
- **LPS**: DICOM standard (X: L→R, Y: P→A, Z: I→S)
- **RAS**: NIfTI/Neuroimaging (X: R→L, Y: A→P, Z: I→S)
- **RPI**: Alternative (X: R→L, Y: P→A, Z: S→I)

## Key Methods

### Query Orientation

```python
# Get current anatomical orientation code
orientation = img.getOrientationCode()  # Returns 'LPS', 'RAS', 'RPI', etc.

# Get direction cosine matrix
direction = img.getDirectionCosines()   # Returns 9-tuple for 3D
# Example: (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0) = identity = LPS
```

### Reorient Image Data (Changes Numpy Array!)

```python
# General method - specify any 3-letter code
img.dicomOrient('LPS')    # Reorient to LPS
img.dicomOrient('RAS')    # Reorient to RAS
img.dicomOrient('RPI')    # Reorient to RPI

# Convenience methods
img.reorientToLPS()       # Same as dicomOrient('LPS')
img.reorientToRAS()       # Same as dicomOrient('RAS')
img.reorientToRPI()       # Same as dicomOrient('RPI')
```

**What these methods do:**
- ✅ **Physically rearrange voxel data** (numpy array changes!)
- ✅ Permute and/or flip axes to match requested orientation
- ✅ Update direction matrix (usually to near-identity)
- ⚡ **No interpolation** (just permutes/flips existing voxels)
- ✅ Best for axis-aligned images that need reorienting

### Resample Oblique Images to Axis-Aligned Grid

```python
# For oblique/rotated acquisitions - resample to standard grid
img.resampleToAxisAligned()  # Sets direction to identity (1,0,0, 0,1,0, 0,0,1)

# Or specify custom target direction
img.changeImageDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

# Check if already axis-aligned
if not img.isAxisAligned():
    img.resampleToAxisAligned()
```

**What these methods do:**
- ✅ **Resamples voxel data** (numpy array changes!)
- ⚠️  **Applies interpolation** (may slightly change voxel values)
- ✅ Converts oblique acquisitions to axis-aligned grid
- ✅ Direction matrix becomes identity: (1,0,0, 0,1,0, 0,0,1)
- ✅ Preserves physical coordinates (mm)
- ✅ Best for oblique acquisitions (e.g., oblique MRI scans)
- ✅ Preserve physical world coordinates (same anatomy in same mm locations)
- ✅ May change image size if axes are permuted

### Change Metadata Only (Does NOT Change Numpy Array!)

```python
# Only update direction matrix, don't move voxels
img.setDirectionCosines((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
img.setImageDirection((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
```

**What these methods do:**
- ✅ Update direction matrix metadata
- ❌ **Do NOT change numpy array** (voxels stay in same positions)
- ❌ Do NOT permute or flip axes
- ⚠️  Changes interpretation of voxel→world mapping only

## Critical Distinction: Data Reorientation vs Metadata Update

### Example 1: Reorienting Data (dicomOrient)

```python
import numpy as np
from pyable_eros_montin.imaginable import Imaginable

# Create image with marker at corner
array = np.zeros((10, 20, 30))
array[0, 0, 0:5] = 100  # Marker in corner
img = Imaginable()
img.setImageFromNumpy(array)

print(f"Original orientation: {img.getOrientationCode()}")  # 'LPS'
print(f"Marker at: array[0,0,0:5] = {img.getImageAsNumpy()[0,0,0:5]}")

# Physically reorient to RAS
img.reorientToRAS()
print(f"New orientation: {img.getOrientationCode()}")  # 'RAS'
new_array = img.getImageAsNumpy()
print(f"Marker at: array[0,0,0:5] = {new_array[0,0,0:5]}")  # [0. 0. 0. 0. 0.]
# Marker moved! Now at different array indices

# The physical location is preserved, but array arrangement changed
```

### Example 2: Metadata Only (setDirectionCosines)

```python
# Create same image
array = np.zeros((10, 20, 30))
array[0, 0, 0:5] = 100
img = Imaginable()
img.setImageFromNumpy(array)

before = img.getImageAsNumpy().copy()

# Only change metadata
img.setDirectionCosines((-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0))

after = img.getImageAsNumpy()
print(f"Arrays equal: {np.array_equal(before, after)}")  # True!
# Array didn't change at all!
```

### Example 3: Handling Oblique Acquisitions

```python
# Load oblique MRI scan (e.g., cardiac oblique, tilted brain scan)
img = Imaginable(imagepath='oblique_cardiac.nii.gz')

print(f"Direction: {img.getDirectionCosines()}")
# Output: (0.866, 0.5, 0.0, -0.5, 0.866, 0.0, 0.0, 0.0, 1.0)  # Oblique!

print(f"Is axis-aligned? {img.isAxisAligned()}")  # False

# Resample to axis-aligned grid
img.resampleToAxisAligned()

print(f"New direction: {img.getDirectionCosines()}")
# Output: (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)  # Identity!

print(f"Is axis-aligned? {img.isAxisAligned()}")  # True

# Now array axes align with anatomical axes
arr = img.getImageAsNumpy()  # Clean (Z,Y,X) with identity direction
```

## Common Use Cases

### For Deep Learning / PyTorch

When preparing images for neural networks, you typically want **consistent array layouts**:

```python
# Load diverse medical images (some may be oblique!)
images = []
for filepath in image_files:
    img = Imaginable(imagepath=filepath)
    
    # Handle oblique acquisitions first
    if not img.isAxisAligned():
        img.resampleToAxisAligned()  # Resample oblique to axis-aligned
    
    # Then standardize anatomical orientation
    img.reorientToLPS()
    
    # Now all images have same array axis meanings:
    # array[k,j,i] where k=inferior→superior, j=posterior→anterior, i=left→right
    
    arr = img.getImageAsNumpy()  # (Z,Y,X) in LPS orientation
    images.append(arr)

# All arrays now have consistent anatomical meaning AND axis-aligned grids
```

### For Visualization

```python
# Display axial slices in standard radiological view
img.reorientToLPS()
axial_slice = img.getImageAsNumpy()[slice_idx, :, :]  # (Y,X) = (P→A, L→R)

# For neuroimaging convention
img.reorientToRAS()
axial_slice = img.getImageAsNumpy()[slice_idx, :, :]  # (Y,X) = (A→P, R→L)
```

### For DICOM Export

```python
# Ensure DICOM standard orientation before export
img.reorientToLPS()
img.writeImageAs('output.dcm')
```

### Working with Multiple Images

```python
# Register two images - ensure same orientation first
img1.reorientToLPS()
img2.reorientToLPS()

# Now perform registration
# Arrays have consistent anatomical meaning
```

## Direction Matrix Details

The direction matrix is a 3x3 matrix (stored as 9-tuple) that defines how voxel indices map to physical coordinates:

```
[dir_x0, dir_x1, dir_x2]   <- Direction cosines for i-axis (columns)
[dir_y0, dir_y1, dir_y2]   <- Direction cosines for j-axis (rows)  
[dir_z0, dir_z1, dir_z2]   <- Direction cosines for k-axis (slices)
```

**Identity matrix = LPS orientation:**
```python
(1.0, 0.0, 0.0,   # i-axis points in +X direction (Left→Right)
 0.0, 1.0, 0.0,   # j-axis points in +Y direction (Posterior→Anterior)
 0.0, 0.0, 1.0)   # k-axis points in +Z direction (Inferior→Superior)
```

**RAS orientation:**
```python
(-1.0, 0.0, 0.0,  # i-axis points in -X direction (Right→Left)
  0.0,-1.0, 0.0,  # j-axis points in -Y direction (Anterior→Posterior)
  0.0, 0.0, 1.0)  # k-axis points in +Z direction (Inferior→Superior)
```

## Best Practices

1. **For Deep Learning**: Always reorient to a standard orientation (LPS or RAS) before extracting numpy arrays
2. **For Consistency**: Use same orientation across all images in a dataset
3. **For DICOM**: Use LPS (DICOM standard)
4. **For NIfTI/Neuroimaging**: Use RAS (common convention)
5. **Check First**: Use `getOrientationCode()` to verify current orientation
6. **For Oblique Acquisitions**: Use `resampleToAxisAligned()` to convert to standard grid
7. **Physical Coordinates**: Remember that physical mm coordinates are preserved regardless of reorientation

## Common Pitfalls

❌ **Wrong**: Assuming all medical images have same array layout
```python
arr = img.getImageAsNumpy()  # Could be any orientation!
# Axes might mean different anatomical directions
```

✅ **Right**: Standardize orientation first
```python
img.reorientToLPS()
arr = img.getImageAsNumpy()  # Now guaranteed LPS layout
```

❌ **Wrong**: Using setDirectionCosines when you need data reorientation
```python
img.setDirectionCosines((-1, 0, 0, 0, -1, 0, 0, 0, 1))  # Only metadata!
arr = img.getImageAsNumpy()  # Array unchanged, just interpretation different
```

✅ **Right**: Use dicomOrient for actual data reorientation
```python
img.reorientToRAS()  # Physically rearranges array
arr = img.getImageAsNumpy()  # Array is actually different
```

❌ **Wrong**: Using dicomOrient on oblique images
```python
# Oblique image with direction like (0.866, 0.5, 0, -0.5, 0.866, 0, 0, 0, 1)
img.dicomOrient('LPS')  # Won't make it axis-aligned!
# Direction stays oblique, just changes anatomical interpretation
```

✅ **Right**: Use resampleToAxisAligned for oblique images
```python
# Check if oblique
if not img.isAxisAligned():
    img.resampleToAxisAligned()  # Resamples to identity direction
# Now direction is (1, 0, 0, 0, 1, 0, 0, 0, 1)
```

## Summary

| Method | Changes Numpy Array? | Changes Direction Matrix? | Interpolation? | Use When |
|--------|---------------------|--------------------------|----------------|----------|
| `dicomOrient('LPS')` | ✅ Yes | ✅ Yes | ❌ No | Need consistent array layout (axis-aligned) |
| `reorientToLPS()` | ✅ Yes | ✅ Yes | ❌ No | Same as above (convenience) |
| `reorientToRAS()` | ✅ Yes | ✅ Yes | ❌ No | Neuroimaging convention (axis-aligned) |
| `reorientToRPI()` | ✅ Yes | ✅ Yes | ❌ No | Alternative orientation (axis-aligned) |
| `resampleToAxisAligned()` | ✅ Yes | ✅ Yes | ⚠️ Yes | **Oblique acquisitions** → identity matrix |
| `changeImageDirection()` | ✅ Yes | ✅ Yes | ⚠️ Yes | Resample to custom direction matrix |
| `isAxisAligned()` | N/A | N/A | N/A | Check if direction is identity |
| `setDirectionCosines()` | ❌ No | ✅ Yes | N/A | Only change interpretation |
| `setImageDirection()` | ❌ No | ✅ Yes | N/A | Only change interpretation |
| `getOrientationCode()` | N/A | N/A | N/A | Check current orientation |
| `getDirectionCosines()` | N/A | N/A | N/A | Get direction matrix |

## Related Documentation

- See `PYABLE_V3_BREAKING_CHANGES.md` for numpy convention changes
- See `test_orientation.py` for working examples
- See `COORDINATE_SYSTEM_REFACTORING.md` for technical details
