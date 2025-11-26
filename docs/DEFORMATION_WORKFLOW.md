# Deformation & Registration Workflow Guide

Complete guide to applying registration transforms and displacement fields to Imaginable objects in pyable.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Transform Types](#transform-types)
3. [Basic Usage](#basic-usage)
4. [Advanced Workflows](#advanced-workflows)
5. [Geometry Management](#geometry-management)
6. [Batch Operations](#batch-operations)
7. [Performance Considerations](#performance-considerations)
8. [Troubleshooting](#troubleshooting)

---

## Quick Start

### Apply a Transform to an Image

```python
from pyable_eros_montin import SITKImaginable

# Load image and apply registration transform
image = SITKImaginable('moving_image.nii.gz')
image.applyTransform('transform.tfm', interpolator='linear')
image.write('warped_image.nii.gz')
```

### Apply a Displacement Field

```python
# Apply displacement field (e.g., from ANTs, elastix)
image = SITKImaginable('moving_image.nii.gz')
image.applyDisplacementField('deformation.mha', target_image='fixed_image.nii.gz')
image.write('warped_image.nii.gz')
```

### Warp a Segmentation/ROI

```python
from pyable_eros_montin import Roiable, LabelMapable

# Warp ROI with label preservation (nearest-neighbor)
roi = Roiable('segmentation.nii.gz')
roi.warpROI('deformation.mha')
roi.write('warped_roi.nii.gz')

# Warp multi-label segmentation
labels = LabelMapable('label_map.nii.gz')
labels.warpLabelMap('deformation.mha')
labels.write('warped_labels.nii.gz')
```

---

## Transform Types

### 1. Affine Transforms (.tfm, .h5)

**Source**: Rigid or affine registration algorithms

**Characteristics**:
- 12 parameters (3x4 matrix)
- Fast to apply
- Suitable for inter-patient alignment

**Usage**:
```python
image = SITKImaginable('moving.nii.gz')
image.applyTransform('affine.tfm', interpolator='linear')
```

### 2. B-spline Transforms (.tfm, .h5)

**Source**: Deformable registration (elastix, ANTs)

**Characteristics**:
- Grid-based deformation
- Smooth, C2-continuous
- Higher computational cost
- Better for organ-level alignment

**Usage**:
```python
image = SITKImaginable('moving.nii.gz')
image.applyTransform('bspline.tfm', interpolator='linear')
```

### 3. Displacement Fields (.mha, .nii.gz)

**Source**: ANTs, elastix, or other registration tools

**Characteristics**:
- Dense vector field
- Defines displacement at each voxel
- Can be non-invertible
- Direct warping without interpolation

**Usage**:
```python
image = SITKImaginable('moving.nii.gz')
image.applyDisplacementField('deformation.mha', target_image='fixed.nii.gz')
```

### 4. Composite Transforms

**Source**: Multi-step registration (e.g., rigid → affine → deformable)

**Characteristics**:
- Multiple transforms applied sequentially
- Preserves registration history
- More accurate than single step

**Usage**:
```python
from pyable_eros_montin.deformations import apply_multi_step_transform

result = apply_multi_step_transform(
    'moving.nii.gz',
    ['rigid.tfm', 'affine.tfm', 'bspline.tfm'],
    target_image='fixed.nii.gz'
)
```

---

## Basic Usage

### Method Chaining

All deformation methods return `self`, enabling method chaining:

```python
image = SITKImaginable('moving.nii.gz')
image.applyTransform('transform.tfm', interpolator='linear')\
      .alignGeometry('fixed.nii.gz')\
      .cast('uint8')\
      .write('result.nii.gz')
```

### Interpolation Methods

| Method | Use Case | Speed | Smoothness |
|--------|----------|-------|-----------|
| `linear` | General (default) | Fast | Good |
| `nearest` | Labels/segmentations | Fastest | Poor |
| `gaussian` | Smooth data | Slow | Best |
| `bspline` | High-quality | Slower | Excellent |

**Example**:
```python
# For continuous intensity images
image.applyTransform('transform.tfm', interpolator='linear')

# For segmentation maps
roi.warpROI('deformation.mha')  # Automatically uses nearest-neighbor
```

### Target Geometry

By default, the output has the same geometry as the input. To resample to a different geometry:

```python
# Output will match 'fixed_image.nii.gz' geometry
image.applyDisplacementField('deformation.mha', target_image='fixed_image.nii.gz')
```

### Default Pixel Value

For regions outside the original image domain:

```python
# Use 0 for background (default)
image.applyTransform('transform.tfm', default_value=0)

# Use -1 to mark out-of-bounds regions
image.applyTransform('transform.tfm', default_value=-1)
```

---

## Advanced Workflows

### Multi-Step Registration

Apply rigid → affine → deformable in sequence:

```python
from pyable_eros_montin.deformations import (
    apply_multi_step_transform,
    create_composite_transform
)

# Method 1: Sequential application
image = SITKImaginable('moving.nii.gz')
image.applyTransform('rigid.tfm', target_image='fixed.nii.gz')
image.applyTransform('affine.tfm', target_image='fixed.nii.gz')
image.applyDisplacementField('bspline.mha', target_image='fixed.nii.gz')
image.write('result.nii.gz')

# Method 2: Composite (single resampling - more accurate)
result = apply_multi_step_transform(
    'moving.nii.gz',
    ['rigid.tfm', 'affine.tfm', 'bspline.tfm'],
    target_image='fixed.nii.gz',
    interpolator='linear'
)
```

### Forward-Backward Consistency

Verify registration quality by checking forward-backward consistency:

```python
from pyable_eros_montin.deformations import apply_inverted_deformation_field

# Apply forward transformation
moving = SITKImaginable('moving.nii.gz')
moving.applyDisplacementField('forward.mha')

# Apply backward transformation
backward = apply_inverted_deformation_field(
    moving.getImage(),
    'forward.mha',
    target_image='moving.nii.gz'
)

# Compare: should be close to original moving image
```

### Invert Displacement Field

For reverse warping (going from template to patient):

```python
df = SITKImaginable('forward_deform.mha')
df.invertDisplacementField()
df.write('backward_deform.mha')

# Now use backward field
moving.applyDisplacementField('backward_deform.mha')
```

### Refine B-spline Grid

For multi-resolution registration refinement:

```python
from pyable_eros_montin.deformations import refine_bspline_grid
import SimpleITK as sitk

# Load coarse B-spline transform
coarse_tfm = sitk.ReadTransform('coarse_bspline.tfm')

# Refine to finer grid
fine_tfm = refine_bspline_grid(coarse_tfm, new_mesh_size=(7, 7, 7))

# Save refined transform
sitk.WriteTransform(fine_tfm, 'fine_bspline.tfm')
```

---

## Geometry Management

### Align Geometry

Fix displacement fields or images with lost/incorrect metadata:

```python
# Problem: Displacement field lost geometry during processing
df = SITKImaginable('deform.mha')
fixed = SITKImaginable('fixed.nii.gz')

# Solution: Align to reference
df.alignGeometry(fixed.getImage())
df.write('deform_aligned.mha')
```

### Check Geometry

Verify geometry before applying transforms:

```python
image = SITKImaginable('moving.nii.gz')
fixed = SITKImaginable('fixed.nii.gz')

print(f"Moving: size={image.getImageSize()}, spacing={image.getImageSpacing()}")
print(f"Fixed: size={fixed.getImageSize()}, spacing={fixed.getImageSpacing()}")
```

---

## Batch Operations

### Warp Multiple Images with Same Transform

```python
import glob
from pyable_eros_montin import SITKImaginable

# Register all images in a dataset
moving_images = glob.glob('data/moving/*.nii.gz')

for moving_path in moving_images:
    img = SITKImaginable(moving_path)
    img.applyTransform('template_to_patient.tfm')
    img.write(f'results/{moving_path.stem}_warped.nii.gz')
```

### Warp All Segmentations

```python
from pyable_eros_montin import LabelMapable

# Warp all label maps to patient space
label_maps = glob.glob('templates/organs/*.nii.gz')
displacement_field = 'registration/template_to_patient.mha'

for label_path in label_maps:
    labels = LabelMapable(label_path)
    labels.warpLabelMap(displacement_field, target_image='patient_image.nii.gz')
    labels.write(f'results/{Path(label_path).stem}_warped.nii.gz')
```

---

## Performance Considerations

### Memory vs Speed Trade-off

```python
# Faster but uses more memory
image.applyTransform('transform.tfm', interpolator='linear')

# Slower but memory-efficient - convert to smaller dtype first
image.cast('uint8')
image.applyTransform('transform.tfm', interpolator='linear')
```

### Large Images

For 3D images >512×512×512, consider:

```python
import SimpleITK as sitk

# Method 1: Stream processing (if supported by resampler)
# Default SimpleITK doesn't support streaming

# Method 2: Process in blocks (manual)
def warp_large_image_in_blocks(image_path, transform_path, block_size=256):
    # Not recommended - SimpleITK should handle large images
    pass

# Method 3: Simplest - just apply normally
image = SITKImaginable(large_image_path)
image.applyTransform(transform_path)
```

### Parallel Processing

```python
from multiprocessing import Pool
import glob

def warp_image(img_path):
    img = SITKImaginable(img_path)
    img.applyTransform('transform.tfm')
    img.write(f'results/{Path(img_path).stem}_warped.nii.gz')

if __name__ == '__main__':
    images = glob.glob('data/*.nii.gz')
    with Pool(4) as p:  # 4 processes
        p.map(warp_image, images)
```

---

## Troubleshooting

### Issue: Transform not found

```
Exception: Error reading file: 'transform.tfm'
```

**Solution**:
```python
from pathlib import Path

# Check file exists
tfm_path = Path('transform.tfm')
if not tfm_path.exists():
    raise FileNotFoundError(f"Transform not found: {tfm_path.absolute()}")

image.applyTransform(str(tfm_path))
```

### Issue: Geometry mismatch

```
Error: Output geometry parameters don't match reference image
```

**Solution**:
```python
# Option 1: Specify target geometry
image.applyDisplacementField('deform.mha', target_image='fixed.nii.gz')

# Option 2: Align first
image.alignGeometry('fixed.nii.gz')
image.applyDisplacementField('deform.mha')
```

### Issue: Labels become gradients after warping

```python
# Problem: Using linear interpolation on label map
labels.applyDisplacementField('deform.mha', interpolator='linear')  # ❌

# Solution: Use nearest-neighbor for labels
labels.warpLabelMap('deform.mha')  # ✓ Automatically uses nearest-neighbor
```

### Issue: Displacement field has wrong geometry

```python
from pyable_eros_montin.deformations import align_geometry

# Fix geometry from reference
df = SITKImaginable('deform.mha')
df.alignGeometry('fixed.nii.gz')
df.write('deform_fixed.mha')

# Now apply normally
image.applyDisplacementField('deform_fixed.mha')
```

### Issue: Transform inverts image coordinates

```python
# Some tools define transforms in different coordinate systems
# If result appears invipped/rotated, try inverse transform

import SimpleITK as sitk

tfm = sitk.ReadTransform('transform.tfm')
# Swap to inverse (if applicable)
# Note: Not all transforms are invertible
```

---

## Complete Example: Atlas-Based Segmentation

```python
from pyable_eros_montin import SITKImaginable, LabelMapable
from pathlib import Path

def atlas_based_segmentation(patient_image, atlas_image, atlas_labels):
    """
    Segment patient using atlas labels via deformable registration.
    
    Parameters
    ----------
    patient_image : str
        Path to patient image
    atlas_image : str
        Path to atlas template image
    atlas_labels : str
        Path to atlas label map
    
    Returns
    -------
    segmentation : LabelMapable
        Patient segmentation from warped atlas
    """
    # Step 1: Load images
    patient = SITKImaginable(patient_image)
    atlas = SITKImaginable(atlas_image)
    atlas_seg = LabelMapable(atlas_labels)
    
    # Step 2: Rigid registration (alignment)
    # In practice, use ANTs or elastix for this
    # For demo, assume we have the transform
    rigid_tfm = 'registration/rigid_atlas_to_patient.tfm'
    
    # Step 3: Deformable registration
    # In practice, output from ANTs or elastix
    deform_field = 'registration/deform_atlas_to_patient.mha'
    
    # Step 4: Warp atlas segmentation to patient space
    atlas_seg.warpLabelMap(deform_field, target_image=patient_image)
    
    # Step 5: Optional refinement (morphological operations)
    atlas_seg.fillHoles()
    atlas_seg.keepLargestObject()
    
    return atlas_seg

# Usage
seg = atlas_based_segmentation(
    'patient.nii.gz',
    'atlas.nii.gz',
    'atlas_labels.nii.gz'
)
seg.write('patient_segmentation.nii.gz')
```

---

## API Reference

### Imaginable Class Methods

| Method | Purpose | Returns |
|--------|---------|---------|
| `applyTransform(transform, ...)` | Apply registration transform | self |
| `applyDisplacementField(field, ...)` | Apply displacement field | self |
| `warpImage(field, ...)` | Alias for applyDisplacementField | self |
| `alignGeometry(reference_image)` | Align geometry to reference | self |
| `invertDisplacementField(...)` | Invert displacement field | self |

### Roiable Class Methods

| Method | Purpose | Returns |
|--------|---------|---------|
| `applyTransform(transform, ...)` | Warp ROI with label-preserving behavior (nearest-neighbor to preserve labels) | self |
| `warpROI(field, ...)` | Warp ROI with displacement field (label-preserving) | self |

### LabelMapable Class Methods

| Method | Purpose | Returns |
|--------|---------|---------|
| `applyTransform(transform, ...)` | Warp multi-label map with label-preserving behavior | self |
| `warpLabelMap(field, ...)` | Warp multi-label with displacement field (label-preserving) | self |

### Deformations Module Functions

| Function | Purpose |
|----------|---------|
| `initialize_deformation_field(ref)` | Create empty deformation field |
| `transform_to_displacement_field(tfm, ...)` | Convert transform to displacement field |
| `apply_deformation_field(image, field, ...)` | Apply displacement field |
| `apply_transform(image, tfm, ...)` | Apply transform |
| `invert_displacement_field(field, ...)` | Invert displacement field |
| `apply_inverted_deformation_field(...)` | Apply inverted field |
| `align_geometry(moving, reference)` | Align geometry |
| `create_composite_transform(tfms, ...)` | Create composite transform |
| `apply_multi_step_transform(image, tfms, ...)` | Apply multiple transforms |

---

## See Also

- [SimpleITK Registration Transforms](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1Transform.html)
- [ANTs (Advanced Normalization Tools)](http://stnava.github.io/ANTs/)
- [Elastix - Image Registration](https://elastix.lumc.nl/)
- [pyable Imaginable API](./QUICKREF.md)

---

**Last Updated**: November 21, 2025  
**Version**: 1.0  
**Status**: Active Development
