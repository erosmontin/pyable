# Vector Fields and Time Series Guide

Complete guide to working with vector fields and time-series images in pyable using the new `Vectorable` and `TimeSeriesable` classes.

## Table of Contents

1. [Overview](#overview)
2. [Vectorable - Vector Field Operations](#vectorable---vector-field-operations)
3. [TimeSeriesable - Time Series Operations](#timeseriesable---time-series-operations)
4. [Common Operations](#common-operations)
5. [Advanced Workflows](#advanced-workflows)
6. [Performance Tips](#performance-tips)
7. [API Reference](#api-reference)

---

## Overview

### Design Philosophy

Both `Vectorable` and `TimeSeriesable` **inherit from `Imaginable`**, providing:

✅ **All inherited capabilities:**
- Geometry operations (rotate, scale, translate)
- Transformations (apply transforms, displacement fields)
- I/O (read/write multiple formats)
- Filtering (Gaussian, median, bilateral)
- NumPy integration
- Undo/redo stack

✅ **Specialized for their data types:**
- **Vectorable:** Vector field specific operations (magnitude, components, normalization)
- **TimeSeriesable:** Temporal operations (frame extraction, temporal statistics)

### When to Use Each Class

| Use Case | Class |
|----------|-------|
| Displacement field warping | `Vectorable` |
| Velocity field smoothing | `Vectorable` |
| Vector magnitude computation | `Vectorable` |
| Cardiac cine sequences | `TimeSeriesable` |
| Dynamic MRI (DCE, DSC) | `TimeSeriesable` |
| Functional MRI time series | `TimeSeriesable` |
| Regular 3D images | `Imaginable` |
| Segmentation masks | `Roiable` |
| Multi-label maps | `LabelMapable` |

---

## Vectorable - Vector Field Operations

### Basic Usage

```python
from pyable import Vectorable

# Load vector field (displacement field, velocity field, etc.)
vf = Vectorable('displacement_field.mha')

# Check properties
print(f"Components: {vf.getNumberOfComponents()}")  # Usually 2 or 3
print(f"Size: {vf.getImageSize()}")
print(f"Dimension: {vf.getImageDimension()}")
```

### Vector-Specific Operations

#### Get Magnitude

```python
vf = Vectorable('displacement.mha')

# Get magnitude at each voxel (scalar image)
magnitude = vf.getMagnitude()
magnitude.write('displacement_magnitude.nii.gz')

# Check statistics
stats = vf.getVectorStatistics()
print(f"Mean displacement: {stats['mean']:.3f} mm")
print(f"Max displacement: {stats['max']:.3f} mm")
```

#### Extract/Set Components

```python
# Extract individual components as scalar images
x_component = vf.getComponent(0)  # X displacement
y_component = vf.getComponent(1)  # Y displacement
z_component = vf.getComponent(2)  # Z displacement

# Process individual component
x_component.applyGaussianSmoothing(sigma=1.0)

# Update vector field with processed component
vf.setComponent(0, x_component)

vf.write('updated_displacement.mha')
```

#### Scale Vectors

```python
vf = Vectorable('deformation.mha')

# Scale all vectors 2x
vf.scaleVector(2.0)

# Scale different components differently
vf.scaleVector([2.0, 2.0, 1.0])  # In-plane 2x, through-plane 1x

vf.write('scaled_deformation.mha')
```

#### Normalize Vectors

```python
# Unit vectors
vf.normalize(target_magnitude=1.0)

# Arbitrary magnitude
vf.normalize(target_magnitude=5.0)

vf.write('normalized_vectors.mha')
```

### Vector Field Filtering

```python
vf = Vectorable('noisy_displacement.mha')

# Gaussian smoothing
vf.applyGaussianSmoothing(sigma=1.5)
vf.write('smooth_displacement.mha')

# Can also use inherited methods
vf.cast('float')
vf.multiply(0.5)  # Scale by 0.5
vf.applyGaussianSmoothing(sigma=2.0)
```

### Vector Field Statistics

```python
vf = Vectorable('displacement.mha')

# Get comprehensive statistics
stats = vf.getVectorStatistics()
print(f"""
Displacement Statistics:
  Mean:  {stats['mean']:.3f} mm
  Std:   {stats['std']:.3f} mm
  Min:   {stats['min']:.3f} mm
  Max:   {stats['max']:.3f} mm
""")

# Get mean vector
mean_vec = vf.getMeanVector()
print(f"Mean vector: ({mean_vec[0]:.1f}, {mean_vec[1]:.1f}, {mean_vec[2]:.1f})")
```

### Vector Field Transformations

```python
# Rotate vector field
vf.rotateImage([10, 0, 0])  # Rotate 10° around X axis

# Scale vector field geometry
vf.scaleImage([2, 2, 2])  # 2x larger domain

# Transform using registration
vf.applyTransform('transform.tfm')

vf.write('transformed_vf.mha')
```

---

## TimeSeriesable - Time Series Operations

### Basic Usage

```python
from pyable import TimeSeriesable

# Load 4D image (e.g., cardiac cine, dynamic MRI)
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Check properties
print(f"Number of frames: {ts.getNumberOfFrames()}")
print(f"Spatial size: {ts.getImageSize()[:3]}")
print(f"Temporal resolution: {ts.getImageSpacing()[3]} s")
```

### Frame Access

```python
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Get specific frame
frame_0 = ts.getFrame(0)  # First frame
frame_5 = ts.getFrame(5)  # 6th frame
diastole = ts.extractPhase(0)  # Semantic alias

# Get frame range
systole_range = ts.getFrameRange(3, 7)  # Frames 3-6
systole_range.write('systole_frames.nii.gz')

# Modify frame
frame_0.applyGaussianSmoothing(sigma=1.0)
ts.setFrame(0, frame_0)
```

### Temporal Statistics

```python
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Temporal mean (average across all frames)
mean_vol = ts.getTemporalMean()
mean_vol.write('mean_cardiac_frame.nii.gz')

# Temporal variance
variance = ts.getTemporalVariance()
variance.write('temporal_variance.nii.gz')

# Temporal standard deviation
std_dev = ts.getTemporalStandardDeviation()
std_dev.write('temporal_std.nii.gz')
```

### Apply Filters to All Frames

```python
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Gaussian smoothing on all frames
ts.applyFilterToAllFrames('gaussian', sigma=1.0)

# Median filtering
ts.applyFilterToAllFrames('median', radius=2)

# Bilateral filtering (edge-preserving)
ts.applyFilterToAllFrames('bilateral', 
                          domain_sigma=1.0,
                          range_sigma=20.0)

ts.write('filtered_cardiac_4d.nii.gz')
```

### Transform All Frames

```python
# Apply same transformation to all frames
ts.transformAllFrames('rigid_transform.tfm', 
                      interpolator='linear')

# Or with displacement field
ts.transformAllFrames('deformation.mha')

ts.write('motion_corrected_4d.nii.gz')
```

### Method Chaining

```python
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Chain operations
ts.applyFilterToAllFrames('gaussian', sigma=1.0)\
  .scaleImage([2, 2, 2])\  # Upsample 2x
  .rotateImage([0, 0, 10])\  # Rotate 10°
  .write('processed_cardiac.nii.gz')
```

---

## Common Operations

### Operation: Extract and Process Single Frame

```python
from pyable import TimeSeriesable, Imaginable

ts = TimeSeriesable('dynamic_mri.nii.gz')

# Extract diastolic frame
diastole = ts.getFrame(0)  # 3D Imaginable

# Process as normal Imaginable
diastole.cast('float')
diastole.divide(diastole.getMaximumValue())  # Normalize
diastole.applyGaussianSmoothing(sigma=2.0)

# Put back in series
ts.setFrame(0, diastole)
```

### Operation: Compare Frames

```python
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Extract two phases
frame_0 = ts.getFrame(0)
frame_5 = ts.getFrame(5)

# Compute difference
diff = frame_0.getDuplicate()
diff_array = sitk.GetArrayFromImage(frame_0.getImage()) - \
             sitk.GetArrayFromImage(frame_5.getImage())
diff.setImageFromNumpy(diff_array)

diff.write('frame_difference.nii.gz')
```

### Operation: Create Maximum/Minimum Projections

```python
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Get temporal mean
array_4d = sitk.GetArrayFromImage(ts.getImage())

# Maximum intensity projection over time
mip = np.max(array_4d, axis=0)
mip_image = sitk.GetImageFromArray(mip)

# Minimum intensity projection
min_ip = np.min(array_4d, axis=0)
min_image = sitk.GetImageFromArray(min_ip)
```

---

## Advanced Workflows

### Workflow 1: Motion Correction in Time Series

```python
from pyable import TimeSeriesable

# Load dynamic series
ts = TimeSeriesable('uncorrected_4d.nii.gz')
ref_frame = ts.getFrame(0)  # Reference frame

# For each frame, compute transform to reference and apply
# (In practice, use registration tool like ANTs)
# ts.transformAllFrames('frame_to_reference.tfm')

# Save motion-corrected
ts.write('motion_corrected_4d.nii.gz')
```

### Workflow 2: Vector Field Manipulation

```python
from pyable import Vectorable

# Load displacement field
vf = Vectorable('forward_displacement.mha')

# Invert (for reverse warping)
vf.invertDisplacementField()

# Smooth to regularize
vf.applyGaussianSmoothing(sigma=2.0)

# Scale
vf.scaleVector(1.5)

# Get magnitude for quality assessment
magnitude = vf.getMagnitude()
mag_stats = magnitude.getImageStatistics()
print(f"Max displacement magnitude: {mag_stats['max']:.2f}")

vf.write('processed_displacement.mha')
```

### Workflow 3: Extract Dynamic Contrast Enhanced (DCE) Kinetics

```python
ts = TimeSeriesable('dce_mri_4d.nii.gz')

# Get temporal statistics per voxel
mean_intensity = ts.getTemporalMean()
temporal_std = ts.getTemporalStandardDeviation()

# Compute enhancement maps
enhancement = temporal_std.getDuplicate()
enhancement_array = sitk.GetArrayFromImage(temporal_std.getImage()) / \
                    sitk.GetArrayFromImage(mean_intensity.getImage())
enhancement.setImageFromNumpy(enhancement_array)

enhancement.write('dce_enhancement_map.nii.gz')
```

### Workflow 4: Cardiac Function Analysis

```python
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Get key phases
ed = ts.extractPhase(0)  # End-diastole
es = ts.extractPhase(5)  # End-systole

# Segment ventricle at ED and ES
# (Use external segmentation tool)

# Calculate strain via temporal statistics
temporal_var = ts.getTemporalVariance()
temporal_var.write('cardiac_strain_map.nii.gz')

# Average across all phases
mean_frame = ts.getTemporalMean()
mean_frame.write('mean_cardiac_frame.nii.gz')
```

---

## Performance Tips

### Memory Efficiency

```python
# For large time series, process frame-by-frame
ts = TimeSeriesable('large_4d_file.nii.gz')

for i in range(ts.getNumberOfFrames()):
    frame = ts.getFrame(i)
    frame.cast('uint8')  # Reduce precision
    ts.setFrame(i, frame)
```

### Batch Processing

```python
import glob
from pyable import TimeSeriesable

# Process all 4D files in directory
for filepath in glob.glob('cardiac_data/*.nii.gz'):
    ts = TimeSeriesable(filepath)
    ts.applyFilterToAllFrames('gaussian', sigma=1.0)
    ts.write(f'processed/{Path(filepath).name}')
```

### Parallel Frame Processing

```python
from multiprocessing import Pool

def process_frame(args):
    ts_path, frame_idx = args
    ts = TimeSeriesable(ts_path)
    frame = ts.getFrame(frame_idx)
    frame.cast('float')
    return frame.getImage()

# Process frames in parallel
if __name__ == '__main__':
    ts_file = 'cardiac_4d.nii.gz'
    ts = TimeSeriesable(ts_file)
    
    with Pool(4) as p:
        frames = p.map(process_frame, 
                      [(ts_file, i) for i in range(ts.getNumberOfFrames())])
```

---

## API Reference

### Vectorable Class

| Method | Purpose | Returns |
|--------|---------|---------|
| `getNumberOfComponents()` | Get # of vector components | int |
| `getMagnitude()` | Get vector magnitudes | Imaginable |
| `getComponent(i)` | Extract component i | Imaginable |
| `setComponent(i, img)` | Set component i | self |
| `scaleVector(factors)` | Scale vectors | self |
| `normalize(magnitude)` | Normalize vectors | self |
| `applyGaussianSmoothing(sigma)` | Smooth vectors | self |
| `getVectorStatistics()` | Statistics on magnitudes | dict |
| `getMeanVector()` | Mean of all vectors | np.ndarray |
| `getDuplicate()` | Copy vector field | Vectorable |

**Inherited from Imaginable:**
- Geometry: `rotateImage()`, `scaleImage()`, `translateImage()`
- Transformations: `applyTransform()`, `applyDisplacementField()`
- Filtering: `cast()`, `multiply()`, `add()`, `divide()`
- I/O: `read()`, `write()`, `getImage()`, `setImage()`

### TimeSeriesable Class

| Method | Purpose | Returns |
|--------|---------|---------|
| `getNumberOfFrames()` | Get # of frames | int |
| `getFrame(i)` | Extract frame i | Imaginable |
| `setFrame(i, img)` | Replace frame i | self |
| `getFrameRange(s, e)` | Extract frames s:e | TimeSeriesable |
| `getTemporalMean()` | Average across time | Imaginable |
| `getTemporalVariance()` | Variance across time | Imaginable |
| `getTemporalStandardDeviation()` | Std dev across time | Imaginable |
| `applyFilterToAllFrames(name, **kw)` | Filter all frames | self |
| `transformAllFrames(tfm, ...)` | Transform all frames | self |
| `extractPhase(n)` | Get frame (semantic alias) | Imaginable |
| `getDuplicate()` | Copy time series | TimeSeriesable |

**Inherited from Imaginable:**
- Geometry: `rotateImage()`, `scaleImage()`, `translateImage()`
- Transformations: `applyTransform()`, `applyDisplacementField()`
- Filtering: All scalar filters work on full 4D
- I/O: `read()`, `write()`, `getImage()`, `setImage()`

---

## See Also

- [Imaginable API](./QUICKREF.md)
- [Deformation Workflow](./DEFORMATION_WORKFLOW.md)
- [SimpleITK Vector Images](https://simpleitk.org/doxygen/latest/html/classitk_1_1simple_1_1Image.html)
- [4D Image Processing](https://simpleitk.org/doxygen/latest/html/structitk_1_1simple_1_1Image.html)

---

**Last Updated**: November 21, 2025  
**Version**: 1.0  
**Status**: Beta
