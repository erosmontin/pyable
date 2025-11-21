# Interactive Plotting and Visualization Guide

Complete guide to the interactive visualization system for pyable's image classes with overlay support, component/frame selection, and image exporting.

## Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Scalar Image Plotting](#scalar-image-plotting)
4. [Vector Field Plotting](#vector-field-plotting)
5. [Time Series Plotting](#time-series-plotting)
6. [Overlay Functionality](#overlay-functionality)
7. [Advanced Features](#advanced-features)
8. [API Reference](#api-reference)

---

## Overview

### Design Philosophy

The plotting system provides **unified visualization** for all data types:

| Data Type | Viewer | Features |
|-----------|--------|----------|
| **Imaginable** (2D/3D scalar) | ScalarPlotter | Slice navigation, overlay support |
| **Vectorable** (vector fields) | VectorPlotter | Component selection, magnitude view |
| **TimeSeriesable** (4D sequences) | TimeSeriesPlotter | Frame navigation, temporal statistics |

### Key Features

✅ **Interactive GUIs** with sliders for navigation  
✅ **Overlay blending** with alpha control  
✅ **Automatic resampling** to match geometries  
✅ **Component/Frame selection** via radio buttons/sliders  
✅ **Method chaining** on all classes  
✅ **Image saving** in any format  
✅ **Matplotlib-based** for seamless integration  

---

## Quick Start

### Basic Scalar Image Viewing

```python
from pyable_eros_montin import Imaginable

# Load image
img = Imaginable('brain.nii.gz')

# View with single line
img.plotOverlay()
```

### Vector Field Viewing

```python
from pyable_eros_montin import Vectorable

# Load displacement field
vf = Vectorable('deformation.mha')

# Interactive component selector GUI
vf.plotOverlay()  

# Or view specific component
vf.plotOverlay(component=0)  # X component
```

### Time Series Viewing

```python
from pyable_eros_montin import TimeSeriesable

# Load 4D cardiac sequence
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Interactive frame selector
ts.plotOverlay()

# Or specific frame
ts.plotOverlay(frame=5)
```

---

## Scalar Image Plotting

### Display Single Image

```python
from pyable_eros_montin import Imaginable

img = Imaginable('image.nii.gz')

# Show full GUI with slice slider (3D)
img.plotOverlay()

# With custom title
img.plotOverlay(title="My Brain Image")

# Start at specific slice
img.plotOverlay(slice_idx=100)
```

### Display with Overlay

```python
img = Imaginable('T1_image.nii.gz')
seg = Imaginable('segmentation.nii.gz')

# Overlay segmentation on image
img.plotOverlay(overlay=seg, alpha=0.5)

# Different opacities
img.plotOverlay(overlay=seg, alpha=0.3)  # More transparent
img.plotOverlay(overlay=seg, alpha=0.7)  # More opaque
```

### Custom Colormaps

```python
from pyable_eros_montin import ScalarPlotter

img = Imaginable('image.nii.gz')
overlay = Imaginable('heatmap.nii.gz')

# Use different colormaps
plotter = ScalarPlotter(img.getImage(), 
                       overlay.getImage(),
                       cmap='gray',           # Primary image
                       cmap_overlay='hot')    # Overlay
plotter.show()
```

---

## Vector Field Plotting

### Interactive Component Selection

```python
from pyable_eros_montin import Vectorable

vf = Vectorable('displacement_field.mha')

# Show GUI with radio button component selector
# User can select X, Y, Z components or Magnitude
vf.plotOverlay()
```

**GUI Features:**
- X, Y, Z component buttons
- Magnitude button
- Slice slider for 3D
- Live update on selection

### View Specific Component

```python
vf = Vectorable('velocity_field.mha')

# Show X component
vf.plotOverlay(component=0)

# Show Y component  
vf.plotOverlay(component=1)

# Show Z component
vf.plotOverlay(component=2)
```

### Vector Statistics

```python
vf = Vectorable('displacement.mha')

# Get statistics
stats = vf.getVectorStatistics()
print(f"Mean displacement: {stats['mean']:.2f} mm")
print(f"Max displacement: {stats['max']:.2f} mm")

# Extract magnitude component
magnitude = vf.getMagnitude()
magnitude.plotOverlay(title="Displacement Magnitude")
```

### Overlay on Vector Component

```python
vf = Vectorable('deformation.mha')
template = Imaginable('template.nii.gz')

# Show X component with anatomical overlay
vf.plotOverlay(overlay=template, component=0, alpha=0.6)
```

---

## Time Series Plotting

### Interactive Frame Selection

```python
from pyable_eros_montin import TimeSeriesable

ts = TimeSeriesable('cardiac_cine.nii.gz')

# Show GUI with frame slider
# User can navigate through frames
ts.plotOverlay()
```

**GUI Features:**
- Frame slider (0 to N-1)
- Frame/total display (e.g., "Frame 5/30")
- Slice slider for 3D frames
- Live update on frame change

### View Specific Frame

```python
ts = TimeSeriesable('dynamic_mri.nii.gz')

# Show frame 0 (first)
ts.plotOverlay(frame=0)

# Show middle frame
ts.plotOverlay(frame=ts.getNumberOfFrames() // 2)

# Show last frame
ts.plotOverlay(frame=ts.getNumberOfFrames() - 1)
```

### Temporal Statistics

```python
ts = TimeSeriesable('cardiac_4d.nii.gz')

# Get mean across all frames
mean_frame = ts.getTemporalMean()
mean_frame.plotOverlay(title="Average Cardiac Frame")

# Get variance (motion map)
variance = ts.getTemporalVariance()
variance.plotOverlay(title="Motion Variance")

# Get frame range
systole_frames = ts.getFrameRange(10, 15)  # Frames 10-14
systole_frames.write('systole_cine.nii.gz')
```

### Per-Frame Processing

```python
ts = TimeSeriesable('4d_sequence.nii.gz')

# Get individual frame
frame_0 = ts.getFrame(0)
frame_0.cast('float')
frame_0.applyGaussianSmoothing(sigma=1.0)

# Put processed frame back
ts.setFrame(0, frame_0)

# Save modified sequence
ts.write('processed_4d.nii.gz')
```

---

## Overlay Functionality

### Automatic Geometry Handling

Overlays are automatically resampled to match the primary image:

```python
img = Imaginable('image_50x50x50.nii.gz')
overlay = Imaginable('segmentation_100x100x100.nii.gz')

# Overlay automatically resampled to 50x50x50
img.plotOverlay(overlay=overlay)
```

### Multiple Overlays (Sequential)

```python
img = Imaginable('image.nii.gz')
seg1 = Imaginable('structure1.nii.gz')
seg2 = Imaginable('structure2.nii.gz')

# View first overlay
img.plotOverlay(overlay=seg1, title="Structure 1", alpha=0.6)

# View second overlay
img.plotOverlay(overlay=seg2, title="Structure 2", alpha=0.6)
```

### Combine Overlays (Manual)

```python
from pyable_eros_montin import Imaginable
import numpy as np

img = Imaginable('image.nii.gz')
mask1 = Imaginable('mask1.nii.gz')
mask2 = Imaginable('mask2.nii.gz')

# Combine masks
combined = mask1.getDuplicate()
combined_array = (sitk.GetArrayFromImage(mask1.getImage()) + 
                 sitk.GetArrayFromImage(mask2.getImage()))
combined.setImageFromNumpy(combined_array)

# Show combined
img.plotOverlay(overlay=combined, alpha=0.5)
```

---

## Advanced Features

### Saving Figures

```python
from pyable_eros_montin import ScalarPlotter

img = Imaginable('image.nii.gz')
overlay = Imaginable('segmentation.nii.gz')

plotter = ScalarPlotter(img.getImage(), overlay.getImage())
plotter.show(slice_idx=50, alpha=0.6)

# Save current view
plotter.saveFigure('visualization.png', dpi=300)
plotter.saveFigure('visualization.pdf')
```

### Convenience Function

```python
from pyable_eros_montin import plotOverlay, Imaginable

# Direct function call (not on object)
img = Imaginable('image.nii.gz')
overlay = Imaginable('segmentation.nii.gz')

viewer = plotOverlay(img, overlay=overlay, alpha=0.6)
```

### Custom Figure Sizes

```python
from pyable_eros_montin import ScalarPlotter

img = Imaginable('image.nii.gz')

plotter = ScalarPlotter(img.getImage(), 
                       figsize=(15, 10))  # Width, height in inches
plotter.show()
```

### Chaining with Processing

```python
from pyable_eros_montin import Imaginable

img = Imaginable('raw_image.nii.gz')

# Chain processing then view
img.cast('float')\
   .applyGaussianSmoothing(sigma=2.0)\
   .plotOverlay(title="Smoothed Image")
```

---

## API Reference

### PlotViewer (Base Class)

Base class for all viewers. Handles 2D/3D image display, overlays, and slicing.

**Initialization**
```python
PlotViewer(image, overlay=None, title="Image Viewer", 
          cmap='gray', cmap_overlay='hot', figsize=(10, 8))
```

**Methods**
| Method | Purpose |
|--------|---------|
| `show(slice_idx=None, alpha=0.5)` | Display viewer with matplotlib |
| `saveFigure(path, dpi=150)` | Save current figure to file |

---

### ScalarPlotter

Specialized for scalar (Imaginable) images.

**Initialization**
```python
ScalarPlotter(image, overlay=None, title="Scalar Image Viewer", **kwargs)
```

**Parameters**
- `image`: sitk.Image or Imaginable
- `overlay`: sitk.Image or Imaginable (optional)
- `title`: str - figure title
- `cmap`: str - colormap for primary image
- `cmap_overlay`: str - colormap for overlay

---

### VectorPlotter

Specialized for vector fields (Vectorable) with component selection.

**Initialization**
```python
VectorPlotter(vector_image, overlay=None, 
             title="Vector Field Viewer", 
             show_vectors=False, **kwargs)
```

**Methods**
| Method | Purpose |
|--------|---------|
| `show(component=0, alpha=0.5, slice_idx=None)` | Display with component selector |
| `_update_component(idx)` | Change displayed component |
| `_update_magnitude()` | Switch to magnitude view |

**Example**
```python
vf = Vectorable('displacement.mha')
plotter = VectorPlotter(vf)
plotter.show(component=0)  # Start with X
```

---

### TimeSeriesPlotter

Specialized for 4D time-series (TimeSeriesable) with frame navigation.

**Initialization**
```python
TimeSeriesPlotter(timeseries_image, overlay=None,
                 title="Time Series Viewer", **kwargs)
```

**Methods**
| Method | Purpose |
|--------|---------|
| `show(frame=0, alpha=0.5, slice_idx=None)` | Display with frame slider |
| `_update_frame(idx)` | Change displayed frame |

**Example**
```python
ts = TimeSeriesable('cardiac_4d.nii.gz')
plotter = TimeSeriesPlotter(ts)
plotter.show(frame=0)  # Start at first frame
```

---

### plotOverlay() Function

Convenience function that auto-detects image type and creates appropriate viewer.

**Signature**
```python
plotOverlay(image, overlay=None, alpha=0.5, title=None,
           component=None, frame=None, slice_idx=None,
           show_vectors=False, **kwargs) -> PlotViewer
```

**Parameters**
- `image`: sitk.Image, Imaginable, Vectorable, or TimeSeriesable
- `overlay`: Same types (optional)
- `alpha`: float (0-1) - overlay opacity
- `title`: str - custom figure title
- `component`: int - vector component (0, 1, 2) or None for selector
- `frame`: int - time frame index for time series
- `slice_idx`: int - slice index for 3D images
- `show_vectors`: bool - display vector arrows

**Returns**
- PlotViewer instance (ScalarPlotter, VectorPlotter, or TimeSeriesPlotter)

**Examples**
```python
# Scalar
from pyable_eros_montin import Imaginable, plotOverlay
img = Imaginable('image.nii.gz')
viewer = plotOverlay(img, title="My Image")

# Vector with component
from pyable_eros_montin import Vectorable
vf = Vectorable('displacement.mha')
viewer = plotOverlay(vf, component=1)  # Y component

# Time series with frame
from pyable_eros_montin import TimeSeriesable
ts = TimeSeriesable('cardiac.nii.gz')
viewer = plotOverlay(ts, frame=10)  # Frame 10
```

---

### Method Shortcuts on Classes

All able classes have `plotOverlay()` methods:

```python
# Imaginable.plotOverlay()
img.plotOverlay(overlay=seg, alpha=0.6)

# Vectorable.plotOverlay()
vf.plotOverlay(component=0, overlay=background)

# TimeSeriesable.plotOverlay()
ts.plotOverlay(frame=5, overlay=mask, alpha=0.5)
```

---

## Troubleshooting

### Issue: "FigureCanvasAgg is non-interactive"

This is expected in headless/terminal environments. The plots still display correctly in Jupyter notebooks and interactive Python environments.

### Issue: Overlay not showing

Check that the overlay has the same dimension as primary image:
```python
print(f"Image: {img.getImage().GetSize()}")
print(f"Overlay: {overlay.getImage().GetSize()}")
```

The system will automatically resample, but same dimension is simpler.

### Issue: Component selector not appearing

Vector fields with fewer than 3 components will show fewer radio buttons:
```python
vf = Vectorable('2d_vectors.mha')
print(f"Components: {vf.getNumberOfComponents()}")  # 2
# Radio buttons: X, Y, Magnitude
```

---

## Performance Tips

- **Large images**: Extract specific slice before processing
- **Memory**: Use frame-by-frame processing for large time series
- **Rendering**: Reduce figure size for faster updates

```python
# Extract slice before overlay
ts = TimeSeriesable('large_4d.nii.gz')
frame = ts.getFrame(0)  # 3D
frame_2d = frame.cropImage([0,0,50], [512,512,50])
frame_2d.plotOverlay()
```

---

## See Also

- [Imaginable API](./QUICKREF.md)
- [Vector Fields and Time Series](./VECTOR_AND_TIMESERIES_GUIDE.md)
- [Deformation Workflows](./DEFORMATION_WORKFLOW.md)
- [Matplotlib Documentation](https://matplotlib.org/)

---

**Last Updated**: November 21, 2025  
**Version**: 1.0  
**Status**: Stable
