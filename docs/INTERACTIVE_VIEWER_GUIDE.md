# Interactive GUI Viewer - Complete Guide

## Overview

The `InteractiveViewer` provides a comprehensive GUI for medical image visualization with advanced controls for all image types (scalar, vector, time-series). It supports:

- **Orientation Selection**: Axial, Sagittal, Coronal views
- **Slice Navigation**: Interactive slider with automatic center detection
- **Multiple Overlays**: Stack multiple images with individual opacity control
- **Vector Components**: X, Y, Z, and Magnitude display for vector fields
- **Time Frame Selection**: Navigate through 4D time-series
- **Real-time Updates**: Smooth interactive experience

---

## Quick Start

### Basic Usage - Scalar Image

```python
from pyable import Imaginable

# Load image
img = Imaginable('image.nii.gz')

# Open interactive viewer
img.viewInteractive(orientation=2)  # 0=Axial, 1=Sagittal, 2=Coronal
```

### With Single Overlay

```python
img = Imaginable('image.nii.gz')
seg = Imaginable('segmentation.nii.gz')

# View with segmentation overlay
img.viewInteractive(overlays=seg, orientation=2, slice_idx=50)
```

### With Multiple Overlays

```python
img = Imaginable('image.nii.gz')
seg1 = Imaginable('organ.nii.gz')
seg2 = Imaginable('lesion.nii.gz')

# Stack multiple overlays with opacity control
img.viewInteractive(overlays=[seg1, seg2], orientation=0)
```

### Vector Field with Component Selection

```python
from pyable import Vectorable

vf = Vectorable('displacement_field.mha')
img = Imaginable('reference.nii.gz')

# View with component selector and overlay
vf.viewInteractive(overlays=img, component=0, orientation=2)
```

### Time Series with Frame Navigation

```python
from pyable import TimeSeriesable

ts = TimeSeriesable('cardiac_4d.nii.gz')
roi = Imaginable('roi_mask.nii.gz')

# View with frame slider and overlay
ts.viewInteractive(overlays=roi, frame=5, orientation=1)
```

---

## API Reference

### Imaginable.viewInteractive()

```python
img.viewInteractive(
    overlays=None,           # Single overlay or list
    orientation=2,           # 0=Axial, 1=Sagittal, 2=Coronal
    slice_idx=None,         # Slice index (auto-center if None)
    title=None,             # Window title
    figsize=(14, 10),       # Figure size
    cmap='gray'             # Colormap
)
```

**Parameters:**
- `overlays`: `Imaginable`, `np.ndarray`, or `list` - Overlay image(s)
- `orientation`: `int` - Initial viewing plane
- `slice_idx`: `int` - Initial slice (None = center)
- `title`: `str` - Custom window title
- `figsize`: `tuple` - Figure dimensions (width, height)
- `cmap`: `str` - Colormap name ('gray', 'hot', 'viridis', etc.)

**Returns:** `InteractiveViewer` instance

**Examples:**

```python
# ROI with specific slice
roi = Imaginable('roi.nii.gz')
seg = Imaginable('segmentation.nii.gz')
roi.viewInteractive(overlays=seg, slice_idx=100, title="ROI Viewer")

# Label map with multiple overlays
labels = Imaginable('labels.nii.gz')
mask1 = Imaginable('mask1.nii.gz')
mask2 = Imaginable('mask2.nii.gz')
labels.viewInteractive(overlays=[mask1, mask2], orientation=0)
```

---

### Roiable.viewInteractive()

```python
roi = Roiable('roi_mask.nii.gz')
img = Imaginable('image.nii.gz')

# Identical to Imaginable (inherits all functionality)
roi.viewInteractive(overlays=img, orientation=1)
```

---

### LabelMapable.viewInteractive()

```python
labels = LabelMapable('anatomy.nii.gz')
ref = Imaginable('reference.nii.gz')

# View label map with reference image
labels.viewInteractive(overlays=ref, slice_idx=75)
```

---

### Vectorable.viewInteractive()

```python
vf.viewInteractive(
    overlays=None,           # Optional overlay
    orientation=2,
    slice_idx=None,
    component=None,         # 0=X, 1=Y, 2=Z (new!)
    title=None,
    figsize=(14, 10),
    cmap='gray'
)
```

**New Parameters:**
- `component`: `int` - Initial vector component (0=X, 1=Y, 2=Z)

**Features:**
- Radio buttons for component selection
- Real-time component updates
- Magnitude computation

**Example:**

```python
vf = Vectorable('displacement.mha')
ref = Imaginable('reference.nii.gz')

# View Y component with overlay
vf.viewInteractive(overlays=ref, component=1, orientation=0)
```

---

### TimeSeriesable.viewInteractive()

```python
ts.viewInteractive(
    overlays=None,           # Optional overlay
    orientation=2,
    slice_idx=None,
    frame=None,             # Initial frame (new!)
    title=None,
    figsize=(14, 10),
    cmap='gray'
)
```

**New Parameters:**
- `frame`: `int` - Initial time frame to display

**Features:**
- Slider for frame navigation
- Real-time frame updates
- Shows current frame in title

**Example:**

```python
ts = TimeSeriesable('cardiac_4d.nii.gz')
roi = Imaginable('roi.nii.gz')

# Start at frame 10
ts.viewInteractive(overlays=roi, frame=10, orientation=2)
```

---

## Interactive Controls

### Main Display Window

**Left Panel:** Main image display
- Shows current slice in selected orientation
- All overlays rendered on top

**Right Panel:** Control elements

### Orientation Selection

```
Radio buttons: [○] Axial  [○] Sagittal  [○] Coronal
```
- Click to switch viewing plane
- Automatically updates slice to center

### Slice Navigation

```
Slice: |————●————————|  (0 - max_slice)
```
- Click or drag slider to change slice
- Real-time display update

### Frame Navigation (Time Series Only)

```
Frame: |————●————————|  (0 - num_frames-1)
```
- Navigate through time steps
- Current frame shown in title

### Component Selection (Vector Fields Only)

```
Radio buttons: [○] X  [○] Y  [○] Z  [○] Magnitude
```
- Choose vector component to display
- Magnitude automatically computed

### Overlay Controls

**Visibility Checkboxes:**
```
Overlays:
 ☑ Overlay_0
 ☑ Overlay_1
 ☐ Overlay_2
```
- Toggle overlay visibility

**Opacity Sliders:**
```
Overlay_0 opacity: |————●————————|  (0.0 - 1.0)
Overlay_1 opacity: |——●──————────|  (0.0 - 1.0)
```
- Individual opacity control (0 = transparent, 1 = opaque)

---

## Advanced Usage

### Multiple Overlays with Different Opacities

```python
img = Imaginable('image.nii.gz')
seg1 = Imaginable('organ.nii.gz')
seg2 = Imaginable('lesion.nii.gz')
seg3 = Imaginable('abnormality.nii.gz')

# Stack overlays
viewer = img.viewInteractive(
    overlays=[seg1, seg2, seg3],
    title="Multi-Overlay Viewer"
)

# Use GUI to adjust individual opacities
```

### Numpy Array Overlays

```python
import numpy as np

img = Imaginable('image.nii.gz')

# Create custom mask
mask = np.zeros((256, 256, 256))
mask[100:150, 100:150, 100:150] = 1

# View with numpy overlay
img.viewInteractive(overlays=mask, slice_idx=125)
```

### Interactive Comparison

```python
# Compare two segmentations
seg1 = Imaginable('algorithm_v1.nii.gz')
seg2 = Imaginable('algorithm_v2.nii.gz')

# View both against reference
ref = Imaginable('reference.nii.gz')
ref.viewInteractive(overlays=[seg1, seg2], orientation=0)
```

### Vector Field with Reference

```python
# Displacement field visualization
displacement = Vectorable('deformation.mha')
reference = Imaginable('template.nii.gz')

# View each component with reference
displacement.viewInteractive(
    overlays=reference,
    component=0,  # X component
    orientation=2
)
```

### 4D Time Series Analysis

```python
# Cardiac imaging workflow
cardiac = TimeSeriesable('cardiac_4d.nii.gz')
myocardium = Imaginable('myocardium_mask.nii.gz')

# Navigate through cardiac phases
viewer = cardiac.viewInteractive(
    overlays=myocardium,
    frame=0,
    orientation=1,
    title="Cardiac Function Analysis"
)
```

---

## InteractiveViewer Class (Advanced)

For advanced users, the `InteractiveViewer` class can be used directly:

### Basic Usage

```python
from pyable import InteractiveViewer
import SimpleITK as sitk

# Create from SimpleITK image
img_sitk = sitk.ReadImage('image.nii.gz')
viewer = InteractiveViewer(img_sitk, title="Direct Viewer")

# Add overlays
viewer.add_overlay(sitk.ReadImage('overlay.nii.gz'), name='Overlay_1')
viewer.add_overlay(np.random.rand(256, 256, 256), name='Mask')

# Display
viewer.show()
```

### OverlayManager (Direct Control)

```python
from pyable import OverlayManager

manager = OverlayManager()

# Add overlays
manager.add('organ', organ_image, cmap='hot', opacity=0.5)
manager.add('lesion', lesion_image, cmap='cool', opacity=0.3)

# Control visibility
manager.set_visible('organ', True)
manager.set_visible('lesion', False)

# Adjust opacity
manager.set_opacity('organ', 0.7)

# Query
count = manager.get_count()
opacity = manager.get_opacity('organ')
is_visible = manager.is_visible('lesion')
```

---

## Use Cases

### Clinical Review

```python
# QC of segmentation
img = Imaginable('patient_mri.nii.gz')
auto_seg = Imaginable('auto_segmentation.nii.gz')
manual_seg = Imaginable('manual_segmentation.nii.gz')

img.viewInteractive(
    overlays=[auto_seg, manual_seg],
    orientation=0,
    title="Segmentation QC"
)
```

### Research Analysis

```python
# Compare registration algorithms
ref = Imaginable('reference.nii.gz')
moving = Imaginable('moving.nii.gz')
reg_v1 = Imaginable('registration_v1.nii.gz')
reg_v2 = Imaginable('registration_v2.nii.gz')

ref.viewInteractive(
    overlays=[moving, reg_v1, reg_v2],
    orientation=2
)
```

### Vector Field Inspection

```python
# Deformation analysis
deformation = Vectorable('deformation_field.mha')
template = Imaginable('template.nii.gz')

# Inspect each component
for component in [0, 1, 2]:
    comp_names = ['X', 'Y', 'Z']
    deformation.viewInteractive(
        overlays=template,
        component=component,
        title=f"Deformation {comp_names[component]} Component"
    )
```

### Temporal Analysis

```python
# Cardiac motion
cardiac = TimeSeriesable('cardiac_4d.nii.gz')
myocardium = Imaginable('myocardium_mask.nii.gz')

cardiac.viewInteractive(
    overlays=myocardium,
    orientation=0,
    title="Cardiac Motion Analysis"
)
```

---

## Tips & Tricks

### Performance

- **Large Images**: Use `slice_idx` to avoid rendering entire volume
- **Many Overlays**: Keep to 2-3 overlays for clarity
- **Low-Spec Hardware**: Use smaller `figsize`

### Visualization

- **Contrast**: Use different colormaps for overlays ('hot', 'cool', 'viridis')
- **Transparency**: Set first overlay to 0.5, others to 0.3
- **Focus**: Hide unnecessary overlays using checkboxes

### Workflow

- Start with **Axial** view for overview
- Switch to **Sagittal/Coronal** for 3D understanding
- Use **multiple overlays** for comparison
- Adjust **opacity dynamically** while viewing

---

## Troubleshooting

### Overlays not showing

```python
# Make sure overlay is visible in GUI
# Check checkboxes in "Overlays" panel
# Verify overlay dimensions match primary image
```

### Slow performance

```python
# Reduce image size before viewing
img.resampleImage([2, 2, 2])  # 2x downsampling
img.viewInteractive()

# Or use smaller figure size
img.viewInteractive(figsize=(10, 8))
```

### Wrong orientation

```python
# Change orientation parameter
img.viewInteractive(orientation=0)  # Axial
img.viewInteractive(orientation=1)  # Sagittal
img.viewInteractive(orientation=2)  # Coronal
```

### Vector component not showing

```python
# Use Vectorable, not Imaginable
vf = Vectorable('displacement.mha')  # Not Imaginable!
vf.viewInteractive(component=0)
```

---

## API Summary

| Class | Method | New Features |
|-------|--------|--------------|
| `Imaginable` | `viewInteractive()` | Multi-overlay, orientation, slice nav |
| `Roiable` | `viewInteractive()` | Inherited from Imaginable |
| `LabelMapable` | `viewInteractive()` | Inherited from Imaginable |
| `Vectorable` | `viewInteractive()` | Component selection, magnitude |
| `TimeSeriesable` | `viewInteractive()` | Frame navigation, time slider |
| `InteractiveViewer` | Direct usage | Advanced overlay control |
| `OverlayManager` | Direct usage | Programmatic overlay management |

---

## Examples Repository

See complete working examples in:
- `docs/examples/interactive_viewer_*.py`
- `tests/test_interactive_viewer.py`

---

**Version**: 3.1.0  
**Last Updated**: November 2025
