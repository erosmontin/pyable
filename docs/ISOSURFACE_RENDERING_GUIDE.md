# VTK Isosurface Rendering Guide

## Overview

The `renderIsosurface()` method provides interactive 3D visualization of image data using VTK's marching cubes algorithm. This allows you to render 3D isosurfaces directly from medical images, segmentations, and vector fields.

## Features

- **Continuous Images**: Render isosurfaces at any intensity value
- **ROI/Segmentations**: Automatically render ROI boundaries
- **Vector Fields**: Render magnitude or component isosurfaces
- **4D Images**: Extract and render specific time frames
- **Multi-component Images**: Extract and render specific components
- **Interactive Visualization**: Full VTK interactivity (rotate, zoom, pan)
- **Batch Processing**: Can render without displaying for automated workflows

## Method Signature

```python
def renderIsosurface(self, isosurface_value=None, component_index=0, time_index=0, 
                    color=(1.0, 0.0, 0.0), opacity=1.0, show=True, title=None)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `isosurface_value` | float or None | None | The isovalue for the isosurface. If None, uses mean intensity for continuous images or 0.5 for ROIs |
| `component_index` | int | 0 | For multi-component images, which component to render (0=magnitude for vector fields) |
| `time_index` | int | 0 | For 4D images, which time frame to render |
| `color` | tuple | (1.0, 0.0, 0.0) | RGB color (0-1 range). Default is red |
| `opacity` | float | 1.0 | Surface opacity (0-1), where 1.0 is fully opaque |
| `show` | bool | True | If True, display interactive window; if False, return actor without showing |
| `title` | str | None | Window title. Auto-generated if None |

## Returns

- If `show=True`: Returns the VTK actor after displaying
- If `show=False`: Returns a tuple `(actor, renderer, window)` for further manipulation

## Usage Examples

### Example 1: Continuous Image with Default Isosurface

Render an MRI scan at its mean intensity:

```python
from pyable import SITKImaginable

# Load image
img = SITKImaginable('mri_scan.nii.gz')

# Render at mean intensity (automatic)
img.renderIsosurface()
```

### Example 2: Continuous Image with Custom Value

Render at a specific intensity value:

```python
# Render at intensity 100
img.renderIsosurface(isosurface_value=100, color=(0.0, 1.0, 0.0))
```

### Example 3: ROI/Segmentation Boundary

Render the boundary of a binary segmentation:

```python
from pyable import Roiable

# Load ROI
roi = Roiable('segmentation.nii.gz')

# Render boundary (automatically uses 0.5)
roi.renderIsosurface(color=(0.0, 1.0, 0.0), opacity=0.9)
```

### Example 4: Multiple Isosurfaces with Different Colors

```python
img = SITKImaginable('image.nii.gz')

# Get some statistics for reasonable isovalues
min_val = img.getMinimumValue()
max_val = img.getMaximumValue()
mean_val = img.getMeanValue()

# Render three isosurfaces
actor1, _, _ = img.renderIsosurface(isosurface_value=min_val*0.5, show=False, 
                                     color=(1.0, 0.0, 0.0))
actor2, _, _ = img.renderIsosurface(isosurface_value=mean_val, show=False,
                                     color=(0.0, 1.0, 0.0))
actor3, _, _ = img.renderIsosurface(isosurface_value=max_val*0.8, show=False,
                                     color=(0.0, 0.0, 1.0))
```

### Example 5: Vector Field/Displacement Field

Render the magnitude of a displacement field:

```python
from pyable import Fieldable

# Load displacement field
df = SITKImaginable('displacement_field.nii.gz')  # Contains magnitude

# Render at a displacement magnitude threshold
df.renderIsosurface(isosurface_value=5.0, color=(1.0, 1.0, 0.0))
```

### Example 6: 4D Image - Specific Time Frame

Render a specific time frame from a 4D (3D+time) image:

```python
# Load 4D time series
img_4d = SITKImaginable('cardiac_series.nii.gz')  # 4D image

# Render isosurface at time frame 10
img_4d.renderIsosurface(time_index=10, isosurface_value=100, show=True)
```

### Example 7: Multi-component Image - Specific Component

Extract and render a specific component from a multi-component image:

```python
# Load multi-component image
img_multi = SITKImaginable('multi_component.nii.gz')

# Render component 2
img_multi.renderIsosurface(component_index=2, isosurface_value=50)
```

### Example 8: Batch Processing (No Display)

Generate isosurfaces for processing without interactive display:

```python
img = SITKImaginable('image.nii.gz')

# Create isosurface without displaying
actor, renderer, window = img.renderIsosurface(show=False, 
                                                isosurface_value=100,
                                                color=(1.0, 0.0, 0.0))

# Access VTK objects for further processing
polydata = actor.GetMapper().GetInput()
num_polygons = polydata.GetNumberOfCells()
num_points = polydata.GetNumberOfPoints()

print(f"Generated isosurface with {num_polygons} polygons, {num_points} points")

# Can modify properties
actor.GetProperty().SetOpacity(0.5)
actor.GetProperty().EdgeVisibilityOn()  # Show edges

# Could add to custom renderer for multi-object visualization
```

### Example 9: Combining Multiple Objects in One Visualization

```python
import vtk

# Create base image and ROI
img = SITKImaginable('image.nii.gz')
roi = Roiable('segmentation.nii.gz')

# Create isosurfaces without showing
actor_img, _, _ = img.renderIsosurface(isosurface_value=100, show=False,
                                        color=(1.0, 0.5, 0.5), opacity=0.3)
actor_roi, _, _ = roi.renderIsosurface(show=False, 
                                        color=(0.0, 1.0, 0.0), opacity=0.9)

# Create custom visualization
renderer = vtk.vtkRenderer()
renderer.AddActor(actor_img)
renderer.AddActor(actor_roi)
renderer.SetBackground(0.1, 0.1, 0.1)
renderer.ResetCamera()

# Render to window
window = vtk.vtkRenderWindow()
window.AddRenderer(renderer)
window.SetSize(800, 600)
window.SetWindowName("Image + ROI Visualization")

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)
style = vtk.vtkInteractorStyleTrackballCamera()
interactor.SetInteractorStyle(style)

interactor.Initialize()
window.Render()
interactor.Start()
```

### Example 10: Transparency and Edge Highlighting

```python
img = SITKImaginable('image.nii.gz')

actor, _, _ = img.renderIsosurface(show=False, isosurface_value=100,
                                    color=(1.0, 0.0, 0.0), opacity=0.7)

# Add edge visualization
actor.GetProperty().EdgeVisibilityOn()
actor.GetProperty().SetEdgeColor(0.0, 0.0, 0.0)
actor.GetProperty().SetLineWidth(0.5)

# Render with modified properties
# (would need to manually set up window and renderer)
```

## Interactive Controls

When rendering with `show=True`, the VTK window provides standard interaction:

- **Rotate**: Left mouse button + drag
- **Zoom**: Right mouse button + drag, or scroll wheel
- **Pan**: Middle mouse button + drag
- **Reset Camera**: Press 'r'
- **Wireframe Toggle**: Press 'w'
- **Surface Toggle**: Press 's'

## Color Reference

Common RGB colors (0-1 scale):

| Color | RGB |
|-------|-----|
| Red | (1.0, 0.0, 0.0) |
| Green | (0.0, 1.0, 0.0) |
| Blue | (0.0, 0.0, 1.0) |
| Yellow | (1.0, 1.0, 0.0) |
| Magenta | (1.0, 0.0, 1.0) |
| Cyan | (0.0, 1.0, 1.0) |
| White | (1.0, 1.0, 1.0) |
| Black | (0.0, 0.0, 0.0) |
| Gray | (0.5, 0.5, 0.5) |

## Automatic Isosurface Value Selection

If `isosurface_value=None` (default):

- **Continuous Images**: Uses mean intensity
  ```
  isosurface_value = mean(image_intensity)
  ```

- **ROIs (Roiable)**: Uses 0.5 (boundary between foreground/background)
  ```
  isosurface_value = 0.5
  ```

- **Vector Fields**: Uses mean magnitude
  ```
  isosurface_value = mean(displacement_magnitude)
  ```

## Performance Notes

- **Polygon Count**: Dependent on isosurface value and image size
  - Example: 100×100×100 image at mean intensity ≈ 40K-50K polygons
  - Can be reduced by choosing higher isosurface values
  
- **Memory Usage**: 
  - Main memory: image data (stored in SITK)
  - GPU memory: polydata sent to VTK mapper
  - Typically manageable for 512×512×512 and smaller

- **Computation Time**:
  - Marching cubes: ~0.1-1s depending on image size
  - Rendering: Real-time with modern GPUs

## Limitations

1. **No Interactive Editing**: Isosurface is fixed once generated; create multiple if needed
2. **No Transparency Culling**: Transparent isosurfaces may be slow with complex geometry
3. **Single Isovalue**: Each call renders one isovalue; use multiple calls for multiple surfaces
4. **4D/Multi-component**: Must manually extract frame/component before rendering

## Related Methods

- `plotOverlay()`: 2D slice-based overlay visualization
- `viewInteractive()`: Interactive 2D viewer with multiple overlays
- `extractRepresentativeSlices()`: Extract representative 2D slices
- `getVTKImage()`: Get underlying VTK image object
- `getImageAsNumpy()`: Export as numpy for external visualization

## Requirements

- **VTK**: Must be installed (`pip install vtk`)
- **SimpleITK**: Already required by pyable
- **Display**: X11 or equivalent (not required for `show=False` batch mode)

## Troubleshooting

### "No module named 'vtk'"
```bash
pip install vtk
```

### Window doesn't appear
- Ensure you have X11 display available
- Use `show=False` to generate geometry without displaying
- Check for VTK version compatibility (works with VTK 8.2+)

### Isosurface doesn't look right
- Try different isosurface values
- Use `renderIsosurface(show=False)` to get polygon count and verify geometry was generated
- Check image is loaded correctly with `img.printImageInfo()`

### Very slow rendering
- Reduce polygon count by using a higher isosurface value
- Reduce image size with resampling before rendering
- Check available GPU memory

## See Also

- VTK Marching Cubes: https://vtk.org/doc/nightly/html/classvtkMarchingCubes.html
- VTK Rendering: https://vtk.org/doc/nightly/html/group__Rendering.html
