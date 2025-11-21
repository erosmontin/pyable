# D2 Utilities Integration - Quick Start Guide

Successfully integrated three powerful utilities from the `d2/` directory into pyable core!

---

## 1. Extract Representative Slices

**Extract 2D slices from 3D volumes around center-of-gravity.**

### Basic Usage

```python
from pyable_eros_montin import Imaginable

# Load image
img = Imaginable('mri_scan.nii.gz')

# Extract slices from all planes
result = img.extractRepresentativeSlices(
    planes='all',           # 'all' or [0, 1, 2]
    offsets=[-10, 0, 10],   # mm offsets from center
    verbose=True
)

# Access results
slices = result['slices']  # List of 2D numpy arrays
labels = result['plane_names']  # [('sagittal', -10), ('sagittal', 0), ...]
cog = result['center_of_gravity']  # Physical coordinates
```

### Use Cases

- **Quick Preview**: Visualize key slices without loading full viewer
- **Batch Processing**: Extract representative slices from many images
- **ML Input**: Prepare 2D slices for machine learning models
- **QC Reports**: Generate automated slice reports

### Advanced Examples

```python
# Single plane extraction
result = img.extractRepresentativeSlices(planes=[0], offsets=[0])

# Fine-grained offsets
result = img.extractRepresentativeSlices(
    planes='all',
    offsets=[-20, -10, -5, 0, 5, 10, 20]
)

# Combine with grid visualization
result = img.extractRepresentativeSlices(planes='all', offsets=[-5, 0, 5])
labels = [f"{name} ({offset}mm)" for name, offset in result['plane_names']]

# Plot in grid
from pyable_eros_montin import GridPlotter
grid = GridPlotter()
grid.show_grid(result['slices'], titles=labels)
```

---

## 2. Batch Process Image Directories

**Process multiple images with a custom function and aggregate results.**

### Basic Usage

```python
from pyable_eros_montin import processImageDirectory

def analyze_image(img):
    """Processor function that returns analysis dict."""
    result = img.extractRepresentativeSlices(planes='all')
    return {
        'size': img.getImageSize(),
        'spacing': img.getImageSpacing(),
        'num_slices': len(result['slices']),
        'center_of_gravity': result['center_of_gravity']
    }

# Process directory
df = processImageDirectory(
    directory='/data/patient_images',
    processor_func=analyze_image,
    file_pattern='*.nii.gz',
    output_csv='/results/image_analysis.csv',
    recursive=True,
    verbose=True
)

# DataFrame has one row per image
print(df.head())
# Output:
#                           filepath      size     spacing  num_slices
# 0  /data/patient_images/img_001.nii.gz  (50, 50, 50)  (1.0, 1.0, 1.0)  9
# 1  /data/patient_images/img_002.nii.gz  (50, 50, 50)  (1.0, 1.0, 1.0)  9
```

### Use Cases

- **Quality Control**: Validate properties across image datasets
- **Preprocessing**: Extract features from all images
- **Reporting**: Generate CSV reports for statistical analysis
- **Pipeline Monitoring**: Track image processing progress

### Advanced Examples

```python
# Complex processor that extracts multiple metrics
def detailed_processor(img):
    from pyable_eros_montin import plotable
    slices = img.extractRepresentativeSlices(planes=[2])
    arr = img.getImageAsNumpy()
    
    return {
        'mean_intensity': arr.mean(),
        'std_intensity': arr.std(),
        'min_intensity': arr.min(),
        'max_intensity': arr.max(),
        'volume_mm3': img.getImageSize()[0] * img.getImageSize()[1] * 
                       img.getImageSize()[2] * 
                       img.getImageSpacing()[0] * img.getImageSpacing()[1] * 
                       img.getImageSpacing()[2]
    }

# Process with detailed metrics
df = processImageDirectory(
    '/data/images',
    detailed_processor,
    output_csv='/results/metrics.csv'
)

# Filter/sort results
high_variance = df[df['std_intensity'] > 50]
large_volumes = df.sort_values('volume_mm3', ascending=False).head(10)
```

---

## 3. Grid Plotter - Multi-Slice Visualization

**Display multiple 2D slices in a grid layout with custom titles and overlays.**

### Basic Usage

```python
from pyable_eros_montin import GridPlotter
import numpy as np

# Create grid plotter
grid = GridPlotter(figsize=(12, 10))

# Create some test slices
slices = [np.random.rand(50, 50) for _ in range(9)]

# Display in 3×3 grid
grid.show_grid(
    slices,
    rows=3,
    cols=3,
    titles=['Slice A', 'Slice B', 'Slice C', ...],
    cmap='gray',
    vmin=0.2,
    vmax=0.8
)
```

### Use Cases

- **Slice Review**: Quick visual inspection of extracted slices
- **QC Interface**: Side-by-side comparison of multiple images
- **Reports**: Generate figures for publications
- **Overlay Comparison**: Display segmentations on images

### Advanced Examples

```python
# Using with extractRepresentativeSlices
img = Imaginable('brain.nii.gz')
result = img.extractRepresentativeSlices(planes='all', offsets=[-5, 0, 5])

labels = [f"{name} ({offset}mm)" for name, offset in result['plane_names']]

grid = GridPlotter(figsize=(15, 12))
grid.show_grid(
    result['slices'],
    titles=labels,
    cmap='viridis',
    vmin=0.1,
    vmax=0.9
)

# With overlays (segmentation on image)
image_slices = [...]  # 9 2D slices
mask_slices = [...]   # 9 corresponding masks

grid.show_grid(
    image_slices,
    overlays=mask_slices,
    titles=labels,
    cmap_overlay='hot',
    alpha_overlay=0.5
)

# Custom grid size
grid.show_grid(slices, rows=2, cols=4)  # 2×4 grid with 8 slices

# Auto grid calculation
grid.show_grid(slices)  # Auto-calculates best grid size
```

---

## Complete Integration Example

Combining all three utilities in a real-world workflow:

```python
from pyable_eros_montin import Imaginable, processImageDirectory, GridPlotter

# Step 1: Define analysis function using slice extraction
def analyze_and_visualize(img):
    result = img.extractRepresentativeSlices(
        planes='all',
        offsets=[-10, 0, 10],
        verbose=False
    )
    
    return {
        'size': img.getImageSize(),
        'spacing': img.getImageSpacing(),
        'num_valid_slices': len(result['slices']),
        'cog_x': result['center_of_gravity'][0],
        'cog_y': result['center_of_gravity'][1],
        'cog_z': result['center_of_gravity'][2]
    }

# Step 2: Process entire directory
df = processImageDirectory(
    directory='/data/patient_mri',
    processor_func=analyze_and_visualize,
    file_pattern='*.nii.gz',
    output_csv='/results/analysis_report.csv',
    recursive=True,
    verbose=True
)

print(f"Processed {len(df)} images")
print(df.describe())

# Step 3: Visualize selected images
selected_img = Imaginable('/data/patient_mri/patient_001.nii.gz')
result = selected_img.extractRepresentativeSlices(planes='all', offsets=[-10, 0, 10])

labels = [f"{name} ({off}mm)" for name, off in result['plane_names']]

grid = GridPlotter(figsize=(14, 12))
grid.show_grid(
    result['slices'],
    titles=labels,
    cmap='gray'
)

# Export summary
print(df.to_string())
df.to_csv('/results/summary.csv', index=False)
```

---

## API Reference

### `Imaginable.extractRepresentativeSlices()`

```python
result = img.extractRepresentativeSlices(
    planes='all',           # str or list - 'all' or [0, 1, 2]
    offsets=[-10, 0, 10],   # list of int - mm from center
    verbose=False           # bool - print progress
)

# Returns dict with keys:
# - 'slices': list of 2D numpy arrays
# - 'plane_names': list of (plane_name, offset) tuples
# - 'offsets': the offsets parameter
# - 'center_of_gravity': (x, y, z) physical coordinates
# - 'center_of_gravity_index': (i, j, k) array indices
```

### `processImageDirectory()`

```python
df = processImageDirectory(
    directory='/path',                  # str - directory path
    processor_func=my_function,         # callable(Imaginable) -> dict
    file_pattern='*.nii.gz',            # str - glob pattern
    output_csv=None,                    # str or None - CSV path
    recursive=True,                     # bool - recursive search
    verbose=False                       # bool - print progress
)

# Returns: pandas.DataFrame with one row per image
```

### `GridPlotter.show_grid()`

```python
grid = GridPlotter(figsize=(10, 8))
grid.show_grid(
    slices,                 # list of 2D numpy arrays
    rows=None,              # int or None - auto-calculate
    cols=None,              # int or None - auto-calculate
    titles=None,            # list of str or None
    cmap='gray',            # str - colormap name
    vmin=None,              # float or None - auto
    vmax=None,              # float or None - auto
    overlays=None,          # list of 2D arrays or None
    cmap_overlay='hot',     # str - overlay colormap
    alpha_overlay=0.5       # float - overlay opacity (0-1)
)
```

---

## Performance Tips

1. **Batch Processing**: Use `processImageDirectory()` with vectorized processors
   ```python
   # Good - process all images efficiently
   df = processImageDirectory('/data', processor, output_csv='/results.csv')
   
   # Avoid - loading/processing images multiple times
   for file in glob.glob('/data/*.nii.gz'):
       # ...
   ```

2. **Slice Extraction**: Cache results if using multiple times
   ```python
   result = img.extractRepresentativeSlices(planes='all')
   # Reuse result['slices'] instead of extracting again
   ```

3. **Grid Visualization**: Use appropriate figure size
   ```python
   # Many slices
   grid = GridPlotter(figsize=(16, 14))
   
   # Few slices  
   grid = GridPlotter(figsize=(8, 6))
   ```

---

## Troubleshooting

### No slices extracted?
- Check image values - need sufficient intensity variation
- Use `verbose=True` to see what's happening
- Try different `offsets` parameters

### Batch processing slow?
- Use `file_pattern` to limit to specific files
- Consider `recursive=False` for shallow directories
- Processor function should be efficient

### Grid visualization looks wrong?
- Adjust `vmin`, `vmax` for better contrast
- Try different `cmap` options
- Check `figsize` for resolution

---

## Related Documentation

- [PLOTTING_GUIDE.md](PLOTTING_GUIDE.md) - Full plotting system documentation
- [D2_UTILITIES_ANALYSIS.md](D2_UTILITIES_ANALYSIS.md) - Technical analysis of d2 utilities
- [VECTOR_AND_TIMESERIES_GUIDE.md](docs/VECTOR_AND_TIMESERIES_GUIDE.md) - Vector field and time-series operations

---

**Status**: ✅ Fully implemented and tested (63 tests passing)  
**Requires**: pandas, SimpleITK, numpy, matplotlib  
**Branch**: v3
