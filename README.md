# pyable

`pyable` is a SimpleITK-first toolkit for day-to-day medical image work. It wraps scalar images, binary ROIs, multi-label segmentations, vector fields, and 4D time series in chainable Python classes so common tasks stay short and readable.

## What it gives you

- A mutable, chainable API around `SimpleITK.Image`
- Explicit NumPy conversion helpers for `(z, y, x)` and `(x, y, z)` layouts
- ROI- and label-preserving transforms for registration/deformation workflows
- Segmentation cleanup, refinement, overlap metrics, and morphometrics
- Plotting helpers and interactive viewers for scalar, vector, and time-series data
- Utility modules for deformations, segmentation, metrics, VTK conversion, and overlays

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install git+https://github.com/erosmontin/pyable.git@v3
```

Requirements:

- Python `>=3.9`
- `numpy`
- `SimpleITK`
- `matplotlib`
- `scikit-image`
- `vtk`
- `pandas`
- `pynico`

## Quick start

```python
from pyable import SITKImaginable, Roiable, LabelMapable

# Scalar image
img = SITKImaginable("image.nii.gz")
img.rotateImage(angle=10).changeImageSpacing([1.0, 1.0, 1.0])
arr = img.getImageAsNumpy()  # (z, y, x)

# Binary ROI
roi = Roiable("mask.nii.gz")
roi.fillBinaryHoles().keepBiggestObj().warpROI("deformation.mha")
metrics = roi.compareTo("reference_mask.nii.gz")

# Multi-label map
labels = LabelMapable("labels.nii.gz")
bone = labels.extractLabel(1)
priors, class_order = labels.buildPriors(tau=0.8)
```

## Main classes

| Class | Role |
| --- | --- |
| `Imaginable` | Base wrapper for scalar `SimpleITK` images |
| `SITKImaginable` | Thin alias/subclass of `Imaginable` |
| `Roiable` | Binary mask / ROI operations |
| `LabelMapable` | Multi-label segmentation workflows |
| `LabelMapableROI` | Legacy label-map helper built from `Roiable` objects |
| `Fieldable` | Vector/displacement field image wrapper |
| `Vectorable` | Explicit vector-field analysis and visualization |
| `TimeSeriesable` | 4D temporal image workflows |
| `PlotViewer`, `ScalarPlotter`, `VectorPlotter`, `TimeSeriesPlotter`, `GridPlotter` | Static plotting helpers |
| `InteractiveViewer`, `OverlayManager` | Interactive browsing and overlay management |
| `RoiComparison` | Classic ROI overlap/similarity metrics |

## Important behavior

- Most mutating methods return `self`, so chaining is the normal usage pattern.
- `getImageAsNumpy()` in v3 returns arrays in `(z, y, x)` order. Use `getImageAsNumpyXYZ()` if you need the legacy `(x, y, z)` view.
- `Roiable` and `LabelMapable` default to nearest-neighbor resampling so labels stay discrete.
- Geometry-sensitive methods assume images are in the same physical space unless they explicitly resample first.

## Documentation

- [Class and method reference](docs/CLASS_REFERENCE.md)
- [Deformation workflow](docs/DEFORMATION_WORKFLOW.md)
- [Plotting guide](docs/PLOTTING_GUIDE.md)
- [Interactive viewer guide](docs/INTERACTIVE_VIEWER_GUIDE.md)
- [Isosurface rendering guide](docs/ISOSURFACE_RENDERING_GUIDE.md)
- [Vector and time-series guide](docs/VECTOR_AND_TIMESERIES_GUIDE.md)

## Testing

Run the workspace copy of the package, not a globally installed version:

```bash
PYTHONPATH=. pytest -q tests/test_phase2_unit_tests.py
PYTHONPATH=. pytest -q tests/test_phase3_integration_tests.py
PYTHONPATH=. pytest -q tests/test_phase5_deformations.py
PYTHONPATH=. pytest -q tests/test_plotable.py
```

## Citation

If `pyable` supports your work, please cite:

1. Montin, E. et al. "A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology." Medical & Biological Engineering & Computing, 2020. https://doi.org/10.1007/s11517-019-02109-4
2. Cavatorta, C. et al. "Retrospective study of late radiation-induced damages after focal radiotherapy for childhood brain tumors." PLOS ONE, 2021. https://doi.org/10.1371/journal.pone.0247748
