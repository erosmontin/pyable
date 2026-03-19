# pyable

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19119236.svg)](https://doi.org/10.5281/zenodo.19119236)

**pyable** is a modern, SimpleITK-first Python toolkit for medical image analysis. It provides a clean, chainable API for working with scalar images, binary ROIs, multi-label segmentations, vector fields, and 4D time series, making complex workflows simple and reproducible.

---

## 🚀 Features

- Chainable, mutable API around `SimpleITK.Image`
- Easy NumPy conversion: `(z, y, x)` and legacy `(x, y, z)` layouts
- ROI- and label-preserving transforms for registration and deformation
- Segmentation cleanup, overlap metrics, morphometrics
- Interactive and static plotting for scalar, vector, and time-series data
- Utilities for deformations, segmentation, metrics, VTK conversion, overlays

---

## 📦 Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install git+https://github.com/erosmontin/pyable.git@v3
```

**Requirements:**

- Python >=3.9
- numpy
- SimpleITK
- matplotlib
- scikit-image
- vtk
- pandas
- pynico

---

## 🏁 Quick Start

```python
from pyable import SITKImaginable, Roiable, LabelMapable

img = SITKImaginable("image.nii.gz")
img.rotateImage(angle=10).changeImageSpacing([1.0, 1.0, 1.0])
arr = img.getImageAsNumpy()  # (z, y, x)

roi = Roiable("mask.nii.gz")
roi.fillBinaryHoles().keepBiggestObj().warpROI("deformation.mha")
metrics = roi.compareTo("reference_mask.nii.gz")

labels = LabelMapable("labels.nii.gz")
bone = labels.extractLabel(1)
priors, class_order = labels.buildPriors(tau=0.8)
```


---

## 🧩 Main Classes

| Class | Description |
| --- | --- |
| `Imaginable` | Base wrapper for scalar `SimpleITK` images |
| `SITKImaginable` | Thin alias/subclass of `Imaginable` |
| `Roiable` | Binary mask / ROI operations |
| `LabelMapable` | Multi-label segmentation workflows |
| `LabelMapableROI` | Legacy label-map helper built from `Roiable` objects |
| `Fieldable` | Vector/displacement field image wrapper |
| `Vectorable` | Vector-field analysis and visualization |
| `TimeSeriesable` | 4D temporal image workflows |
| `PlotViewer`, `ScalarPlotter`, `VectorPlotter`, `TimeSeriesPlotter`, `GridPlotter` | Static plotting helpers |
| `InteractiveViewer`, `OverlayManager` | Interactive browsing and overlay management |
| `RoiComparison` | Classic ROI overlap/similarity metrics |


---

## ⚡️ Usage Notes

- Most mutating methods return `self` for easy chaining.
- `getImageAsNumpy()` returns arrays in `(z, y, x)` order (v3+). Use `getImageAsNumpyXYZ()` for legacy `(x, y, z)`.
- `Roiable` and `LabelMapable` use nearest-neighbor resampling by default to preserve labels.
- Geometry-sensitive methods assume images are in the same physical space unless resampled.


---

## 📚 Documentation

- [Class and method reference](docs/CLASS_REFERENCE.md)
- [Deformation workflow](docs/DEFORMATION_WORKFLOW.md)
- [Plotting guide](docs/PLOTTING_GUIDE.md)
- [Interactive viewer guide](docs/INTERACTIVE_VIEWER_GUIDE.md)
- [Isosurface rendering guide](docs/ISOSURFACE_RENDERING_GUIDE.md)
- [Vector and time-series guide](docs/VECTOR_AND_TIMESERIES_GUIDE.md)


---

## 🧪 Testing

Run the workspace copy of the package (not a globally installed version):

```bash
PYTHONPATH=. pytest -q tests/test_phase2_unit_tests.py
PYTHONPATH=. pytest -q tests/test_phase3_integration_tests.py
PYTHONPATH=. pytest -q tests/test_phase5_deformations.py
PYTHONPATH=. pytest -q tests/test_plotable.py
```


---

## 📖 Citation

If `pyable` supports your work, please cite:

1. Montin, E. et al. "A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology." Medical & Biological Engineering & Computing, 2020. https://doi.org/10.1007/s11517-019-02109-4
2. Cavatorta, C. et al. "Retrospective study of late radiation-induced damages after focal radiotherapy for childhood brain tumors." PLOS ONE, 2021. https://doi.org/10.1371/journal.pone.0247748
